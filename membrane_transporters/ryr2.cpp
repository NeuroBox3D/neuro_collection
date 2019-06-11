/*
 *	Discretization for the RyR calcium channel in the ER membrane
 *
 *  Created on: Nov 12, 2015
 *      Author: marcuskessler, mbreit
 */

#include "ryr2.h"

#include "lib_grid/algorithms/debug_util.h"   // for ElementDebugInfo
#include "lib_grid/grid/grid_base_objects.h"  // for VERTEX ...
#include "lib_grid/tools/surface_view.h"      // for MG_ALL

namespace ug{
namespace neuro_collection{

template<typename TDomain>
RyR2<TDomain>::
RyR2
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
KAminus(28.8), KBminus(385.9), KCminus(0.1),
MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
m_time(0.0), m_initiated(false)
{
	construct(subsets, approx);
}

template<typename TDomain>
RyR2<TDomain>::
RyR2
(
	const char* fcts,
	const char* subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
KAminus(28.8), KBminus(385.9), KCminus(0.1),
MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
m_time(0.0), m_initiated(false)
{
	construct(TokenizeString(subsets), approx);
}


template<typename TDomain>
void RyR2<TDomain>::construct
(
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
{
	// save underlying domain
	m_dom = approx->domain();

	// save underlying multigrid
	m_mg = m_dom->grid();

	// save underlying surface dof distribution
	// TODO: think about: do we also need level DDs in the multi-grid context?
	//       after all, the derivatives in the level matrices depend on pOpen!
	//       the surface attachments should be copied down
	m_dd = approx->dof_distribution(GridLevel(), true);


// process subsets
	std::vector<std::string> vsSubset(subsets);

	//	remove white space
	for (size_t i = 0; i < vsSubset.size(); ++i)
		RemoveWhitespaceFromString(vsSubset[i]);

	//	if no subset passed, clear subsets
	if (vsSubset.size() == 1 && vsSubset[0].empty())
		vsSubset.clear();

	//	if subsets passed with separator, but not all tokens filled, throw error
	for (size_t i = 0; i < vsSubset.size(); ++i)
	{
		if (vsSubset[i].empty())
		{
			UG_THROW("Error while setting subsets in " << name() << ": passed "
					 "subset string lacks a subset specification at position "
					 << i << "(of " << vsSubset.size()-1 << ")");
		}
	}

	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(m_dom->subset_handler(), vsSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	for (std::size_t si = 0; si < ssGrp.size(); si++)
		m_vSubset.push_back(ssGrp[si]);


// manage attachments
	if (m_mg->template has_attachment<Vertex>(this->m_aO2))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be created, since it already exists.");
	m_mg->template attach_to<Vertex>(this->m_aO2);

	if (m_mg->template has_attachment<Vertex>(this->m_aC1))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be created, since it already exists.");
	m_mg->template attach_to<Vertex>(this->m_aC1);

	if (m_mg->template has_attachment<Vertex>(this->m_aC2))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be created, since it already exists.");
	m_mg->template attach_to<Vertex>(this->m_aC2);

	if (m_mg->template has_attachment<Vertex>(this->m_aOavg))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be created, since it already exists.");
	m_mg->template attach_to<Vertex>(this->m_aOavg);

	if (m_mg->template has_attachment<Vertex>(this->m_aCaOld))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be created, since it already exists.");
	m_mg->template attach_to<Vertex>(this->m_aCaOld);

	m_aaO2 = Grid::AttachmentAccessor<Vertex, ADouble>(*m_mg, m_aO2);
	m_aaC1 = Grid::AttachmentAccessor<Vertex, ADouble>(*m_mg, m_aC1);
	m_aaC2 = Grid::AttachmentAccessor<Vertex, ADouble>(*m_mg, m_aC2);
	m_aaOavg = Grid::AttachmentAccessor<Vertex, ANumber>(*m_mg, m_aOavg);
	m_aaCaOld = Grid::AttachmentAccessor<Vertex, ANumber>(*m_mg, m_aCaOld);
}


template<typename TDomain>
RyR2<TDomain>::~RyR2()
{
	m_mg->template detach_from<Vertex>(this->m_aO2);
	m_mg->template detach_from<Vertex>(this->m_aC1);
	m_mg->template detach_from<Vertex>(this->m_aC2);
	m_mg->template detach_from<Vertex>(this->m_aOavg);
	m_mg->template detach_from<Vertex>(this->m_aCaOld);
};


template <typename TDomain>
void RyR2<TDomain>::prepare_timestep(number future_time, const number time, VectorProxyBase* upb)
{
	// before the first step: initiate to equilibrium (or init again; stationary case)
	if (!m_initiated || future_time == m_initTime)
		init(time, upb);

    // get global fct index for ccyt function
    FunctionGroup fctGrp(m_dd->dof_distribution_info());
    fctGrp.add(this->m_vFct);
    size_t ind_ccyt = fctGrp.unique_id(_CCYT_);

	// for DoF index storage
	std::vector<DoFIndex> dofIndex;

	// update time
	m_time = future_time;
	number dt = m_time - time;

	// loop sides and update potential and then gatings
	typedef typename DoFDistribution::traits<Vertex>::const_iterator it_type;

	size_t si_sz = m_vSubset.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = m_dd->begin<Vertex>(m_vSubset[si], SurfaceView::MG_ALL); // shadow rim copy are required!
		it_type it_end = m_dd->end<Vertex>(m_vSubset[si], SurfaceView::MG_ALL);

		for (; it != it_end; ++it)
		{
			number& pO2 = m_aaO2[*it];
			number& pC2 = m_aaC2[*it];
			number& pC1 = m_aaC1[*it];
			number& pOavg = m_aaOavg[*it];
			number pO1 = 1.0 - (pO2 + pC1 + pC2);
			number& pCaOld = m_aaCaOld[*it];

			// get ca_cyt
			// we suppose our approx space to be 1st order Lagrange (linear, DoFs in the vertices)
			// and interpolate value at the center of the element
			number ca_cyt = 0.0;
            if (!this->has_constant_value(_CCYT_, ca_cyt))
            {
                // we suppose our approx space to be 1st order Lagrange (linear, DoFs in the vertices)
                m_dd->dof_indices(*it, ind_ccyt, dofIndex, true, true);
                UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << ind_ccyt
                	<< " in vertex " << ElementDebugInfo(*m_mg, *it));
                ca_cyt = upb->evaluate(dofIndex[0]);
            }
            // else the constant value has been written to ca_cyt by has_constant_value()

            // scale by appropriate factor for correct unit
            ca_cyt *= this->scale_input(_CCYT_);

            // approximate the Ca derivative and retain the calcium value for future time step
            const number caDeriv = (ca_cyt - pCaOld) / dt;
            pCaOld = ca_cyt;


			// We use backwards Euler here to evolve o1, o2, c1, c2:
			// the following relation must hold:
			//
			//     u_new = u_old + dt * Au_new
			//
			// where u_new, u_old are three-component vectors belonging to o2, c1, c2,
			// A is a 3x4 matrix defined by the Markov model and u is the vector
			// (1 o2_new c1_new c2_new)^T.
			// Additionally, we always need o1+o2+c1+c2 = 1, which becomes the first equation.
			// This is equivalent to solving the system:
			//
			//     (  1    1    1    1   )  (o1_new)   (   1  )
			//     ( -b1  1+b2  0    0   )  (o2_new)   (o2_old)
			//     ( -a2   0   1+a1  0   )  (c1_new) = (c1_old)
			//     ( -c1   0    0   1+c2 )  (c2_new)   (c2_old)
			//
			// with coefficients as defined in the code directly below.
			// We solve this by transformation into a lower-left triangular matrix
			// (i.e., solving for o1_new) and then inverting the rest iteratively.
            number inner_dt = dt;
            size_t nSteps = 1;
            const number thresh = 1e-6;
            if (fabs(inner_dt) > thresh)
            {
            	nSteps = (size_t) ceil(fabs(dt) / thresh);
            	inner_dt = dt / nSteps;
            }

            // forward step (implicit)
            if (inner_dt > 0)
            {
            	pOavg = 0.0;
				for (size_t i = 0; i < nSteps; ++i)
				{
					// estimate current calcium
					ca_cyt += caDeriv*inner_dt;

					const number a1 = inner_dt * KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt;
					const number a2 = inner_dt * KAminus;
					const number b1 = inner_dt * KBplus * ca_cyt*ca_cyt*ca_cyt;
					const number b2 = inner_dt * KBminus;
					const number c1 = inner_dt * KCplus;
					const number c2 = inner_dt * KCminus;

					pO1 = (1.0 - pO2/(1.0+b2) - pC1/(1.0+a1) - pC2/(1.0+c2))
						/ (1.0 + b1/(1.0+b2)  + a2/(1.0+a1)  + c1/(1.0+c2));

					pO2 = (pO2 + b1*pO1) / (1.0 + b2);
					pC1 = (pC1 + a2*pO1) / (1.0 + a1);
					pC2 = 1.0 - (pO1 + pO2 + pC1); // make sure sum is 1;

					pOavg += inner_dt * (pO1 + pO2);
				}
				pOavg /= dt;
            }

            // backward step (explicit)
            else if (inner_dt < 0)
            {
            	for (size_t i = 0; i < nSteps; ++i)
				{
            		// TODO: this might be very wrong if the previous forward step was large
					// estimate current calcium
					ca_cyt -= caDeriv*inner_dt;

					const number a1 = inner_dt * KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt;
					const number a2 = inner_dt * KAminus;
					const number b1 = inner_dt * KBplus * ca_cyt*ca_cyt*ca_cyt;
					const number b2 = inner_dt * KBminus;
					const number c1 = inner_dt * KCplus;
					const number c2 = inner_dt * KCminus;

					pC2 += c2 * pC2 - c1 * pO1;
					pC1 += a1 * pC1 - a2 * pO1;
					pO2 += b2 * pO2 - b1 * pO1;
					pO1 = 1.0 - (pO2 + pC1 + pC2);
				}

            	// TODO: How is pOavg to be treated? Atm, take constant pO1+pO2.
            	pOavg = pO1 + pO2;
            }

            // if time step is zero, do nothing
		}
	}
}


template<typename TDomain>
void RyR2<TDomain>::init(number time, VectorProxyBase* upb)
{
	this->m_time = time;
	this->m_initTime = time;

	/*
	// get solution u with which to prepare time step (this code only accepts CPUAlgebra type)
	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;
	*/

	// get global fct index for ccyt function
	FunctionGroup fctGrp(m_dd->dof_distribution_info());
	fctGrp.add(this->m_vFct);
	size_t ind_ccyt = fctGrp.unique_id(_CCYT_);

	// for DoF index storage
	std::vector<DoFIndex> dofIndex;

	typedef typename DoFDistribution::traits<Vertex>::const_iterator it_type;

	size_t si_sz = m_vSubset.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = m_dd->begin<Vertex>(m_vSubset[si], SurfaceView::MG_ALL); // shadow rim copy are required!
		it_type it_end = m_dd->end<Vertex>(m_vSubset[si], SurfaceView::MG_ALL);

		for (; it != it_end; ++it)
		{
			number& pO2 = m_aaO2[*it];
			number& pC2 = m_aaC2[*it];
			number& pC1 = m_aaC1[*it];
			number& pOavg = m_aaOavg[*it];
			number& pCaOld = m_aaCaOld[*it];

			// get ca_cyt
			number ca_cyt = 0.0;
			if (!this->has_constant_value(_CCYT_, ca_cyt))
			{
				// we suppose our approx space to be 1st order Lagrange (linear, DoFs in the vertices)
				m_dd->dof_indices(*it, ind_ccyt, dofIndex, true, true);
                UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << ind_ccyt
                	<< " in vertex " << ElementDebugInfo(*m_mg, *it));
                ca_cyt = upb->evaluate(dofIndex[0]);
			}
			// else the constant value has been written to ca_cyt by has_constant_value()

			// scale by appropriate factor for correct unit
			ca_cyt *= this->scale_input(_CCYT_);
			pCaOld = ca_cyt;

			// calculate equilibrium
			number KA = KAplus/KAminus * ca_cyt*ca_cyt*ca_cyt*ca_cyt;
			number KB = KBplus/KBminus * ca_cyt*ca_cyt*ca_cyt;
			number KC = KCplus/KCminus;

			number denom_inv = 1.0 / (1.0 + KC + 1.0/KA + KB);

			pO2 = KB * denom_inv;
			pC1 = denom_inv / KA;
			pC2 = KC * denom_inv;
			pOavg = 1.0 - (pC1 + pC2);
		}
	}

	m_initiated = true;
}



template<typename TDomain>
template<typename TBaseElem>
number RyR2<TDomain>::open_prob(GridObject* o) const
{
	TBaseElem* e = static_cast<TBaseElem*>(o);
	number pOpen = 0.0;
	const size_t nVrt = e->num_vertices();
	for (size_t v = 0; v < nVrt; ++v)
	{
		Vertex* vrt = e->vertex(v);
//if (m_aaC1[vrt] == 0.0)
//	UG_LOGN("accessing vrt " << ElementDebugInfo(*m_mg, vrt));
		pOpen += m_aaOavg[vrt];
	}
	return pOpen /= nVrt;
}


template<typename TDomain>
number RyR2<TDomain>::open_prob_for_grid_object(GridObject* o) const
{
	switch (o->base_object_id())
	{
		case VERTEX:
		{
			Vertex* vrt = static_cast<Vertex*>(o);
			return m_aaOavg[vrt];
		}
		case EDGE:
			return open_prob<Edge>(o);
		case FACE:
			return open_prob<Face>(o);
		default:
		{
			UG_THROW("Base object id must be VERTEX, EDGE or FACE, but is "
				<< o->base_object_id() << ".");
		}
	}
}


template<typename TDomain>
void RyR2<TDomain>::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration

	// membrane current corresponding to diffusion pressure
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);

	// open probability
	number pOpen = open_prob_for_grid_object(e);

	flux[0] = pOpen * current;

	//UG_COND_THROW(pOpen != pOpen || current != current,
	//	"RyR NaN: pOpen = " << pOpen << ", current = " << current);

	/*
	static size_t cnt = 0;
	if (!cnt)
	{
		UG_LOGN("RyR2 single channel flux: " << flux[0] << ",  pOpen = " << pOpen);
		++cnt;
	}
	*/
}


template<typename TDomain>
void RyR2<TDomain>::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	number pOpen = open_prob_for_grid_object(e);
	number deriv_value = pOpen * R*T/(4*F*F) * MU_RYR/REF_CA_ER;

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = -deriv_value;
		++i;
	}
	if (!has_constant_value(_CER_))
	{
		flux_derivs[0][i].first = local_fct_index(_CER_);
		flux_derivs[0][i].second = deriv_value;
		++i;
	}
}


template<typename TDomain>
size_t RyR2<TDomain>::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_CCYT_))
		n--;
	if (has_constant_value(_CER_))
		n--;

	return n;
}


template<typename TDomain>
size_t RyR2<TDomain>::n_fluxes() const
{
	return 1;
};


template<typename TDomain>
const std::pair<size_t,size_t> RyR2<TDomain>::flux_from_to(size_t flux_i) const
{
    size_t from, to;
    if (allows_flux(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
    if (allows_flux(_CER_)) from = local_fct_index(_CER_); else from = InnerBoundaryConstants::_IGNORE_;

    return std::pair<size_t, size_t>(from, to);
}


template<typename TDomain>
const std::string RyR2<TDomain>::name() const
{
	return std::string("RyR2");
};


template<typename TDomain>
void RyR2<TDomain>::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CER_))
	{
		UG_THROW("Supplying neither cytosolic nor endoplasmic calcium concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}
}


template<typename TDomain>
void RyR2<TDomain>::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_er]   mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


// explicit template specializations
#ifdef UG_DIM_1
	template class RyR2<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class RyR2<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class RyR2<Domain3d>;
#endif


} // namespace neuro_collection
} // namespace ug



