/*
 *	Discretization for the RyR calcium channel in the ER membrane
 *
 *  Created on: Nov 12, 2015
 *      Author: marcuskessler
 */

#include "ryr2.h"


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
KAplus(1500.0e12), KBplus(1500.0e12), KCplus(1.75),
KAminus(28.8), KBminus(385.9), KCminus(0.1),
MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
m_time(0.0), m_oldTime(0.0), m_initiated(false)
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
KAplus(1500.0e12), KBplus(1500.0e12), KCplus(1.75),
KAminus(28.8), KBminus(385.9), KCminus(0.1),
MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
m_time(0.0), m_oldTime(0.0), m_initiated(false)
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
	if (m_mg->template has_attachment<side_t>(this->m_aO1))
	UG_THROW("Attachment necessary for RyR channel dynamics "
			 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_aO1);

	if (m_mg->template has_attachment<side_t>(this->m_aO2))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_aO2);

	if (m_mg->template has_attachment<side_t>(this->m_aC1))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_aC1);

	if (m_mg->template has_attachment<side_t>(this->m_aC2))
		UG_THROW("Attachment necessary for RyR channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_aC2);

	m_aaO1 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_aO1);
	m_aaO2 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_aO2);
	m_aaC1 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_aC1);
	m_aaC2 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_aC2);
}


template<typename TDomain>
RyR2<TDomain>::~RyR2()
{
	m_mg->template detach_from<side_t>(this->m_aO1);
	m_mg->template detach_from<side_t>(this->m_aO2);
	m_mg->template detach_from<side_t>(this->m_aC1);
	m_mg->template detach_from<side_t>(this->m_aC2);
};


template <typename TDomain>
void RyR2<TDomain>::prep_timestep(const number time, VectorProxyBase* upb)
{
	// before the first step: initiate to equilibrium
	if (!m_initiated)
	{
		init(time, upb);
		return;
	}

	// get solution u with which to prepare time step (this code only accepts CPUAlgebra type)
	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;

	// for DoF index storage
	std::vector<DoFIndex> dofIndex;

	// update time
	if (time != m_time)
	{
		m_oldTime = m_time;
		m_time = time;
	}
	number dt = m_time - m_oldTime;

	// loop sides and update potential and then gatings
	typedef typename DoFDistribution::traits<side_t>::const_iterator it_type;

	size_t si_sz = m_vSubset.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = m_dd->begin<side_t>(m_vSubset[si]);
		it_type it_end = m_dd->end<side_t>(m_vSubset[si]);

		for (; it != it_end; ++it)
		{
			number& pO2 = m_aaO2[*it];
			number& pO1 = m_aaO1[*it];
			number& pC2 = m_aaC2[*it];
			number& pC1 = m_aaC1[*it];

			// get ca_cyt and ca_er;
			// we suppose our approx space to be 1st order Lagrange (linear, DoFs in the vertices)
			// and interpolate value at the center of the element
			number ca_cyt = 0.0;
			m_dd->dof_indices(*it, _CCYT_, dofIndex, false, true);
			UG_ASSERT(dofIndex.size() > 0, "No DoF found for function " << _CCYT_ << "on element.");
			for (size_t i = 0; i < dofIndex.size(); ++i) ca_cyt += DoFRef(u, dofIndex[i]);
			ca_cyt /= dofIndex.size();

			// I am not sure whether this is correct ...
			// It is definitely incorrect, e.g.: exponents for ca_cyt mixed up
			//number pO1num = (1 + dt * KBminus) * ((-pC2 + 1 + dt * KCminus) * (1 + dt * KAplus * pow(ca_cyt, 3)) - pC1 * (1 + dt * KCminus)) - pO2 * (1 + dt * KCminus) * (1 + dt * KAplus * pow(ca_cyt, 3));
			//number pO1den = (1 + dt * KBminus) * ((1 + dt * KCminus) * (dt * KAminus) + (dt * KCplus + 1 + dt * KCminus) * (1 + dt * KAplus * pow(ca_cyt, 3))) + (dt * KBplus * pow(ca_cyt, 4)) * (1 + dt * KCminus) * (1 + dt * KAplus * pow(ca_cyt, 3));
			//pO1 = pO1num / pO1den;
			//pO2 = (pO2 + dt * KBplus * pow(ca_cyt, 4) * pO1) / (1 + dt * KBminus);
			//pC1 = (pC1 + dt * KAminus * pO1) / (1 + dt * KAplus * pow(ca_cyt, 3));
			//pC2 = (pC2 + dt * KCplus * pO1) / (1 + dt * KCminus);


			// We use backwards Euler here to evolve o1, o2, c1, c2:
			// the following relation must hold:
			//
			//     u_new = u_old + dt * Au
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

			number a1 = dt * KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt;
			number a2 = dt * KAminus;
			number b1 = dt * KBplus * ca_cyt*ca_cyt*ca_cyt;
			number b2 = dt * KBminus;
			number c1 = dt * KCplus;
			number c2 = dt * KCminus;

			pO1 = (1.0 - pO2/(1.0+b2) - pC1/(1.0+a1) - pC2/(1.0+c2))
				/ (1.0 + b1/(1.0+b2)  + a2/(1.0/a1)  + c1/(1.0+c2) );

			pO2 = (pO2 + b1*pO1) / (1.0 + b2);
			pC1 = (pC1 + a2*pO1) / (1.0 + a1);
			pC2 = (pC2 + c1*pO1) / (1.0 + c2);
		}
	}
}


template<typename TDomain>
void RyR2<TDomain>::init(number time, VectorProxyBase* upb)
{
	this->m_time = time;

	// get solution u with which to prepare time step (this code only accepts CPUAlgebra type)
	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;

	// for DoF index storage
	std::vector<DoFIndex> dofIndex;

	typedef typename DoFDistribution::traits<side_t>::const_iterator it_type;

	size_t si_sz = m_vSubset.size();
	for (size_t si = 0; si < si_sz; ++si)
	{
		it_type it = m_dd->begin<side_t>(m_vSubset[si]);
		it_type it_end = m_dd->end<side_t>(m_vSubset[si]);

		for (; it != it_end; ++it)
		{
			number& pO2 = m_aaO2[*it];
			number& pO1 = m_aaO1[*it];
			number& pC2 = m_aaC2[*it];
			number& pC1 = m_aaC1[*it];

			// get ca_cyt and ca_er;
			// we suppose our approx space to be 1st order Lagrange (linear, DoFs in the vertices)
			// and interpolate value at the center of the element
			number ca_cyt = 0.0;
			m_dd->dof_indices(*it, _CCYT_, dofIndex, false, true);
			UG_ASSERT(dofIndex.size() > 0, "No DoF found for function " << _CCYT_ << " on element.");
			for (size_t i = 0; i < dofIndex.size(); ++i) ca_cyt += DoFRef(u, dofIndex[i]);
			ca_cyt /= dofIndex.size();

			// calculate equilibrium
			//pO1 = (1.0 / (1.0 + KCplus/KCminus + (1/KAplus/KAminus*pow(ca_cyt,3)) + KBplus/KBminus*pow(ca_cyt,4)));
			//pO2 = (KBplus/KBminus*pow(ca_cyt,3) / (1.0 + KCplus/KCminus + (1/KAplus/KAminus*pow(ca_cyt,3)) + KBplus/KBminus*pow(ca_cyt,4)));
			//pC1 = ((1/KAplus/KAminus*pow(ca_cyt,3)) / (1.0 + KCplus/KCminus + (1/KAplus/KAminus*pow(ca_cyt,3)) + KBplus/KBminus*pow(ca_cyt,4)));
			//pC2 = (KCplus/KCminus / (1.0 + KCplus/KCminus + (1/KAplus/KAminus*pow(ca_cyt,3)) + KBplus/KBminus*pow(ca_cyt,4)));

			number KA = KAplus/KAminus * ca_cyt*ca_cyt*ca_cyt*ca_cyt;
			number KB = KBplus/KBminus * ca_cyt*ca_cyt*ca_cyt;
			number KC = KCplus/KCminus;

			number denom_inv = 1.0 / (1.0 + KC + 1.0/KA + KB);

			pO1 = denom_inv;
			pO2 = KB * denom_inv;
			pC1 = denom_inv / KA;
			pC2 = KC * denom_inv;
		}
	}

	this->m_initiated = true;
}

template<typename TDomain>
void RyR2<TDomain>::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	// cast to side_type
	side_t* elem = dynamic_cast<side_t*>(e);
	if (!elem) UG_THROW("RyR2 fluxDensityFunction called with the wrong type of element.");

	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration

	// membrane current corresponding to diffusion pressure
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);

	// open probability
	number pOpen = m_aaO1[elem] + m_aaO2[elem];
	//number pOpen =  (1.0 + KB*pow(caCyt,3)) / (1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3));

	flux[0] = pOpen * current;
}


template<typename TDomain>
void RyR2<TDomain>::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// cast to side_type
	side_t* elem = dynamic_cast<side_t*>(e);
	if (!elem) UG_THROW("RyR2 fluxDensityFunction called with the wrong type of element.");

	number deriv_value = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (m_aaO1[elem] + m_aaO2[elem]);

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



