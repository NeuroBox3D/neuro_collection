/*
 * ryr_implicit.cpp
 * Fully implicit RyR implementation.
 *
 * Date:   2017-08-16
 * Author: mbreit
 */

#include "ryr_implicit.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"  // for GeomProvider


namespace ug {
namespace neuro_collection {

template <typename TDomain>
RyRImplicit<TDomain>::
RyRImplicit(const std::vector<std::string>& fcts, const std::vector<std::string>& subsets)
: IMembraneTransporter(fcts),
  IElemDisc<TDomain>(fcts, subsets),
  R(8.314), T(310.0), F(96485.0),
  KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
  KAminus(28.8), KBminus(385.9), KCminus(0.1),
  MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
#if 0
  m_dd(SPNULL),
  m_oldSol(NULL),
  m_nTSteps(1),
#endif
  m_bNonRegularGrid(false),
  m_bCurrElemIsHSlave(false)
{}

template <typename TDomain>
RyRImplicit<TDomain>::RyRImplicit(const char* fcts, const char* subsets)
: IMembraneTransporter(fcts),
  IElemDisc<TDomain>(fcts, subsets),
  R(8.314), T(310.0), F(96485.0),
  KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
  KAminus(28.8), KBminus(385.9), KCminus(0.1),
  MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
#if 0
  m_dd(SPNULL),
  m_oldSol(NULL),
  m_nTSteps(1),
#endif
  m_bNonRegularGrid(false),
  m_bCurrElemIsHSlave(false)
{}


template<typename TDomain>
RyRImplicit<TDomain>::~RyRImplicit()
{}


#if 0
template <typename TDomain>
void RyRImplicit<TDomain>::prep_timestep(number future_time, const number time, VectorProxyBase* upb)
{
	// get solution u with which to prepare time step (this code only accepts CPUAlgebra type)
	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	m_oldSol = &up->m_v;

	number dt = future_time - time;
	m_nTSteps = 1;
	const number thresh = 1e-6;
	if (fabs(dt) > thresh)
		m_nTSteps = (size_t) ceil(fabs(dt) / thresh);
}
#endif

template <typename TDomain>
void RyRImplicit<TDomain>::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration
	number c1 = u[_C1_];
	number c2 = u[_C2_];

	// membrane current corresponding to diffusion pressure
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);

	// open probability
	number pOpen = 1.0 - (c1 + c2);

	// TODO: properly integrate over pOpen (changes faster than time step)
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


template <typename TDomain>
void RyRImplicit<TDomain>::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration
	number c1 = u[_C1_];
	number c2 = u[_C2_];

	number constFactor = R*T/(4*F*F) * MU_RYR/REF_CA_ER;
	number pOpenPart = 1.0 - (c1 + c2);
	number calciumPart = caER - caCyt;
	number deriv_value_ca = pOpenPart * constFactor;
	number deriv_value_states = calciumPart * constFactor;

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = -deriv_value_ca;
		++i;
	}
	if (!has_constant_value(_CER_))
	{
		flux_derivs[0][i].first = local_fct_index(_CER_);
		flux_derivs[0][i].second = deriv_value_ca;
		++i;
	}
	flux_derivs[0][i].first = local_fct_index(_C1_);
	flux_derivs[0][i].second = -deriv_value_states;
	++i;
	flux_derivs[0][i].first = local_fct_index(_C2_);
	flux_derivs[0][i].second = -deriv_value_states;
	++i;
}


template <typename TDomain>
size_t RyRImplicit<TDomain>::n_dependencies() const
{
	size_t n = 4;
	if (has_constant_value(_CCYT_))
		n--;
	if (has_constant_value(_CER_))
		n--;

	return n;
}


template <typename TDomain>
size_t RyRImplicit<TDomain>::n_fluxes() const
{
	return 1;
};


template<typename TDomain>
const std::pair<size_t,size_t> RyRImplicit<TDomain>::flux_from_to(size_t flux_i) const
{
    size_t from, to;
    if (allows_flux(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
    if (allows_flux(_CER_)) from = local_fct_index(_CER_); else from = InnerBoundaryConstants::_IGNORE_;

    return std::pair<size_t, size_t>(from, to);
}


template <typename TDomain>
const std::string RyRImplicit<TDomain>::name() const
{
	return std::string("RyRImplicit");
};


template <typename TDomain>
void RyRImplicit<TDomain>::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CER_))
	{
		UG_THROW("Supplying neither cytosolic nor endoplasmic calcium concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}

	// check that neither state variable is set const
	if (has_constant_value(_O2_) || has_constant_value(_C1_) || has_constant_value(_C2_))
	{
		UG_THROW("None of the RyR channel state variables (O2, C1, C2) can be set const.");
	}
}


template <typename TDomain>
void RyRImplicit<TDomain>::print_units() const
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
	UG_LOG("|      O2        1                                                             |"<< std::endl);
	UG_LOG("|      C1        1                                                             |"<< std::endl);
	UG_LOG("|      C2        1                                                             |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}



template <typename TDomain>
void RyRImplicit<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check that Lagrange 1st order
	for (size_t i = 0; i < vLfeID.size(); ++i)
		if (vLfeID[i].type() != LFEID::LAGRANGE || vLfeID[i].order() != 1)
			UG_THROW("RyRImplicit: 1st order Lagrange functions expected.");

	// update assemble functions
	m_bNonRegularGrid = bNonRegularGrid;
	register_all_fv1_funcs();
}


template <typename TDomain>
bool RyRImplicit<TDomain>::use_hanging() const
{
	return true;
}

#if 0
template <typename TDomain>
void RyRImplicit<TDomain>::approximation_space_changed()
{
	m_dd = this->approx_space()->dof_distribution(GridLevel(), false);

    // get global fct index for ccyt function
    FunctionGroup fctGrp(m_dd->dof_distribution_info());
    fctGrp.add(this->IMembraneTransporter::m_vFct);
    m_globInd[_CCYT_] = fctGrp.unique_id(_CCYT_);
    m_globInd[_CER_] = fctGrp.unique_id(_CER_);
    m_globInd[_O2_] = fctGrp.unique_id(_O2_);
    m_globInd[_C1_] = fctGrp.unique_id(_C1_);
    m_globInd[_C2_] = fctGrp.unique_id(_C2_);
}
#endif

template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::fsh_elem_loop()
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
#ifdef UG_PARALLEL
	DistributedGridManager& dgm = *this->approx_space()->domain()->grid()->distributed_grid_manager();
	m_bCurrElemIsHSlave = dgm.get_status(elem) & ES_H_SLAVE;
#endif

	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("RyRImplicit::prep_elem: Cannot update finite volume geometry.");
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::add_def_A_elem
(
	LocalVector& d,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// strictly speaking, we only need ODE assemblings here,
	// but it does not hurt to integrate over the boxes either
	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom::BF& bf = fvgeom.bf(i);

		// get associated node
		const int co = bf.node_id();
#if 0
		// get associated vertex
		side_t* side = dynamic_cast<side_t*>(elem);
		UG_COND_THROW(!side, "Wrong side type.")
		Vertex* vrt = side->vertex(co);

		// get old solution values
		std::vector<DoFIndex> dofIndex;

		number ca_cyt_old = 0.0;
		if (!this->has_constant_value(_CCYT_, ca_cyt_old))
		{
			m_dd->dof_indices(vrt, m_globInd[_CCYT_], dofIndex, true, true);
	        UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << m_globInd[_CCYT_]
	        	<< " in vertex " << ElementDebugInfo(*this->approx_space()->domain()->grid(), vrt));
	        ca_cyt_old = DoFRef(*m_oldSol, dofIndex[0]);
		} // else the constant value has been written to ca_cyt by has_constant_value()
		ca_cyt_old *= this->scale_input(_CCYT_);

		number o2_old = 0.0;
		m_dd->dof_indices(vrt, m_globInd[_O2_], dofIndex, true, true);
		UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << m_globInd[_O2_]
	       	<< " in vertex " << ElementDebugInfo(*this->approx_space()->domain()->grid(), vrt));
		o2_old = DoFRef(*m_oldSol, dofIndex[0]);

		number c1_old = 0.0;
		m_dd->dof_indices(vrt, m_globInd[_C1_], dofIndex, true, true);
		UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << m_globInd[_C1_]
			<< " in vertex " << ElementDebugInfo(*this->approx_space()->domain()->grid(), vrt));
		c1_old = DoFRef(*m_oldSol, dofIndex[0]);

		number c2_old = 0.0;
		m_dd->dof_indices(vrt, m_globInd[_C2_], dofIndex, true, true);
		UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << m_globInd[_C2_]
			<< " in vertex " << ElementDebugInfo(*this->approx_space()->domain()->grid(), vrt));
		c2_old = DoFRef(*m_oldSol, dofIndex[0]);


		number ca_cyt_new = u(_CCYT_, co) * scale_input(_CCYT_);
		number fac = 1.0 / m_nTSteps;

		for (size_t step = 1; step <= m_nTSteps; ++step)
		{
			number ca_cyt = ca_cyt_old + step * (ca_cyt_new - ca_cyt_old) * fac;
			number o2 = o2_old + step * (u(_O2_, co) - o2_old) * fac;
			number c1 = c1_old + step * (u(_C1_, co) - c1_old) * fac;
			number c2 = c2_old + step * (u(_C2_, co) - c2_old) * fac;
			number o1 = 1.0 - (o2 + c1 + c2);

			d(_O2_, co) -= (KBplus * ca_cyt*ca_cyt*ca_cyt * o1 - KBminus * o2) * fac * bf.volume();
			d(_C1_, co) -= (KAminus * o1 - KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt * c1) * fac * bf.volume();
			d(_C2_, co) -= (KCplus * o1 - KCminus * c2) * fac * bf.volume();
		}
#endif
		number ca_cyt = u(_CCYT_, co) * scale_input(_CCYT_);
		number o1 = 1.0 - (u(_O2_, co) + u(_C1_, co) + u(_C2_, co));

		d(_O2_, co) -= (KBplus * ca_cyt*ca_cyt*ca_cyt * o1 - KBminus * u(_O2_, co)) * bf.volume();
		d(_C1_, co) -= (KAminus * o1 - KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt * u(_C1_, co)) * bf.volume();
		d(_C2_, co) -= (KCplus * o1 - KCminus * u(_C2_, co)) * bf.volume();
	}
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::add_def_M_elem
(
	LocalVector& d,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// strictly speaking, we only need ODE assemblings here,
	// but it does not hurt to integrate over the boxes either
	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom::BF& bf = fvgeom.bf(i);

		// get associated node
		const int co = bf.node_id();

		d(_O2_, co) += u(_O2_, co) * bf.volume();
		d(_C1_, co) += u(_C1_, co) * bf.volume();
		d(_C2_, co) += u(_C2_, co) * bf.volume();
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::add_jac_A_elem
(
	LocalMatrix& J,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// strictly speaking, we only need ODE assemblings here,
	// but it does not hurt to integrate over the boxes either
	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom::BF& bf = fvgeom.bf(i);

		// get associated node
		const int co = bf.node_id();

#if 0
		// get associated vertex
		side_t* side = dynamic_cast<side_t*>(elem);
		UG_COND_THROW(!side, "Wrong side type.")
		Vertex* vrt = side->vertex(co);

		// FIXME: This is simply wrong if level matrices are assembled.
		// get old solution values
		std::vector<DoFIndex> dofIndex;

		number ca_cyt_old = 0.0;
		if (!this->has_constant_value(_CCYT_, ca_cyt_old))
		{
			m_dd->dof_indices(vrt, m_globInd[_CCYT_], dofIndex, true, true);
	        UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << m_globInd[_CCYT_]
	        	<< " in vertex " << ElementDebugInfo(*this->approx_space()->domain()->grid(), vrt));
	        ca_cyt_old = DoFRef(*m_oldSol, dofIndex[0]);
		} // else the constant value has been written to ca_cyt by has_constant_value()
		ca_cyt_old *= this->scale_input(_CCYT_);

		number o2_old = 0.0;
		m_dd->dof_indices(vrt, m_globInd[_O2_], dofIndex, true, true);
		UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << m_globInd[_O2_]
	       	<< " in vertex " << ElementDebugInfo(*this->approx_space()->domain()->grid(), vrt));
		o2_old = DoFRef(*m_oldSol, dofIndex[0]);

		number c1_old = 0.0;
		m_dd->dof_indices(vrt, m_globInd[_C1_], dofIndex, true, true);
		UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << m_globInd[_C1_]
			<< " in vertex " << ElementDebugInfo(*this->approx_space()->domain()->grid(), vrt));
		c1_old = DoFRef(*m_oldSol, dofIndex[0]);

		number c2_old = 0.0;
		m_dd->dof_indices(vrt, m_globInd[_C2_], dofIndex, true, true);
		UG_ASSERT(dofIndex.size() == 1, "Not exactly 1 DoF found for function " << m_globInd[_C2_]
			<< " in vertex " << ElementDebugInfo(*this->approx_space()->domain()->grid(), vrt));
		c2_old = DoFRef(*m_oldSol, dofIndex[0]);


		number ca_cyt_new = u(_CCYT_, co) * scale_input(_CCYT_);
		number fac = 1.0 / m_nTSteps;

		for (size_t step = 1; step <= m_nTSteps; ++step)
		{
			number ca_cyt = ca_cyt_old + step * (ca_cyt_new - ca_cyt_old) * fac;
			number o2 = o2_old + step * (u(_O2_, co) - o2_old) * fac;
			number c1 = c1_old + step * (u(_C1_, co) - c1_old) * fac;
			number c2 = c2_old + step * (u(_C2_, co) - c2_old) * fac;
			number o1 = 1.0 - (o2 + c1 + c2);

			number fac1 = step * fac*fac * bf.volume();

			J(_O2_, co, _CCYT_, co) -= KBplus * 3.0 * ca_cyt*ca_cyt * scale_input(_CCYT_) * o1 * fac1;
			J(_O2_, co, _O2_, co) += (KBminus + KBplus * ca_cyt*ca_cyt*ca_cyt) * fac1;
			J(_O2_, co, _C1_, co) += KBplus * ca_cyt*ca_cyt*ca_cyt * fac1;
			J(_O2_, co, _C2_, co) += KBplus * ca_cyt*ca_cyt*ca_cyt * fac1;

			J(_C1_, co, _CCYT_, co) += KAplus * 4.0 * ca_cyt*ca_cyt*ca_cyt * scale_input(_CCYT_) * c1 * fac1;
			J(_C1_, co, _O2_, co) += KAminus * fac1;
			J(_C1_, co, _C1_, co) += (KAminus + KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt) * fac1;
			J(_C1_, co, _C2_, co) += KAminus * fac1;

			J(_C2_, co, _O2_, co) += KCplus * fac1;
			J(_C2_, co, _C1_, co) += KCplus * fac1;
			J(_C2_, co, _C2_, co) += (KCplus + KCminus) * fac1;
		}
#endif
		number ca_cyt = u(_CCYT_, co) * scale_input(_CCYT_);
		number o1 = 1.0 - (u(_O2_, co) + u(_C1_, co) + u(_C2_, co));

		J(_O2_, co, _CCYT_, co) -= KBplus * 3.0 * ca_cyt*ca_cyt * scale_input(_CCYT_) * o1 * bf.volume();
		J(_O2_, co, _O2_, co) += (KBminus + KBplus * ca_cyt*ca_cyt*ca_cyt) * bf.volume();
		J(_O2_, co, _C1_, co) += KBplus * ca_cyt*ca_cyt*ca_cyt * bf.volume();
		J(_O2_, co, _C2_, co) += KBplus * ca_cyt*ca_cyt*ca_cyt * bf.volume();

		J(_C1_, co, _CCYT_, co) += KAplus * 4.0 * ca_cyt*ca_cyt*ca_cyt * scale_input(_CCYT_) * u(_C1_, co) * bf.volume();
		J(_C1_, co, _O2_, co) += KAminus * bf.volume();
		J(_C1_, co, _C1_, co) += (KAminus + KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt) * bf.volume();
		J(_C1_, co, _C2_, co) += KAminus * bf.volume();

		J(_C2_, co, _O2_, co) += KCplus * bf.volume();
		J(_C2_, co, _C1_, co) += KCplus * bf.volume();
		J(_C2_, co, _C2_, co) += (KCplus + KCminus) * bf.volume();
	}
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::add_jac_M_elem
(
	LocalMatrix& J,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// on horizontal interfaces: only treat hmasters
	if (m_bCurrElemIsHSlave) return;

	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// strictly speaking, we only need ODE assemblings here,
	// but it does not hurt to integrate over the boxes either
	for (size_t i = 0; i < fvgeom.num_bf(); ++i)
	{
		// get current BF
		const typename TFVGeom::BF& bf = fvgeom.bf(i);

		// get associated node
		const int co = bf.node_id();

		J(_O2_, co, _O2_, co) += bf.volume();
		J(_C1_, co, _C1_, co) += bf.volume();
		J(_C2_, co, _C2_, co) += bf.volume();
	}
}



// register for 2D and 3d
template <typename TDomain>
void RyRImplicit<TDomain>::
register_all_fv1_funcs()
{
//	get all grid element types in the dimension below
	typedef typename domain_traits<dim>::ManifoldElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList>(RegisterFV1(this));
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef RyRImplicit<TDomain> T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(		id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	 	id, &T::template add_rhs_elem<TElem, TFVGeom>);
}






template <typename TDomain>
RyRImplicit_1drotsym<TDomain>::
RyRImplicit_1drotsym(const std::vector<std::string>& fcts, const std::vector<std::string>& subsets)
: IElemDisc<TDomain>(fcts, subsets),
  R(8.314), T(310.0), F(96485.0),
  KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
  KAminus(28.8), KBminus(385.9), KCminus(0.1),
  MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
  m_scale_cc(1.0),
  m_bNonRegularGrid(false)
{}

template <typename TDomain>
RyRImplicit_1drotsym<TDomain>::RyRImplicit_1drotsym(const char* fcts, const char* subsets)
: IElemDisc<TDomain>(fcts, subsets),
  R(8.314), T(310.0), F(96485.0),
  KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
  KAminus(28.8), KBminus(385.9), KCminus(0.1),
  MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
  m_scale_cc(1.0),
  m_bNonRegularGrid(false)
{}


template<typename TDomain>
RyRImplicit_1drotsym<TDomain>::~RyRImplicit_1drotsym()
{}


template<typename TDomain>
void RyRImplicit_1drotsym<TDomain>::set_calcium_scale(number scale_cc)
{
	m_scale_cc = scale_cc;
}


template <typename TDomain>
void RyRImplicit_1drotsym<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check that Lagrange 1st order
	for (size_t i = 0; i < vLfeID.size(); ++i)
		if (vLfeID[i].type() != LFEID::LAGRANGE || vLfeID[i].order() != 1)
			UG_THROW("RyRImplicit_1drotsym: 1st order Lagrange functions expected.");

	// update assemble functions
	m_bNonRegularGrid = bNonRegularGrid;
	register_all_fv1_funcs();
}


template <typename TDomain>
bool RyRImplicit_1drotsym<TDomain>::use_hanging() const
{
	return true;
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::fsh_elem_loop()
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::
prep_elem
(
	const LocalVector& u,
	GridObject* elem,
	const ReferenceObjectID roid,
	const MathVector<dim> vCornerCoords[]
)
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("RyRImplicit_1drotsym::prep_elem: Cannot update finite volume geometry.");
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::add_def_A_elem
(
	LocalVector& d,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// strictly speaking, we only need ODE assemblings here,
	// but it does not hurt to integrate over the boxes either
	for (size_t i = 0; i < fvgeom.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(i);

		// get associated node
		const int co = scv.node_id();

		number ca_cyt = u(_CCYT_, co) * m_scale_cc;
		number o1 = 1.0 - (u(_O2_, co) + u(_C1_, co) + u(_C2_, co));

		d(_O2_, co) -= (KBplus * ca_cyt*ca_cyt*ca_cyt * o1 - KBminus * u(_O2_, co)) * scv.volume();
		d(_C1_, co) -= (KAminus * o1 - KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt * u(_C1_, co)) * scv.volume();
		d(_C2_, co) -= (KCplus * o1 - KCminus * u(_C2_, co)) * scv.volume();
	}
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::add_def_M_elem
(
	LocalVector& d,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// strictly speaking, we only need ODE assemblings here,
	// but it does not hurt to integrate over the boxes either
	for (size_t i = 0; i < fvgeom.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(i);

		// get associated node
		const int co = scv.node_id();

		d(_O2_, co) += u(_O2_, co) * scv.volume();
		d(_C1_, co) += u(_C1_, co) * scv.volume();
		d(_C2_, co) += u(_C2_, co) * scv.volume();
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::add_jac_A_elem
(
	LocalMatrix& J,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// strictly speaking, we only need ODE assemblings here,
	// but it does not hurt to integrate over the boxes either
	for (size_t i = 0; i < fvgeom.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(i);

		// get associated node
		const int co = scv.node_id();

		number ca_cyt = u(_CCYT_, co) * m_scale_cc;
		number o1 = 1.0 - (u(_O2_, co) + u(_C1_, co) + u(_C2_, co));

		J(_O2_, co, _CCYT_, co) -= KBplus * 3.0 * ca_cyt*ca_cyt * m_scale_cc * o1 * scv.volume();
		J(_O2_, co, _O2_, co) += (KBminus + KBplus * ca_cyt*ca_cyt*ca_cyt) * scv.volume();
		J(_O2_, co, _C1_, co) += KBplus * ca_cyt*ca_cyt*ca_cyt * scv.volume();
		J(_O2_, co, _C2_, co) += KBplus * ca_cyt*ca_cyt*ca_cyt * scv.volume();

		J(_C1_, co, _CCYT_, co) += KAplus * 4.0 * ca_cyt*ca_cyt*ca_cyt * m_scale_cc * u(_C1_, co) * scv.volume();
		J(_C1_, co, _O2_, co) += KAminus * scv.volume();
		J(_C1_, co, _C1_, co) += (KAminus + KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt) * scv.volume();
		J(_C1_, co, _C2_, co) += KAminus * scv.volume();

		J(_C2_, co, _O2_, co) += KCplus * scv.volume();
		J(_C2_, co, _C1_, co) += KCplus * scv.volume();
		J(_C2_, co, _C2_, co) += (KCplus + KCminus) * scv.volume();
	}
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::add_jac_M_elem
(
	LocalMatrix& J,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// strictly speaking, we only need ODE assemblings here,
	// but it does not hurt to integrate over the boxes either
	for (size_t i = 0; i < fvgeom.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(i);

		// get associated node
		const int co = scv.node_id();

		J(_O2_, co, _O2_, co) += scv.volume();
		J(_C1_, co, _C1_, co) += scv.volume();
		J(_C2_, co, _C2_, co) += scv.volume();
	}
}


template <typename TDomain>
void RyRImplicit_1drotsym<TDomain>::
register_all_fv1_funcs()
{
	if (m_bNonRegularGrid)
		register_fv1_func<RegularEdge, HFV1Geometry<RegularEdge, dim> >();
	else
		register_fv1_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRImplicit_1drotsym<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef RyRImplicit_1drotsym<TDomain> T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(		id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	 	id, &T::template add_rhs_elem<TElem, TFVGeom>);
}




// explicit template specializations
#ifdef UG_DIM_1
	template class RyRImplicit<Domain1d>;
	template class RyRImplicit_1drotsym<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class RyRImplicit<Domain2d>;
	template class RyRImplicit_1drotsym<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class RyRImplicit<Domain3d>;
	template class RyRImplicit_1drotsym<Domain3d>;
#endif


} // namespace neuro_collection
} // namespace ug



