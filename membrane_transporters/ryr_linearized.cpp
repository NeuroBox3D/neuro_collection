/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-08-16
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "ryr_linearized.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"  // for GeomProvider


namespace ug {
namespace neuro_collection {

template <typename TDomain>
RyRLinearized<TDomain>::
RyRLinearized(const std::vector<std::string>& fcts, const std::vector<std::string>& subsets)
: IMembraneTransporter(fcts),
  IElemDisc<TDomain>(fcts, subsets),
  R(8.314), T(310.0), F(96485.0),
  KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
  KAminus(28.8), KBminus(385.9), KCminus(0.1),
  MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
  m_bNonRegularGrid(false),
  m_bCurrElemIsHSlave(false)
{}

template <typename TDomain>
RyRLinearized<TDomain>::RyRLinearized(const char* fcts, const char* subsets)
: IMembraneTransporter(fcts),
  IElemDisc<TDomain>(fcts, subsets),
  R(8.314), T(310.0), F(96485.0),
  KAplus(1500.0e12), KBplus(1500.0e9), KCplus(1.75),
  KAminus(28.8), KBminus(385.9), KCminus(0.1),
  MU_RYR(5.0e-11), REF_CA_ER(2.5e-1),
  m_bNonRegularGrid(false),
  m_bCurrElemIsHSlave(false)
{}


template<typename TDomain>
RyRLinearized<TDomain>::~RyRLinearized()
{}


template <typename TDomain>
void RyRLinearized<TDomain>::calc_flux
(
	const std::vector<number>& u,
	const std::vector<number>& uOld,
	GridObject* e,
	std::vector<number>& flux
) const
{
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration
	number c1 = uOld[_C1_];
	number c2 = uOld[_C2_];

	// membrane current corresponding to diffusion pressure
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);

	// open probability
	number pOpen = 1.0 - (c1 + c2);

	flux[0] = pOpen * current;

	//UG_COND_THROW(pOpen != pOpen || current != current,
	//	"RyR NaN: pOpen = " << pOpen << ", current = " << current);

	/*
	static size_t cnt = 0;
	if (!cnt)
	{
		UG_LOGN("RyRLinearized single channel flux: " << flux[0] << ",  pOpen = " << pOpen);
		++cnt;
	}
	*/
}


template <typename TDomain>
void RyRLinearized<TDomain>::calc_flux_deriv
(
	const std::vector<number>& u,
	const std::vector<number>& uOld,
	GridObject* e,
	std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs
) const
{
	number c1 = uOld[_C1_];
	number c2 = uOld[_C2_];

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = - (1.0 - (c1 + c2)) * R*T/(4*F*F) * MU_RYR/REF_CA_ER;
		++i;
	}
	if (!has_constant_value(_CER_))
	{
		flux_derivs[0][i].first = local_fct_index(_CER_);
		flux_derivs[0][i].second = (1.0 - (c1 + c2)) * R*T/(4*F*F) * MU_RYR/REF_CA_ER;
		++i;
	}
	flux_derivs[0][i].first = local_fct_index(_C1_);
	flux_derivs[0][i].second = 0.0;
	++i;
	flux_derivs[0][i].first = local_fct_index(_C2_);
	flux_derivs[0][i].second = 0.0;
	++i;
}


template <typename TDomain>
size_t RyRLinearized<TDomain>::n_dependencies() const
{
	size_t n = 4;
	if (has_constant_value(_CCYT_))
		n--;
	if (has_constant_value(_CER_))
		n--;

	return n;
}


template <typename TDomain>
size_t RyRLinearized<TDomain>::n_fluxes() const
{
	return 1;
};


template<typename TDomain>
const std::pair<size_t,size_t> RyRLinearized<TDomain>::flux_from_to(size_t flux_i) const
{
    size_t from, to;
    if (is_supplied(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
    if (is_supplied(_CER_)) from = local_fct_index(_CER_); else from = InnerBoundaryConstants::_IGNORE_;

    return std::pair<size_t, size_t>(from, to);
}


template <typename TDomain>
const std::string RyRLinearized<TDomain>::name() const
{
	return std::string("RyRLinearized");
};


template <typename TDomain>
void RyRLinearized<TDomain>::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!is_supplied(_CCYT_) && !is_supplied(_CER_))
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
void RyRLinearized<TDomain>::print_units() const
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
bool RyRLinearized<TDomain>::needsPreviousSolution() const
{
	return true;
}



template <typename TDomain>
void RyRLinearized<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check that Lagrange 1st order
	for (size_t i = 0; i < vLfeID.size(); ++i)
		if (vLfeID[i].type() != LFEID::LAGRANGE || vLfeID[i].order() != 1)
			UG_THROW("RyRLinearized: 1st order Lagrange functions expected.");

	// update assemble functions
	m_bNonRegularGrid = bNonRegularGrid;
	register_all_fv1_funcs();
}


template <typename TDomain>
bool RyRLinearized<TDomain>::use_hanging() const
{
	return true;
}


template <typename TDomain>
template <typename TAlgebra>
void RyRLinearized<TDomain>::prep_timestep(number future_time, number time, VectorProxyBase* upb)
{
	// the proxy provided here will not survive, but the vector it wraps will,
	// so we re-wrap it in a new proxy for ourselves
	typedef VectorProxy<typename TAlgebra::vector_type> tProxy;
	tProxy* proxy = static_cast<tProxy*>(upb);
	m_spOldSolutionProxy = make_sp(new tProxy(proxy->m_v));
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRLinearized<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRLinearized<TDomain>::fsh_elem_loop()
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRLinearized<TDomain>::
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
	UG_CATCH_THROW("RyRLinearized::prep_elem: Cannot update finite volume geometry.");

	// copy local vector of current solution to get indices and map access right,
	// then copy correct values
	m_locUOld = u;
	const LocalIndices& ind = m_locUOld.get_indices();
	for (size_t fct = 0; fct < m_locUOld.num_all_fct(); ++fct)
	{
		for(size_t dof = 0; dof < m_locUOld.num_all_dof(fct); ++dof)
		{
			const DoFIndex di = ind.multi_index(fct, dof);
			m_locUOld.value(fct, dof) = m_spOldSolutionProxy->evaluate(di);
		}
	}
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRLinearized<TDomain>::add_def_A_elem
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

		number ca_cyt = m_locUOld(_CCYT_, co) * scale_input(_CCYT_);
		number o1_o2 = 1.0 - (u(_O2_, co) + m_locUOld(_C1_, co) + m_locUOld(_C2_, co));
		number o1_c1 = 1.0 - (m_locUOld(_O2_, co) + u(_C1_, co) + m_locUOld(_C2_, co));
		number o1_c2 = 1.0 - (m_locUOld(_O2_, co) + m_locUOld(_C1_, co) + u(_C2_, co));

		d(_O2_, co) -= (KBplus * ca_cyt*ca_cyt*ca_cyt * o1_o2 - KBminus * u(_O2_, co)) * bf.volume();
		d(_C1_, co) -= (KAminus * o1_c1 - KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt * u(_C1_, co)) * bf.volume();
		d(_C2_, co) -= (KCplus * o1_c2 - KCminus * u(_C2_, co)) * bf.volume();
	}
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRLinearized<TDomain>::add_def_M_elem
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
void RyRLinearized<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRLinearized<TDomain>::add_jac_A_elem
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

		number ca_cyt = m_locUOld(_CCYT_, co) * scale_input(_CCYT_);

		J(_O2_, co, _O2_, co) += (KBminus + KBplus * ca_cyt*ca_cyt*ca_cyt) * bf.volume();
		J(_C1_, co, _C1_, co) += (KAminus + KAplus * ca_cyt*ca_cyt*ca_cyt*ca_cyt) * bf.volume();
		J(_C2_, co, _C2_, co) += (KCplus + KCminus) * bf.volume();
	}
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRLinearized<TDomain>::add_jac_M_elem
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
void RyRLinearized<TDomain>::
register_all_fv1_funcs()
{
	// register assembling functions
	typedef typename domain_traits<dim>::ManifoldElemList ElemList;
	boost::mpl::for_each<ElemList>(RegisterFV1(this));

	// register prepare_timestep functions
	RegisterPrepTimestep<bridge::CompileAlgebraList>(this);
}


template <typename TDomain>
template <typename TElem, typename TFVGeom>
void RyRLinearized<TDomain>::
register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef RyRLinearized<TDomain> T;

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
	template class RyRLinearized<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class RyRLinearized<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class RyRLinearized<Domain3d>;
#endif


} // namespace neuro_collection
} // namespace ug



