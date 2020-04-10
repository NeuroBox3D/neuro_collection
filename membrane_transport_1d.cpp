/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-08-15
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

#include "membrane_transport_1d.h"
#include "bindings/lua/lua_user_data.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_grid/global_attachments.h"  // for GlobalAttachments


namespace ug {
namespace neuro_collection {



template<typename TDomain>
MembraneTransport1d<TDomain>::MembraneTransport1d(const char* subsets, SmartPtr<IMembraneTransporter> mt)
: IElemDisc<TDomain>("", ""), m_radiusFactor(1.0), m_constRadius(1e-6), m_bConstRadiusSet(false),
  m_aDiameter(GlobalAttachments::attachment<ANumber>("diameter")),
  m_spDensityFct(SPNULL), m_spMembraneTransporter(mt), m_currSI(-1)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	this->IElemDisc<TDomain>::set_subsets(subsets);
	this->IElemDisc<TDomain>::set_functions(mt->symb_fcts());
}


template<typename TDomain>
MembraneTransport1d<TDomain>::MembraneTransport1d(const std::vector<std::string>& subsets, SmartPtr<IMembraneTransporter> mt)
: IElemDisc<TDomain>("", ""), m_radiusFactor(1.0), m_constRadius(1e-6), m_bConstRadiusSet(false),
  m_aDiameter(GlobalAttachments::attachment<ANumber>("diameter")),
  m_spDensityFct(SPNULL), m_spMembraneTransporter(mt), m_currSI(-1)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	this->IElemDisc<TDomain>::set_subsets(subsets);
	this->IElemDisc<TDomain>::set_functions(mt->symb_fcts());
}


template<typename TDomain>
MembraneTransport1d<TDomain>::~MembraneTransport1d()
{}


template<typename TDomain>
void MembraneTransport1d<TDomain>::set_density_function(SmartPtr<CplUserData<number,dim> > densityFct)
{
	this->m_spDensityFct = densityFct;
}

template<typename TDomain>
void MembraneTransport1d<TDomain>::set_density_function(const number dens)
{
	set_density_function(make_sp(new ConstUserNumber<dim>(dens)));
}

template<typename TDomain>
void MembraneTransport1d<TDomain>::set_density_function(const char* name)
{
	// name must be a valid lua function name conforming to LuaUserNumber specs
	if (LuaUserData<number, dim>::check_callback_returns(name))
	{
		set_density_function(LuaUserDataFactory<number, dim>::create(name));
		return;
	}

	// no match found
	if (!CheckLuaCallbackName(name))
		UG_THROW("Lua-Callback with name '" << name << "' does not exist.");

	// name exists, but wrong signature
	UG_THROW("Cannot find matching callback signature. Use:\n"
			"Number - Callback\n" << (LuaUserData<number, dim>::signature()) << "\n");
}

template<typename TDomain>
void MembraneTransport1d<TDomain>::set_radius(number r)
{
	m_constRadius = r;
	m_bConstRadiusSet = true;
}

template<typename TDomain>
void MembraneTransport1d<TDomain>::set_radius_factor(number r)
{
	m_radiusFactor = r;
}


template<typename TDomain>
bool MembraneTransport1d<TDomain>::fluxDensityFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxCond& fc
)
{
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// resize flux structure
	fc.flux.resize(n_flux);
	fc.from.resize(n_flux);
	fc.to.resize(n_flux);

	// calculate single-channel flux
	m_spMembraneTransporter->flux(u, e, fc.flux);

	// get density in membrane
	if (!this->m_spDensityFct.valid())
	{
		UG_THROW("No density information available for " << m_spMembraneTransporter->name()
				 << " membrane transport mechanism. Please set using set_density_function().");
	}
	number density;
	(*this->m_spDensityFct)(density, coords, this->time(), si);

	// save currents
	for (size_t i = 0; i < n_flux; ++i)
	{
		fc.flux[i] *= density;
		fc.from[i] = m_spMembraneTransporter->flux_from_to(i).first;
		fc.to[i] = m_spMembraneTransporter->flux_from_to(i).second;
	}

	return true;
}


template<typename TDomain>
bool MembraneTransport1d<TDomain>::fluxDensityDerivFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxDerivCond& fdc
)
{
	const size_t n_dep = m_spMembraneTransporter->n_dependencies();
	const size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// resize flux derivs structure
	fdc.fluxDeriv.resize(n_flux);
	fdc.from.resize(n_flux);
	fdc.to.resize(n_flux);
	for (size_t i = 0; i < n_flux; ++i)
		fdc.fluxDeriv[i].resize(n_dep);

	// calculate single-channel flux derivs
	m_spMembraneTransporter->flux_deriv(u, e, fdc.fluxDeriv);

	// get density in membrane
	if (!this->m_spDensityFct.valid())
	{
		UG_THROW("No density information available for " << m_spMembraneTransporter->name()
				 << " membrane transport mechanism. Please set using set_density_function().");
	}
	number density;
	(*this->m_spDensityFct)(density, coords, this->time(), si);

	// save current derivs
	for (size_t i = 0; i < n_flux; ++i)
	{
		for (size_t j = 0; j < n_dep; ++j)
			fdc.fluxDeriv[i][j].second *= density;
		fdc.from[i] = m_spMembraneTransporter->flux_from_to(i).first;
		fdc.to[i] = m_spMembraneTransporter->flux_from_to(i).second;
	}

	return true;
}



template<typename TDomain>
void MembraneTransport1d<TDomain>::approximation_space_changed()
{
	SmartPtr<MultiGrid> grid = this->approx_space()->domain()->grid();

	// handle diameter attachment
	if (!grid->has_attachment<Vertex>(m_aDiameter))
		grid->attach_to_vertices_dv(m_aDiameter, 2.0*m_constRadius);
	else
	{
		if (m_bConstRadiusSet)
		{
			UG_LOG("HINT: Even though you have explicitly set a constant diameter to the domain\n"
				   "      MembraneTransport1d will use the diameter information attached to the grid\n"
				   "      you specified.\n");
		}
	}

	// this will distribute the attachment values to the whole grid
	m_dah.set_attachment(m_aDiameter);
	m_dah.set_grid(grid);

    m_aaDiameter = Grid::VertexAttachmentAccessor<ANumber>(*grid, m_aDiameter);
}


template<typename TDomain>
void MembraneTransport1d<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// check that Lagrange 1st order
	for (size_t i = 0; i < vLfeID.size(); ++i)
		if (vLfeID[i].type() != LFEID::LAGRANGE || vLfeID[i].order() != 1)
			UG_THROW("MembraneTransport1d: 1st order Lagrange expected.");

	// update assemble functions
	register_assembling_funcs();
}

template<typename TDomain>
bool MembraneTransport1d<TDomain>::
use_hanging() const
{
	return false;
}


template<typename TDomain>
//template <typename TAlgebra>
void MembraneTransport1d<TDomain>::prep_timestep
(
    number future_time,
    number time,
    VectorProxyBase* upb
)
{
	/*
	// TODO: should be (with a TAlgebra-templated version of this method):
	typedef typename TAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;
	m_spMembraneTransporter->prep_timestep(future_time, time, u);
	*/
	m_spMembraneTransporter->prepare_timestep(future_time, time, upb);
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	m_currSI = si;
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
	// update geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try {geo.update(elem, vCornerCoords, &(this->subset_handler()));}
	UG_CATCH_THROW("MembraneTransport1d::prep_elem: "
						"Cannot update Finite Volume Geometry.");
}

// assemble stiffness part of Jacobian
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	const static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type (in order to allow access to attachments)
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	FluxDerivCond fdc;
	const size_t nFct = u.num_fct();
	std::vector<LocalVector::value_type> uAtCorner(nFct);
	for (size_t i = 0; i < fvgeom.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(i);

		// get associated node
		const int co = scv.node_id();

		// get solution at the corner of the scv
		for (size_t fct = 0; fct < nFct; ++fct)
			uAtCorner[fct] = u(fct, co);

		// get corner coordinates
		const MathVector<dim>& cc = scv.global_corner(0);

		if (!fluxDensityDerivFct(uAtCorner, elem, cc, m_currSI, fdc))
			UG_THROW("MembraneTransport1d::add_jac_A_elem:"
				" Call to fluxDensityDerivFct resulted did not succeed.");

		// scale with volume of SCV
		const number scale = scv.volume() * PI * m_aaDiameter[pElem->vertex(co)] * m_radiusFactor;

		// add to Jacobian
		for (size_t j = 0; j < fdc.fluxDeriv.size(); ++j)
		{
			for (size_t k = 0; k < fdc.fluxDeriv[j].size(); ++k)
			{
				if (fdc.from[j] != InnerBoundaryConstants::_IGNORE_)
					J(fdc.from[j], co, fdc.fluxDeriv[j][k].first, co) += fdc.fluxDeriv[j][k].second * scale;
				if (fdc.to[j] != InnerBoundaryConstants::_IGNORE_)
					J(fdc.to[j], co, fdc.fluxDeriv[j][k].first, co)	-= fdc.fluxDeriv[j][k].second * scale;
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// get finite volume geometry
	static TFVGeom& fvgeom = GeomProvider<TFVGeom>::get();

	// cast elem to appropriate type (in order to allow access to attachments)
	TElem* pElem = dynamic_cast<TElem*>(elem);
	if (!pElem) {UG_THROW("Wrong element type.");}

	FluxCond fc;
	const size_t nFct = u.num_fct();
	std::vector<LocalVector::value_type> uAtCorner(nFct);

	// loop boundary Faces
	for (size_t i = 0; i < fvgeom.num_scv(); ++i)
	{
		// get current SCV
		const typename TFVGeom::SCV& scv = fvgeom.scv(i);

		// get associated node
		const int co = scv.node_id();

		// get solution at the corner of the scv
		for (size_t fct = 0; fct < nFct; ++fct)
			uAtCorner[fct] = u(fct, co);

		// get corner coordinates
		const MathVector<dim>& cc = scv.global_corner(0);

		// get flux densities in that node
		if (!fluxDensityFct(uAtCorner, elem, cc, m_currSI, fc))
		{
			UG_THROW("MembraneTransport1d::add_def_A_elem:"
				" Call to fluxDensityFct did not succeed.");
		}

		// scale with volume of SCV
		const number scale = scv.volume() * PI * m_aaDiameter[pElem->vertex(co)] * m_radiusFactor;

		// add to defect
		for (size_t j = 0; j < fc.flux.size(); ++j)
		{
			if (fc.from[j] != InnerBoundaryConstants::_IGNORE_)
				d(fc.from[j], co) += fc.flux[j] * scale;
			if (fc.to[j] != InnerBoundaryConstants::_IGNORE_)
				d(fc.to[j], co) -= fc.flux[j] * scale;
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void MembraneTransport1d<TDomain>::
add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}




template<typename TDomain>
void MembraneTransport1d<TDomain>::register_assembling_funcs()
{
	// register prep_timestep function for all known algebra types
	RegisterPrepTimestepFct<bridge::CompileAlgebraList>(this);

	// register assembling functionality
	typedef RegularEdge TElem;
	typedef FV1Geometry<TElem, dim> TFVGeom;
	typedef MembraneTransport1d<TDomain> T;
	const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFVGeom>);
}


// explicit template specializations
#ifdef UG_DIM_1
	template class MembraneTransport1d<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class MembraneTransport1d<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class MembraneTransport1d<Domain3d>;
#endif


} // end namespace neuro_collection
} // end namespace ug

