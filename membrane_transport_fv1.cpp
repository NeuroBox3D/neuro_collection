/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2011-12-20
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

#include "membrane_transport_fv1.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary_impl.h"
#include "bindings/lua/lua_user_data.h"


namespace ug {
namespace neuro_collection {


template<typename TDomain>
MembraneTransportFV1<TDomain>::MembraneTransportFV1(const char* subsets, SmartPtr<IMembraneTransporter> mt)
: base_type(),
  R(8.314), T(310.0), F(96485.0), m_spMembraneTransporter(mt)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	// elem discs (subsets and) functions need only be set after the previous check
	this->IElemDisc<TDomain>::set_subsets(subsets);
	this->IElemDisc<TDomain>::set_functions(mt->symb_fcts());
}

template<typename TDomain>
MembraneTransportFV1<TDomain>::MembraneTransportFV1(const std::vector<std::string>& subsets, SmartPtr<IMembraneTransporter> mt)
: base_type(),
  R(8.314), T(310.0), F(96485.0), m_spMembraneTransporter(mt)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	// elem discs (subsets and) functions need only be set after the previous check
	this->IElemDisc<TDomain>::set_subsets(subsets);
	this->IElemDisc<TDomain>::set_functions(mt->symb_fcts());
}

template<typename TDomain>
MembraneTransportFV1<TDomain>::~MembraneTransportFV1()
{
	// nothing to do
}


template<typename TDomain>
void MembraneTransportFV1<TDomain>::set_density_function(SmartPtr<CplUserData<number,dim> > densityFct)
{
	this->m_spDensityFct = densityFct;
}

template<typename TDomain>
void MembraneTransportFV1<TDomain>::set_density_function(const number dens)
{
	set_density_function(make_sp(new ConstUserNumber<dim>(dens)));
}

template<typename TDomain>
void MembraneTransportFV1<TDomain>::set_density_function(const char* name)
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
void MembraneTransportFV1<TDomain>::set_membrane_transporter(SmartPtr<IMembraneTransporter> mt)
{
	m_spMembraneTransporter = mt;
}


template<typename TDomain>
bool MembraneTransportFV1<TDomain>::fluxDensityFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxCond& fc
)
{
	return fluxDensityFctImpl(
		[&u, this](GridObject* e, FluxCond& fc) {m_spMembraneTransporter->flux(u, e, fc.flux);},
		e, coords, si, fc);
}


template <typename TDomain>
bool MembraneTransportFV1<TDomain>::fluxDensityFct
(
	const std::vector<LocalVector::value_type>& u,
	const std::vector<LocalVector::value_type>& uOld,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxCond& fc
)
{
	return fluxDensityFctImpl(
		[&u, &uOld, this](GridObject* e, FluxCond& fc) {m_spMembraneTransporter->flux(u, uOld, e, fc.flux);},
		e, coords, si, fc);
}


template <typename TDomain>
template <typename CallToMemTransporter>
bool MembraneTransportFV1<TDomain>::fluxDensityFctImpl
(
	const CallToMemTransporter& call,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxCond& fc
)
{
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// calculate single-channel flux
	fc.flux.resize(n_flux);
	fc.from.resize(n_flux);
	fc.to.resize(n_flux);

	// get solution at the corner of the bf
	call(e, fc);

	// get density in membrane
	if (!this->m_spDensityFct.valid())
	{
		UG_THROW("No density information available for " << m_spMembraneTransporter->name()
				 << " membrane transport mechanism. Please set using set_density_function().");
	}
	number density;
	(*this->m_spDensityFct)(density, coords, this->time(), si);

	for (size_t i = 0; i < n_flux; i++)
	{
		fc.flux[i] *= density;
		fc.from[i] = m_spMembraneTransporter->flux_from_to(i).first;
		fc.to[i] = m_spMembraneTransporter->flux_from_to(i).second;
	}

	return true;
}


template<typename TDomain>
bool MembraneTransportFV1<TDomain>::fluxDensityDerivFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxDerivCond& fdc
)
{
	return fluxDensityDerivFctImpl(
		[&u, this](GridObject* e, FluxDerivCond& fdc) {m_spMembraneTransporter->flux_deriv(u, e, fdc.fluxDeriv);},
		e, coords, si, fdc);
}


template<typename TDomain>
bool MembraneTransportFV1<TDomain>::fluxDensityDerivFct
(
	const std::vector<LocalVector::value_type>& u,
	const std::vector<LocalVector::value_type>& uOld,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxDerivCond& fdc
)
{
	return fluxDensityDerivFctImpl(
		[&u, &uOld, this](GridObject* e, FluxDerivCond& fdc) {m_spMembraneTransporter->flux_deriv(u, uOld, e, fdc.fluxDeriv);},
		e, coords, si, fdc);
}


template <typename TDomain>
template <typename CallToMemTransporter>
bool MembraneTransportFV1<TDomain>::fluxDensityDerivFctImpl
(
	const CallToMemTransporter& call,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxDerivCond& fdc
)
{
	size_t n_dep = m_spMembraneTransporter->n_dependencies();
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// calculate single-channel flux
	fdc.fluxDeriv.resize(n_flux);
	fdc.from.resize(n_flux);
	fdc.to.resize(n_flux);
	for (size_t i = 0; i < n_flux; i++)
		fdc.fluxDeriv[i].resize(n_dep);

	call(e, fdc);

	number density;
	if (this->m_spDensityFct.valid())
		(*this->m_spDensityFct)(density, coords, this->time(), si);
	else
	{
		UG_THROW("No density information available for " << m_spMembraneTransporter->name()
				<< " membrane transport mechanism. Please set using set_density_function().");
	}

	for (size_t i = 0; i < n_flux; i++)
	{
		for (size_t j = 0; j < n_dep; j++)
			fdc.fluxDeriv[i][j].second *= density;
		fdc.from[i] = m_spMembraneTransporter->flux_from_to(i).first;
		fdc.to[i] = m_spMembraneTransporter->flux_from_to(i).second;
	}

	return true;
}


template<typename TDomain>
void MembraneTransportFV1<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// set assembling functions from base class first
	this->base_type::prepare_setting(vLfeID, bNonRegularGrid);

	// update assemble functions
	register_all_fv1_funcs();
}


template <typename TDomain>
template <typename TAlgebra>
void MembraneTransportFV1<TDomain>::prep_timestep
(
    number future_time,
    number time,
    VectorProxyBase* upb
)
{
	this->previous_solution_required(m_spMembraneTransporter->needsPreviousSolution());
	this->base_type::template prep_timestep<TAlgebra>(future_time, time, upb);

	m_spMembraneTransporter->prepare_timestep(future_time, time, upb);
}


template<typename TDomain>
void MembraneTransportFV1<TDomain>::register_all_fv1_funcs()
{
	// register prep_timestep function for all known algebra types
	RegisterPrepTimestep<bridge::CompileAlgebraList>(this);
}



// explicit template specializations
#ifdef UG_DIM_1
	template class MembraneTransportFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class MembraneTransportFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class MembraneTransportFV1<Domain3d>;
#endif


} // end namespace neuro_collection
} // end namespace ug

