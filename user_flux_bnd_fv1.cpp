/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2014-01-16
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

#include "user_flux_bnd_fv1.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary_impl.h"


namespace ug {
namespace neuro_collection {

template <typename TDomain>
UserFluxBoundaryFV1<TDomain>::UserFluxBoundaryFV1(const char* functions, const char* subsets)
: base_type(functions, subsets)
{}


template <typename TDomain>
UserFluxBoundaryFV1<TDomain>::UserFluxBoundaryFV1
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
: base_type(functions, subsets)
{}


template <typename TDomain>
void UserFluxBoundaryFV1<TDomain>::set_flux_function(SmartPtr<CplUserData<number, dim> > fluxFct)
{
	m_fluxFct = fluxFct;
}


template <typename TDomain>
void UserFluxBoundaryFV1<TDomain>::set_flux_function(number constValue)
{
	SmartPtr<CplUserData<number, dim> > sp = make_sp(new ConstUserNumber<dim>(constValue));
	set_flux_function(sp);
}


template <typename TDomain>
void UserFluxBoundaryFV1<TDomain>::set_flux_function(const char* name)
{
	// name must be a valid lua function name conforming to LuaUserNumber specs
	if (LuaUserData<number, dim>::check_callback_returns(name))
	{
		set_flux_function(LuaUserDataFactory<number, dim>::create(name));
		return;
	}

	// no match found
	if (!CheckLuaCallbackName(name))
		UG_THROW("Lua-Callback with name '" << name << "' does not exist.");

	// name exists, but wrong signature
	UG_THROW("Cannot find matching callback signature. Use:\n"
			"Number - Callback\n" << (LuaUserData<number, dim>::signature()) << "\n");
}


template <typename TDomain>
bool UserFluxBoundaryFV1<TDomain>::fluxDensityFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxCond& fc
)
{
	number fluxDensity;
	if (m_fluxFct.valid())
		(*m_fluxFct)(fluxDensity, coords, this->time(), si);
	else fluxDensity = 0.0;

	fc.flux.resize(1, 0.0);
	fc.flux[0] = fluxDensity;

	fc.to.resize(1);
	fc.to[0] = 0;

	// if a source is given, then use it; otherwise don't
	fc.from.resize(1);
	if (this->m_vFct.size() > 1)
		fc.from[0] = 1;
	else
		fc.from[0] = InnerBoundaryConstants::_IGNORE_;

	return true;
}


template <typename TDomain>
bool UserFluxBoundaryFV1<TDomain>::fluxDensityDerivFct
(
	const std::vector<LocalVector::value_type>& u,
	GridObject* e,
	const MathVector<dim>& coords,
	int si,
	FluxDerivCond& fdc
)
{
	// nothing to compute here, since the flux does not depend on any variables
	return true;
}


// explicit template specializations
#ifdef UG_DIM_1
	template class UserFluxBoundaryFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class UserFluxBoundaryFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class UserFluxBoundaryFV1<Domain3d>;
#endif



} // namespace neuro_collection
} // namespace ug


