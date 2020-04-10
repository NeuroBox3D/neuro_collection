/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2013-02-05
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

#include "vdcc_bg_userdata.h"


namespace ug{
namespace neuro_collection{


template <typename TDomain>
VDCC_BG_UserData<TDomain>::VDCC_BG_UserData
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: VDCC_BG<TDomain>(fcts, subsets, approx), m_bIsConstData(false)
{
	// nothing to do
}

template <typename TDomain>
VDCC_BG_UserData<TDomain>::VDCC_BG_UserData
(
	const char* fcts,
	const char* subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: VDCC_BG<TDomain>(fcts, subsets, approx), m_bIsConstData(false)
{
	// nothing to do
}

template <typename TDomain>
VDCC_BG_UserData<TDomain>::~VDCC_BG_UserData()
{
	// nothing to do
}

template<typename TDomain>
void VDCC_BG_UserData<TDomain>::set_potential_function(const char* name)
{
	// name must be a valid lua function name conforming to LuaUserNumber specs
	if (LuaUserData<number, dim>::check_callback_returns(name))
	{
		set_potential_function(LuaUserDataFactory<number, dim>::create(name));
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
void VDCC_BG_UserData<TDomain>::set_potential_function(const number value)
{
	// if value is null, smart pointer will point to null
	if (value == 0.0) set_potential_function(SmartPtr<CplUserData<number, dim> >());
	else set_potential_function(make_sp(new ConstUserNumber<dim>(value)));

	m_bIsConstData = true;
}

template<typename TDomain>
void VDCC_BG_UserData<TDomain>::set_potential_function(SmartPtr<CplUserData<number, dim> > spPotFct)
{
	m_spPotential = spPotFct;
	m_bIsConstData = false;
}


template<typename TDomain>
void VDCC_BG_UserData<TDomain>::update_potential(vm_grid_object* elem)
{
	// only work if really necessary
	if (m_bIsConstData && this->m_initiated) return;

	// fill attachments with renewed values
	const typename TDomain::position_type& coords = CalculateCenter(elem, this->m_aaPos);
	number vm = -0.065;
	if (this->m_spPotential.valid())
		(*this->m_spPotential)(vm, coords, this->m_time, this->m_sh->get_subset_index(elem));

	// set membrane potential value
	this->m_aaVm[elem] = vm;
}



// explicit template specializations
#ifdef UG_DIM_1
	template class VDCC_BG_UserData<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class VDCC_BG_UserData<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class VDCC_BG_UserData<Domain3d>;
#endif

} // neuro_collection
} // namespace ug
