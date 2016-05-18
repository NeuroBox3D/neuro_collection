/*
 * vdcc_bg_userdata.cpp
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
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
void VDCC_BG_UserData<TDomain>::update_potential(side_t* elem)
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
