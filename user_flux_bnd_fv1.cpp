#include "user_flux_bnd_fv1.h"


namespace ug {
namespace neuro_collection {

template <typename TDomain>
UserFluxBoundaryFV1<TDomain>::UserFluxBoundaryFV1(const char* functions, const char* subsets)
: FV1InnerBoundaryElemDisc<TDomain>(functions, subsets)
{}


template <typename TDomain>
UserFluxBoundaryFV1<TDomain>::UserFluxBoundaryFV1
(
	const std::vector<std::string>& functions,
	const std::vector<std::string>& subsets
)
: FV1InnerBoundaryElemDisc<TDomain>(functions, subsets)
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


