/*
 *  membrane_transport_fv1.cpp
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#include "membrane_transport_fv1.h"

namespace ug {
namespace neuro_collection {


template<typename TDomain>
MembraneTransportFV1<TDomain>::MembraneTransportFV1(const char* functions, const char* subsets)
: FV1InnerBoundaryElemDisc<TDomain>(functions, subsets),
  R(8.314), T(310.0), F(96485.0), m_bNonRegularGrid(false)
{
	// nothing else to do
}

template<typename TDomain>
MembraneTransportFV1<TDomain>::MembraneTransportFV1(const char* subsets, SmartPtr<IMembraneTransporter> mt)
: FV1InnerBoundaryElemDisc<TDomain>(),
  R(8.314), T(310.0), F(96485.0), m_spMembraneTransporter(mt), m_bNonRegularGrid(false)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	static_cast<IElemDisc<TDomain>*>(this)->set_subsets(subsets);
	static_cast<IElemDisc<TDomain>*>(this)->set_functions(mt->symb_fcts());
}

template<typename TDomain>
MembraneTransportFV1<TDomain>::MembraneTransportFV1(const std::vector<std::string>& subsets, SmartPtr<IMembraneTransporter> mt)
: FV1InnerBoundaryElemDisc<TDomain>(),
  R(8.314), T(310.0), F(96485.0), m_spMembraneTransporter(mt), m_bNonRegularGrid(false)
{
	// check validity of transporter setup and then lock
	mt->check_and_lock();

	static_cast<IElemDisc<TDomain>*>(this)->set_subsets(subsets);
	static_cast<IElemDisc<TDomain>*>(this)->set_functions(mt->symb_fcts());
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
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// calculate single-channel flux
	fc.flux.resize(n_flux);
	fc.from.resize(n_flux);
	fc.to.resize(n_flux);

	m_spMembraneTransporter->flux(u, e, fc.flux);

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
	size_t n_dep = m_spMembraneTransporter->n_dependencies();
	size_t n_flux = m_spMembraneTransporter->n_fluxes();

	// calculate single-channel flux
	fdc.fluxDeriv.resize(n_flux);
	fdc.from.resize(n_flux);
	fdc.to.resize(n_flux);
	for (size_t i = 0; i < n_flux; i++)
		fdc.fluxDeriv[i].resize(n_dep);

	m_spMembraneTransporter->flux_deriv(u, e, fdc.fluxDeriv);

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
	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// set assembling functions from base class first
	this->FV1InnerBoundaryElemDisc<TDomain>::prepare_setting(vLfeID, bNonRegularGrid);

	// update assemble functions
	register_all_fv1_funcs();
}


template<typename TDomain>
void MembraneTransportFV1<TDomain>::prep_timestep
(
	number time,
	VectorProxyBase* upb
)
{
	/*
	typedef CPUAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	UG_COND_THROW(!up, "Wrong algebra type!");
	const v_type& u = up->m_v;
	*/
	m_spMembraneTransporter->prep_timestep(time);
}


template<typename TDomain>
void MembraneTransportFV1<TDomain>::register_all_fv1_funcs()
{
	// register prep_timestep function for all known algebra types
	Register<bridge::CompileAlgebraList>(this);
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

