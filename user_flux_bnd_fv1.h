/*
 *  FV1UserFluxBoundary.h
 *
 *  Created on: 16.01.2014
 *      Author: mbreit
 */

#ifndef __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__FV1_USER_FLUX_BOUNDARY__
#define __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__FV1_USER_FLUX_BOUNDARY__


#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"
#include "bindings/lua/lua_user_data.h"
#include "common/util/smart_pointer.h"




namespace ug{
namespace neuro_collection{

/// \addtogroup plugin_neuro_collection
/// \{

/// Finite Volume element discretization for a Neumann Boundary
///	that only needs a boundary subset
/**
 * This class implements the inner_boundary interface to provide a normal
 * Neumann boundary the flux over which is defined by an instance of UserNumberData.
 */

template<typename TDomain>
class UserFluxBoundaryFV1
: public FV1InnerBoundaryElemDisc<TDomain>
{
	private:
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxCond FluxCond;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxDerivCond FluxDerivCond;

		///	world dimension
		static const int dim = TDomain::dim;

	public:
		/// constructor with c-strings
		UserFluxBoundaryFV1(const char* functions, const char* subsets)
			: FV1InnerBoundaryElemDisc<TDomain>(functions, subsets) {};

		/// constructor with vectors
		UserFluxBoundaryFV1(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
			: FV1InnerBoundaryElemDisc<TDomain>(functions, subsets) {};

	public:
		/// setting flux information
		void set_flux_function(SmartPtr<CplUserData<number, dim> > fluxFct)
		{
			this->m_fluxFct = fluxFct;
		}

		/// setting flux information
		void set_flux_function(number constValue)
		{
			SmartPtr<CplUserData<number, dim> > sp = make_sp(new ConstUserNumber<dim>(constValue));
			set_flux_function(sp);
		}

		/// setting flux information
		void set_flux_function(const char* name)
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


		// virtual functions inherited from FV1InnerBoundaryElemDisc
		/// calculates the flux density
		bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, GridObject* e, const MathVector<dim>& coords, int si, FluxCond& fc)
		{
			number fluxDensity;
			if (this->m_fluxFct.valid())
				(*this->m_fluxFct)(fluxDensity, coords, this->time(), si);
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

		/// calculates the flux density derivatives
		bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, GridObject* e, const MathVector<dim>& coords, int si, FluxDerivCond& fdc)
		{
			// nothing to compute here, since the flux does not depend on any variables
			return true;
		}

	protected:
		SmartPtr<UserData<number,dim> > m_fluxFct;
};

/// \}

} // namespace neuro_collection
} // end namespace ug


#endif //__H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__FV1_USER_FLUX_BOUNDARY__
