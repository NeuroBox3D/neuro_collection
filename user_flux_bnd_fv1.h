/*
 *  FV1UserFluxBoundary.h
 *
 *  Created on: 16.01.2014
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__USER_FLUX_BND_FV1_H
#define UG__PLUGINS__NEURO_COLLECTION__USER_FLUX_BND_FV1_H


#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"
#include "bindings/lua/lua_user_data.h"
#include "common/util/smart_pointer.h"


namespace ug {
namespace neuro_collection {

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
		UserFluxBoundaryFV1(const char* functions, const char* subsets);

		/// constructor with vectors
		UserFluxBoundaryFV1
		(
			const std::vector<std::string>& functions,
			const std::vector<std::string>& subsets
		);

	public:
		/// setting flux information
		void set_flux_function(SmartPtr<CplUserData<number, dim> > fluxFct);

		/// setting flux information
		void set_flux_function(number constValue);

		/// setting flux information
		void set_flux_function(const char* name);


		// virtual functions inherited from FV1InnerBoundaryElemDisc
		/// calculates the flux density
		bool fluxDensityFct
		(
			const std::vector<LocalVector::value_type>& u,
			GridObject* e,
			const MathVector<dim>& coords,
			int si,
			FluxCond& fc
		);


		/// calculates the flux density derivatives
		bool fluxDensityDerivFct
		(
			const std::vector<LocalVector::value_type>& u,
			GridObject* e,
			const MathVector<dim>& coords,
			int si,
			FluxDerivCond& fdc
		);

	protected:
		SmartPtr<UserData<number,dim> > m_fluxFct;
};

/// \}

} // namespace neuro_collection
} // end namespace ug


#endif  // UG__PLUGINS__NEURO_COLLECTION__USER_FLUX_BND_FV1_H
