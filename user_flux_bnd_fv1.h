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
