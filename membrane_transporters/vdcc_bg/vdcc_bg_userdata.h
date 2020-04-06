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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_USERDATA_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_USERDATA_H

#include "vdcc_bg.h"

namespace ug{
namespace neuro_collection{

///@addtogroup plugin_neuro_collection
///@{


/// Borg Graham type VGCCs with UserData membrane potential supply.
/** This class is a specialization of the Borg-Graham interface.
 *	It supplies the channel with the necessary membrane potential values by a UserData object,
 *	i.e. constant UserData or UserData provided by a lua function.
**/
template<typename TDomain>
class VDCC_BG_UserData : public VDCC_BG<TDomain>
{
	public:
		static const int dim = TDomain::dim;	//!< world dimension

	protected:
		using typename VDCC_BG<TDomain>::vm_grid_object;
		using VDCC_BG<TDomain>::R;			//!< universal gas constant
		using VDCC_BG<TDomain>::T;			//!< temperature (310K)
		using VDCC_BG<TDomain>::F;			//!< Faraday constant
		using VDCC_BG<TDomain>::has_hGate;	//!< Faraday constant


	public:
		/**
		 * @brief constructor with vectors
		 *
		 * @param fcts		functions as vector of string
		 * @param subsets	subsets as vector of string
		 * @param approx	approximation space
		 */
		VDCC_BG_UserData
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		/**
		 * @brief constructor with c-strings
		 *
		 * @param fcts		functions as comma-separated c-string
		 * @param subsets	subsets as comma-separated c-string
		 * @param approx	approximation space
		 */
		VDCC_BG_UserData
		(
			const char* fcts,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		/// destructor
		virtual ~VDCC_BG_UserData();

		/// adding potential information for pumps/channels in membrane
		void set_potential_function(const char* name);
		void set_potential_function(const number value);
		void set_potential_function(SmartPtr<CplUserData<number, dim> > spPotFct);

		/// @copydoc VDCC_BG<TDomain>::update_potential()
		virtual void update_potential(vm_grid_object* elem);

	private:
		SmartPtr<CplUserData<number,dim> > m_spPotential;		//!< the UserData for potential
		bool m_bIsConstData;
};

///@}


} // namespace neuro_collection
} // namespace ug


#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_USERDATA_H
