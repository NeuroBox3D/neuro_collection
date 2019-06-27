/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2014-12-01
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_NEURON_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_NEURON_H

#include "../../plugins/MembranePotentialMapping/vm2ug_mpm.h"
#include "../../plugins/MembranePotentialMapping/neuron_mpm.h"

#include "vdcc_bg.h"


namespace ug{
namespace neuro_collection{

///@addtogroup plugin_neuro_collection
///@{

/// Borg Graham type VGCCs with Vm2uG membrane potential supply by NEURON.
/** This class is a specialization of the Borg-Graham interface.
 *	It supplies the channel with the necessary membrane potential values by a Vm2uG object,
 *	that has to be fed file names for the files containing the V_m values produced by Neuron.
**/
template<typename TDomain>
class VDCC_BG_VM2UG_NEURON : public VDCC_BG<TDomain>
{
private:
	SmartPtr<membrane_potential_mapping::Transformator> m_NrnInterpreter;
	SmartPtr<membrane_potential_mapping::Mapper<TDomain::dim, number> > m_vmProvider;
	SmartPtr<membrane_potential_mapping::NeuronMPM> m_mapper;
	//SmartPtr<Vm2uG<std::string> > m_vmProvider;

	protected:
		using VDCC_BG<TDomain>::R;		//!< universal gas constant
		using VDCC_BG<TDomain>::T;		//!< temperature (310K)
		using VDCC_BG<TDomain>::F;		//!< Faraday constant

		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		//typedef Vm2uG<std::string> vmProvType;
		typedef membrane_potential_mapping::Mapper<TDomain::dim, number> vmProvType;

	public:
		/**
		 * @brief constructor with vectors
		 *
		 * @param fcts			functions as vector of string
		 * @param subsets		subsets as vector of string
		 * @param approx		approximation space
		 * @param transformator	transformator object (default: Transformator)
		 * @param baseName		base file name for V_m file input (default: "timesteps/timestep_")
		 * @param timeFmt		format of time in file names (default: "%.4f")
		 * @param ext			file extension of input files (default: ".dat")
		 * @param posCanChange	whether coordinates in input files can change (default: false)
		 */
		VDCC_BG_VM2UG_NEURON
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
			SmartPtr<membrane_potential_mapping::Transformator> transformator = make_sp(new membrane_potential_mapping::Transformator),
			const std::string baseName = "timesteps/timestep_",
			const char* timeFmt = "%.4f",
			const std::string ext = ".dat",
			const bool posCanChange = false
		);

		/**
		 * @brief constructor with c-strings
		 *
		 * @param fcts			functions as comma-separated c-string
		 * @param subsets		subsets as comma-separated c-string
		 * @param approx		approximation space
		 * @param transformator	transformator object (default: Transformator)
		 * @param baseName		base file name for V_m file input (default: "timesteps/timestep_")
		 * @param timeFmt		format of time in file names (default: "%.4f")
		 * @param ext			file extension of input files (default: ".dat")
		 * @param posCanChange	whether coordinates in input files can change (default: false)
		 */
		VDCC_BG_VM2UG_NEURON
		(
			const char* fcts,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
			SmartPtr<membrane_potential_mapping::Transformator> transformator = make_sp(new membrane_potential_mapping::Transformator),
			const std::string baseName = "timesteps/timestep_",
			const char* timeFmt = "%.4f",
			const std::string ext = ".dat",
			const bool posCanChange = false
		);

		/// destructor
		virtual ~VDCC_BG_VM2UG_NEURON();

		// inherited from BorgGraham
		virtual void init(number time);

		// update internal time if necessary
		virtual void update_time(number newTime);

		// update membrane potential
		virtual void update_potential(side_t* elem);

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const;

		/// transformator and mapper soon obsolete!

		// set the transformator
		inline void set_transformator(SmartPtr<membrane_potential_mapping::Transformator> transformator) {
			this->m_NrnInterpreter = transformator;
		}

		// set the provider
		inline void set_provider(SmartPtr<membrane_potential_mapping::Mapper<TDomain::dim, number> > provider) {
			this->m_vmProvider = provider;
		}

		// set the mapper
		inline void set_mapper(SmartPtr<membrane_potential_mapping::NeuronMPM> mapper) {
			this->m_mapper = mapper;
		}

	private:
		/// whether this channel has an inactivating gate
		bool has_hGate() {return this->m_channelType == VDCC_BG<TDomain>::BG_Ntype
								|| this->m_channelType == VDCC_BG<TDomain>::BG_Ttype;}

	private:
		std::string m_tFmt;				//!< time format for the membrane potential files
		number m_vmTime;
};


///@}


} // namespace neuro_collection
} // namespace ug


#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_NEURON_H
