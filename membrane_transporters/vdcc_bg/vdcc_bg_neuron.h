/*
 * vdcc_bg_neuron.h
 *
 *  Created on: 05.02.2013
 *      Author: sgrein
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_NEURON_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_NEURON_H__

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


#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_NEURON_H__
