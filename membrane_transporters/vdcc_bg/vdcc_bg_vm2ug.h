/*
 * vdcc_bg.h
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_VM2UG_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_VM2UG_H__

#include "../../plugins/MembranePotentialMapping/vm2ug_mpm.h"
#include "../../plugins/MembranePotentialMapping/neuron_mpm.h"

#include "vdcc_bg.h"


namespace ug{
namespace neuro_collection{

///@addtogroup plugin_neuro_collection
///@{

/// Borg Graham type VGCCs with Vm2uG membrane potential supply.
/** This class is a specialization of the Borg-Graham interface.
 *	It supplies the channel with the necessary membrane potential values by a Vm2uG object,
 *	that has to be fed file names for the files containing the V_m values produced by Neuron.
**/

template<typename TDomain>
class VDCC_BG_VM2UG : public VDCC_BG<TDomain>
{

	protected:
		using VDCC_BG<TDomain>::R;			//!< universal gas constant
		using VDCC_BG<TDomain>::T;			//!< temperature (310K)
		using VDCC_BG<TDomain>::F;			//!< Faraday constant
		using VDCC_BG<TDomain>::has_hGate;
		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		///typedef Vm2uG<std::string> vmProvType;

	public:
		/**
		 * @brief constructor with vectors
		 *
		 * @param fcts			functions as vector of string
		 * @param subsets		subsets as vector of string
		 * @param approx		approximation space
		 * @param baseName		base file name for V_m file input (default: "timesteps/timestep_")
		 * @param timeFmt		format of time in file names (default: "%.4f")
		 * @param ext			file extension of input files (default: ".dat")
		 * @param posCanChange	whether coordinates in input files can change (default: false)
		 */
		VDCC_BG_VM2UG
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
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
		 * @param baseName		base file name for V_m file input (default: "timesteps/timestep_")
		 * @param timeFmt		format of time in file names (default: "%.4f")
		 * @param ext			file extension of input files (default: ".dat")
		 * @param posCanChange	whether coordinates in input files can change (default: false)
		 */
		VDCC_BG_VM2UG
		(
			const char* fcts,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx,
			const std::string baseName = "timesteps/timestep_",
			const char* timeFmt = "%.4f",
			const std::string ext = ".dat",
			const bool posCanChange = false
		);

		/// destructor
		virtual ~VDCC_BG_VM2UG();

		/// @copydoc VDCC_BG<TDomain>::init()
		virtual void init(number time);

		/// @copydoc VDCC_BG<TDomain>::update_potential()
		virtual void update_potential(side_t* elem);

		/// @copydoc VDCC_BG<TDomain>::update_time()
		virtual void update_time(number newTime);

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const;

		/**
		 * @brief setting the times for which files are present
		 * Using this method, it is possible to tell the VDCC implementation at which
		 * points in time there exists a file with the corresponding membrane potential
		 * values. Points in time must be separated by a constant interval and may be
		 * shifted by an offset.
		 *
		 * @param fileInterval	constant interval separating two points in time
		 * @param fileOffset	offset (the first point in time, default: 0.0)
		 */
		void set_file_times(const number fileInterval, const number fileOffset = 0.0)
		{
			m_fileInterval = fileInterval;
			m_fileOffset = fileOffset;
		}

	private:
		Vm2uGMPM<TDomain::dim> m_vmProvider;		    //!< the Vm2uG object
		std::string m_tFmt;				//!< time format for the membrane potential files
		number m_fileInterval;			//!< intervals in which voltage files are available
		number m_fileOffset;			//!< offset of time intervals for which voltage files are available

		std::string m_timeAsString;
		std::string m_baseName;
		std::string m_ext;
};

///@}

} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_VM2UG_H__
