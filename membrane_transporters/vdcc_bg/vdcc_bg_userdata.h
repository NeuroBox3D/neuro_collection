/*
 * vdcc_bg.h
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_USERDATA_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_USERDATA_H__

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
		typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

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
		virtual void update_potential(side_t* elem);

	private:
		SmartPtr<CplUserData<number,dim> > m_spPotential;		//!< the UserData for potential
		bool m_bIsConstData;
};

///@}


} // namespace neuro_collection
} // namespace ug


#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__VDCC_BG_USERDATA_H__
