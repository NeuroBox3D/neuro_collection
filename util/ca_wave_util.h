/*
 * ca_wave_util.h
 *
 *  Created on: 2017-09-18
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__CA_WAVE_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__CA_WAVE_UTIL_H

#include "common/types.h"  // number
#include "common/util/smart_pointer.h"
#include "../membrane_transporters/ryr_implicit.h"


namespace ug {
namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{

/**
 * @brief Find maximal RyR flux density through a membrane subset.
 *
 * This function assumes that all of the following functions are
 * defined on the membrane subset:
 * - c_c, c_e (cytosolic and ER calcium concentrations)
 * - c1, c2 (RyR state variables from the Keizer & Levine model)
 * Their names must be provided in the fctNames argument in this order.
 *
 * The function is designed to work with the fully implicit RyR
 * implementation.
 *
 * @param u             solution grid function
 * @param fctNames      names for functions ca_cyt, ca_er, c1, c2 (in this order)
 * @param subsetNames   names of all ER membrane subsets with RyR channels
 * @param ryr           RyR channel object
 * @return              maximal flux density through RyR channel (mol/(m^2*s))
 */
template <typename TGridFunction>
number maxRyRFluxDensity
(
	ConstSmartPtr<TGridFunction> u,
	const char* fctNames,
	const char* subsetNames,
	ConstSmartPtr<RyRImplicit<typename TGridFunction::domain_type> > ryr
);


/**
 * @brief Find rightmost point of RyR activity.
 *
 * This function is useful in determining the wave front position
 * of a RyR-induced calcium wave from the left to the right.
 * The wave front is determined by a user-specified threshold on the
 * open probability of the RyR channel.
 *
 * The function is designed to work with the fully implicit RyR
 * implementation.
 *
 * @param u             solution grid function
 * @param fctNames      names for functions c1, c2 (in this order)
 * @param subsetNames   names of all ER membrane subsets with RyR channels
 * @param thresh        threshold value of (1-(c1+c2))
 * @return              rightmost vertex where threshold value is exceeded
 */
template <typename TGridFunction>
number waveFrontX
(
	ConstSmartPtr<TGridFunction> u,
	const char* fctNames,
	const char* subsetNames,
	number thresh
);


///@}

} // namespace ug
} // namespace neuro_collection

#include "ca_wave_util_impl.h"

#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__CA_WAVE_UTIL_H
