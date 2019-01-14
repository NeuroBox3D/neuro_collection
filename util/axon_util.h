/*
 * axon_util.h
 *
 *  Created on: 2018-10-29
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__AXON_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__AXON_UTIL_H

#include <string>

#include "common/util/smart_pointer.h"  // for SmartPtr, ConstSmartPtr
#include "lib_grid/refinement/refiner_interface.h"  // for IRefiner (ptr only), Refinement...


namespace ug {

// forward declarations
class IRefiner;
template <typename TDomain> class ApproximationSpace;

namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{


template <typename TDomain>
void unmark_ranvier_areas
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const std::string& ranvierSubsets,
	bool doUnmark = true
);



///@}

} // namespace neuro_collection
} // namespace ug


#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__AXON_UTIL_H
