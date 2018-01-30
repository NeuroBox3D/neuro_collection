/*
 * misc_util.h
 *
 *  Created on: 2018-01-30
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__MISC_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__MISC_UTIL_H

#include <vector>                       // for vector

#include "common/types.h"               // for number
#include "common/util/smart_pointer.h"  // for SmartPtr, ConstSmartPtr


namespace ug {
namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{


template <typename TGridFunction>
void scale_dimless_vector
(
	SmartPtr<TGridFunction> scaledVecOut,
	ConstSmartPtr<TGridFunction> dimlessVecIn,
	const std::vector<number>& scalingFactors
);


///@}

} // namespace ug
} // namespace neuro_collection

#include "misc_util_impl.h"


#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__MISC_UTIL_H
