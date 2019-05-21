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
#include "lib_disc/domain.h"            // for Domain3d, ...
#include "lib_grid/refinement/refiner_interface.h"  // for IRefiner (ptr only), Refinement...


namespace ug {

// forward declarations
class IRefiner;
template <typename TDomain> class ApproximationSpace;

namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{

/**
 * \brief Scales all functions contained in a grid function.
 *
 *	Each function has a separate scaling factor.
 *
 * \param scaledVecOut    the scaled grid function (output)
 * \param dimlessVecIn	  the original grid function (input)
 * \param scalingFactors  vector of scales for each of the composite functions
 */
template <typename TGridFunction>
void scale_dimless_vector
(
	SmartPtr<TGridFunction> scaledVecOut,
	ConstSmartPtr<TGridFunction> dimlessVecIn,
	const std::vector<number>& scalingFactors
);



template <typename TDomain>
void mark_global(SmartPtr<IRefiner> refiner, SmartPtr<TDomain> domain);

template <typename TDomain>
void mark_anisotropic
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	number thresholdRatio
);

template <typename TDomain>
void mark_anisotropic_onlyX
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	number thresholdRatio
);


void MarkNeuriteForAxialRefinement(SmartPtr<IRefiner> refiner, SmartPtr<Domain3d> domain);


template <typename TDomain>
void RemoveAllNonDefaultRefinementProjectors(SmartPtr<TDomain> dom);



///@}

} // namespace ug
} // namespace neuro_collection

#include "misc_util_impl.h"


#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__MISC_UTIL_H
