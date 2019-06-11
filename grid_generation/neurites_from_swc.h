/*
 * neurites_from_swc.h
 *
 *  Created on: 2016-12-27
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__NEURITES_FROM_SWC_H
#define UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__NEURITES_FROM_SWC_H

#include "common/types.h"  // number
#include "common/math/ugmath_types.h"  // vector3
#include "lib_grid/grid/grid.h"  // Grid
#include "lib_grid/refinement/projectors/neurite_projector.h"
#include "lib_grid/tools/subset_handler_grid.h"  // SubsetHandler

#include <cstdlib>
#include <vector>
#include <string>
#include <utility>

namespace ug {
namespace neuro_collection {
namespace neurites_from_swc {


void import_neurites_from_swc
(
	const std::string& fileName,
	number anisotropy = 2.0,
	size_t numRefs = 0
);
void import_er_neurites_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number erScaleFactor,
	number anisotropy = 2.0,
	size_t numRefs = 0
);
void import_1d_neurites_from_swc
(
	const std::string& fileName,
	number anisotropy = 2.0,
	size_t numRefs = 0,
	number scale = 1e-6
);


} // namespace neurites_from_swc
} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__NEURITES_FROM_SWC_H
