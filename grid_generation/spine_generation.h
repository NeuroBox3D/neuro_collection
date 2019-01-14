/*
 * spine_generation.h
 *
 *  Created on: 16.09.2014
 *      Author: Marcus Kessler, de-ProMesh-ified by mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__SPINE_GENERATION_H
#define UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__SPINE_GENERATION_H

// system includes
#include <string>
#include <vector>

#include "common/types.h"  // for number


namespace ug {

/**
 * \brief Builds a spine with part of the connected dendrite and ER
 *
 * \param paramVector  Contains geometry parameters in the following order:
 *                       - cytosol radius (um)
 *                       - ER radius (um)
 *                       - dendrite length (um)
 *                       - spine position (um)
 *                       - spine ER neck radius (um)
 *                       - spine ER neck length (um)
 *                       - spine ER head radius (in addition to spine ER neck radius) (um)
 *                       - spine ER head length (um)
 *                       - spine neck radius (um)
 *                       - spine neck length (um)
 *                       - spine head radius (in addition to spine neck radius) (um)
 *                       - spine head length (um)
 * \param boolVector  Contains switches for optional configuration:
 *                      - build a synapse?  -- unused
 *                      - build ER?
 *                      - build spine ER?
 *                      - synapse at different location?  -- unused
 *                      - build spine ER head?
 * \param fileName  File name for the produced grid file. Needs to be *.ugx.
 */
void BuildSpine
(
	const std::vector<number>& paramVector,
	const std::vector<bool>& boolVector,
	const std::string& fileName
);


} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__SPINE_GENERATION_H

