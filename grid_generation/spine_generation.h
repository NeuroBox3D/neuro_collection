/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Marcus Kessler
 * Creation date: 2014-09-16
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__SPINE_GENERATION_H
#define UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__SPINE_GENERATION_H

// system includes
#include <string>
#include <vector>

#include "common/types.h"  // for number


namespace ug {

/**
 * \brief Builds a 3D spine geometry with part of the connected dendrite and ER
 *
 * The spine is mushroom-like with defined head and neck and contains a piece
 * of ER growing from the base of the spine up to a specified position in the spine.
 * A PSD zone is marked by a separate subset at the tip of the spine head.
 *
 * The generated geometries are intended to be used for simulations of calcium
 * dynamics following synaptic release.
 *
 * This function has been used for the creation of the geomtries in:
 * M. Breit et al.: "Spine-to-dendrite calcium modeling discloses relevance for
 * precise positioning of ryanodine receptor-containing spine endoplasmic
 * reticulum", Scientific Reports (2018).
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

