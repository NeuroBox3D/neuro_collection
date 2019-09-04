/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-12-27
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


/**
 * @brief Creates a 3D grid representation of a cable geometry given as SWC file.
 *
 * The position and radius information of the SWC file is interpolated using
 * cubic splines. These define a kind of backbone of the neurites (dendrites or
 * axons), around which the volume of the cytosol is inflated with the appropriate
 * radius.
 *
 * To make the grid as coarse as possible, the mesh consists of anisotropic
 * hexahedral elements that each connect two locations sampled from the backbone
 * splines. Rather isotropic hexahedra are used to mesh branching points.
 *
 * In a simulation, the geometry will be continually refined anisotropically
 * first, so as to reduce the anisotropy of the elements until all of them are
 * relatively isotropic.
 * The grid is created with a NeuriteProjector that can be used during refinement
 * to move vertices created in the refinement to appropriate positions along the
 * smooth model spline.
 *
 * For testing purposes, one can directly have the geometry refined (but only
 * isotropically) here.
 *
 * @param fileName    input file (SWC)
 * @param anisotropy  target anisotropy (ratio longest edge / shortest edge) of
 *                    the elements in the coarse grid
 * @param numRefs     number of isotropic refinements to perform
 */
void import_neurites_from_swc
(
	const std::string& fileName,
	number anisotropy = 2.0,
	size_t numRefs = 0
);

/**
 * @brief Creates a 3D grid representation of a cable geometry given as SWC file.
 *
 * The position and radius information of the SWC file is interpolated using
 * cubic splines. These define a kind of backbone of the neurites (dendrites or
 * axons), around which the volume of the cytosol is inflated with the appropriate
 * radius.
 * The mesh contains a centrally positioned ER, which is a scaled version of the
 * cytosol.
 *
 * To make the grid as coarse as possible, the mesh consists of anisotropic
 * hexahedral elements that each connect two locations sampled from the backbone
 * splines. Rather isotropic hexahedra are used to mesh branching points.
 *
 * In a simulation, the geometry will be continually refined anisotropically
 * first, so as to reduce the anisotropy of the elements until all of them are
 * relatively isotropic.
 * The grid is created with a NeuriteProjector that can be used during refinement
 * to move vertices created in the refinement to appropriate positions along the
 * smooth model spline.
 *
 * For testing purposes, one can directly have the geometry refined (but only
 * isotropically) here.
 *
 * @param fileNameIn     input file (SWC)
 * @param fileNameOut    output file (UGX)
 * @param erScaleFactor  ratio ER radius / cytosol radius (must be < 1)
 * @param anisotropy     target anisotropy (ratio longest edge / shortest edge) of
 *                       the elements in the coarse grid
 * @param numRefs        number of isotropic refinements to perform
 */
void import_er_neurites_from_swc
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number erScaleFactor,
	number anisotropy = 2.0,
	size_t numRefs = 0
);

/**
 * @brief Creates a 1D grid representation of a cable geometry given as SWC file.
 *
 * The position and radius information of the SWC file is interpolated using
 * cubic splines. These define a kind of backbone of the neurites (dendrites or
 * axons), around which the volume of the cytosol is inflated with the appropriate
 * radius.
 *
 * The spline backbone is sampled at regular intervals and these are connected by
 * edges. The sampling frequency is determined by the desired ratio between the
 * length of the edges and the radius of the neurites at their locations.
 *
 * The grid is created with a NeuriteProjector that can be used during refinement
 * to move vertices created in the refinement to appropriate positions along the
 * smooth model spline.
 *
 * For testing purposes, one can directly have the geometry refined (but only
 * isotropically) here.
 *
 * As the resulting geometries are usually meant for simulation with the
 * cable_neuron plugin, the geometry is scaled (default from um scale to m scale)
 * in the end.
 *
 * @param fileName    input file (SWC)
 * @param anisotropy  target anisotropy (ratio edge length / radius) of the
 *                    elements in the coarse grid
 * @param numRefs     number of isotropic refinements to perform
 * @param scale       factor for final scaling of the geometry
 */
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
