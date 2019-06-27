/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Martin Stepniewski
 * Creation date: 2013-09-03
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


#ifndef UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__BOUTON_GENERATOR_H
#define UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__BOUTON_GENERATOR_H

/* system includes */
#include <stddef.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "lib_grid/lib_grid.h"
#include "lib_grid/algorithms/remeshing/delaunay_triangulation.h"



using namespace std;


namespace ug {
namespace neuro_collection{

///@addtogroup plugin_neuro_collection
///@{

/// Function for spherical point distribution
/** This function distributes N points on a sphere
 * 	of a given radius evenly.
**/
void GetNEvenlyDistributedSphereCoords(vector<vector3>& coords, int N, double radius);


/// Function for creating a bouton
/** This function creates a presynaptic bouton of the drosophila
 * 	neuromuscular junction
 *
 * 	@param bExtSpace			flag for creation of extracellur space
 * 	@param radius				bouton radius in um
 * 	@param numRefinements		number of refinements of the bouton icosphere (smoothness)
 * 	@param numReleaseSites		number of synapses
 * 	@param TbarHeight			height of synaptic T-bar
 * 	@param TbarLegRadius		leg radius of synaptic T-bar
 * 	@param TbarTopRadius		plate radius of synaptic T-bar
 * 	@param TbarTopHeight		plate height of synaptic T-bar
 * 	@param fileName				output filename of ugx grid
**/
void BuildBouton(	bool bExtSpace, number radius, int numRefinements, int numReleaseSites,
					number TbarHeight,
					number TbarLegRadius,
					number TbarTopRadius,
					number TbarTopHeight,
					string fileName);


/// Function for creating a synaptic T-bar
/** This function creates a synaptic T-bar inside a presynaptic bouton
 * 	of the drosophila neuromuscular junction
 *
 * 	@param grid				grid reference
 * 	@param sh_orig			subset handler reference
 * 	@param vrt				vertex reference on top of which the T-bar will be created
 * 	@param aaPos			position accessor
 * 	@param aaNorm			norm accessor
 * 	@param TbarHeight		height of synaptic T-bar
 * 	@param TbarLegRadius	leg radius of synaptic T-bar
 * 	@param TbarTopRadius	plate radius of synaptic T-bar
 * 	@param TbarTopHeight	plate height of synaptic T-bar

**/
void BuildTbar(	Grid& grid, SubsetHandler& sh_orig, Vertex* vrt,
				Grid::VertexAttachmentAccessor<APosition>& aaPos,
				Grid::VertexAttachmentAccessor<ANormal>& aaNorm,
				int si,
				number TbarHeight,
				number TbarLegRadius,
				number TbarTopRadius,
				number TbarTopHeight);


////////////////////////////////////////////////////////////////////////////////////////////
//	SaveSelectionStatesToFile
////////////////////////////////////////////////////////////////////////////////////////////
void SaveSelectionStatesToFile(Grid& mg, Selector& msel, const char* filename);

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__BOUTON_GENERATOR_H

