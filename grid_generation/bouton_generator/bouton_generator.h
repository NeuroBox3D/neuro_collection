/*
 * bouton_generator.h
 *
 *  Created on: 03.09.2013
 *      Author: Martin Stepniewski
 */


#ifndef __BOUTON_GENERATOR_H__
#define __BOUTON_GENERATOR_H__

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

#endif // __BOUTON_GENERATOR_H__

