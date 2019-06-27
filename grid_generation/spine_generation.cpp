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

#include "spine_generation.h"

#include <cmath>

#include "common/math/math_vector_matrix/math_matrix_vector_functions.h"  // for VecAdd
#include "common/math/ugmath_types.h"  // for vector3
#include "lib_grid/algorithms/extrusion/extrude.h"  // for Extrude
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"  // for CalculateCenter, SplitEdge
#include "lib_grid/algorithms/geom_obj_util/misc_util.h"  // for FindClosestByCoordinate
#include "lib_grid/algorithms/geom_obj_util/volume_util.h"  // for CalculateCenter
#include "lib_grid/algorithms/grid_generation/tetrahedralization.h"  // for Tetrahedralize
#include "lib_grid/algorithms/selection_util.h"  // for SelectAssociatedVertices, ...
#include "lib_grid/algorithms/subset_util.h"  // for SeparateSubsetsByLowerDimSubsets ...
#include "lib_grid/common_attachments.h"  // for aPosition, aNormal
#include "lib_grid/file_io/file_io.h"  // for SaveGridToFile
#include "lib_grid/grid/grid.h"  // for Grid
#include "lib_grid/tools/selector_grid.h"  // for Selector
#include "lib_grid/tools/subset_handler_grid.h"  // for SubsetHandler


namespace ug {

//#define DG_DEBUG 1


static void CreateCircle
(
	Grid& grid,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	const vector3& center,
	number radius,
  	int numRimVertices
)
{
	// create the vertices
	// create one up front, the others in a loop
	Vertex* firstVrt = *grid.create<RegularVertex>();
	aaPos[firstVrt] = vector3(0, radius, 0);
	VecAdd(aaPos[firstVrt], aaPos[firstVrt], center);
	Vertex* lastVrt = firstVrt;
	for (int i = 1; i < numRimVertices; ++i)
	{
		//	create a new vertex
		number ia = (float)i / (float)numRimVertices;
		Vertex* vNew = *grid.create<RegularVertex>();
		aaPos[vNew] = vector3(sin(2. * PI * ia), cos(2. * PI * ia), 0);
		VecScale(aaPos[vNew], aaPos[vNew], radius);
		VecAdd(aaPos[vNew], aaPos[vNew], center);

		//	create a new edge
		grid.create<RegularEdge>(EdgeDescriptor(lastVrt, vNew));

		//	prepare the next iteration
		lastVrt = vNew;
	}

	// one triangle / edge is still missing
	grid.create<RegularEdge>(EdgeDescriptor(lastVrt, firstVrt));
}


static void ScaleAroundCenter
(
	const std::vector<Vertex*>& vrts,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	const vector3& scale
)
{
	vector3 center = CalculateCenter(vrts.begin(), vrts.end(), aaPos);

	std::vector<Vertex*>::const_iterator iter = vrts.begin();
	for (; iter != vrts.end(); ++iter)
	{
		vector3& v = aaPos[*iter];
		VecSubtract(v, v, center);
		v[0] *= scale[0];
		v[1] *= scale[1];
		v[2] *= scale[2];
		VecAdd(v, v, center);
	}
}





void BuildSpine
(
	const std::vector<number>& paramVector,
	const std::vector<bool>& boolVector,
	const std::string& fileName
)
{
	// geometric parameters
	number cyt_radius = paramVector[0];
	number er_radius = paramVector[1];
	number dend_length = paramVector[2];
	number spine_pos = paramVector[3];
	number app_neck_radius = paramVector[4];
	number app_neck_length = paramVector[5];
	number app_head_radius = paramVector[6];
	number app_head_length = paramVector[7];
	number spine_neck_radius = paramVector[8];
	number spine_neck_length = paramVector[9];
	number spine_head_radius = paramVector[10];
	number spine_head_length = paramVector[11];

	// functional switches
	bool opt_ER = boolVector[1];
	bool opt_app = boolVector[2];
	bool bMakeAppHead = true;
	if (boolVector.size() >= 5)
		bMakeAppHead = boolVector[4];



	// create grid and subset handler
	Grid grid;
	SubsetHandler sh(grid);
	sh.set_default_subset_index(0);

	// setup coordinate attachments and accessors
	grid.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	// prepare a selector
	Selector sel(grid);

// //////////////////////
// create cytosol / ER //
// //////////////////////
	vector3 cyt_start(0.0, 0.0, 0.0);
	vector3 cyt_end(0.0, 0.0, dend_length);
	vector3 er_start(0.0, 0.0, 1.0);
	vector3 er_end(0.0, 0.0, dend_length - 1.0);

	const int numRimVertices = 8;

	// create left-hand circle for cytosol
	sel.enable_autoselection(true);
	CreateCircle(grid, aaPos, cyt_start, cyt_radius, numRimVertices);
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), 0);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), 0);

	// extrude to build dendrite
	std::vector<Vertex*> vrts;
	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	std::vector<Edge*> edges;
	edges.assign(sel.edges_begin(), sel.edges_end());

	vector3 cyt_ExtrDir(0.0, 0.0, (spine_pos - spine_neck_radius) / 10.0);
	for (size_t i = 0; i < 10; ++i)
		Extrude(grid, &vrts, &edges, NULL, cyt_ExtrDir, aaPos, EO_CREATE_FACES, NULL);

	cyt_ExtrDir[2] = spine_neck_radius;
	Extrude(grid, &vrts, &edges, NULL, cyt_ExtrDir, aaPos, EO_CREATE_FACES, NULL);
	Extrude(grid, &vrts, &edges, NULL, cyt_ExtrDir, aaPos, EO_CREATE_FACES, NULL);

	cyt_ExtrDir[2] = (dend_length - spine_neck_radius - spine_pos) / 10.0;
	for (size_t i = 0; i < 10; ++i)
		Extrude(grid, &vrts, &edges, NULL, cyt_ExtrDir, aaPos, EO_CREATE_FACES, NULL);

	// create left-hand circle for ER
	sel.clear();
	CreateCircle(grid, aaPos, er_start, er_radius, numRimVertices);
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), 1);
	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), 1);

	// extrude to build ER
	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	edges.assign(sel.edges_begin(), sel.edges_end());

	vector3 er_ExtrDir(0.0, 0.0, (spine_pos - spine_neck_radius - 1.0) / 9.0);
	for (size_t i = 0; i < 9; ++i)
		Extrude(grid, &vrts, &edges, NULL, er_ExtrDir, aaPos, EO_CREATE_FACES, NULL);

	er_ExtrDir[2] = spine_neck_radius;
	Extrude(grid, &vrts, &edges, NULL, er_ExtrDir, aaPos, EO_CREATE_FACES, NULL);
	Extrude(grid, &vrts, &edges, NULL, er_ExtrDir, aaPos, EO_CREATE_FACES, NULL);

	er_ExtrDir[2] = (dend_length - spine_neck_radius - spine_pos - 1.0) / 9.0;
	for (size_t i = 0; i < 9; ++i)
		Extrude(grid, &vrts, &edges, NULL, er_ExtrDir, aaPos, EO_CREATE_FACES, NULL);


	// generate faces for start of dendrite and for measurement zone border
	RegularVertex* vrt_cytStart = *grid.create<RegularVertex>();
	aaPos[vrt_cytStart] = cyt_start;
	sh.assign_subset(vrt_cytStart, 0);
	RegularVertex* vrt_cytEnd = *grid.create<RegularVertex>();
	aaPos[vrt_cytEnd] = cyt_end;
	sh.assign_subset(vrt_cytEnd, 0);

	RegularVertex* vrt_erStart = *grid.create<RegularVertex>();
	aaPos[vrt_erStart] = er_start;
	sh.assign_subset(vrt_erStart, 1);
	RegularVertex* vrt_erEnd = *grid.create<RegularVertex>();
	aaPos[vrt_erEnd] = er_end;
	sh.assign_subset(vrt_erEnd, 1);

	for (int i = 0; i < numRimVertices; ++i)
	{
		// left-hand side cytosolic boundary
		sh.set_default_subset_index(0);
		vector3 searchPos(sin(2.0*PI*i/numRimVertices)*cyt_radius, cos(2.0*PI*i/numRimVertices)*cyt_radius, 0.0);
		Vertex* vrt1 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
		searchPos = vector3(sin(2.0*PI*(i+1)/numRimVertices)*cyt_radius, cos(2.0*PI*(i+1)/numRimVertices)*cyt_radius, 0.0);
		Vertex* vrt2 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
		Triangle* tri = *grid.create<Triangle>(TriangleDescriptor(vrt_cytStart, vrt1, vrt2));
		sh.assign_subset(tri, 0);

		// left-hand side ER boundary
		sh.set_default_subset_index(1);
		searchPos = vector3(sin(2.0*PI*i/numRimVertices)*er_radius, cos(2.0*PI*i/numRimVertices)*er_radius, 1.0);
		vrt1 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		searchPos = vector3(sin(2.0*PI*(i+1)/numRimVertices)*er_radius, cos(2.0*PI*(i+1)/numRimVertices)*er_radius, 1.0);
		vrt2 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		tri = *grid.create<Triangle>(TriangleDescriptor(vrt_erStart, vrt1, vrt2));
		sh.assign_subset(tri, 1);

		// left measurement zone boundary
		sh.set_default_subset_index(4);
		searchPos = vector3(sin(2.0*PI*i/numRimVertices)*cyt_radius, cos(2.0*PI*i/numRimVertices)*cyt_radius, spine_pos - spine_neck_radius);
		vrt1 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
		searchPos = vector3(sin(2.0*PI*(i+1)/numRimVertices)*cyt_radius, cos(2.0*PI*(i+1)/numRimVertices)*cyt_radius, spine_pos - spine_neck_radius);
		vrt2 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
		searchPos = vector3(sin(2.0*PI*(i+1)/numRimVertices)*er_radius, cos(2.0*PI*(i+1)/numRimVertices)*er_radius, spine_pos - spine_neck_radius);
		Vertex* vrt3 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		searchPos = vector3(sin(2.0*PI*i/numRimVertices)*er_radius, cos(2.0*PI*i/numRimVertices)*er_radius, spine_pos - spine_neck_radius);
		Vertex* vrt4 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		Quadrilateral* quad = *grid.create<Quadrilateral>(QuadrilateralDescriptor(vrt1, vrt2, vrt3, vrt4));
		sh.assign_subset(quad, 4);

		// right measurement zone boundary
		sh.set_default_subset_index(4);
		searchPos = vector3(sin(2.0*PI*i/numRimVertices)*cyt_radius, cos(2.0*PI*i/numRimVertices)*cyt_radius, spine_pos + spine_neck_radius);
		vrt1 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
		searchPos = vector3(sin(2.0*PI*(i-1)/numRimVertices)*cyt_radius, cos(2.0*PI*(i-1)/numRimVertices)*cyt_radius, spine_pos + spine_neck_radius);
		vrt2 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
		searchPos = vector3(sin(2.0*PI*(i-1)/numRimVertices)*er_radius, cos(2.0*PI*(i-1)/numRimVertices)*er_radius, spine_pos + spine_neck_radius);
		vrt3 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		searchPos = vector3(sin(2.0*PI*i/numRimVertices)*er_radius, cos(2.0*PI*i/numRimVertices)*er_radius, spine_pos + spine_neck_radius);
		vrt4 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		quad = *grid.create<Quadrilateral>(QuadrilateralDescriptor(vrt1, vrt2, vrt3, vrt4));
		sh.assign_subset(quad, 4);

		// right-hand side cytosolic boundary
		sh.set_default_subset_index(0);
		searchPos = vector3(sin(2.0*PI*i/numRimVertices)*cyt_radius, cos(2.0*PI*i/numRimVertices)*cyt_radius, dend_length);
		vrt1 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
		searchPos = vector3(sin(2.0*PI*(i-1)/numRimVertices)*cyt_radius, cos(2.0*PI*(i-1)/numRimVertices)*cyt_radius, dend_length);
		vrt2 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
		tri = *grid.create<Triangle>(TriangleDescriptor(vrt_cytEnd, vrt1, vrt2));
		sh.assign_subset(tri, 0);

		// right-hand side ER boundary
		sh.set_default_subset_index(1);
		searchPos = vector3(sin(2.0*PI*i/numRimVertices)*er_radius, cos(2.0*PI*i/numRimVertices)*er_radius, dend_length - 1.0);
		vrt1 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		searchPos = vector3(sin(2.0*PI*(i-1)/numRimVertices)*er_radius, cos(2.0*PI*(i-1)/numRimVertices)*er_radius, dend_length - 1.0);
		vrt2 = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		tri = *grid.create<Triangle>(TriangleDescriptor(vrt_erEnd, vrt1, vrt2));
		sh.assign_subset(tri, 1);
	}

#ifdef DG_DEBUG
	if (!SaveGridToFile(grid, sh, (fileName+"_shaft.ugx").c_str()))
		UG_THROW("Error while saving dendrite to file '" << fileName+"_shaft.ugx" << "'.");
#endif


// //////////////////
// create spine ER //
// //////////////////
	if (opt_app)
	{
		sh.set_default_subset_index(1);

		number oct_side_len = 2 * sin(PI/8) * er_radius;

		vector3 searchPos(sin(3*PI/8)*oct_side_len/2, er_radius - cos(3*PI/8)*oct_side_len/2, spine_pos);
		Edge* edge = FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos);
		Vertex* leftVrt = SplitEdge<RegularVertex>(grid, edge);
		aaPos[leftVrt] = vector3(app_neck_radius, er_radius, spine_pos);


		searchPos = vector3(-sin(3*PI/8)*oct_side_len/2, er_radius - cos(3*PI/8)*oct_side_len/2, spine_pos);
		edge = FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos);
		Vertex* rightVrt = SplitEdge<RegularVertex>(grid, edge);
		aaPos[rightVrt] = vector3(-app_neck_radius, er_radius, spine_pos);


		searchPos = vector3(0.0, er_radius, spine_pos + spine_neck_radius/2);
		edge = FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos);
		Vertex* upperVrt = SplitEdge<RegularVertex>(grid, edge);
		aaPos[upperVrt] = vector3(0.0, er_radius, spine_pos + app_neck_radius);


		searchPos = vector3(0.0, er_radius, spine_pos - spine_neck_radius/2);
		edge = FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos);
		Vertex* lowerVrt = SplitEdge<RegularVertex>(grid, edge);
		aaPos[lowerVrt] = vector3(0.0, er_radius, spine_pos - app_neck_radius);


		vector3 upperLeftVrtPos;
		VecScaleAdd(upperLeftVrtPos, 0.5, aaPos[upperVrt], 0.5, aaPos[leftVrt]);
		edge = FindClosestByCoordinate<Edge>(upperLeftVrtPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos);
		Vertex* upperLeftVrt = SplitEdge<RegularVertex>(grid, edge);
		aaPos[upperLeftVrt] = vector3(app_neck_radius/sqrt(2), er_radius, spine_pos + app_neck_radius/sqrt(2));


		vector3 upperRightVrtPos;
		VecScaleAdd(upperRightVrtPos, 0.5, aaPos[upperVrt], 0.5, aaPos[rightVrt]);
		edge = FindClosestByCoordinate<Edge>(upperRightVrtPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos);
		Vertex* upperRightVrt = SplitEdge<RegularVertex>(grid, edge);
		aaPos[upperRightVrt] = vector3(- app_neck_radius/sqrt(2), er_radius, spine_pos + app_neck_radius/sqrt(2));


		vector3 lowerLeftVrtPos;
		VecScaleAdd(lowerLeftVrtPos, 0.5, aaPos[lowerVrt], 0.5, aaPos[leftVrt]);
		edge = FindClosestByCoordinate<Edge>(lowerLeftVrtPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos);
		Vertex* lowerLeftVrt = SplitEdge<RegularVertex>(grid, edge);
		aaPos[lowerLeftVrt] = vector3(app_neck_radius/sqrt(2), er_radius, spine_pos - app_neck_radius/sqrt(2));


		vector3 lowerRightVrtPos;
		VecScaleAdd(lowerRightVrtPos, 0.5, aaPos[lowerVrt], 0.5, aaPos[rightVrt]);
		edge = FindClosestByCoordinate<Edge>(lowerRightVrtPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos);
		Vertex* lowerRightVrt = SplitEdge<RegularVertex>(grid, edge);
		aaPos[lowerRightVrt] = vector3(-app_neck_radius/sqrt(2), er_radius, spine_pos - app_neck_radius/sqrt(2));


		// pull outer vertices to the er_radius y plane
		// although it seems negligible, this is absolutely necessary for the tetrahedralization to work!
		vector3 offset = vector3(0, cos(3*PI/8)*oct_side_len, 0);
		searchPos = vector3(sin(3*PI/8)*oct_side_len, er_radius - cos(3*PI/8)*oct_side_len, spine_pos);
		Vertex* vrt = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		VecAdd(aaPos[vrt], aaPos[vrt], offset);
		searchPos = vector3(-sin(3*PI/8)*oct_side_len, er_radius - cos(3*PI/8)*oct_side_len, spine_pos);
		vrt = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		VecAdd(aaPos[vrt], aaPos[vrt], offset);
		searchPos = vector3(sin(3*PI/8)*oct_side_len, er_radius - cos(3*PI/8)*oct_side_len, spine_pos + spine_neck_radius);
		vrt = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		VecAdd(aaPos[vrt], aaPos[vrt], offset);
		searchPos = vector3(-sin(3*PI/8)*oct_side_len, er_radius - cos(3*PI/8)*oct_side_len, spine_pos + spine_neck_radius);
		vrt = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		VecAdd(aaPos[vrt], aaPos[vrt], offset);
		searchPos = vector3(sin(3*PI/8)*oct_side_len, er_radius - cos(3*PI/8)*oct_side_len, spine_pos - spine_neck_radius);
		vrt = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		VecAdd(aaPos[vrt], aaPos[vrt], offset);
		searchPos = vector3(-sin(3*PI/8)*oct_side_len, er_radius - cos(3*PI/8)*oct_side_len, spine_pos - spine_neck_radius);
		vrt = FindClosestByCoordinate<Vertex>(searchPos, sh.begin<Vertex>(1), sh.end<Vertex>(1), aaPos);
		VecAdd(aaPos[vrt], aaPos[vrt], offset);


		// extrude to make the actual neck and head
		sel.clear();
		VecScaleAdd(searchPos, 0.5, aaPos[upperVrt], 0.5, upperLeftVrtPos);
		sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos));
		VecScaleAdd(searchPos, 0.5, aaPos[upperVrt], 0.5, upperRightVrtPos);
		sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos));
		VecScaleAdd(searchPos, 0.5, aaPos[lowerVrt], 0.5, lowerLeftVrtPos);
		sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos));
		VecScaleAdd(searchPos, 0.5, aaPos[lowerVrt], 0.5, lowerRightVrtPos);
		sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos));
		VecScaleAdd(searchPos, 0.5, aaPos[leftVrt], 0.5, upperLeftVrtPos);
		sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos));
		VecScaleAdd(searchPos, 0.5, aaPos[rightVrt], 0.5, upperRightVrtPos);
		sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos));
		VecScaleAdd(searchPos, 0.5, aaPos[leftVrt], 0.5, lowerLeftVrtPos);
		sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos));
		VecScaleAdd(searchPos, 0.5, aaPos[rightVrt], 0.5, lowerRightVrtPos);
		sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(1), sh.end<Edge>(1), aaPos));

		SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());

		sh.assign_subset(sel.edges_begin(), sel.edges_end(), 3);
		sh.assign_subset(sel.vertices_begin(), sel.vertices_end(), 3);

		vector3 app_dir(0.0, app_neck_length / 10.0, 0.0);
		if (app_head_length + app_neck_length + er_radius == cyt_radius)
			app_dir = vector3(0.0, (app_neck_length - 0.001) / 10.0, 0.0);

		vrts.assign(sel.vertices_begin(), sel.vertices_end());
		edges.assign(sel.edges_begin(), sel.edges_end());

		for (size_t i = 0; i < 10; ++i)
			Extrude(grid, &vrts, &edges, NULL, app_dir, aaPos, EO_CREATE_FACES, NULL);

		if (bMakeAppHead)
		{
			vector3 app_head_dir(0.0, app_head_length/3.0, 0.0);
			vector3 growth_vector((app_neck_radius+app_head_radius)/app_neck_radius, 1, (app_neck_radius+app_head_radius)/app_neck_radius);
			vector3 shrink_vector(app_neck_radius/(app_neck_radius+app_head_radius), 1, app_neck_radius/(app_neck_radius+app_head_radius));

			Extrude(grid, &vrts, &edges, NULL, app_head_dir, aaPos, EO_CREATE_FACES, NULL);
			ScaleAroundCenter(vrts, aaPos, growth_vector);
			Extrude(grid, &vrts, &edges, NULL, app_head_dir, aaPos, EO_CREATE_FACES, NULL);
			Extrude(grid, &vrts, &edges, NULL, app_head_dir, aaPos, EO_CREATE_FACES, NULL);
			ScaleAroundCenter(vrts, aaPos, shrink_vector);
		}

		Extrude(grid, &vrts, &edges, NULL, vector3(0.0, 0.0, 0.0), aaPos, EO_CREATE_FACES, NULL);
		ScaleAroundCenter(vrts, aaPos, vector3(0.0, 0.0, 0.0));

		sel.clear();
		for (size_t i = 0; i < 8; ++i)
		{
			searchPos = vector3(sin(PI*(i+0.5)/4.0) * app_neck_radius/2.0, er_radius, spine_pos + cos(PI*(i+0.5)/4.0) * app_neck_radius/2.0);
			Face* f = FindClosestByCoordinate<Face>(searchPos, sh.begin<Face>(1), sh.end<Face>(1), aaPos);
			sel.select(f);
		}
		SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
		SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());

		sh.assign_subset(sel.faces_begin(), sel.faces_end(), 2);
		sh.assign_subset(sel.edges_begin(), sel.edges_end(), 2);
		sh.assign_subset(sel.vertices_begin(), sel.vertices_end(), 2);
	}

#ifdef DG_DEBUG
	if (!SaveGridToFile(grid, sh, (fileName+"_spineER.ugx").c_str()))
		UG_THROW("Error while saving dendrite to file '" << fileName+"_spineER.ugx" << "'.");
#endif



// /////////////////
// generate spine //
// /////////////////

	sh.set_default_subset_index(0);
	number oct_side_len = 2 * sin(PI/8) * cyt_radius;

	vector3 searchPos(sin(3*PI/8)*oct_side_len/2, cyt_radius - cos(3*PI/8)*oct_side_len/2, spine_pos);
	Edge* edge = FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos);
	Vertex* leftVrt = SplitEdge<RegularVertex>(grid, edge);
	aaPos[leftVrt] = vector3(spine_neck_radius, cyt_radius, spine_pos);


	searchPos = vector3(-sin(3*PI/8)*oct_side_len/2, cyt_radius - cos(3*PI/8)*oct_side_len/2, spine_pos);
	edge = FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos);
	Vertex* rightVrt = SplitEdge<RegularVertex>(grid, edge);
	aaPos[rightVrt] = vector3(-spine_neck_radius, cyt_radius, spine_pos);


	vector3 upperVrtPos = vector3(0, cyt_radius, spine_pos + spine_neck_radius);
	vector3 lowerVrtPos = vector3(0, cyt_radius, spine_pos - spine_neck_radius);


	vector3 upperLeftVrtPos;
	VecScaleAdd(upperLeftVrtPos, 0.5, upperVrtPos, 0.5, aaPos[leftVrt]);
	edge = FindClosestByCoordinate<Edge>(upperLeftVrtPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos);
	Vertex* upperLeftVrt = SplitEdge<RegularVertex>(grid, edge);
	aaPos[upperLeftVrt] = vector3(spine_neck_radius/sqrt(2), cyt_radius, spine_pos + spine_neck_radius/sqrt(2));


	vector3 upperRightVrtPos;
	VecScaleAdd(upperRightVrtPos, 0.5, upperVrtPos, 0.5, aaPos[rightVrt]);
	edge = FindClosestByCoordinate<Edge>(upperRightVrtPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos);
	Vertex* upperRightVrt = SplitEdge<RegularVertex>(grid, edge);
	aaPos[upperRightVrt] = vector3(-spine_neck_radius/sqrt(2), cyt_radius, spine_pos + spine_neck_radius/sqrt(2));


	vector3 lowerLeftVrtPos;
	VecScaleAdd(lowerLeftVrtPos, 0.5, lowerVrtPos, 0.5, aaPos[leftVrt]);
	edge = FindClosestByCoordinate<Edge>(lowerLeftVrtPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos);
	Vertex* lowerLeftVrt = SplitEdge<RegularVertex>(grid, edge);
	aaPos[lowerLeftVrt] = vector3(spine_neck_radius/sqrt(2), cyt_radius, spine_pos - spine_neck_radius/sqrt(2));


	vector3 lowerRightVrtPos;
	VecScaleAdd(lowerRightVrtPos, 0.5, lowerVrtPos, 0.5, aaPos[rightVrt]);
	edge = FindClosestByCoordinate<Edge>(lowerRightVrtPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos);
	Vertex* lowerRightVrt = SplitEdge<RegularVertex>(grid, edge);
	aaPos[lowerRightVrt] = vector3(-spine_neck_radius/sqrt(2), cyt_radius, spine_pos - spine_neck_radius/sqrt(2));


	sel.clear();
	vector3 centerVertexPos(0.0, cyt_radius, spine_pos);
	Vertex* centerVrt = FindClosestByCoordinate<Vertex>(centerVertexPos, sh.begin<Vertex>(0), sh.end<Vertex>(0), aaPos);
	sel.select(centerVrt);
	ExtendSelection(sel, 1);
	SelectAssociatedEdges(sel, sel.begin<Vertex>(), sel.end<Vertex>());
	SelectAssociatedFaces(sel, sel.begin<Edge>(), sel.end<Edge>());
	SelectInnerSelectionFaces(sel);
	SelectInnerSelectionEdges(sel);
	SelectInnerSelectionVertices(sel);
	grid.erase(sel.begin<Face>(), sel.end<Face>());
	grid.erase(sel.begin<Edge>(), sel.end<Edge>());
	grid.erase(sel.begin<Vertex>(), sel.end<Vertex>());


	// extrude to make the actual neck and head
	sel.clear();
	VecScaleAdd(searchPos, 0.5, upperVrtPos, 0.5, upperLeftVrtPos);
	sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos));
	VecScaleAdd(searchPos, 0.5, upperVrtPos, 0.5, upperRightVrtPos);
	sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos));
	VecScaleAdd(searchPos, 0.5, lowerVrtPos, 0.5, lowerLeftVrtPos);
	sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos));
	VecScaleAdd(searchPos, 0.5, lowerVrtPos, 0.5, lowerRightVrtPos);
	sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos));
	VecScaleAdd(searchPos, 0.5, aaPos[leftVrt], 0.5, upperLeftVrtPos);
	sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos));
	VecScaleAdd(searchPos, 0.5, aaPos[rightVrt], 0.5, upperRightVrtPos);
	sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos));
	VecScaleAdd(searchPos, 0.5, aaPos[leftVrt], 0.5, lowerLeftVrtPos);
	sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos));
	VecScaleAdd(searchPos, 0.5, aaPos[rightVrt], 0.5, lowerRightVrtPos);
	sel.select(FindClosestByCoordinate<Edge>(searchPos, sh.begin<Edge>(0), sh.end<Edge>(0), aaPos));

	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>());

	sh.assign_subset(sel.edges_begin(), sel.edges_end(), 5);
	sh.assign_subset(sel.vertices_begin(), sel.vertices_end(), 5);

	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	edges.assign(sel.edges_begin(), sel.edges_end());

	vector3 spine_dir(0.0, spine_neck_length/2.0, 0.0);
	Extrude(grid, &vrts, &edges, NULL, spine_dir, aaPos, EO_CREATE_FACES, NULL);
	Extrude(grid, &vrts, &edges, NULL, spine_dir, aaPos, EO_CREATE_FACES, NULL);
	vector3 spine_head_dir(0, spine_head_length/3.0, 0);
	Extrude(grid, &vrts, &edges, NULL, spine_head_dir, aaPos, EO_CREATE_FACES, NULL);
	vector3 spine_growth((spine_neck_radius+spine_head_radius)/spine_neck_radius, 1, (spine_neck_radius+spine_head_radius)/spine_neck_radius);
	ScaleAroundCenter(vrts, aaPos, spine_growth);
	VecScale(spine_head_dir, spine_head_dir, 1.0/3.0);
	for (size_t i = 0; i < 3; ++i)
		Extrude(grid, &vrts, &edges, NULL, spine_head_dir, aaPos, EO_CREATE_FACES, NULL);
	VecScale(spine_head_dir, spine_head_dir, 3.0);
	Extrude(grid, &vrts, &edges, NULL, spine_head_dir, aaPos, EO_CREATE_FACES, NULL);
	vector3 spine_shrink(spine_neck_radius/(spine_neck_radius+spine_head_radius), 1, spine_neck_radius/(spine_neck_radius+spine_head_radius));
	ScaleAroundCenter(vrts, aaPos, spine_shrink);

	const size_t num_faces = grid.num_faces();

	Extrude(grid, &vrts, &edges, NULL, vector3(0.0, 0.0, 0.0), aaPos, EO_CREATE_FACES, NULL);
	ScaleAroundCenter(vrts, aaPos, vector3(0.0, 0.0, 0.0));

	sel.clear();
	size_t counter = 0;
	FaceIterator iter = grid.begin<Face>();
	for (; iter != grid.end<Face>(); ++iter, ++counter)
		if (counter == num_faces)
			break;
	UG_COND_THROW(iter == grid.end<Face>(), "No new faces present.");
	for (; iter != grid.end<Face>(); ++iter)
		sel.select(*iter);
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>());
	SelectAssociatedEdges(sel, sel.begin<Face>(), sel.end<Face>());

	sh.assign_subset(sel.faces_begin(), sel.faces_end(), 6);
	sh.assign_subset(sel.edges_begin(), sel.edges_end(), 6);
	sh.assign_subset(sel.vertices_begin(), sel.vertices_end(), 6);

#ifdef DG_DEBUG
	if (!SaveGridToFile(grid, sh, (fileName+"_surfaceComplete.ugx").c_str()))
		UG_THROW("Error while saving dendrite to file '" << fileName+"_surfaceComplete.ugx" << "'.");
#endif


// ///////////////////
// generate volumes //
// ///////////////////

	sh.set_default_subset_index(-1);
	UG_COND_THROW(!Tetrahedralize(grid, sh, 18.0, false, false, aPosition, 0),
		"Tetrahedralization failed.");

	SeparateSubsetsByLowerDimSubsets<Volume>(grid, sh, true);
	sh.join_subsets(0, 0, 5, true);

#ifdef DG_DEBUG
	if (!SaveGridToFile(grid, sh, (fileName+"_tets.ugx").c_str()))
		UG_THROW("Error while saving dendrite to file '" << fileName+"_tets.ugx" << "'.");
#endif


// //////////////
// fix subsets //
// //////////////

	// get the subset indices for cyt, er and app
	int cyt_index = 0;
	int er_index = 0;
	int app_index = 0;

	Volume* vol = FindClosestByCoordinate<Volume>(vector3(0.0, 0.0, 0.0), grid.begin<Volume>(), grid.end<Volume>(), aaPos);
	cyt_index = sh.get_subset_index(vol);
	sh.move_subset(cyt_index, 0);

	vol = FindClosestByCoordinate<Volume>(vector3(0.0, 0.0, dend_length), grid.begin<Volume>(), grid.end<Volume>(), aaPos);
	cyt_index = sh.get_subset_index(vol);
	sh.join_subsets(0, 0, cyt_index, true);
	sh.join_subsets(0, 0, 5, true);

	vol = FindClosestByCoordinate<Volume>(vector3(0.0, 0.0, 2.0), grid.begin<Volume>(), grid.end<Volume>(), aaPos);
	er_index = sh.get_subset_index(vol);
	sh.move_subset(er_index, 1);
	sh.join_subsets(1, 1, 4, true);

	if (opt_app)
	{
		vol = FindClosestByCoordinate<Volume>(vector3(0.0, er_radius + app_neck_length / 2.0, spine_pos), grid.begin<Volume>(), grid.end<Volume>(), aaPos);
		UG_COND_THROW(sh.get_subset_index(vol) == 0 || sh.get_subset_index(vol) == 1,
			"Spine ER is not separated from rest volume.");
		app_index = sh.get_subset_index(vol);
		sh.move_subset(app_index, 2);
	}

	// split measure zones
	int numSs = sh.num_subsets();
	for (VolumeIterator vIter = grid.volumes_begin(); vIter != grid.volumes_end(); ++vIter)
	{
		Volume* vol = *vIter;
		vector3 center = CalculateCenter(vol, aaPos);
		if (center[2] > spine_pos - spine_neck_radius && center[2] < spine_pos + spine_neck_radius && center[1] < cyt_radius && sh.get_subset_index(vol) == numSs - 1)
			sh.assign_subset(vol, 11);
		if (center[2] > spine_pos - spine_neck_radius && center[2] < spine_pos + spine_neck_radius && center[1] <= cyt_radius + spine_neck_length && sh.get_subset_index(vol) == numSs - 1)
			sh.assign_subset(vol, 12);
		else if (center[1] > cyt_radius + spine_neck_length && sh.get_subset_index(vol) == numSs - 1)
			sh.assign_subset(vol, 13);
	}


	// erase empty subsets
	int i = 0;
	while (i < sh.num_subsets())
	{
		if (sh.empty(i))
			sh.erase_subset(i);
		else
			++i;
	}

	// except those for spine ER
	if (!opt_app)
	{
		sh.insert_subset(2);
		sh.insert_subset(5);
	}


	// make sure all faces, edges and vertices belong to a subset
	sh.swap_subsets(0, 4);  // somehow, these are necessary
	sh.assign_subset(grid.vertices_begin(), grid.vertices_end(), -1);
	sh.assign_subset(grid.edges_begin(), grid.edges_end(), -1);

	for (int i = 0; i < sh.num_subsets(); ++i)
	{
		CopySubsetIndicesToSides(sh, sh.begin<Face>(i), sh.end<Face>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Edge>(i), sh.end<Edge>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Vertex>(i), sh.end<Vertex>(i), true);
	}

	for (int i = 0; i < sh.num_subsets(); ++i)
	{
		CopySubsetIndicesToSides(sh, sh.begin<Volume>(i), sh.end<Volume>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Face>(i), sh.end<Face>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Edge>(i), sh.end<Edge>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Vertex>(i), sh.end<Vertex>(i), true);
	}
	sh.swap_subsets(0, 4);  //  somehow, these are necessary

	sh.set_subset_name("cyt", 0);
	sh.set_subset_name("er", 1);
	sh.set_subset_name("app", 2);
	sh.set_subset_name("mem_cyt", 3);
	sh.set_subset_name("mem_er", 4);
	sh.set_subset_name("mem_app", 5);
	sh.set_subset_name("syn", 6);
	sh.set_subset_name("measZone_dend", 7);
	sh.set_subset_name("measZone_neck", 8);
	sh.set_subset_name("measZone_head", 9);

	if (!opt_ER)
	{
		sh.join_subsets(0, 0, 1, false);
		sh.join_subsets(0, 0, 4, false);
	}

	// colorize the subsets
	for (int i = 0; i < sh.num_subsets(); ++i)
	{
		SubsetInfo& si = sh.subset_info(i);
		vector3 col = GetColorFromStandardPalette(i);
		si.color.x() = col.x();
		si.color.y() = col.y();
		si.color.z() = col.z();
		si.color.w() = 1.f;
	}

	// write file
	if (!SaveGridToFile(grid, sh, fileName.c_str()))
		UG_THROW("Error while saving dendrite to file '" << fileName << "'");
}

} // namespace ug
