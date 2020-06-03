/*
 * Copyright (c) 2009-2020: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2020-05-19
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

#include "polygonal_mesh_from_txt.h"
#include "lib_grid/lib_grid.h"
#include <fstream>
#include "common/math/ugmath_types.h"
#include "lib_grid/refinement/regular_refinement.h"
#include "lib_grid/algorithms/remeshing/delaunay_triangulation.h"

namespace ug {
	namespace neuro_collection {
		////////////////////////////////////////////////////////////////////////
		/// polygonal_mesh_from_txt
		////////////////////////////////////////////////////////////////////////
		void polygonal_mesh_from_txt
		(
			const std::string& fileName,
			const ug::vector3& p1,
			const ug::vector3& p2,
			const ug::vector3& p3,
			const ug::vector3& p4
		) {
			/// setup grid
			Grid g;
			SubsetHandler sh(g);
			sh.set_default_subset_index(0);
			g.attach_to_vertices(aPosition);
			Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
			AInt aInt;
			g.attach_to_vertices(aInt);
			const size_t numRefs = 2;
			const number minAngle = 20;

			/// create vertices
			double x, y;
			const double z = 0;
			std::vector<Vertex*> vertices;
			std::ifstream infile(fileName.c_str());
			while (infile >> x >> y)
			{
				Vertex* v = *g.create<RegularVertex>();
				aaPos[v] = vector3(x, y, z);
				vertices.push_back(v);
			}

			/// create edges
			std::vector<Vertex*>::const_iterator vit = vertices.begin();
			for (size_t i = 0; i < vertices.size()-1; i++) {
				*g.create<RegularEdge>(EdgeDescriptor(vertices[i], vertices[i+1]));
			}
			*g.create<RegularEdge>(EdgeDescriptor(vertices.front(), vertices.back()));

			/// rectangle
			ug::Vertex* v1 = *g.create<RegularVertex>();
			aaPos[v1] = p1; // 6
			ug::Vertex* v2 = *g.create<RegularVertex>();
			aaPos[v2] = p2; // 7
			ug::Vertex* v3 = *g.create<RegularVertex>();
			aaPos[v3] = p3; // 8
			ug::Vertex* v4= *g.create<RegularVertex>();
			aaPos[v4] = p4; // 9
			ug::RegularEdge* e1 = *g.create<RegularEdge>(EdgeDescriptor(v1, v2));
			sh.assign_subset(e1, 3);
			ug::RegularEdge* e2 = *g.create<RegularEdge>(EdgeDescriptor(v1, v3));
			sh.assign_subset(e2, 2);
			ug::RegularEdge* e3 = *g.create<RegularEdge>(EdgeDescriptor(v2, v4));
			sh.assign_subset(e3, 5);
			ug::RegularEdge* e4 = *g.create<RegularEdge>(EdgeDescriptor(v3, v4));
			sh.assign_subset(e4, 4);
			sh.assign_subset(v1, 2);
			sh.assign_subset(v2, 3);
			sh.assign_subset(v3, 4);
			sh.assign_subset(v4, 5);

			/// refinement
			Selector sel(g);
			for (int i = 0; i < sh.num_subsets(); i++) {
				SelectSubset(sel, sh, i, true);
			}

			for (size_t i = 0; i < numRefs; i++) {
				Refine(g, sel);
			}
			sel.clear();

			// triangulation
			sh.set_default_subset_index(7);
			TriangleFill_SweepLine(g, g.edges_begin(), g.edges_end(), aPosition, aInt, &sh, 7);
			QualityGridGeneration(g, g.faces_begin(), g.faces_end(), aaPos, minAngle);

			/// Separate faces
			SelectSubset(sel, sh, 0, true);
			SeparateSubsetsByLowerDimSelection<Face>(g, sh, sel);

			/// Subset assignment
			sel.clear();
			SelectSubset(sel, sh, 7, true);
			AssignSelectionToSubset(sel, sh, 0);
			sel.clear();
			SelectSubset(sel, sh, 1, true);
			CloseSelection(sel);
			AssignSelectionToSubset(sel, sh, 1);
			sel.clear();
			SelectAreaBoundary(sel, sh.begin<Edge>(1), sh.end<Edge>(1));
			SelectAreaBoundary(sel, sh.begin<Face>(1), sh.end<Face>(1));
			CloseSelection(sel);
			AssignSelectionToSubset(sel, sh, 8);
			sel.clear();

			/// Colors and names
			EraseEmptySubsets(sh);
			AssignSubsetColors(sh);
			sh.subset_info(0).name = "vol";
			sh.subset_info(1).name = "tower";
			sh.subset_info(2).name = "left";
			sh.subset_info(3).name = "bottom";
			sh.subset_info(4).name = "top";
			sh.subset_info(5).name = "right";
			sh.subset_info(6).name = "tower bnd";

			/// save grid
			std::string outFileName = FilenameWithoutExtension(fileName) + ".ugx";
			SaveGridToFile(g, sh, outFileName.c_str());
		}
	} // neuro_collection
} // ug
