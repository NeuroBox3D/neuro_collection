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

namespace ug {
	namespace neuro_collection {
		////////////////////////////////////////////////////////////////////////
		/// polygonal_mesh_from_txt
		////////////////////////////////////////////////////////////////////////
		void polygonal_mesh_from_txt(const std::string& fileName) {
			/// setup grid
			Grid g;
			SubsetHandler sh(g);
			sh.subset_info(0).name = "polygon";
			sh.set_default_subset_index(0);
			g.attach_to_vertices(aPosition);
			Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

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

			/// save grid
			AssignSubsetColors(sh);
			std::string outFileName = FilenameWithoutExtension(fileName) + ".ugx";
			SaveGridToFile(g, sh, outFileName.c_str());
		}
	} // neuro_collection
} // ug
