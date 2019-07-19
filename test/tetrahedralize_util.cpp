/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-07-19
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

#include <vector>
#include <sstream>
#include "tetrahedralize_util.h"
#include <lib_grid/algorithms/grid_generation/tetrahedralization.h>

#ifdef UG_TETGEN
	#include "tetgen.h"
#endif

using namespace std;

namespace ug {
	namespace neuro_collection {
        ////////////////////////////////////////////////////////////////////////
        /// Tetrahedralize
        ////////////////////////////////////////////////////////////////////////
		bool Tetrahedralize
		(
			Selector& sel,
		    Grid& grid,
		    ISubsetHandler* pSH,
			number quality,
			bool preserveBnds,
			bool preserveAll,
			APosition& aPos,
			int verbosity
		)
		{
			#ifdef UG_TETGEN
				UG_COND_THROW(!grid.has_vertex_attachment(aPos), "Grid has no position attachment.");
				Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

				/// TODO: store old quadrilaterals - will be needed later on to be added to the grid again
				if (sel.num<Quadrilateral>() > 0) {
					Triangulate(grid, sel.begin<Quadrilateral>(), sel.end<Quadrilateral>(), &aaPos);
				}

				if(sel.num<Edge>() > 0 || sel.num<Face>() > 0) {
					size_t numVrtsRemoved = 0;
					for (VertexIterator iter = sel.begin<Vertex>();
							iter != sel.end<Vertex>();)
					{
						Vertex* v = *iter;
						++iter;
						if ((NumAssociatedEdges(grid, v) == 0) ||
							(NumAssociatedFaces(grid, v) == 0))
						{
							grid.erase(v);
							++numVrtsRemoved;
						}
					}

					UG_COND_LOGN(numVrtsRemoved, "WARNING in Tetrahedralize: Removed " <<
							numVrtsRemoved << " vertices which were not connected to any "
							<< "edge or face.");
				}

				RemoveDuplicates(grid, sel.begin<Face>(), sel.end<Face>());

				//	attach an index to the vertices
				AInt aInd;
				grid.attach_to_vertices(aInd);
				Grid::VertexAttachmentAccessor<AInt> aaInd(grid, aInd);

				//	datastructures to communicate with tetgenio
				tetgenio in, out;

				//	setup points
				{
					in.numberofpoints = sel.num<Vertex>();
					in.pointlist = new REAL[in.numberofpoints*3];

					//	copy position data
					int counter = 0;

					for (VertexIterator iter = sel.vertices_begin();
							iter != sel.vertices_end(); ++iter, ++counter)
					{
						aaInd[*iter] = counter;
						vector3& v = aaPos[*iter];
						// position types are casted to float, since this circumvents an
						// error that occurs on some geometries. Somehow tetgen constructs
						// selfintersecting facets otherwise (sometimes). I didn't really understand
						// this behaviour yet.
						in.pointlist[counter * 3] = (float)v.x();
						in.pointlist[counter * 3 + 1] = (float)v.y();
						in.pointlist[counter * 3 + 2] = (float)v.z();
					}
				}

				/// setup facets
				{
					if(sel.num<Face>() > 0) {
						in.numberoffacets = sel.num<Face>();
						in.facetlist = new tetgenio::facet[in.numberoffacets];
						in.facetmarkerlist = new int[in.numberoffacets];

						int counter = 0;
						for (FaceIterator iter = sel.faces_begin();
								iter != sel.faces_end(); ++iter, ++counter)
						{
							Face* f = *iter;

							//	copy the face to the facetlist
							tetgenio::facet* tf = &in.facetlist[counter];
							tf->numberofpolygons = 1;
							tf->polygonlist = new tetgenio::polygon[tf->numberofpolygons];
							tf->numberofholes = 0;
							tf->holelist = NULL;
							tetgenio::polygon* p = &tf->polygonlist[0];
							p->numberofvertices = f->num_vertices();
							p->vertexlist = new int[p->numberofvertices];
							for (int i = 0; i < p->numberofvertices; ++i)
								p->vertexlist[i] = aaInd[f->vertex(i)];

							//	set the face mark
							if(pSH)
								in.facetmarkerlist[counter] = pSH->get_subset_index(f);
							else
								in.facetmarkerlist[counter] = 0;
						}
					}
				}

				//	the aInd attachment is no longer required
				grid.detach_from_vertices(aInd);
				aaInd.invalidate();

				//	call tetrahedralization
				try {
					stringstream ss;
					ss << VerbosityToTetgenParam(verbosity);
					if(grid.num_faces() > 0){
						ss << "p";
						if(quality > SMALL)
							ss << "qq" << quality;
						if(preserveBnds || preserveAll)
							ss << "Y";
						if(preserveAll)
							ss << "Y";	// if inner bnds shall be preserved "YY" has to be passed to tetgen
					}
					ss << "Q";//"Q";
					tetrahedralize(const_cast<char*>(ss.str().c_str()), &in, &out);
				}
				catch (int errCode) {
					UG_LOGN("  aborting tetrahedralization. Received error: " << errCode);
					return false;
				}

				//	add new vertices to the grid. store all vertices in a vector.
				vector<Vertex*> vVrts(out.numberofpoints);
				{
					int counter = 0;
					//	add the old ones to the vector
					for (VertexIterator iter = sel.begin<Vertex>();
							iter != sel.end<Vertex>(); ++iter, ++counter)
					{
						aaPos[*iter].x() = out.pointlist[counter*3];
						aaPos[*iter].y() = out.pointlist[counter*3+1];
						aaPos[*iter].z() = out.pointlist[counter*3+2];
						vVrts[counter] = *iter;
						if (counter == out.numberofpoints -1){
							UG_LOGN("	WARNING: Unused points may remain!");
							break;
						}
					}

					//	create new ones and add them to the vector
					for(; counter < out.numberofpoints; ++counter)
					{
						RegularVertex* v = *grid.create<RegularVertex>();
						aaPos[v].x() = out.pointlist[counter*3];
						aaPos[v].y() = out.pointlist[counter*3+1];
						aaPos[v].z() = out.pointlist[counter*3+2];
						vVrts[counter] = v;
					}
				}

				//	erase edges if boundary segments were not preserved
				if(!preserveAll){
					grid.erase(sel.begin<Edge>(), sel.end<Edge>());
				}

				//	add new faces
				grid.erase(sel.begin<Face>(), sel.end<Face>());
				for (int i = 0; i < out.numberoftrifaces; ++i)
				{
					Triangle* tri = *grid.create<Triangle>(TriangleDescriptor(vVrts[out.trifacelist[i*3]],
												vVrts[out.trifacelist[i*3 + 1]],
												vVrts[out.trifacelist[i*3 + 2]]));

					if(pSH && out.trifacemarkerlist)
						pSH->assign_subset(tri, out.trifacemarkerlist[i]);
				}

				if(out.numberoftetrahedra < 1)
					return false;

				//	add new volumes
				for (int i = 0; i < out.numberoftetrahedra; ++i)
				{
					Tetrahedron* tet = *grid.create<Tetrahedron>(
								TetrahedronDescriptor(vVrts[out.tetrahedronlist[i*4]],
														vVrts[out.tetrahedronlist[i*4 + 1]],
														vVrts[out.tetrahedronlist[i*4 + 2]],
														vVrts[out.tetrahedronlist[i*4 + 3]]));
					if(pSH)
						pSH->assign_subset(tet, 0);
				}
				return true;
		#else
				UG_THROW("\nPerformTetrahedralization: Tetgen is not available in the "
						"current build.\nRecompile with Tetgen support to use tetrahedralization.\n");
				return false;
		#endif
		}
	}
}
