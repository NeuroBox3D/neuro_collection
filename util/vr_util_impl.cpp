/*
 * Copyright (c) 2009-2021: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2021-04-20
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

#include <lib_grid/grid/grid.h>
#include <lib_disc/domain_util.h>
#include <lib_grid/global_attachments.h>
#include <lib_grid/refinement/projectors/neurite_projector.h>
#include <lib_grid/file_io/file_io.h>
#include <map>

namespace ug {
    namespace neuro_collection {
        ///////////////////////////////////////////////////////////////////////
        /// Write3dMeshTo1D
        ///////////////////////////////////////////////////////////////////////
        void Write3dMeshTo1d(const std::string& filename) {
            /// Load 3d mesh
	        Domain3d dom;
        	try {LoadDomain(dom, filename.c_str());}
        	UG_CATCH_THROW("Failed loading domain from '" << filename << "'.");

            /// Get 3d mesh's mapping and surface parameters attachments
            Attachment<NeuriteProjector::Mapping> aMapping = GlobalAttachments::attachment<Attachment<NeuriteProjector::Mapping> >("npMapping");
            UG_COND_THROW(!dom.grid().get()->has_attachment<Vertex>(aMapping), "Grid does not have a 'npMapping' attachment.");
            Attachment<NeuriteProjector::SurfaceParams> aSurfaceParams = GlobalAttachments::attachment<Attachment<NeuriteProjector::SurfaceParams> > ("npSurfParams");
            UG_COND_THROW(!dom.grid().get()->has_attachment<Vertex>(aSurfaceParams), "Grid does not have a 'npSurfParams' attachment.");
	    	Grid::AttachmentAccessor<Vertex, Attachment<NeuriteProjector::Mapping> > aaMapping(*dom.grid(), aMapping);
	    	Grid::AttachmentAccessor<Vertex, Attachment<NeuriteProjector::SurfaceParams> > aaSurfaceParams(*dom.grid(), aSurfaceParams);
            /// Set 1d mesh's diameter attachment
            ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter"); 

            /// Iterate over 3d mesh's edges to generate a 1d mesh
            ConstEdgeIterator eit = dom.grid().get()->begin<Edge>();
		    ConstEdgeIterator eit_end = dom.grid().get()->end<Edge>();

            dom.grid().get()->attach_to_vertices(aPosition);
            Grid::VertexAttachmentAccessor<APosition> aaPos2(*dom.grid().get(), aPosition);

            /// Some room for optimization below
            MGSubsetHandler* sh = dom.subset_handler().get();
            std::map<vector3, std::vector<vector3> > edgePairs;
            std::map<vector3, std::vector<int> > fromSI;
            std::map<vector3, std::vector<int> > toSI;
            std::map<vector3, std::vector<number> > fromDiam;
            std::map<vector3, std::vector<number> > toDiam;
	    	for (; eit != eit_end; ++eit) {
                const Edge* e = *eit;
                const NeuriteProjector::Mapping& m1 = aaMapping[e->vertex(0)];
                const NeuriteProjector::Mapping& m2 = aaMapping[e->vertex(1)];

                /// At branching regions, we generate still one wrong edge
                if (aaSurfaceParams[e->vertex(0)].neuriteID != aaSurfaceParams[e->vertex(1)].neuriteID) {
                    /// TODO: FIXME
                    continue;
                }

                /// If v1 and v2 are not on an edge of a polygon (radial)
                if (! ((m1.v1 == m2.v1) || (m1.v2 == m2.v2)) || (m1.v1 == m2.v2) || (m1.v2 == m2.v1) ) {
                    edgePairs[m1.v1].push_back(m2.v1);
                    fromSI[m1.v1].push_back(sh->get_subset_index(e->vertex(0)));
                    toSI[m1.v1].push_back(sh->get_subset_index(e->vertex(1)));
                    fromDiam[m1.v1].push_back(aaSurfaceParams[e->vertex(0)].radial);
                    toDiam[m1.v1].push_back(aaSurfaceParams[e->vertex(1)].radial);
                }
	    	}

            /// Produce 1d mesh from 3d mesh
            Grid g;
            SubsetHandler sh2(g);
            g.attach_to_vertices(aPosition);
            g.attach_to_vertices(aDiam);
            Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
            Grid::VertexAttachmentAccessor<ANumber> aaDiam(g, aDiam);

            for (auto const& edges : edgePairs)
            {
                for (size_t i = 0; i < edges.second.size(); i++) {
                    RegularVertex* v1 = *g.create<RegularVertex>();
                    RegularVertex* v2 = *g.create<RegularVertex>();
                    RegularEdge* e = *g.create<RegularEdge>(EdgeDescriptor(v1, v2));

                    aaPos[v1] = edges.first;
                    aaPos[v2] = edges.second[i];
                    aaDiam[v1] = fromDiam[edges.first][i];
                    aaDiam[v2] = toDiam[edges.first][i];
                    /*
                    sh2.assign_subset(v1, fromSI[edges.first]);
                    sh2.assign_subset(v2, toSI[edges.second][i]);
                    sh2.assign_subset(e, toSI[edges.second][i]);
                    */
                }
            }

            /// Save the mesh
            SaveGridToFile(g, sh2, "1dmesh.ugx");
        }
    }
}