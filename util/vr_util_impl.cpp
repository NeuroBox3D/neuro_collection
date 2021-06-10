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
#include <lib_grid/refinement/projectors/projection_handler.h>
#include <map>
#include "vr_util.h"

namespace ug {
    namespace neuro_collection {
        void PostProcessMesh
        (
            const std::string& fileName
        )
        {
           	Domain3d dom;
		    try {LoadDomain(dom, fileName.c_str());}
	    	UG_CATCH_THROW("Failed loading domain from '" << fileName << "'.");
         
            auto sh = dom.subset_handler();
            auto grid = dom.grid();
            Selector sel(*grid.get());
            /// Remove all volumes and ER from mesh and clear subsets
            SelectSubset(sel, *sh.get(), 3, true);
            SelectAssociatedEdges(sel, sel.begin<Vertex>(), sel.end<Vertex>());
            SelectAssociatedFaces(sel, sel.begin<Edge>(), sel.end<Edge>());             
            SelectAssociatedVolumes(sel, sel.begin<Face>(), sel.end<Face>());
            EraseSelectedObjects(sel);
            EraseEmptySubsets(*sh.get());

            typedef NeuriteProjector::Mapping NPMapping;
	        UG_COND_THROW(!GlobalAttachments::is_declared("npMapping"),
	    		"GlobalAttachment 'npMapping' was not declared.");
	        Attachment<NPMapping> aNPMapping = GlobalAttachments::attachment<Attachment<NPMapping> >("npMapping");
	        if (!grid->has_vertex_attachment(aNPMapping)) {
	          	grid->attach_to_vertices(aNPMapping);
	        } 

            Grid::VertexAttachmentAccessor<Attachment<NPMapping> > aaMapping;
	        aaMapping.access(*grid.get(), aNPMapping);

            /// Assign all neurite caps to neurite set
            for (size_t i = 2; i < sh->num_subsets(); i++) {
                SelectSubset(sel, *sh.get(), i, true);
            }
            AssignSelectionToSubset(sel, *sh.get(), 0);
            EraseEmptySubsets(*sh.get());

            /// Save regular (quad surface mesh) and triangulated surface mesh
            SaveGridToFile(*grid.get(), *sh.get(), "test_new.ugx");
		    Triangulate(*grid.get(), grid->begin<ug::Quadrilateral>(), grid->end<ug::Quadrilateral>());
            SaveGridToFile(*grid.get(), *sh.get(), "test_new_tri.ugx");
        }

        ///////////////////////////////////////////////////////////////////////
        /// Write3dMeshTo1d
        ///////////////////////////////////////////////////////////////////////
        void Write3dMeshTo1d
        (
            SmartPtr<Domain3d> dom,
            const size_t gridLevel
        ) 
        {
            /// Get 3d mesh's mapping and surface parameters attachments
            auto aMapping = GlobalAttachments::attachment<Attachment<NeuriteProjector::Mapping> >("npMapping");
            UG_COND_THROW(!dom->grid()->has_attachment<Vertex>(aMapping), "Grid does not have a 'npMapping' attachment.");
            auto aSurfaceParams = GlobalAttachments::attachment<Attachment<NeuriteProjector::SurfaceParams> > ("npSurfParams");
            UG_COND_THROW(!dom->grid()->has_attachment<Vertex>(aSurfaceParams), "Grid does not have a 'npSurfParams' attachment.");
	    	Grid::AttachmentAccessor<Vertex, Attachment<NeuriteProjector::Mapping> > aaMapping(*dom->grid().get(), aMapping);
	    	Grid::AttachmentAccessor<Vertex, Attachment<NeuriteProjector::SurfaceParams> > aaSurfaceParams(*dom->grid().get(), aSurfaceParams);
            /// Set 1d mesh's diameter attachment
            auto aDiam = GlobalAttachments::attachment<ANumber>("diameter"); 

            /// Iterate over 3d mesh's edges to generate a 1d mesh
            ConstEdgeIterator eit = dom->grid()->begin<Edge>(gridLevel);
		    ConstEdgeIterator eit_end = dom->grid()->end<Edge>(gridLevel);

            /// Some room for optimization below
            auto sh = dom->subset_handler();

            std::map<vector3, std::vector<vector3> > edgePairs;
            std::map<vector3, std::vector<number> > fromDiam;
            std::map<vector3, std::vector<number> > toDiam;
            std::map<vector3, std::vector<int> > fromSI;
            std::map<vector3, std::vector<int> > toSI;
	    	for (; eit != eit_end; ++eit) {
                const Edge* e = *eit;
                const auto& m1 = aaMapping[e->vertex(0)];
                const auto& m2 = aaMapping[e->vertex(1)];

                /// At branching regions do not create edges between child neurites
                if (aaSurfaceParams[e->vertex(0)].neuriteID != aaSurfaceParams[e->vertex(1)].neuriteID) {
                    continue;
                }

                /// Only if v1 and v2 are not an edge of a (radial) polygon we create an axial edge along the neurite
                if (! ((m1.v1 == m2.v1) || (m1.v2 == m2.v2)) || (m1.v1 == m2.v2) || (m1.v2 == m2.v1) ) {
                    auto itFrom = std::find (edgePairs[m1.v1].begin(), edgePairs[m1.v1].end(), m1.v1);
                    auto itTo = std::find (edgePairs[m2.v1].begin(), edgePairs[m2.v1].end(), m2.v1);
                    auto itFrom2 = std::find (edgePairs[m1.v1].begin(), edgePairs[m1.v1].end(), m2.v1);
                    auto itTo2 = std::find (edgePairs[m2.v1].begin(), edgePairs[m2.v1].end(), m1.v1);
                    /// Only if v1 and v2 have not yet been created yet 
                    if ((itFrom  == edgePairs[m1.v1].end() && itTo  == edgePairs[m2.v1].end()) &&
                        (itFrom2 == edgePairs[m1.v1].end() && itTo2 == edgePairs[m2.v1].end())) {
                       edgePairs[m1.v1].push_back(m2.v1);
                       // fromDiam[m1.v1].push_back(aaSurfaceParams[e->vertex(0)].radial);
                       // toDiam[m1.v1].push_back(aaSurfaceParams[e->vertex(1)].radial);
                       fromDiam[m1.v1].push_back(GetRadius(e->vertex(0), dom));
                       toDiam[m1.v1].push_back(GetRadius(e->vertex(1), dom));
                       toSI[m1.v1].push_back(sh->get_subset_index(e->vertex(1)));
                       fromSI[m1.v1].push_back(sh->get_subset_index(e->vertex(0)));
                    }
                }
	    	}


            /// Generate 1d mesh from 3d mesh
            Grid g;
            SubsetHandler sh2(g);
            g.attach_to_vertices(aPosition);
            g.attach_to_vertices(aDiam);
            Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
            Grid::VertexAttachmentAccessor<ANumber> aaDiam(g, aDiam);

            for (auto const& edges : edgePairs)
            {
                for (size_t i = 0; i < edges.second.size(); i++) {

                    Vertex* v1;
                    Vertex* v2;
                    bool create;

                    Vertex* vv = FindClosestByCoordinate<Vertex>(edges.first, g.vertices_begin(), g.vertices_end(), aaPos);
                    if (!vv) {
                        v1 = *g.create<RegularVertex>();
                        create = true;
                    } else {
                        if (VecDistance(aaPos[vv], edges.first) > SMALL) {
                            v1 = *g.create<RegularVertex>();
                            create = true;
                        } else {
                            v1 = vv;
                        }
                    }

                    Vertex* vv2 = FindClosestByCoordinate<Vertex>(edges.second[i], g.vertices_begin(), g.vertices_end(), aaPos);
                    if (!vv2) {
                        v2 = *g.create<RegularVertex>();
                        create = true;
                    } else {
                        if (VecDistance(aaPos[vv2], edges.second[i]) > SMALL) {
                            v2 = *g.create<RegularVertex>();
                            create = true;
                        } else {
                            v2 = vv2;
                        }
                    }

                    if (create) {
                        auto* e = *g.create<RegularEdge>(EdgeDescriptor(v1, v2));
                        aaPos[v1] = edges.first;
                        aaPos[v2] = edges.second[i];
                        aaDiam[v1] = fromDiam[edges.first][i];
                        aaDiam[v2] = toDiam[edges.first][i];
                        sh2.assign_subset(v1, fromSI[edges.first][i]);
                        sh2.assign_subset(v2, toSI[edges.first][i]);
                        sh2.assign_subset(e, toSI[edges.first][i]);
                    }
                }
            }

            /// Save the mesh 
            EraseEmptySubsets(sh2);
            AssignSubsetColors(sh2);
            sh2.subset_info(0).name = "Neurite";
            sh2.subset_info(1).name = "Soma";
            std::stringstream ss;
            ss << "1dmesh" << "_" << gridLevel << ".ugx";
            SaveGridToFile(g, sh2, ss.str().c_str());
        }

        ///////////////////////////////////////////////////////////////////////
        /// Write3dMeshTo1d
        ///////////////////////////////////////////////////////////////////////
        void LoadAndWrite3dMeshTo1d
        (
            const std::string& fileName,
            const size_t gridLevel
        ) {
		    Domain3d dom;
	    	dom.create_additional_subset_handler("projSH");
    		try {LoadDomain(dom, fileName.c_str());}
    		UG_CATCH_THROW("Failed loading domain from '" << fileName << "'.");
            Write3dMeshTo1d(make_sp(&dom), gridLevel);
        }

        ///////////////////////////////////////////////////////////////////////
        /// GetRadius
        ///////////////////////////////////////////////////////////////////////
        number GetRadius
        (
            const Vertex* const vertex, 
            SmartPtr<Domain3d> dom
        ) {
            /// Get the NeuriteProjector from the domain
            auto GetNeuriteProjector = [&]
            {
                auto ph = dynamic_cast<ProjectionHandler*>(dom->refinement_projector().get());
                UG_COND_THROW(!ph, "No projection handler available in the provided domain.")
                /// Assume neurite projector is the default projector in the domain
                auto np = (NeuriteProjector*) dynamic_cast<NeuriteProjector*>(ph->default_projector().get());
                if (!np) {
                   /// If neurite projector not the default projector, try all other projectors in domain
                   for (size_t i = 0; i < ph->num_projectors(); i++) {
                       np = dynamic_cast<NeuriteProjector*>(ph->projector(i).get());
                       if (np) {
                           return np;
                       }
                   }
                }
                /// Not any neurite projector found in the list of possible projectors
                UG_COND_THROW(!np, "Neurite projector not available in the provided domain.");
                return np;
            };

            /// Get the NeuriteProjector from the domain and the surface parameters
            auto np = GetNeuriteProjector();

            /// Surface parameters
            auto aSurfParams = GlobalAttachments::attachment<Attachment<NeuriteProjector::SurfaceParams> > ("npSurfParams");
          	auto aaSurfParams = np->surface_params_accessor();
            aaSurfParams.access(*dom->grid().get(), aSurfParams);

            /// Parameters for vertex
            auto neuriteID = aaSurfParams[vertex].neuriteID;
            auto t = aaSurfParams[vertex].axial;

            const auto plainNID = (neuriteID << 12) >> 12;
            NeuriteProjector::Section cmpSec(t);
            const auto& neurite = np->neurites()[plainNID];
            auto secIt = std::lower_bound(neurite.vSec.begin(), neurite.vSec.end(), cmpSec, NeuriteProjector::CompareSections());
            UG_COND_THROW(secIt == neurite.vSec.end(), "Could not find section for parameter t = " << t << " in neurite " << neuriteID << ".");

            /// Get correct radius for the associated 3d vertex
            auto radius = 0;
            auto te = secIt->endParam;
            const auto* s = &secIt->splineParamsX[0];
            s = &secIt->splineParamsR[0];
            radius = s[0]*(te-t) + s[1];
            radius = radius*(te-t) + s[2];
            radius = radius*(te-t) + s[3];
            return radius;
        }

        ///////////////////////////////////////////////////////////////////////
        /// GetCenter
        ///////////////////////////////////////////////////////////////////////
        vector3 GetCenter
        (
            const Vertex* const vertex,
            SmartPtr<Domain3d> dom
        ) 
        {
            /// Get the NeuriteProjector from the domain
            auto GetNeuriteProjector = [&] 
            {
                auto ph = dynamic_cast<ProjectionHandler*>(dom->refinement_projector().get());
                UG_COND_THROW(!ph, "No projection handler available in the provided domain.")
                /// Assume neurite projector is the default projector in the domain
                auto np = (NeuriteProjector*) dynamic_cast<NeuriteProjector*>(ph->default_projector().get());
                if (!np) {
                   /// If neurite projector not the default projector, try all other projectors in domain
                   for (size_t i = 0; i < ph->num_projectors(); i++) {
                       np = dynamic_cast<NeuriteProjector*>(ph->projector(i).get());
                       if (np) {
                           return np;
                       }
                   }
                }
                /// Not any neurite projector found in the list of possible projectors
                UG_COND_THROW(!np, "Neurite projector not available in the provided domain.");
                return np;
            };

            /// Get the NeuriteProjector from the domain and the surface parameters
            auto np = GetNeuriteProjector();

            /// Surface parameters
            auto aSurfParams = GlobalAttachments::attachment<Attachment<NeuriteProjector::SurfaceParams> > ("npSurfParams");
            auto aaSurfParams = np->surface_params_accessor();
            aaSurfParams.access(*dom->grid().get(), aSurfParams);

            /// Parameters for vertex
            auto neuriteID = aaSurfParams[vertex].neuriteID;
            auto t = aaSurfParams[vertex].axial;

            const auto plainNID = (neuriteID << 12) >> 12;
            NeuriteProjector::Section cmpSec(t);
            const auto& neurite = np->neurites()[plainNID];
            auto secIt = std::lower_bound(neurite.vSec.begin(), neurite.vSec.end(), cmpSec, NeuriteProjector::CompareSections());
            UG_COND_THROW(secIt == neurite.vSec.end(), "Could not find section for parameter t = " << t << " in neurite " << neuriteID << ".");

            /// Get velocity of vertex
         	vector3 vel;
			const auto h = neurite.vSec[0].endParam;
			vel[0] = -3.0 * neurite.vSec[0].splineParamsX[0] * h * h
					- 2.0 * neurite.vSec[0].splineParamsX[1] * h - neurite.vSec[0].splineParamsX[2];
			vel[1] = -3.0 * neurite.vSec[0].splineParamsY[0] * h * h
					- 2.0 * neurite.vSec[0].splineParamsY[1] * h - neurite.vSec[0].splineParamsY[2];
			vel[2] = -3.0 * neurite.vSec[0].splineParamsZ[0] * h * h
					- 2.0 * neurite.vSec[0].splineParamsZ[1] * h - neurite.vSec[0].splineParamsZ[2];
			VecNormalize(vel, vel);

           	number segAxPos = t;
			for (size_t curSec = 0; curSec < neurite.vSec.size(); ++curSec) {
				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				if (sec.endParam >= segAxPos) {
                    break;
                }
            }
        
            /// Get center of vertex
            vector3 curPos;
            const auto monom = secIt->endParam - segAxPos;
            const auto* sp = &secIt->splineParamsX[0];
            auto& p0 = curPos[0];
            auto& v0 = vel[0];
            p0 = sp[0] * monom + sp[1];
            p0 = p0 * monom + sp[2];
            p0 = p0 * monom + sp[3];
            v0 = -3.0 * sp[0] * monom - 2.0 * sp[1];
            v0 = v0 * monom - sp[2];

            sp = &secIt->splineParamsY[0];
            auto& p1 = curPos[1];
            auto& v1 = vel[1];
            p1 = sp[0] * monom + sp[1];
            p1 = p1 * monom + sp[2];
            p1 = p1 * monom + sp[3];
            v1 = -3.0 * sp[0] * monom - 2.0 * sp[1];
            v1 = v1 * monom - sp[2];

            sp = &secIt->splineParamsZ[0];
            auto& p2 = curPos[2];
            auto& v2 = vel[2];
            p2 = sp[0] * monom + sp[1];
            p2 = p2 * monom + sp[2];
            p2 = p2 * monom + sp[3];
            v2 = -3.0 * sp[0] * monom - 2.0 * sp[1];
            v2 = v2 * monom - sp[2];
            
            return curPos;
        }
    } // end namespace neuro_collection
} // end namespace ug