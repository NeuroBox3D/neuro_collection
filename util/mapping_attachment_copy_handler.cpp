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

#include "mapping_attachment_copy_handler.h"
#include "lib_grid/global_attachments.h" 
#include "vr_util.h"

extern ug::DebugID NC_TNP;
namespace ug
{
        namespace neuro_collection
        {
                ///////////////////////////////////////////////////////////////
                /// copy_from_other_elem_type                               
                ///////////////////////////////////////////////////////////////
                void MappingAttachmentHandler::copy_from_other_elem_type(GridObject* parent, Vertex* child)
                {
                        // ensure that parent is an edge
                        Edge* parentEdge = dynamic_cast<Edge*>(parent);
                        if (!parentEdge) {
                            return;
                        }

                        // v_1^' ..... v_1 (parentEdge(0)) --------(child vertex)---------- v_2^' (parentEdge(1)) ..... v_2
                        auto& mFrom = m_aa[(*parentEdge)[0]];
                        auto& mTo = m_aa[(*parentEdge)[1]];

                        /// Refine edge in axial direction: New vertex must be in the center
                        UG_COND_THROW(!spDom.valid(), "Domain not set up for MappingAttachmentHandler.");

                        /// Calculate (refined) child vertex center
                        auto center = GetCenter(child, spDom);

                        /// parentEdge's mapping attachment to the first vertex (parentEdge(0))
                        auto v_1_p = mFrom.v1; // start vertex of parent edge
                        auto center_0 = center; // new vertex in the middle of the parent edge (creating two new edges by split edge)

                        /// child attachment 
                        auto v_1 = center; // center of the refined coarse edge: start vertex of refined edge
                        auto v_2_p = mTo.v1; // end vertex of the refined edge: end vertex of (coarse) parent edge

                        /// populate attachments
                        NeuriteProjector::Mapping mCenter;
                        mCenter.v1 = v_1;
                        mCenter.v2 = v_2_p;

                        /// Orientation of edge not known a-priori, override the best guess
                        if (std::fabs(VecDistance(mCenter.v1, mCenter.v2)) < SMALL) {
                           mCenter.v2 = mTo.v2;
                        }

                        /// Debug coordinates
                        IF_DEBUG(NC_TNP, 0) {
                           UG_LOGN("v1: " << mCenter.v1); UG_LOGN("v2: " << mCenter.v2)
                           UG_LOGN("v_1: " << v_1); UG_LOGN("center_0:" << center_0)
                           UG_LOGN("v_1_p:" << v_1_p); UG_LOGN("v_2_p: " << v_2_p)
                           UG_LOGN("center: " << center)
                        }

                        /// Set child mapping attachment at center (between parent edge's vertices)
                        mCenter.lambda = 0.5;
                        m_aa[child] = mCenter;
                }

                ///////////////////////////////////////////////////////////////
                /// copy
                ///////////////////////////////////////////////////////////////
                void MappingAttachmentHandler::copy(Vertex* parent, Vertex* child)
		{
                        auto& mappingParent = m_aa[parent];
                        auto& mappingChild = m_aa[child];
                        auto center = GetCenter(parent, spDom);
                        /// Orientation of edge not known a-priori, make a best guess
                        mappingChild.v2 = mappingParent.v1;
                        mappingChild.v1 = center;
		}

                ///////////////////////////////////////////////////////////////
                /// AddMappingAttachmentHandlerToGrid                      
                ///////////////////////////////////////////////////////////////
                void AddMappingAttachmentHandlerToGrid(SmartPtr<Domain3d> dom)
                {
                        Attachment<NeuriteProjector::Mapping> aMapping = GlobalAttachments::attachment<Attachment<NeuriteProjector::Mapping> >("npMapping");
                        UG_COND_THROW(!dom->grid()->has_attachment<Vertex>(aMapping),
                                      "Grid does not have a 'npMapping' attachment.");
                        SmartPtr<MappingAttachmentHandler> spMah(new MappingAttachmentHandler(dom));
                        spMah->set_attachment(aMapping);
                        spMah->set_grid(dom->grid());
                }
        } // end namespace neuro_collection
} // end namespace ug
