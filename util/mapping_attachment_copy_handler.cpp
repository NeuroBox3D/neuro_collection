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

namespace ug
{
        namespace neuro_collection
        {
                ///////////////////////////////////////////////////////////////
                /// copy_from_other_elem_type                               
                ///////////////////////////////////////////////////////////////
                void MappingAttachmentHandler::copy_from_other_elem_type(GridObject *parent, Vertex *child)
                {
                        // ensure that parent is an edge
                        auto* parentEdge = dynamic_cast<Edge*>(parent);
                        if (!parentEdge) {
                            return;
                        }

                        auto& mFrom = m_aa[(*parentEdge)[0]];
                        auto& mTo = m_aa[(*parentEdge)[1]];

                        /// Refine edge in axial direction: New vertex must be in the center
                        auto center = GetCenter(child, spDom);
                        NeuriteProjector::Mapping mCenter;
                        mFrom.v2 = center;
                        mTo.v1 = center;
                        mCenter.v1 = mTo.v1;
                        mCenter.v2 = mFrom.v2;
                        mCenter.lambda = 0.5;
                        m_aa[child] = mCenter;
                }

                ///////////////////////////////////////////////////////////////
                /// AddMappingAttachmentHandlerToGrid                      
                ///////////////////////////////////////////////////////////////
                void AddMappingAttachmentHandlerToGrid(SmartPtr<Domain3d> dom)
                {
                        Attachment<NeuriteProjector::Mapping> aMapping = GlobalAttachments::attachment<Attachment<NeuriteProjector::Mapping> >("npMapping");
                        UG_COND_THROW(!dom->grid()->has_attachment<Vertex>(aMapping),
                                      "Grid does not have a 'npMapping' attachment.");
                        SmartPtr<MappingAttachmentHandler> spMah(new MappingAttachmentHandler());
                        spMah->set_attachment(aMapping);
                        spMah->set_grid(dom->grid());
                        spMah->set_domain(dom);
                }
        } // end namespace neuro_collection
} // end namespace ug
