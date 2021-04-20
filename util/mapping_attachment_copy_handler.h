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

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__MAPPING_ATTACHMENT_COPY_HANDLER_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__MAPPING_ATTACHMENT_COPY_HANDLER_H

#include "lib_grid/grid/grid_base_objects.h"
#include "lib_grid/tools/copy_attachment_handler.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"

namespace ug
{
    namespace neuro_collection
    {
        /**
         * @brief copy attachment handler for the mapping attachment
         *
         * The mapping vertex attachment needs to be propagated not only by copying from
         * base vertices to their child vertices and grand-child vertices (and so on);
         * it is also required on vertices that appear as children of refined edges.
         * 
         * Note that it is assumed, the mapping attachment handler is used only in
         * conjunction with the NeuriteAxialRefinement marker and may fail otherwise!
         */
        class MappingAttachmentHandler
            : public CopyAttachmentHandler<Vertex, Attachment<NeuriteProjector::Mapping> >
        {
        public:
            /// Ctor
            MappingAttachmentHandler(){};
            /// Dtor
            virtual ~MappingAttachmentHandler(){};

        protected:
            /*!
             * \brief Copy from other elem type to child 
             * \param[in] parent 
             * \param[out] child
             */ 
            virtual void copy_from_other_elem_type(GridObject *parent, Vertex *child);
        };

        /**
         * \brief Add the mapping attachment handler to the grid
         * \param[in] grid
         */
        void AddMappingAttachmentHandlerToGrid(SmartPtr<MultiGrid> grid);
    } // end namespace neuro_collection
} // end namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__MAPPING_ATTACHMENT_COPY_HANDLER_H