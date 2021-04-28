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

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__VR_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__VR_UTIL_H

#include <lib_disc/domain.h>
#include <lib_grid/refinement/projectors/neurite_projector.h>

namespace ug {
    namespace neuro_collection {
        /*!
         * \brief Write a 1d mesh from the corresponding 3d mesh in the domain
         * \param[in] dom domain
         * \param[in] gridLevel level of (refined) grid hierarchy
         */
        void Write3dMeshTo1d
        (
            SmartPtr<Domain3d> dom,
            size_t gridLevel
        );

        /*!
         * \brief Write a 1d mesh loaded from a 3d input mesh from storage
         * \param[in] filename name of file on storage
         * \param[in] gridLevel level of (refined) grid hierarchy
         */
        void LoadAndWrite3dMeshTo1d
        (
            const std::string& filename,
            size_t gridLevel
        );

        /*!
         * \brief Get the radius information from a 3d vertex
         * \param[in] vertex desired vertex in domain
         * \param[in] dom domain
         */ 
        number GetRadius
        (
            const Vertex* const vertex, 
            SmartPtr<Domain3d> dom
        );

        /*!
         * \brief Get the center information from a 3d vertex in domain
         * \param[in] vertex desired vertex in domain
         * \param[in] dom domain
         */
        vector3 GetCenter
        (
            const Vertex* const vertex,
            SmartPtr<Domain3d> dom
        );
    } // end namespace neuro_collection
} // end namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__VR_UTIL_H