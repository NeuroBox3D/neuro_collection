/*
 * Copyright (c) 2009-2021: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2021-02-10
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__CONSISTENCY_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__CONSISTENCY_UTIL_H

#include "lib_grid/grid/grid.h"
#include "lib_grid/grid/geometry.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"
#include "common/math/ugmath_types.h"
#include "types.h"

namespace ug {
	namespace neuro_collection {
        /// FRAGMENTS: Pair<Position, Radius>
        typedef std::pair<std::vector<std::vector<vector3> >, std::vector<std::vector<number> > > FRAGMENTS;
        /*!
         * \brief check if diameters of neurite fragments do vary significantly
         * \param[in] fragments a list of neurite fragments
         * \param[in] eps threshold
         * \param[in] running_window number of section to calculate variance
         * \return \c bool
         */
        bool check_diameter_variability
        (
            const FRAGMENTS& fragments, 
            const number eps=0.20, 
            const int running_window=4,
            const bool ignore_first_vertex_after_bp=true
        );

        /*!
         * \brief check if branching points 
         * Branching points are said to be close with respect 
         * to their average radii of the start and end vertex
         * and some added safety margin (eps) in percentage
         * \param[in] fragments a list of neurite fragments
         * \param[in] eps safety margin (default: 0.20)
         * \return \c bool
         */
        bool check_for_close_branching_points
        (
           const FRAGMENTS& fragments, 
           const number eps=0.20
        );
    }
}

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__CONSISTENCY_UTIL_H