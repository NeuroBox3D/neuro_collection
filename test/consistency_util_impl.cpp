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

#include "consistency_util.h"
#include "neurite_runtime_error.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/range/algorithm/for_each.hpp> 
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <sstream>

namespace ug
{
    namespace neuro_collection
    {
   		////////////////////////////////////////////////////////////////////////
		/// CHECK_DIAMETER_VARIABILIY
		////////////////////////////////////////////////////////////////////////
        void check_diameter_variability
        (
            const FRAGMENTS &fragments, 
            const number eps, 
            const int running_window,
            const bool ignore_first_vertex_after_bp
        ) 
        {
            using namespace boost::accumulators;
            using namespace boost;
            const int offset = ignore_first_vertex_after_bp ? 1 : 0;
            const int n = fragments.second.size(); 
            /// for each fragment
            for (int i = 0; i < n; i++) {
                const int m = fragments.second[i].size();
                /// for all fragment's sections
                for (int j = offset; j < m-offset-running_window; j+=running_window) {
                    std::vector<number> temp;
                    for (int k = 0; k < running_window; k++) {
                       temp.push_back(fragments.second[i][j+k]);
                    }
                    accumulator_set<number, stats<tag::mean, tag::variance> > acc;
                    for_each(temp.begin(), temp.end(), bind<void>(ref(acc), _1));
                    if (variance(acc) > (mean(acc) * eps)) {
                        std::stringstream ss;
                        ss << "Variance of the diameters is " << variance(acc) 
                            << "µm along some neurite of the specified geometry "
                            << "(with a running window of #" << running_window 
                            << " sections) " << " exceeds allowed tolerance of " 
                            << eps*100 << " [%].";
                        throw HighDiameterVariability(ss.str());
                    }
                 }
            }
        }

   		////////////////////////////////////////////////////////////////////////
		/// CHECK_FOR_CLOSE_BRANCHING_POINTS
		////////////////////////////////////////////////////////////////////////
        void check_for_close_branching_points
        (
            const FRAGMENTS &fragments, 
            const number eps
        )
        {
            const int n = fragments.first.size();
            for (int i = 0; i < n; i++) {
                const int m = fragments.second[i].size();
                number dist = 0;
                for (int j = 0; j < m-1; j++) {
                    dist += VecDistance(fragments.first[i][j], fragments.first[i][j+1]);
                }
                const number threshold = 2 * (fragments.second[i].front() + fragments.second[i].back());
                if (threshold > dist * (1-eps)) {
                    std::stringstream ss;
                    ss << "Close by branching points detected in specified " <<
                    "geometry with distance " << dist << " µm which is below the allowed " <<
                    "minimum distance threshold of " << threshold << " [µm] including an " <<
                    "added safety margin of " << 100*eps << " [%].";
                    throw BranchingPointClustering(ss.str());
                }
            }
        }

   		////////////////////////////////////////////////////////////////////////
		/// CHECK_FOR_SMALL_RADII
		////////////////////////////////////////////////////////////////////////
        void check_for_small_radii
        (
            const FRAGMENTS &fragments, 
            const number eps
        )
        {
            const int n = fragments.first.size();
            for (int i = 0; i < n; i++) {
                const int m = fragments.second[i].size();
                for (int j = 0; j < m; j++) {
                    const number radius = fragments.second[i][j];
                    if (radius < SMALL) {
                       std::stringstream ss;
                       ss << "Radius of section small or negative: r=" 
                          << radius << " for point p=" << fragments.first[i][j] << ".";
                       throw SmallOrNegativeRadius(ss.str());
                    }
                }
            }
        }
    }
}