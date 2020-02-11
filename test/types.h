/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-11-13
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__TYPES_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__TYPES_H

#include <utility>
#include <lib_grid/grid/grid_base_objects.h>

namespace ug {
	namespace neuro_collection {
		/*!
		 * \brief Functor to compare edges by length
		 */
		struct EdgeLengthCompare
		{
			/*!
			 * \param[in] e1
			 * \param[in] e2
			 *
			 * \return \c bool
			 */
			bool operator()
			(
				const std::pair<Edge*, number> e1,
				const std::pair<Edge*, number> e2
			)
			{
				return e1.second > e2.second;
			}
		};

		/// SWC types
		enum swc_type
		{
			SWC_UNDF = 0, ///< UNDEFINED
			SWC_SOMA = 1, ///< SOMA
			SWC_AXON = 2, ///< AXON
			SWC_DEND = 3, ///< DENDRITE
			SWC_APIC = 4, ///< APICAL DENDRITE
			SWC_FORK = 5, ///< FORK POINT
			SWC_END = 6, ///< END POINT
			SWC_CUSTOM = 7 ///< CUSTOM
		};

		/*!
		 * \brief SWCPoint
		 * A struct representing an SWC point
		 */
		struct SWCPoint
		{
			vector3 coords; ///< coordinates
			number radius; ///< radius
			swc_type type; ///< type
			std::vector<size_t> conns; ///< connections
		};

		/*!
		 * \brief FindSWCPoint
		 * A functor which finds a point by it's coordinates
		 */
		struct FindSWCPoint
		{
			vector3 coords;

		public:
			FindSWCPoint(const vector3& coords) : coords(coords) {}

		    bool operator()(const vector3& vec)
		    {
		        return vec == coords;
		    }
		};
	}
}

#endif ///  UG__PLUGINS__NEURO_COLLECTION__TEST__TYPES_H
