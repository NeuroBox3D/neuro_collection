/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-08-09
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__DENDRITE_GENERATOR_H
#define UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__DENDRITE_GENERATOR_H

#include <cstddef>                                               // for size_t
#include <string>                                                // for string

#include "common/types.h"                                        // for number


namespace ug {
namespace neuro_collection {


class DendriteGenerator
{
	public:
		DendriteGenerator();

		void set_dendrite_length(number l);
		void set_dendrite_radius(number r);
		void set_er_radius(number r);
		void set_synapse_area(number a);
		void set_num_segments(size_t n);

		number num_segments() const;

		/**
		 * @brief set ER to not be contiguous, but consist of many parts
		 *
		 * The ER is then created as an alternating sequence of ER and holes (cytosol),
		 * therefore, the axial length of ER and holes (in units of segment length)
		 * must be provided.
		 * It is recommended to use set_num_segments() beforehand.
		 *
		 * @note This is not implemented for the generation of 1d geometries or geometries
		 *       with discrete RyR channels.
		 *
		 * @param erSegLength    number of segments making up one block of ER
		 * @param holeSegLength  number of segments making up one block of hole
		 */
		void set_bobbel_er(size_t erSegLength, size_t holeSegLength);

		/// creates a 2d rotationally symmetric dendrite
		void create_dendrite_middle_influx(const std::string& filename);

		/// creates a 2d rotationally symmetric dendrite (with influx to the left)
		void create_dendrite(const std::string& filename);

		/// creates a 1d dendrite (with influx to the left)
		void create_dendrite_1d(const std::string& filename);

		/// creates a 1d dendrite (with influx to the left)
		void create_dendrite_discreteRyR(const std::string& filename, number channelDistance);

	private:
		number m_dendrite_length;
		number m_dendrite_radius;
		number m_er_radius;
		number m_synapse_width;
		size_t m_numSegments;
		bool m_bNumSegSet;

		bool m_bBobbelER;
		size_t m_numERBlockSegments;
		size_t m_numHoleBlockSegments;

		const int CYT_SI;
		const int ER_SI;
		const int PM_SI;
		const int ERM_SI;
		const int RYR_SI;
		const int SYN_SI;
		const int BND_CYT_SI;
		const int BND_ER_SI;
};



} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__GRID_GENERATION__DENDRITE_GENERATOR_H
