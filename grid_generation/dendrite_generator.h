/*
 * dendrite_generator.h
 *
 *  Created on: 09.08.2017
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__DENDRITE_GENERATOR_H
#define UG__PLUGINS__NEURO_COLLECTION__DENDRITE_GENERATOR_H

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

#endif // UG__PLUGINS__NEURO_COLLECTION__DENDRITE_GENERATOR_H
