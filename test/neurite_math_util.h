/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-10-17
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_MATH_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_MAHT_UTIL_H

#include "common/math/ugmath_types.h"

namespace ug {
	namespace neuro_collection {
		/*!
		 * \brief Rotates a vector (OP) around axis (OQ)
		 * \param[in] p end point of vector
		 * \param[in] q end point of axis
		 * \param[in] o origin
		 * \param[in] theta degree in radians
		 * \param[out] v_rot
		 */
		void RotateVectorAroundAxis
		(
			const vector3& p,
			const vector3& q,
			const vector3& o,
			number theta,
			vector3 v_rot
		);

		/*!
		 * \brief Finds the closest rotated vector's end point to point r
		 * \param[in] p end point of vector
		 * \param[in] q end point of axis
		 * \param[in] o origin
		 * \param[in] r point for which a closest rotated vector should be found
		 * \param[in] step resolution of search
		 * \param[ou] v_closest
		 */
		void FindClosestPointToRotatedVector
		(
			const vector3& p,
			const vector3& q,
			const vector3& o,
			const vector3& r,
			number step,
			vector3& v_closest
		);

		/*!
		 * \brief Calculates angle in radians between vectors p and q
		 *
		 * \param[in] p
		 * \param[in] q
		 * \return \c angle in radians [0, 2*pi]
		 */
		number AngleBetweenDirections
		(
			const vector3& p,
			const vector3& q
		);

		/*!
		 * \brief Calculates angle between vectors p relative to reference point
		 *
		 * \param[in] p
		 * \param[in] n normal of plane
		 * \param[in] s reference point in plane
		 * \return \c angle in degree [0, 360]
		 */
		number SignedAngleBetweenDirsInPlane
		(
			const vector3& p,
			const vector3& n,
			const vector3& s
		);

		/*!
		 * \brief Brute force strategy to find permissible render vector
		 * Note, a better strategy and probably more efficient strategy
		 * would be to translate the cartesian coordinates of the direction
		 * vectors to polar coordinates and then search on a sphere with radius
		 * of 1 on the 2d-manifold on the sphere surface a permissible render
		 * vector. Each vector will create a deadzone of minAngle on the surface
		 * and this will result in a complicated 2d-shape (polygon) to manage.
		 *
		 * \param[in] directions
		 * \param[in] minAngle
		 * \param[in] maxIter
		 * \param[out] v_permissible
		 */
		void FindPermissibleRenderVector
		(
			const std::vector<vector3>& directions,
			number minAngle,
			size_t maxIter,
			vector3& v_permissible
		);
	}
}

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_MATH_UTIL_H
