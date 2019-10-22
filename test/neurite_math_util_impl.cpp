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

#include "neurite_math_util.h"
#include "common/error.h"
#include "common/math/misc/math_util.h"
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include <boost/algorithm/clamp.hpp>
#include <boost/math/special_functions/sign.hpp>

namespace ug {
	namespace neuro_collection {
		using namespace std;
		////////////////////////////////////////////////////////////////////////
		/// RotateVectorAroundAxis
		////////////////////////////////////////////////////////////////////////
		void RotateVectorAroundAxis
		(
			const vector3& p,
			const vector3& q,
			const vector3& o,
			number theta,
			vector3 v_rot
		) {
			// axis and vector
			vector3 v, k;
			VecSubtract(v, p, o);
			VecSubtract(k, q, o);
			VecNormalize(k, k);

			// rotated vector
			vector3 vScaleC, vScaleS, vCross, vScaledCneg;
			VecScale(vScaleC, v, cos(theta));
			VecCross(vCross, k, v);
			VecScale(vScaleS, vCross, sin(theta));
			VecScale(vScaledCneg, k, VecDot(k, v));
			VecScale(vScaledCneg, vScaledCneg, 1.0 - cos(theta));
			VecScaleAdd(v_rot, 1, v_rot, 1, vScaleC, 1, vScaleS, 1, vScaledCneg);
		}

		////////////////////////////////////////////////////////////////////////
		/// FindClosestPointToRotatedVector
		////////////////////////////////////////////////////////////////////////
		void FindClosestPointToRotatedVector
		(
			const vector3& p,
			const vector3& q,
			const vector3& o,
			const vector3& r,
			number step,
			vector3& v_closest
		) {
			number dist = numeric_limits<number>::infinity();
			for (size_t theta = 0; theta < 2*PI; theta+=step) {
				vector3 v_rot;
				RotateVectorAroundAxis(p, q, o, theta, v_rot);
				VecScaleAdd(v_rot, 1, o, -1, r);
				if (VecLength(v_rot) < dist) {
					dist = VecLength(v_rot);
					v_closest = v_rot;
				}
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// AngleBetweenDirections
		////////////////////////////////////////////////////////////////////////
		number AngleBetweenDirections
		(
			const vector3& p,
			const vector3& q
		) {
			return acos(boost::algorithm::clamp(VecDot(p, q) / (VecLength(p) * VecLength(q)), -1, 1));
		}

		////////////////////////////////////////////////////////////////////////
		/// SignedAngleBetweenDirectiosnInPlane
		////////////////////////////////////////////////////////////////////////
		number SignedAngleBetweenDirsInPlane
		(
			const vector3& p,
			const vector3& n,
			const vector3& s
		) {

			vector3 cross;
			VecCross(cross, s, p);
			int signum = boost::math::sign(VecDot(n, cross) / (VecLength(n) * VecLength(cross)));

			return signum * fmod(rad_to_deg(AngleBetweenDirections(s, p)), 360);

			/// TODO: fmod(.,.) does not behave as expected - why? Angle not in [0,360]
			/*
			 * if (signum == -1) {
			 * 		return fmod(rad_to_deg(AngleBetweenDirections(s, p)), 360) + 360
			 * } else {
			 * 		return fmod(rad_to_deg(AngleBetweenDirections(s, p)), 360);
			 * }
			 */
		}

		////////////////////////////////////////////////////////////////////////
		/// FindPermissibleRenderVector
		/// TODO: Test in grid generation method test_import_swc_general_var
		////////////////////////////////////////////////////////////////////////
		void FindPermissibleRenderVector
		(
			const std::vector<vector3>& directions,
			number minAngle,
			size_t maxIter,
			vector3& v_permissible
		) {
			number angleMin = minAngle;
			vector<vector3>::const_iterator it = directions.begin();
			bool found = false;
			size_t iter = 0;
			while (!found) {
				bool below = false;
				vector3 renderVec(rand(), rand(), rand());
				for (; it != directions.end(); ++it) {
					number angle = 0;
					if (angle < angleMin) {
						below = true;
					}
				}
				if (below) {
					renderVec = vector3(rand(), rand(), rand());
					below = false;
				} else {
					found = false;
					v_permissible = renderVec;
				}
				if (iter > maxIter) {
					angleMin -= 1.0;
					iter = 1;
				}

				UG_COND_THROW((std::abs(angleMin) <= SMALL * std::abs(angleMin)),
						"No permissible render vector found for minimal angle");
			}
		}
	}
}
