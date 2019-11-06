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
#include <lib_algebra/vector_interface/vec_functions.h>

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
		/// SignedAngleBetweenDirectionsInPlane
		////////////////////////////////////////////////////////////////////////
		number SignedAngleBetweenDirsInPlane
		(
			const vector3& p,
			const vector3& n,
			const vector3& s
		) {
			/// TODO: replace VecCross with VecDot and take sign of this expression
			vector3 cross;
			VecCross(cross, s, p);
			UG_COND_WARNING(VecLength(cross) < SMALL, "Potential null-vector encountered. "
					"Expect bogus output / behaviour from the underlying algorithm relying"
					"on the angle.")

			/// TODO replace and test if vectors are parallel or antiparallel => 180 or -180 degree
			/*
			if (VecLength(cross) < SMALL) {
				return 180;
			}
			*/

			/*
			int signum = boost::math::sign(VecDot(n, cross) / (VecLength(n) * VecLength(cross)));
			number mod = fmod(signum * rad_to_deg(AngleBetweenDirections(s, p)), 360);

			if (signum == -1) {
				return fmod(mod + 360, 360);
			} else {
				return mod;
			}*/

			int signum = boost::math::sign(VecDot(n, cross));
			number mod = fmod(signum * rad_to_deg(AngleBetweenDirections(s, p)), 360);

			if (signum == -1) {
				return fmod(mod + 360, 360);
			} else {
				return mod;
			}

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

		ostream& operator<< (ostream &os, const Cylinder& cyl) {
		    return (os << "Cylinder (" << cyl.c << ", " << cyl.w << " with radius "
		    		<< cyl.r << " and height " << cyl.h << std::endl);
		}


		bool SeparatedByFirstCylindersUnitAxis
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			vector3 cross;
			VecCross(cross, a.w, b.w);
			number dot = VecDot(a.w, b.w);

			vector3 diff;
			VecSubtract(diff, a.c, b.c);
			vector1 awScaled = VecDot(a.w, diff);
			return std::fabs(b.r*VecNorm2(cross)+a.h/2+b.h/2*dot-VecNorm2(awScaled)) < SMALL;
		}


		bool SeparatedBySecondCylindersUnitAxis
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			vector3 cross;
			VecCross(cross, a.w, b.w);
			number dot = VecDot(a.w, b.w);

			vector3 diff;
			VecSubtract(diff, a.c, b.c);
			vector1 awScaled = VecDot(b.w, diff);
			return std::fabs(a.r*VecNorm2(cross)+a.h/2+b.h/2*dot-VecNorm2(awScaled)) < SMALL;
		}

		bool SeparatedByUnitAxesCrossed
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			vector3 cross;
			VecCross(cross, a.w, b.w);
			vector3 diff;
			VecSubtract(diff, a.c, b.c);
			vector1 dot = VecDot(cross, diff);
			return std::fabs((a.r+b.r)*VecNorm2(cross)-VecNorm2(dot)) < SMALL;
		}

		bool SeparatedByFirstCylinderPerpendicular
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			/// TODO: Implement: This will be more involving as an iterative method, bisect is used
			UG_THROW("Not implemented!");
			return false;
		}

		bool SeparatedBySecondCylinderPerpendicular
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			/// TODO: Implement: This will be more involving as an iterative method, bisect is used
			UG_THROW("Not implemented!");
			return false;

		}

		bool SeparatedByOtherDirections
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			/// TODO: Implement: This will be more involving as an iterative method, e.g. CG is used
			UG_THROW("Not implemented!");
			return false;
		}

		bool SeparatedByHeight
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			vector3 diff;
			VecSubtract(diff, a.c, b.c);
			vector1 dot = VecDot(a.w, diff);
		    return std::fabs((a.h/2+b.h/2 - VecNorm2(dot)) < SMALL);
		}

		bool SeparatedRadially
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			vector3 diff;
			VecSubtract(diff, a.c, b.c);
			vector3 diffdot;
			vector3 dotVec;
			VecScale(dotVec, a.w, VecDot(a.w, diff));
			VecSubtract(diffdot, diff, dotVec);
			return std::fabs((a.r+b.r-VecNorm2(diffdot)) < SMALL);
		}



		bool CylinderCylinderSeparationTest
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			/// test functions
			const size_t numPar = 2;
			const size_t numOther = 6;
			bool (*cylTestPar[numPar])(const Cylinder&, const Cylinder&);
			bool (*cylTestOther[numOther])(const Cylinder&, const Cylinder&);
			cylTestPar[0] = SeparatedByHeight;
			cylTestPar[1] = SeparatedRadially;
			cylTestOther[0] = SeparatedByFirstCylindersUnitAxis;
			cylTestOther[1] = SeparatedBySecondCylindersUnitAxis;
			cylTestOther[2] = SeparatedByUnitAxesCrossed;
			cylTestOther[3] = SeparatedByFirstCylinderPerpendicular;
			cylTestOther[4] = SeparatedBySecondCylinderPerpendicular;
			cylTestOther[5] = SeparatedByOtherDirections;

			vector3 cross;
			VecCross(cross, a.w, b.w);
			if (VecNorm2(cross) > 0) {
				// All other directions
				for (size_t i = 0; i < numOther; i++) {
					if ((*cylTestOther[i])(a, b)) return true;
				}
			} else {
				/// Parallel cases
				for (size_t i = 0; i < numPar; i++) {
					if ((*cylTestPar[i])(a, b)) return true;
				}
			}
			return false;
		}
	}
}
