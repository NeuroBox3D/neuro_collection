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

#include "common/error.h"
#include "common/math/misc/math_util.h"
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include <lib_algebra/vector_interface/vec_functions.h>
#include <queue>
#include <boost/algorithm/clamp.hpp>
#include <boost/math/special_functions/sign.hpp>
#include "neurite_math_util.h"
#include <lib_grid/algorithms/space_partitioning/lg_ntree.h>
#include "lib_disc/domain_util.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/algorithms/remeshing/resolve_intersections.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <ctime>

extern ug::DebugID NC_TNP;

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
			const number theta,
			vector3& v_rot
		)
		{
			// define axis and vector
			vector3 v, k;
			VecSubtract(v, p, o);
			VecSubtract(k, q, o);
			VecNormalize(k, k);

			// rotated vector with Rodrigues' formula
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
		/// TODO: Test in grid generation method test_import_swc_general_var
		////////////////////////////////////////////////////////////////////////
		void FindClosestPointToRotatedVector
		(
			const vector3& p,
			const vector3& q,
			const vector3& o,
			const vector3& r,
			const number step,
			vector3& v_closest
		)
		{
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
		)
		{
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
		)
		{
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
		////////////////////////////////////////////////////////////////////////
		void FindPermissibleRenderVector
		(
			const std::vector<vector3>& directions,
			const number minAngle,
			const size_t maxIter,
			vector3& v_permissible
		)
		{
			typedef boost::minstd_rand base_generator_type;
			base_generator_type generator(std::time(NULL));
			boost::uniform_real<> uni_dist(0, 1);
			boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
			number angleMin = minAngle;
			bool found = false;
			size_t iter = 0;
			size_t upperLimit = 0;
			while (!found) {
				bool below = false;
				vector3 renderVec(uni(), uni(), uni());
				for (size_t i = 0; i < directions.size(); i++) {
					number angle = rad_to_deg(AngleBetweenDirections(renderVec, directions[i]));
					UG_DLOGN(NC_TNP, 0, "angle: " << angle);
					if ( (angle < angleMin) || (angle > (180-angleMin))) {
						below = true;
					}
				}

				iter++;

				if (below) {
					renderVec = vector3(uni(), uni(), uni());
					below = false;
				} else {
					found = true;
					v_permissible = renderVec;
					UG_DLOGN(NC_TNP, 0, "Found render vec: " << v_permissible);
					UG_LOGN("Found render vector")
				}
				if (iter > maxIter) {
					angleMin -= 1.0;
					iter = 0;
					upperLimit++;
				}
				if (upperLimit == maxIter) {
					break;
				}
			}

			UG_COND_THROW(!found, "No permissible render vector could be found with "
					"parameters (minAngle/maxIter): " << minAngle << ", " << maxIter);
		}

		////////////////////////////////////////////////////////////////////////
		/// SeparatedByFirstCylindersUnitAxis
		////////////////////////////////////////////////////////////////////////
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


		////////////////////////////////////////////////////////////////////////
		/// SeparatedBySecondCylindersUnitAxis
		////////////////////////////////////////////////////////////////////////
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

		////////////////////////////////////////////////////////////////////////
		/// SeparatedByUnitAxesCrossed
		////////////////////////////////////////////////////////////////////////
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

		////////////////////////////////////////////////////////////////////////
		/// SeparatedByFirstCylinderPerpendicular
		////////////////////////////////////////////////////////////////////////
		bool SeparatedByFirstCylinderPerpendicular
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			UG_WARNING("Not yet implemented!");
			return false;
		}

		////////////////////////////////////////////////////////////////////////
		/// SeparatedBySecondCylinderPerpendicular
		////////////////////////////////////////////////////////////////////////
		bool SeparatedBySecondCylinderPerpendicular
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			UG_WARNING("Not yet implemented!");
			return false;

		}

		////////////////////////////////////////////////////////////////////////
		/// SeparatedByOtherDirections
		////////////////////////////////////////////////////////////////////////
		bool SeparatedByOtherDirections
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			UG_WARNING("Not yet implemented!");
			return false;
		}

		////////////////////////////////////////////////////////////////////////
		/// SeparatedByHeight
		////////////////////////////////////////////////////////////////////////
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

		////////////////////////////////////////////////////////////////////////
		/// SeparatedRadially
		////////////////////////////////////////////////////////////////////////
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

		////////////////////////////////////////////////////////////////////////
		/// CylinderCylinderSeparationTest
		////////////////////////////////////////////////////////////////////////
		bool CylinderCylinderSeparationTest
		(
			const Cylinder& a,
			const Cylinder& b
		) {
			/// test functions
			const size_t numPar = 2;
			const size_t numOther = 6;
			bool (*cylTestPar[numPar])(const Cylinder&, const Cylinder&) =
			{
					SeparatedByHeight,
					SeparatedRadially
			};
			bool (*cylTestOther[numOther])(const Cylinder&, const Cylinder&) =
			{
					SeparatedByFirstCylindersUnitAxis,
					SeparatedBySecondCylindersUnitAxis,
					SeparatedByUnitAxesCrossed,
					SeparatedByUnitAxesCrossed,
					SeparatedBySecondCylinderPerpendicular,
					SeparatedByOtherDirections
			};

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

		////////////////////////////////////////////////////////////////////////
		/// operator<< for Cylinder
		////////////////////////////////////////////////////////////////////////
		ostream& operator<< (ostream &os, const Cylinder& cyl) {
		    return (os << "Cylinder (" << cyl.c << ", " << cyl.w << " with radius "
		    		<< cyl.r << " and height " << cyl.h << std::endl);
		}

		////////////////////////////////////////////////////////////////////////
		/// is_cyclic
		////////////////////////////////////////////////////////////////////////
		bool is_cyclic
		(
			const vector<vector<int> > adj,
			const int s,
            const int V,
            vector<bool>& visited
        ) {
			vector<int> parent(V, -1);
			std::queue<int> q;
			visited[s] = true;
			q.push(s);

			while (!q.empty()) {

			   int u = q.front();
			   q.pop();

			   for (size_t i = 0; i < adj[u].size(); i++) {
			      if (!visited[adj[u][i]]) {
			           visited[adj[u][i]] = true;
			            q.push(adj[u][i]);
			            parent[adj[u][i]] = u;
			       }
			       else if (parent[u] != adj[u][i])
			             return true;
			       }
			    }
			   return false;
		}

		////////////////////////////////////////////////////////////////////////
		/// is_cyclic
		////////////////////////////////////////////////////////////////////////
		bool is_cyclic
		(
			const vector<vector<int> > adj,
			const int V
		)
		{
		    vector<bool> visited(V, false);

		    for (int i = 0; i < V; i++) {
		        if (!visited[i] && is_cyclic(adj, i, V, visited)) {
		            return true;
		        }
		    }

		    return false;
		}

		////////////////////////////////////////////////////////////////////////
		/// makeCombiUtil
		////////////////////////////////////////////////////////////////////////
		void makeCombiUtil
		(
			vector<vector<int> >& ans,
		    vector<int>& tmp,
		    const int n,
		    const int left,
		    const int k
		)
		{
		    if (k == 0)
		    {
		        ans.push_back(tmp);
		        return;
		    }

		    for (int i = left; i <= n; ++i)
		    {
		        tmp.push_back(i);
		        makeCombiUtil(ans, tmp, n, i + 1, k - 1);
		        tmp.pop_back();
		    }
		}

		////////////////////////////////////////////////////////////////////////
		/// makeCombi
		////////////////////////////////////////////////////////////////////////
		vector<vector<int> > makeCombi
		(
			const int n,
			const int k
		)
		{
		    vector<vector<int> > ans;
		    vector<int> tmp;
		    makeCombiUtil(ans, tmp, n, 1, k);
		    return ans;
		}

		////////////////////////////////////////////////////////////////////////
		/// ContainsCycle
		////////////////////////////////////////////////////////////////////////
		bool ContainsCycle
		(
			const vector<SWCPoint>& vPoints
		)
		{
			const size_t nPts = vPoints.size();
			std::vector<std::vector<int> > adj;
      adj.resize(nPts);
			for (size_t i = 0; i < nPts; i++) {
				const vector<size_t>& conns = vPoints[i].conns;
				for (size_t j = 0; j < conns.size(); j++) {
					adj[i].push_back(conns[j]);
					adj[conns[j]].push_back(i);
				}
			}
			return !is_cyclic(adj, nPts);
		}

		////////////////////////////////////////////////////////////////////////
		/// HasUndesiredAngles
		////////////////////////////////////////////////////////////////////////
		bool HasUndesiredAngles(
			const std::vector<SWCPoint>& vPoints,
			const number angleThresholdMin,
			const number angleThresholdMax
		) {
			const size_t nPts = vPoints.size();
			for (size_t i = 0; i < nPts; i++) {
				const vector<size_t>& conns = vPoints[i].conns;
				for (size_t j = 1; j < conns.size(); j++) {
					vector3 v1, v2;
					VecSubtract(v1, vPoints[conns[0]].coords, vPoints[i].coords);
					VecSubtract(v2, vPoints[conns[j]].coords, vPoints[i].coords);
					const number angle = AngleBetweenDirections(v1, v2);
					if ( (angle < angleThresholdMin) || (angle > angleThresholdMax)) {
						UG_LOGN("First offending angle encountered at: "
								<< vPoints[conns[0]].coords << ", "
								<< vPoints[conns[j]].coords << ", "
								<< vPoints[i].coords << std::endl)
						return true;
					}
				}
			}
			return false;
		}

		////////////////////////////////////////////////////////////////////////
		/// CylinderCylinderSomaSeparationTest
		////////////////////////////////////////////////////////////////////////
		bool CylinderCylinderSomaSeparationTest
		(
			const vector<SWCPoint>& vSomaPoints
		) {
			for (size_t i = 0; i < vSomaPoints.size(); i++) {
				for (size_t j = 0; j < vSomaPoints.size(); j++) {
					if (i != j) {
						const number dist = VecDistance(vSomaPoints[i].coords,  vSomaPoints[j].coords);
						if (dist < vSomaPoints[i].radius || dist < vSomaPoints[j].radius) {
							return false;
						}
					}
				}
			}
			return true;
		}


		////////////////////////////////////////////////////////////////////////
		/// CylinderCylinderSeparationTest
		////////////////////////////////////////////////////////////////////////
		bool CylinderCylinderSeparationTest
		(
			const vector<SWCPoint>& vPoints
		)
		{
			typedef vector<vector<int> > combination_t;
			typedef vector<vector<int> >::const_iterator combination_iter;
			std::vector<Cylinder> cylinders;
			size_t nPts = vPoints.size();
			for (size_t i = 0; i < nPts; i++) {
				/// Regular points and root points only, since BPs may overlap usually
				if (vPoints[i].conns.size() == 2) {
					ug::vector3 dir;
					VecSubtract(dir, vPoints[vPoints[i].conns[0]].coords, vPoints[vPoints[i].conns[1]].coords);
					number len = VecLength(dir);
					VecNormalize(dir, dir);
					cylinders.push_back(Cylinder(vPoints[i].coords, dir, vPoints[i].radius, len));
					}
				}
			// create all pairwise cylinders, then test for intersection
			combination_t combinations = makeCombi(cylinders.size(), 2);
			combination_iter it = combinations.begin();
			while (it != combinations.end()) {
				if (!CylinderCylinderSeparationTest(cylinders[(*it)[0]], cylinders[(*it)[1]])) {
					return false;
				}
				++it;
			}
			return true;
		}

		////////////////////////////////////////////////////////////////////////
		/// CheckRootToSomaNeuriteDiameters
		////////////////////////////////////////////////////////////////////////
		void  CheckRootToSomaNeuriteDiameters
		(
			const std::vector<number>& radii,
			const number somaRadius
		)
		{
			for (size_t i = 0; i < radii.size()/12; i++) {
				UG_LOGN("Soma/Neurite ratio for neurite with index #" << i+1 << ": " << somaRadius/radii[i*12]);
				UG_COND_WARNING(somaRadius/ radii[i*12] > 10, "Soma/Neurite ratio highly "
						"anistropic for neurite with index " << i+1 << " thus one might "
						"expect non-optimal meshes at connecting region. Refine "
						"mesh (locally) adaptive around connecting region or choose"
						"less soma refinements. Consider enabling tapering toggle.");
			}
		}
		////////////////////////////////////////////////////////////////////////
		/// CalculateCenter
		////////////////////////////////////////////////////////////////////////
		void CalculateCenter
		(
			const std::vector<Vertex*>& vertices,
			const Grid::VertexAttachmentAccessor<APosition>& aaPos,
			const size_t size,
			ug::vector3& center
		)
		{
			VecSet(center, 0);
			UG_LOGN("0ed vector")
			for (size_t i = 0; i < size; i++) {
				UG_LOGN("i: " << i)
				VecAdd(center, center, aaPos[vertices[i]]);
			}
			UG_LOGN("scale it")
			VecScale(center, center, 1./(number)size);
		}

		////////////////////////////////////////////////////////////////////////
		/// GetNumberOfTriangleIntersections
		////////////////////////////////////////////////////////////////////////
		int GetNumberOfTriangleIntersections
		(
			const std::string& fileName,
			const number snapThreshold,
			const bool verbose
		) {
			Domain3d dom;
			try {LoadDomain(dom, fileName.c_str()); }
			UG_CATCH_THROW("Failed loading domain from '" << fileName << "'.");
			ug::MultiGrid* grid = dom.grid().get();
			Grid::VertexAttachmentAccessor<APosition> aaPos(*grid, aPosition);
			Selector sel(*grid);
			MultiGridSubsetHandler* sh = dom.subset_handler().get();
			Triangulate(*grid, grid->begin<ug::Quadrilateral>(), grid->end<ug::Quadrilateral>());
			RemoveDoubles<3>(*grid, grid->begin<Vertex>(), grid->end<Vertex>(), aPosition, SMALL);
			int numSubsets = sh->num_subsets();
			int numIntersections = 0;
			typedef lg_ntree<3, 3, Face> octree_t;
			typedef octree_t::box_t	box_t;
			size_t triCounter = 0;
			octree_t octree(*grid, aPosition);
			octree.create_tree(grid->begin<Triangle>(), grid->end<Triangle>());
			for(TriangleIterator triIter1 = grid->begin<Triangle>();
				triIter1 != grid->end<Triangle>(); ++triIter1, ++triCounter)
			{
				std::vector<Face*> closeTris;
				Face* t1 = *triIter1;

				box_t bbox;
				bbox.min = bbox.max = aaPos[t1->vertex(0)];
				for(size_t i = 1; i < t1->num_vertices(); ++i)
					bbox = box_t(bbox, aaPos[t1->vertex(i)]);
				bbox.min -= vector3(snapThreshold, snapThreshold, snapThreshold);
				bbox.max += vector3(snapThreshold, snapThreshold, snapThreshold);

				FindElementsInIntersectingNodes(closeTris, octree, bbox);

				for (size_t i_close = 0; i_close < closeTris.size(); ++i_close) {
					Face* t2 = closeTris[i_close];

					// avoid checking t1 with t1 for intersections
					if(grid->get_attachment_data_index(t1) >= grid->get_attachment_data_index(t2)) {
						continue;
					}


					// coplanar triangles
					vector3 n1, n2;
					CalculateNormal(n1, t1, aaPos);
					CalculateNormal(n2, t2, aaPos);
					number d = VecDot(n1, n2);
					if(fabs(d) > 1. - snapThreshold) {
						continue;
					}

					// potential intersection
					vector3 ip[2];
					if(TriangleTriangleIntersection(aaPos[t1->vertex(0)], aaPos[t1->vertex(1)],
						aaPos[t1->vertex(2)], aaPos[t2->vertex(0)], aaPos[t2->vertex(1)], aaPos[t2->vertex(2)],
						&ip[0], &ip[1], SMALL) == 1)
					{
						// shared edge is no intersection
						size_t numShared = NumSharedVertices(*grid, t1, t2);
						if(! (numShared >= 1)) {
							numIntersections++;
							// report some information if intersection occured
							if (verbose) {
								UG_LOGN("Intersection of triangle # " << grid->get_attachment_data_index(t1) <<
									"(" << aaPos[t1->vertex(0)] << ", " << aaPos[t1->vertex(1)] << ", " << aaPos[t1->vertex(2)] << ") " <<
									"with # " << grid->get_attachment_data_index(t2) <<
									"(" << aaPos[t2->vertex(0)] << ", " << aaPos[t2->vertex(1)] << ", " << aaPos[t2->vertex(2)] << ") " <<
									"at (intersection point): " << ip[0] << ", " << ip[1] << std::endl);
							}
							sel.select(t1); sel.select(t2);
							sel.select(t1->vertex(0)); sel.select(t1->vertex(1)); sel.select(t1->vertex(2));
							sel.select(t2->vertex(0)); sel.select(t2->vertex(1)); sel.select(t2->vertex(2));
						}
					}
				}
			}

			AssignSelectionToSubset(sel, *sh, numSubsets);
			sh->subset_info(numSubsets).name = "intersections";
			AssignSubsetColors(*sh);
			EraseEmptySubsets(*sh);
			std::string outFileName = fileName.substr(0, fileName.size()-4) + "_intersections.ugx";
			SaveGridToFile(*grid, *sh, outFileName.c_str());
			return sel.num<Triangle>();
		}
	}
}
