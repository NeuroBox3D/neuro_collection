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
#define UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_MATH_UTIL_H

#include "common/math/ugmath_types.h"
#include "types.h"
#include <lib_grid/grid/grid.h>
#include <lib_grid/common_attachments.h>

#include <lib_grid/algorithms/geom_obj_util/misc_util.h>
#include "lib_disc/domain_util.h"
#include "test_neurite_proj.h"

namespace ug {
	namespace neuro_collection {
		/*!
		 * \brief signum function {-1,0,1}
		 * \tparam[T] type
		 * \param[in] val - value
		 * \return \c integer \f$ i \in {-1,0,1} \f$
		 */
		template <typename T> int sgn(T val) {
			return (T(0) < val) - (val < T(0));
		}

		/*!
		 * \brief Rotate a vector (OP) around axis (OQ) using Rodrigues' rotation
		 * \f$ \vec{v}_{rot} =
		 * 		\vec{v} cos(\theta) + (\vec{v} \times \vec{k}) sin(\theta) +
		 * 		\vec{k} (\vec{k} \cdot \vec{v}) cos(\theta)\$
		 *
		 * \param[in] p end point of vector
		 * \param[in] q end point of axis
		 * \param[in] o origin (from where p and q are emanating)
		 * \param[in] theta degree in radians
		 * \param[out] v_rot the rotated vector
		 */
		void RotateVectorAroundAxis
		(
			const vector3& p,
			const vector3& q,
			const vector3& o,
			number theta,
			vector3& v_rot
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

		/*!
		 * \brief Determines if two cylinders a and b are separated
		 * \param[in] a - first cylinder
		 * \param[in] b - second cylinder
		 *
		 * \return \c true if separated otherwise false
		 */
		struct Cylinder;
		bool CylinderCylinderSeparationTest
		(
			const Cylinder& a,
			const Cylinder& b
		);

		/*!
		 * \brief A simple cylinder-cylinder separation test for soma sphere
		 *
		 * The pairwise neurite start points on the soma surface are checked
		 * for an overlap of their associated cylinders (resepctively radii)
		 * To this end the method calculates the geodesic between any two points
		 * on the soma surface (which is assumed to be a perfect sphere). The
		 * the distance is compared with the sum of the radii under an epsilon
		 * threshold criterion. The method returns true if at least one of
		 * the cylinders intersect - otherwise the method returns false.
		 *
		 * \param[in] vSomaPoints - neurite start points on soma surface
		 * \param[in] somaPoint - center of soma
		 * \param[in] eps - epsilon
		 *
		 * \return \c bool
		 */
		bool CylinderCylinderSomaSeparationTest
		(
			const std::vector<SWCPoint>& vSomaPoints,
			const SWCPoint& somaPoint,
			number eps
		);


		struct Cylinder {
			vector3 c; //<! center (point)
			vector3 w; //<! unit-length axis (vector)
			number r; //<! radius
			number h; //<! height
			Cylinder(const vector3& c, const vector3& w, number r, number h)
			: c(c),
			  w(w),
			  r(r),
			  h(h)
			{ }
		};

		/*!
		 * \brief create a cylinder from SWC point
		 * \param[in] p SWC point
		 * \return Cylinder
		 */
		Cylinder makeCyl
		(
			const SWCPoint& p
		);

		/*!
		 * \brief creates combinations by a DFS algorithm
		 * \param[in] n
		 * \param[in] k
		 * \return \c (n over k) combinations
		 */
		std::vector<std::vector<int> > makeCombi
		(
			int n,
			int k
		);

		/*!
		 * \brief utility function which actually creates the combinations
		 * \param[in] ans
		 * \param[in] tmp
		 * \param[in] n
		 * \param[in] left
		 * \param[in] k
		 */
		void makeCombiUtil(
			std::vector<std::vector<int> >& ans,
		    std::vector<int>& tmp,
		    int n,
		    int left,
		    int k
		);

		/*!
		 * \brief Check for a cycle in a graph given by an adjacency list
		 * \param[in] adj adjacency list
		 * \param[in] V index of vertex to start BFS
		 * The function returns true if at least one cycle found otherwise false
		 * Runtime complexity of this BFS algorithm: O(|V|+|E|)
		 * \return \c \bool
		 */
		/*
		bool is_cyclic
		(
			const std::vector<int> adj[],
			int V
		);
		*/

		/*!
		 * \brief Check for a cycle in the SWC graph
		 * The function returns true if at least one cycle found otherwise false
		 * \param[in] vPoints list of SWC points
		 * \return \c bool
		 */
		bool ContainsCycle
		(
			const std::vector<SWCPoint>& vPoints
		);

		/*!
		 * \brief Check for undesired angles
		 * \param[in] vPoints list of SWC points
		 * \param[in] angleThresholdMin
		 * \param[in] angleThresholdMax
		 */
		bool HasUndesiredAngles(
			const std::vector<SWCPoint>& vPoints,
			number angleThresholdMin,
			number angleThresholdMax
		);

		/*!
		 * \brief Check for any intersection of two SWC cylinders
		 * The function returns true if at least one intersection found else false
		 * \param[in] vPoints list of SWC points
		 * \return \c bool
		 */
		bool CylinderCylinderSeparationTest
		(
			const std::vector<SWCPoint>& vPoints
		);

		/*!
		 * \brief Check ratio of soma and neurites
		 * \param[in] radii diameter of root neurite vertices
		 * \param[in] somaRadius soma's radius
		 */
		void CheckRootToSomaNeuriteDiameters
		(
			const std::vector<number>& radii,
			number somaRadius
		);

		/*!
		 * \brief CalculateCenter
		 * \param[in] vertices
		 * \param[in] aaPos
		 * \param[in] size
		 * \param[out] center
		 */
		void CalculateCenter
		(
			const std::vector<Vertex*>& vertices,
			const Grid::VertexAttachmentAccessor<APosition>& aaPos,
			size_t size,
			ug::vector3& center
		);

		/*
		 * \brief Report the number of triangle intersections in grid
		 * \param[in] fileName name of the grid
		 * \param[in] snapThreshold
		 * \param[in] verbose report triangle indices, coordinates and intersection points
		 */
		int GetNumberOfTriangleIntersections
		(
			const std::string& fileName,
			number snapThreshold,
			bool verbose
		);

		/*!
		 * \brief Calculate the path length between two subsets
		 */
		template <typename TDomain>
		number PathLength1D(
			const std::string& fileName,
			const std::string& fromSubset,
			const std::string& toSubset,
			TDomain& dom
		) {
			/// 1D grid
			Domain3d dom1d;
			try { LoadDomain(dom1d, fileName.c_str()); }
			UG_CATCH_THROW("Failed loading domain from '" << fileName << "'.");

			/// 3D grid
			Grid& grid = *dom.grid();
			const MGSubsetHandler& sh = *dom.subset_handler();
			grid.attach_to_vertices(aPosition);
			Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

			/// Check for subset presence in 3D grid, then calculate center in 3D grids
			const int siFrom = sh.get_subset_index(fromSubset.c_str());
			UG_COND_THROW( (siFrom == -1), "Measuring subset index (#" << siFrom << ") not found for "
					<< "measuring subset name (" << fromSubset << ") in domain");
			const vector3 centerFrom = CalculateCenter(sh.begin<Edge>(siFrom, 0),
					sh.end<Edge>(siFrom, 0), aaPos);

			const int siTo = sh.get_subset_index(toSubset.c_str());
			UG_COND_THROW( (siTo == -1), "Measuring subset index (#" << siTo << ") not found for "
					<< "measuring subset name (" << toSubset << ") in domain");
			const vector3 centerTo = CalculateCenter(sh.begin<Edge>(siTo, 0),
					sh.end<Edge>(siTo, 0), aaPos);

			/// Find the corresponding 1D vertices to the two 3D measurement subsets
			typename TDomain::position_type v;
			typename TDomain::grid_type& g = *dom1d.grid();
			/// Start of 1D path
			ug::Vertex* v1 = FindClosestByCoordinate<Vertex>(centerFrom, g.template begin<Vertex>(),
					g.template end<Vertex>(), dom1d.position_accessor());
			UG_COND_THROW(!v1, "Start vertex of 1D path not found!")

			/// End of 1D path
			ug::Vertex* v2 = FindClosestByCoordinate<Vertex>(centerTo, g.template begin<Vertex>(),
					g.template end<Vertex>(), dom1d.position_accessor());
			UG_COND_THROW(!v2, "End vertex of 1D path not found!")

			UG_COND_WARNING(v1 == v2, "Start and end vertex of 1D path are identical. Intentional?")

			/// Populate graph
			AInt aInt;
			g.attach_to_vertices(aInt);
			g.attach_to_vertices(aPosition);
			Grid::VertexAttachmentAccessor<AInt> aaInt(g, aInt);

			VertexIterator vit = g.template begin<Vertex>();
			size_t size = 0;
			for (; vit != g.template end<Vertex>(); ++vit) {
				aaInt[*vit] = ++size;
			}

			Graph graph(size+1);

			ConstEdgeIterator eit = g.template begin<Edge>();
			for (; eit != g.template end<Edge>(); ++eit) {
				ug::Edge* edge = *eit;
				number dist = VecDistance(dom1d.position_accessor()[edge->vertex(0)],
						dom1d.position_accessor()[edge->vertex(1)]);
				graph.addEdge(aaInt[edge->vertex(0)], aaInt[edge->vertex(1)], dist);
			}

			return graph.shortestPath(aaInt[v1])[aaInt[v2]];
		}
	}
}

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_MATH_UTIL_H
