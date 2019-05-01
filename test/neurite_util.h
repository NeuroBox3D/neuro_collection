/*!
 * \file neurite_util.h
 *
 * TODO: Cleanup and commenting of code - throw old legacy code before volumes
 * TODO: Move all useful code out of testNeuriteProjector.cpp to the util class
 * TODO: Add unit tests for the current code base
 *
 *  Created on: Apr 22, 2019
 *      Author: Stephan Grein
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_UTIL_H
#include "lib_grid/grid/grid.h"
#include "lib_grid/file_io/file_io_ugx.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/grid/geometry.h"
#include "lib_grid/global_attachments.h"
#include "common/math/ugmath_types.h"
#include "lib_grid/algorithms/subset_color_util.h"
#include "lib_grid/algorithms/remeshing/resolve_intersections.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"
#include "test_neurite_proj.h"

namespace ug {
	namespace neuro_collection {
		/*!
		 * \brief calculates the angle enclosed by the points point and origin
		 * The calculation is relative to pos and the points are supposed to be
		 * contained in the plane denoted by pos and the the orthogonal vector n
		 * \param[in] pos
		 * \param[in] origin
		 * \param[in] point
		 * \param[in] n
		 * \return \c angle in interval [-pi, pi]
		 *
		 */
		float calculate_angle
		(
			const vector3& pos,
			const vector3& origin,
			const vector3& point,
			const vector3& n
		);

		/*!
		 * \brief calculates the angle enclosed by the points point and origin
		 * The calculation is relative to pos and the points don't lie in a plane
		 * \param[in] pos
		 * \param[in] origin
		 * \param[in] point
		 * \return \c angle in interval [0, 2*pi]
		 */
		float calculate_angle
		(
			const vector3& pos,
			const vector3& origin,
			const vector3& point
		);

		/*!
		 * \brief converts angles in interval [-pi, pi] to interval [0, 2*pi]
		 * fmod might not be have as expected in C++ (sign of reminder is kept)
		 * \param[in] angle
		 * \return \c angle in interval [0, 2*pi]
		 */
		float deg_to_full_range(float angle);

		/**
		 * \brief calculates the angles according to a center vertex with an implicitly specified reference vertex or specified by index
		 * The calculation is done for all points in the vector points, if no
		 * reference point is chosen, then the first point in the vector is used.
		 * \param[in] originCenter
		 * \param[in] points
		 * \param[out] angles
		 * \param[in] refIndex
		 */
		void calculate_angles
		(
			const vector3& originCenter,
			const std::vector<ug::vector3>& points,
			std::vector<number>& angles,
			size_t refIndex=-1
		);

		/**
		 * \brief calculates the angles according to a center vertex with required reference vertex
		 * The calculation is done for all points in the vector points
		 * \param[in] originCenter
		 * \param[in] points
		 * \param[out] angles
		 * \param[in] normals
		 * \param[in] ref
		 */
		void calculate_angles
		(
			const vector3& originCenter,
			const std::vector<ug::vector3>& points,
			std::vector<number>& angles,
			const std::vector<ug::vector3>& normals,
			const ug::vector3& ref
		);

		/**
		 * \brief calculates the angles according to a center vertex with implicitly specified reference vertex or specified by index
		 * The calculation is done for all points in the vector points
		 * \param[in] originCenter
		 * \param[in] points
		 * \param[out] angles
		 * \param[in] normals
		 * \param[in] ref
		 */
		void calculate_angles
		(
			const vector3& originCenter,
			const std::vector<ug::vector3>& points,
			std::vector<number>& angles,
			std::vector<ug::vector3>& normals,
			size_t refIndex=-1
		);

		/*!
		 * \brief sorts a vector's values and stores the original indices of unsorted vector
		 * Note that the original vector is not changed
		 * \param[in] values
		 * \param[in] indices
		 */
		void sortIndices
		(
			const std::vector<number>& values,
			std::vector<size_t>& indices
		);

		/*!
		 * \brief "comparable" to sort a vector of pair<ug::Vertex*, number> by the number
		 * \param[in] a
		 * \param[in] b
		 * \return \c bool
		 */
		bool sortbysec
		(
			const std::pair<ug::Vertex*, number> &a,
	        const std::pair<ug::Vertex*, number> &b
	    );

		/*!
		 * \brief convert a angle in range [-pi,pi] to [0, 2*pi] reliable
		 * \param[in] a
		 * \return \c angle in range [0, 2*pi]
		 * Note that the fmod C++ function might not work as expected
		 */
		number deg360
		(
			number a
		);


		/*!
		 * \brief this connects the outer soma surface with the root neurites vertices ER (smaller quad) and PM (larger quad)
		 * First the root neurite vertices are projected to the outer sphere / soma surface. Then the smallest angle pairs
		 * are calculated and then the associated unprojected vertices are connected with the small and large quad soma sphere vertices
		 * 1. Find the closest vertex from each of the inner quad on the soma surface to the outer quad of the soma surface
	     * 2. The corresponding inner root neurite vertex (indices will be the same as for the outer root neurite vertices)
	     * 3. Connect the corresponding inner quad vertex with the inner root neurite
	     * Could also call this method two times but the reference point of angle calculation has to be the same
	     * Thus we find for each outer quad soma vertices the closest inner quad soma vertices - These vertices correspond to
	     * the rootNeuritesInner list and have the same index as the vertices used when connecting the outer soma quad vertices
	     * with rootNeurites list. Thus remember for each pair of closest soma inner quad vertex and corresponding outer quad vertex
	     * the index with which the outer quad vertex was connected with the outer rootNeurite vertices. This index will be used to
	     * connect the closest soma inner vertex with the corresponding vertex in rootNeuritesInner vertices with the same index.
		 * \param[in] somaIndex
		 * \param[in] numQuads
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] sh
		 * \param[in] rootNeurites
		 * \param[in] rootNeuritesInner
		 */
		void connect_outer_and_inner_root_neurites_to_outer_soma
		(
			size_t somaIndex,
			size_t numQuads,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			std::vector<ug::Vertex*>& rootNeurites,
			std::vector<ug::Vertex*>& rootNeuritesInner
		);



		/*!
		 * \brief connects the inner neurites (ER) to the inner sphere's (ER) surface quads
		 * Projected vertices onto the inner sphere's quad plane will correspond to an unprojected outer sphere's quad vertex.
		 * The projected vertices are used to find the vertices with smallest angle difference between the inner sphere's quad vertices and the projected vertices onto that plane
		 * Then the unprojected vertices and the corresponding vertex on the inner sphere's quad are connected by edges and faces.
		 * This procedure is also used in connect_outer_and_inner_root_neurites_to_outer_soma which uses angle differences -
		 * this method connect_inner_neurites_to_inner_soma used the smallest difference of positions which should suffice
		 * \param[in] somaIndex
		 * \param[in] numQuads
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] sh
		 */
		void connect_inner_neurites_to_inner_soma
		(
				size_t somaIndex, /// inner soma index: beginning of inner sphere's quads (ER) is somaIndex+1, outer sphere's quads (ER) is somaIndex-numQuads-1
				size_t numQuads, /// number of total surface quads or neurites to connect to
			    Grid& g,
			    Grid::VertexAttachmentAccessor<APosition>& aaPos,
			    SubsetHandler& sh
		);

		/*!
		 * \brief finds the quadrilateral vertices on soma surface
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] oldVertices
		 * \param[in] outRads
		 * \param[in] si
		 * \param[in] sh
		 * \param[in] rimSnaptThresholdFactor
		 * \param[in] numQuads
		 */
		std::vector<ug::vector3> find_quad_verts_on_soma
		(
		   Grid& g,
		   Grid::VertexAttachmentAccessor<APosition>& aaPos,
		   std::vector<ug::Vertex*> verticesOld,
		   std::vector<std::vector<number> > outRads,
		   size_t si,
		   SubsetHandler& sh,
		   number rimSnapThresholdFactor,
		   size_t numQuads
		);

		/*!
		 * \brief connects neurites to soma
		 * \param[in] g
		 * \param[in] aaPos
		 * TODO: Find suitable parameters for tangential smooth and resolve intersection
		 */
		void connect_neurites_with_soma
		(
		   Grid& g,
		   Grid::VertexAttachmentAccessor<APosition>& aaPos,
		   Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		   std::vector<Vertex*> outVerts,
		   std::vector<Vertex*> outVertsInner,
		   std::vector<number> outRads,
		   std::vector<Vertex*>& smallerQuadVerts,
		   size_t si,
		   SubsetHandler& sh,
		   const std::string& fileName,
		   number rimSnapThresholdFactor,
		   std::vector<std::pair<size_t, std::pair<ug::vector3, ug::vector3> > >& axisVectors,
		   std::vector<NeuriteProjector::Neurite>& vNeurites,
		   std::vector<std::vector<ug::Vertex*> >& connectingVertices,
		   std::vector<std::vector<ug::Vertex*> >& connectingVerticesInner,
		   std::vector<std::vector<ug::Edge*> >& connectingEdges,
	       std::vector<std::vector<ug::Edge*> >& connectingEdgesInner,
		   bool createInner=true,
		   number alpha=0.01,
		   int numIterations=10,
		   number resolveThreshold=0.00001,
		   number scale=0.5
		);

		/*!
		 * \brief shrinks a quadrilateral and creates a copy of the smaller
		 * \param[in] vVrt
		 * \param[out] outvVrt
		 * \param[in] oldVertices
		 * \param[out] outvEdge
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] percentage
		 * \param[in] createFacs
		 * \param[in] outSel
		 * \param[in] currentDir
		 */
		void shrink_quadrilateral_copy
		(
			const std::vector<Vertex*>& vVrt,
			std::vector<Vertex*>& outvVrt,
			const std::vector<Vertex*>& oldVertices,
			std::vector<Edge*>& outvEdge,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			bool createFaces=true,
			ISelector* outSel = NULL,
			ug::vector3* currentDir = NULL
		);

		/*!
		 * \brief shrinks quadrilateral and overwrites old quadrilateral's vertices
		 * \param[in] vVrt
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] percentage
		 */
		void shrink_quadrilateral
		(
			std::vector<Vertex*> vVrt,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage
		);

		/*!
		 * \brief creates the soma
		 * \param[in] somaPts
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] sh
		 * \param[in] si
		 * \param[in] numRefs
		 */
		void create_soma
		(
				const std::vector<SWCPoint>& somaPts,
				Grid& g,
				Grid::VertexAttachmentAccessor<APosition>& aaPos,
				SubsetHandler& sh,
				size_t si,
				size_t numRefs = 2
		);

		/*!
		 * \brief creates the soma
		 * \param[in] somaPts
		 * \param[in] g
		 * \param[in] aaPos
		 */
		void create_soma
		(
				const std::vector<SWCPoint>& somaPts,
				Grid& g,
				Grid::VertexAttachmentAccessor<APosition>& aaPos
		);
	}
}


#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_UTIL_H
