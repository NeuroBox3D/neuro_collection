/*!
 * \file neurite_util.h
 *
 * TODO: Cleanup comments, add unit tests and const correctness to current code base
 *
 *  Created on: Apr 22, 2019
 *      Author: Stephan Grein
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_UTIL_H
#include "test_neurite_proj.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/file_io/file_io_ugx.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/grid/geometry.h"
#include "lib_grid/global_attachments.h"
#include "lib_grid/algorithms/subset_color_util.h"
#include "lib_grid/algorithms/remeshing/resolve_intersections.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"
#include "common/math/ugmath_types.h"

namespace ug {
	namespace neuro_collection {
		/*!
		 * \brief generic comparator for SurfaceParams
		 */
		typedef float (NeuriteProjector::SurfaceParams::*membervar);
		template<membervar m> struct CompareBy {
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> > m_aaSurfParams;
			bool operator()(const ug::Vertex* const a, const ug::Vertex* const b) const {
				return m_aaSurfParams[a].*m < m_aaSurfParams[b].*m ;
			}
			CompareBy(const Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams) {
				m_aaSurfParams = aaSurfParams;
			}
		};

		/*!
		 * \brief comparator for elements in vector
		 */
		template <typename TElem>
		struct ExistsInVector
		{
			/*!
			 * \brief initializes vector
			 * \param[in] vec
			 */
			ExistsInVector(const std::vector<TElem>& vec) : m_vec(vec) {
			}

			/*!
			 * \brief operator()
			 * \param[in] elem
			 *
			 * \return \c true if elemenet found otherwise false
			 */
			bool operator() (TElem elem) {
				return (std::find(m_vec.begin(), m_vec.end(), elem) != m_vec.end());
			}
		private:
			const std::vector<TElem>& m_vec;
		};


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
			const std::vector<ug::vector3>& normals,
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
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] sh
		 * \param[in] rootNeurites
		 * \param[in] rootNeuritesInner
		 */
		void connect_outer_and_inner_root_neurites_to_outer_soma
		(
			size_t somaIndex,
			size_t numQuads,
			Grid& grid,
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
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] sh
		 */
		void connect_inner_neurites_to_inner_soma
		(
			size_t somaIndex, /// inner soma index: beginning of inner sphere's quads (ER) is somaIndex+1, outer sphere's quads (ER) is somaIndex-numQuads-1
			size_t numQuads, /// number of total surface quads or neurites to connect to
		    Grid& grid,
		    Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh
		);

		/*!
		 * \brief finds the quadrilateral vertices on soma surface
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] oldVertices
		 * \param[in] outRads
		 * \param[in] si
		 * \param[in] sh
		 * \param[in] rimSnapThresholdFactor
		 * \param[in] numQuads
		 * \param[in] numVerts
		 */
		std::vector<ug::vector3> find_quad_verts_on_soma
		(
		   Grid& grid,
		   Grid::VertexAttachmentAccessor<APosition>& aaPos,
		   std::vector<ug::Vertex*> verticesOld,
		   std::vector<std::vector<number> > outRads,
		   size_t si,
		   SubsetHandler& sh,
		   number rimSnapThresholdFactor,
		   size_t numQuads,
		   size_t numVerts=4
		);

		/*!
		 * \brief connects neurites to soma
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] aaSurfParams
		 * \param[out] outVerts
		 * \param[out] outVertsInner
		 * \param[out] outRads
		 * \param[out] smallerQuadVerts
		 * \param[in] si
		 * \param[in] sh
		 * \param[in] fileName
		 * \param[in] rimSnapThreshold
		 * \param[out] axisVectors
		 * \param[in] vNeurites
		 * \param[in] connectingVertices
		 * \param[in] connectingVerticesInner
		 * \param[in] connectingEdges
		 * \param[in] connectingEdgesInner
		 * \param[in] createInner
		 * \param[in] alpha
		 * \param[in] numIterations
		 * \param[in] resolveThreshold
		 * \param[in] scale
		 * \param[in] numVerts
		 * \param[in] numQuads
		 *
		 * TODO: Find suitable parameters for tangential smooth and resolve intersection
		 */
		void connect_neurites_with_soma
		(
		   Grid& grid,
		   Grid::VertexAttachmentAccessor<APosition>& aaPos,
		   Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		   std::vector<Vertex*> outVerts,
		   std::vector<Vertex*> outVertsInner,
		   std::vector<number> outRads,
		   std::vector<Vertex*>& smallerQuadVerts,
		   int si,
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
		   number scale=0.5,
		   size_t numVerts=4,
		   size_t numQuads=1
		);

		/*!
		 * \brief shrinks a quadrilateral and creates a copy of the smaller
		 * \param[in] vVrt
		 * \param[out] outvVrt
		 * \param[in] oldVertices
		 * \param[out] outvEdge
		 * \param[in,out] grid
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
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			bool createFaces=true,
			ISelector* outSel = NULL,
			ug::vector3* currentDir = NULL
		);

		/*!
		 * \brief shrinks quadrilateral and overwrites old quadrilateral's vertices
		 * \param[in,out] vVrt
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] percentage
		 */
		void shrink_quadrilateral
		(
			std::vector<Vertex*> vVrt,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage
		);

		/*!
		 * \brief creates the soma
		 * \param[in] somaPts
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] sh
		 * \param[in] si
		 * \param[in] numRefs
		 */
		void create_soma
		(
			const std::vector<SWCPoint>& somaPts,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			size_t si,
			size_t numRefs = 2
		);

		/*!
		 * \brief creates the soma
		 * \param[in] somaPts
		 * \param[in,out] grid
		 * \param[in] aaPos
		 */
		void create_soma
		(
			const std::vector<SWCPoint>& somaPts,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos
		);

		/*!
		 * \brief split a quadrilateral along its edges
		 * \param[in] vVrt
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] percentage
		 * \param[in] vecDir
		 * \param[out] vertices
		 * \param[out] edges
		 * \param[in] conservative
		 */
		void split_quadrilateral_along_edges
		(
			std::vector<Vertex*> vVrt,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			ug::vector3 vecDir,
			std::vector<ug::Vertex*>& vertices,
			std::vector<ug::Edge*>& edges,
			bool conservative = true
		);


		/*!
		 * \brief shrink a quadrilateral towards its center
		 * \param[in] vVrt
		 * \param[in,out] grid
		 * \param[in] aaPOs
		 * \param[in] percentage
		 * \param[in] center
		 */
		void shrink_quadrilateral_center
		(
			std::vector<Vertex*>& vVrt,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			ug::vector3& center
		);

		/*!
		 * \brief reorders the vertices accordingly
		 * \param[in] v
		 * \param[in] e
		 */
		void reorder_connecting_elements
		(
			std::vector<ug::Vertex*>& v,
			std::vector<ug::Edge*> e
		);

		/*!
		 * \brief Correcting inner branching points of neurites
		 * Note: In case of very small shrinkage factor might result in intersections
		 * \param[in] verts
		 * \param[in] edges
		 * \param[in] oldVertsSorted
		 * \param[in] aaSurfParams
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] scale
		 */
		void correct_edges
		(
			std::vector<ug::Vertex*>& verts,
			std::vector<ug::Edge*>& edges,
			std::vector<ug::Vertex*>& oldVertsSorted,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number scale
		);

		/*!
		 * \brief helper method to correct one side of the quadrilateral
		 * \param[in] verts
		 * \param[in] vertsOpp,
		 * \param[in] edges
		 * \param[in] edgesOpp
		 * \param[in] aaSurfParams
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] scale
		 */
		void correct_edges_all
		(
			std::vector<ug::Vertex*>& verts,
			std::vector<ug::Vertex*>& vertsOpp,
			std::vector<ug::Edge*>& edges,
			std::vector<ug::Edge*>& edgesOpp,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number scale
		);

		/*!
		 * \brief corrects the axial offset at the inner branching points
		 * This means, we move the points with smaller axial value further down
		 * the current neurite and the larger axial values further back
		 * \param[in] verts
		 * \param[in] aaSurfParams
		 * \param[in] aaPos
		 * \param[in] scale
		 */
		void correct_axial_offset
		(
			std::vector<ug::Vertex*>& verts,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number scale
		);

		/*!
		 * \brief add the new soma surface points to the precondioned swc file.
		 * Note that we assume every dendrite is connected to the ROOT soma point,
		 * and this might, depending on the reconstruction of the SWC file not be true.
		 * \param[in] lines
		 * \param[in] fn_precond
		 * \param[in] fn_precond_with_soma
		 * \param[in] vPointsSomaSurface
		 */
		void add_soma_surface_to_swc
		(
			const size_t& lines,
			const std::string& fn_precond,
			const std::string& fn_precond_with_soma,
			const std::vector<ug::vector3>& vPointsSomaSurface
		);

		/*!
		 * \brief get closest points to soma
		 * \param[in] fn_precond
		 * \param[in] vPos
		 * \param[in] lines
		 */
		void get_closest_points_to_soma
		(
			const std::string& fn_precond,
			std::vector<ug::vector3>& vPos,
			size_t& lines
		);

		/*!
		 * \brief get closest vertices on soma
		 * \param[in] vPos
		 * \param[in] vPointsSomaSurface
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] sh
		 * \param[in] si
		 */
		void get_closest_vertices_on_soma
		(
			const std::vector<ug::vector3>& vPos,
			std::vector<ug::Vertex*>& vPointsSomaSurface,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			size_t si
		);

		/*!
		 * \brief get closest points on soma
		 * \param[in] vPos
		 * \param[in] vPointsSomaSurface
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] sh
		 * \param[in] si
		 */
		void get_closest_points_on_soma
		(
			const std::vector<ug::vector3>& vPos,
			std::vector<ug::vector3>& vPointsSomaSurface,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			size_t si
		);

		/*!
		 * \brief replace first root neurite vertex in swc
		 * \param[in] lines
		 * \param[in] fn_precond
		 * \param[in] fn_precond_with_soma
		 * \param[in] vPointsSomaSurface
		 */
		void replace_first_root_neurite_vertex_in_swc
		(
			const size_t& lines,
			const std::string& fn_precond,
			const std::string& fn_precond_with_soma,
			const std::vector<ug::vector3>& vPointsSomaSurface
		);

		/*!
		 * \brief shrinks a polygon towards the barycenter
		 * \param[in] vVrt
		 * \param[out] outvVrt
		 * \param[in] oldVertices
		 * \param[out] outvEdge
		 * \param[in,out] grid
		 * \param[in] aaPos
		 * \param[in] percentage
		 * \param[in] createFacs
		 * \param[in] outSel

		 */
		void shrink_polygon_copy
		(
			const std::vector<Vertex*>& vVrt,
			std::vector<Vertex*>& outvVrt,
			const std::vector<Vertex*>& oldVertices,
			std::vector<Edge*>& outvEdge,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			bool createFaces,
			ISelector* outSel
		);

		/*!
		 * \brief connect outer and inner root neurites to outer soma variant
		 * \param[in] somaIndex
		 * \param[in] numQuads
		 * \param[in,out] grid
		 * \param[in] aaPos
		 */
		void connect_outer_and_inner_root_neurites_to_outer_soma_variant
		(
			size_t somaIndex,
			size_t numQuads,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			std::vector<ug::Vertex*>& rootNeurites,
			size_t numVerts
		);

		/*!
		 * \brief fill volume with tetrahedra
		 * \param[in, out] grid
		 * \param[in] quality
		 * \param[in] preserveBnds
		 * \param[in] preserveAll
		 * \param[in] aaPos
		 * \param[in] verbosity
		 */
		void tetrahedralize_soma
		(
			Grid& grid,
			number quality,
			bool preserveBnds,
			bool preserveAll,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			int verbosity=1
		);
	}
}


#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_UTIL_H
