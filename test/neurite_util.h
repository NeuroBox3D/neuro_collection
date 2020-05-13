/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-04-22
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
		 * \brief vector exists comparator
		 * Checks if an element is available in the given vector
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
		float deg_to_full_range
		(
			float angle
		);

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
	     * !!! WARNING !!!: Works only if inner and outer root neurite have the same number of connecting vertices
		 * \param[in] somaIndex
		 * \param[in] numQuads
		 * \param[in,out] grid
		 * \param[in, out] aaPos
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
		 * !!! WARNING !!!: Might not work if using distance-based criterion - should use angle based criterion after projection to plane
		 * \param[in] somaIndex
		 * \param[in] numQuads
		 * \param[in,out] grid
		 * \param[in, out] aaPos
		 * \param[in] sh
		 * \param[in, out] aaSurfParams
		 * \param[in] scale
		 */
		void connect_inner_neurites_to_inner_soma
		(
			size_t somaIndex, /// inner soma index: beginning of inner sphere's quads (ER) is somaIndex+1, outer sphere's quads (ER) is somaIndex-numQuads-1
			size_t numQuads, /// number of total surface quads or neurites to connect to
		    Grid& grid,
		    Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
		    Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		    number scale
		);

		/*!
		 * \brief finds the quadrilateral vertices on a soma surface
		 * \param[in,out] grid
		 * \param[in, out] aaPos
		 * \param[in] oldVertices
		 * \param[in] outRads
		 * \param[in] si
		 * \param[in] sh
		 * \param[in] rimSnapThresholdFactor
		 * \param[in] numQuads
		 * \param[in] numVerts
		 */
		std::vector<ug::vector3> FindSomaSurfaceCenters
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
		 * \param[in, out] aaPos
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
		 * FIXME: Find suitable parameters for tangential smooth and resolve intersection
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
		   size_t numQuads,
		   bool createInner=true,
		   number alpha=0.01,
		   int numIterations=10,
		   number resolveThreshold=0.00001,
		   number scale=0.5,
		   size_t numVerts=4
		);

		/*!
		 * \brief shrinks a quadrilateral and creates a copy of the smaller one
		 * \param[in] vVrt
		 * \param[out] outvVrt
		 * \param[in] oldVertices
		 * \param[out] outvEdge
		 * \param[in,out] grid
		 * \param[in, out] aaPos
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
		 * \param[in, out] aaPos
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
		 * \brief creates the soma as icosphere
		 * \param[in] somaPts
		 * \param[in,out] grid
		 * \param[in, out] aaPos
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
		 * \brief sets the somata's axial parameters (ER and PM surfaces)
		 * \param[in,out] grid
		 * \param[in] sh
		 * \param[in,out] aaSurfParams
		 * \param[in] somaIndex
		 * \param[in] erIndex
		 */
		void set_somata_axial_parameters
		(
			Grid& grid,
			SubsetHandler& sh,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			size_t somaIndex,
			size_t erIndex
		);

		/*!
		 * \brief fix somata axial parameters
		 * \param[in,out] grid
		 * \param[in,out] sh
		 * \param[in,out] aaSurfParams
		 * \param[in] aaPos
		 * \param[in] somaIndex
		 * \param[in] erIndex
		 * \param[in] somaPoint
		 * \param[in] scaleER
		 */
		void fix_axial_parameters
		(
			Grid& g,
			SubsetHandler& sh,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		    const Grid::VertexAttachmentAccessor<APosition>& aaPos,
			size_t somaIndex,
			size_t erIndex,
            const SWCPoint& somaPoint,
            number scaleER
		);

		/*!
		 * \brief reassign soma volumes specified by soma indices
		 * \param[in,out] grid
		 * \param[in,out] sh
		 * \param[in] somaIndexInner
		 * \param[in] somaIndexOuter
		 * \param[in] scale
		 * \param[in] savedSomaPoint
		 * \param[in,out] aaPos
		 */
		void reassign_volumes
		(
			Grid& grid,
			SubsetHandler& sh,
			size_t somaIndexInner,
			size_t somaIndexOuter,
			number scale,
			const SWCPoint& point,
			Grid::VertexAttachmentAccessor<APosition>& aaPos
		);

		/*!
		 * \brief creates the soma as icosahedron
		 * \param[in] somaPts
		 * \param[in,out] grid
		 * \param[in, out] aaPos
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
		 * \param[in, out] aaPos
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
		 * \param[in, out] aaPos
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
		 * \param[in, out] aaPos
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
		 * \param[in, out] aaPos
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
		 * \param[in, out] aaPos
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
		 * \param[in, out] aaPos
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
		void ReplaceFirstRootNeuriteVertexInSWC
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
		 * \param[in, out] aaPos
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
		 * \param[in, out] aaPos
		 * \param[in] merge
		 */
		void connect_outer_and_inner_root_neurites_to_outer_soma_variant
		(
			size_t somaIndex,
			size_t numQuads,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			std::vector<ug::Vertex*>& rootNeurites,
			size_t numVerts,
			bool merge
		);

		/*!
		 * \brief fill volume with tetrahedra
		 * \param[in, out] grid
		 * \param[in] sh
		 * \param[in, out] aaPos
		 * \param[in] somaIndex
		 * \param[in] erIndex
		 * \param[in] scale
		 * \param[in] somaPoint
		 */
		void tetrahedralize_soma
		(
			Grid& grid,
			SubsetHandler& sh,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			size_t somaIndex,
			size_t erIndex,
			const std::vector<SWCPoint>& somaPoint,
			number scale = -1.0
		);

		/*!
		 * \brief create pyramid
		 * \param[in,out] grid
		 * \param[in] quad
		 * \param[in,out] aaPos
		 * \param[in] scale
		 * \param[in,out] aaSurfParams
		 * \return \c pointer to new pyramid
		 */
		Pyramid* create_pyramid
		(
			Grid& grid,
			const Quadrilateral* const quad,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number scale = 1.0,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >* aaSurfParams = NULL
		);

		/*!
		 * \brief saves a grid to file and prepares the grid before
		 * Assigns subsets colors and erases empty subsets then saves the grid
		 * \param[in, out] grid
		 * \param[in, out] sh
		 * \param[in] fileName
		 */
		void SavePreparedGridToFile
		(
			Grid& grid,
			ISubsetHandler& sh,
			const char* const fileName
		);

		/*!
		 * \brief finds the quadrilaterals constrained to have a certain axial and radial parameters
		 * \param[in,out] grid
		 * \param[in,out] aaSurfParams
		 * \param[out] quadCont
		 * \param[in] axial
		 * \param[in] scale
		 * \param[in] numVertices
		 */
		void find_quadrilaterals_constrained
		(
			Grid& grid,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid::traits<Quadrilateral>::secure_container& quadCont,
			number axial = 0.0,
			number scale = 1.0,
			size_t numVertices = 2
		);

		/*!
		 * \brief extend ER part of neurites into soma in normal direction
		 * \param[in,out] grid
		 * \param[in,out] sh
		 * \param[in,out] aaPos
		 * \param[in,out] aaSurfParams
		 * \param[in] somaIndex
		 * \param[in] numQuads
		 * \param[in] scale
		 * \param[out] outVertsInner
		 */
		 void extend_ER_within
		 (
			Grid& grid,
			SubsetHandler& sh,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::Mapping> >& aaMapping,
			int somaIndex,
			size_t numQuads,
			number scale,
			std::vector<ug::Vertex*>& outVertsInner,
			const SWCPoint& somaPoint
		);

		/*!
		 * \brief splits a grid into two subgrids based on the provided subset indices
		 * \param[in] gridIn
		 * \param[in] srcSh
		 * \param[out] gridOut
		 * \param[out] destSh
		 * \param[in,out] aaPos
		 * \param[in] vSi
		 */
		void split_grid_based_on_subset_indices
		(
			Grid& gridIn,
			ISubsetHandler& srcSh,
			Grid& gridOut,
			ISubsetHandler& destSh,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			const std::vector<size_t>& vSi
		);

		/*!
		 * \brief splits a grid into two subgrids based on the provided selection
		 * \param[in] gridIn
		 * \param[in] srcSh
		 * \param[out] gridOut
		 * \param[out] destSh
		 * \param[in,out] aaPos
		 * \param[in] somaPoint
		 */
		void split_grid_based_on_selection
		(
			Grid& gridIn,
			ISubsetHandler& srcSh,
			Grid& gridOut,
			ISubsetHandler& destSh,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			const std::vector<SWCPoint>& somaPoint
		);

		/*!
		 * \brief selects elements whose center lies within or on a sphere specified by center and radius
		 * \param[in] grid
		 * \param[out] sel
		 * \param[in] center
		 * \param[in] radius
		 * \param[in,out] aaPos
		 */
		template <class TElem>
		void SelectElementsInSphere
		(
			Grid& grid,
			Selector& sel,
			const ug::vector3& center,
			const number radius,
			Grid::VertexAttachmentAccessor<APosition>& aaPos
		)
		{
			for(typename Grid::traits<TElem>::iterator iter = grid.begin<TElem>();
				iter != grid.end<TElem>(); ++iter)
			{
				vector3 c = CalculateCenter(*iter, aaPos);
				vector3 diff;
				VecSubtract(diff, c, center);
				VecPow(diff, diff, 2.0);
				number s = diff.x() + diff.y() + diff.z();
				if(s <= (sq(radius) + SMALL)) {
					sel.select(*iter);
				}
			}
		}



		/*!
		 * \brief selects elements whose axial surface parameter is smaller than given by axial parameter
		 * \param[in] grid
		 * \param[out] sel
		 * \param[in] axial
		 * \param[in,out] aaPos
		 * \param[in] aaSurfParams
		 */
		template <class TElem>
		void SelectElementsByAxialPosition
		(
			Grid& grid,
			Selector& sel,
			const number axial,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams
		);

		/*!
		 * \brief selects elements whose axial surface parameter is smaller than given by axial parameter
		 * \param[in] grid
		 * \param[out] sel
		 * \param[in] axial
		 * \param[in,out] aaPos
		 * \param[in] aaSurfParams
		 * \param[in] sh
		 * \param[in] si
		 * \param[in] scale
		 */
		template <class TElem>
		void SelectElementsByAxialPositionInSubset
		(
			Grid& grid,
			Selector& sel,
			const number axial,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			SubsetHandler& sh,
			const int si,
			const number scale
		);

		/*!
		 * \brief Deletes inner edges from quadrilaterals
		 * \param[in,out] grid
		 * \param[in,out] sh
		 * \param[in] si
		 */
		void DeleteInnerEdgesFromQuadrilaterals
		(
			Grid& grid,
			SubsetHandler& sh,
			int si = 0
		);

		/*!
		 * \brief Rotates vector around axis TODO: Old and deprecated
		 * The rotated vectors (CW and CCW) are stored as xPrime and xPrime2
		 * \param[in] vector the vector to be rotated around an axis
		 * \param[in] axis the axis a vector will be rotated around
		 * \param[in] origin of vector as well as axis
		 * \param[out] xPrime rotated vector by theta CCW
		 * \param[out] xPrime2 rotated vector by theta CW
		 * \param[in] theta angle in degree
		 */
		void rotate_vector_around_axis
		(
			const vector3& vector,
			const vector3& axis,
			const vector3& origin,
			vector3& xPrime,
			vector3& xPrime2,
			number theta
		);

		/*!
		 * \brief Adapts the surface grid to square
		 * \param[in] i
		 * TODO: Documentation
		 */
		void adapt_surface_grid_to_square
		(
			size_t i
		);

		/*!
		 * \brief new strategy: TODO refactor
		 * TODO: DOcumentation
		 */
		void connect_neurites_with_soma_var
		(
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			std::vector<Vertex*> outVerts,
			std::vector<number> outRads,
			int si,
			SubsetHandler& sh,
			const std::string& fileName,
			number rimSnapThresholdFactor,
			std::vector<NeuriteProjector::Neurite>& vNeurites,
			size_t numVerts,
			size_t numDodecagons
		);

		/*!
		 * \brief TODO refactor
		 * TODO: Documentation
		 */
		void connect_polygon_with_polygon
		(
			const std::vector<Vertex*>& from,
			const std::vector<Vertex*>& to,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			std::vector<std::pair<Vertex*, Vertex*> >& pairs
		);

		/*!
		 * \brief connew new TODO refactor
		 * TODO: Documentation
		 */
		void connect_new
		(
			Grid& g,
			SubsetHandler& sh,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			size_t newSomaIndex,
			size_t numDodecagons,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			SmartPtr<NeuriteProjector>
		);

		/*!
		 * \brief connects the plasma membrane with the soma surface
		 * TODO: Documentation
		 */
		void connect_pm_with_soma
		(
			size_t somaIndex,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			std::vector<std::vector<ug::Vertex*> >& rootNeurites,
			bool merge=false,
			size_t offset=0,
			bool mergeFirst=true
		);


		/*!
		 * \brief connects the er with the inner sphere surface (er)
		 * TODO: Documentation
		 */
		void connect_er_with_er
		(
			size_t somaIndex,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			std::vector<std::vector<ug::Vertex*> >& rootNeuritesInner,
			size_t offset,
			bool merge=false,
			bool mergeFirst=false
		);

		/*!
		 * \brief set somata mapping parameters
		 */
		void set_somata_mapping_parameters
		(
			Grid& g,
			SubsetHandler& sh,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::Mapping> >& aaMapping,
			size_t somaIndex,
			size_t erIndex,
			const SWCPoint& somaPoint
		);

		/*!
		 * \brief check if an edge is a root edge
		 * \param[in] edge
		 * \param[in] vPoints
		 * \return \c true if edge is a root edge and false otherwise
		 */
		bool IsRootEdge
		(
			const ug::Edge& edge,
			const std::vector<SWCPoint>& vPoints,
			const Grid::VertexAttachmentAccessor<APosition>& aaPos
		);
	}
}

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_UTIL_H
