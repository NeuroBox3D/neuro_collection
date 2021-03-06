/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-05-06
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_GRID_GENERATION_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_GRID_GENERATION_H

#include "lib_grid/grid/grid.h"
#include "lib_grid/grid/geometry.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"
#include "common/math/ugmath_types.h"
#include "types.h"

namespace ug {
	namespace neuro_collection {
		/*!
		 * \brief creates neurites with ER
		 *
		 * \param[in] vNeurites
		 * \param[in] vPos
		 * \param[in] vR
		 * \param[in] nid
		 * \param[in] erScaleFactor
		 * \param[in] anisotropy
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] aaSurfParams
		 * \param[in] sh
		 * \param[in] connectingVrts
		 * \param[in] connectingEdges
		 * \param[in] connectingFaces
		 * \param[in] initialOffset
		 * \param[out] outVerts
		 * \param[out] outVertsInner
		 * \param[out] outRads
		 * \param[out] outRadsInner
		 */
		void create_neurite_with_er
		(
			const std::vector<NeuriteProjector::Neurite>& vNeurites,
			const std::vector<std::vector<vector3> >& vPos,
			const std::vector<std::vector<number> >& vR,
			size_t nid,
			number erScaleFactor,
			number anisotropy,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::Mapping> > aaMappingParams,
			SubsetHandler& sh,
			number blowUpFactor = 1.0,
			std::vector<Vertex*>* connectingVrts = NULL,
			std::vector<Edge*>* connectingEdges = NULL,
			std::vector<Face*>* connectingFaces = NULL,
			number initialOffset = 0.0,
			std::vector<Vertex*>* outVerts = NULL,
			std::vector<Vertex*>* outVertsInner = NULL,
			std::vector<number>* outRads = NULL,
			std::vector<number>* outRadsInner = NULL,
			std::vector<SWCPoint>* points = NULL,
			MeasuringSubsetCollection* subsets = NULL,
			int bid = -1,
			const std::string& option = "user",
			number desiredSegLength = -1,
			bool onlyCaps = false
		);

		/*!
		 * \brief creates neurites with ER
		 *
		 * \param[in] vNeurites
		 * \param[in] vPos
		 * \param[in] vR
		 * \param[in] nid
		 * \param[in] erScaleFactor
		 * \param[in] anisotropy
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] aaSurfParams
		 * \param[in] sh
		 * \param[out] outVerts
		 * \param[out] outVertsInner
		 * \param[out] outRads
		 * \param[out] outRadsInner
		 */
		void create_neurite_with_er
		(
			const std::vector<NeuriteProjector::Neurite>& vNeurites,
			const std::vector<std::vector<vector3> >& vPos,
			const std::vector<std::vector<number> >& vR,
			size_t nid,
			number erScaleFactor,
			number anisotropy,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::Mapping> > aaMappingParams,
			SubsetHandler& sh,
			number blowUpFactor,
			std::vector<Vertex*>* outVerts,
			std::vector<Vertex*>* outVertsInner,
			std::vector<number>* outRads,
			std::vector<number>* outRadsInner,
			std::vector<SWCPoint>* points,
			MeasuringSubsetCollection* subsets,
			int bid,
			const std::string& option,
			number segLength,
			bool onlyCaps
		);


		/*!
		 * \brief creates neurite surface
		 *
		 * \param[in] vNeurites
		 * \param[in] vPos
		 * \param[in] vR
		 * \param[in] nid
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] aaSurfParams
		 * \param[in] sh
		 * \param[in] connectingVrts
		 * \param[in] connectingEdges
		 * \param[in] initialOffset
		 */
		void create_neurite_surf(
			const std::vector<NeuriteProjector::Neurite>& vNeurites,
			const std::vector<std::vector<vector3> >& vPos,
			const std::vector<std::vector<number> >& vR,
			size_t nid,
			number anisotropy, Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<
			Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			std::vector<Vertex*>* connectingVrts = NULL,
			std::vector<Edge*>* connectingEdges = NULL,
			number initialOffset = 0.0
		);

		/*!
		 * \brief Calculate and returns length of the spline in units of radii
		 *
		 * \param[in] t_start
		 * \param[in] t_end
		 * \param[in] neurite
		 * \param[in] startSec
		 *
		 * \return \c number
		 */
		number calculate_length_over_radius
		(
			number t_start,
			number t_end,
			const NeuriteProjector::Neurite& neurite,
			size_t startSec
		);

		/*!
		 * \brief Calculates and returns length of the spline in arclength
		 *
		 * \param[in] t_start
		 * \param[in] t_end
		 * \param[in] neurite
		 * \param[in] startSec
		 *
		 * \return \c number
		 */
		number calculate_length_over_radius_variant
		(
			number t_start,
			number t_end,
			const NeuriteProjector::Neurite& neurite,
			size_t startSec
		);

		/*!
		 * \brief Calculate the axial positions to achieved desired anisotropy
		 *
		 * This is the original method to use for achieving the desired
		 * anisotropy in the refinement limit on the surface
		 *
		 * \param[out] segAxPosOut evaluation points
		 * \param[in] t_start start position within any spline section
		 * \param[in] t_end end position within any spline section
		 * \param[in] neurite current neurite
		 * \param[in] startSec first spline section of neurite
		 * \param[in] segLength segment length
		 */
		void calculate_segment_axial_positions
		(
			std::vector<number>& segAxPosOut,
			number t_start,
			number t_end,
			const NeuriteProjector::Neurite& neurite,
			size_t startSec,
			number segLength
		);

		/*!
		 * \brief Force evaluation positions to coincide with support nodes
		 *
		 * The methods forces the evaluation positions along the spline curve
		 * to coincide with the (axial) 1d support nodes of the actual spline
		 * for a given neurite
		 *
		 * \param[out] segAxPosOut evaluation points
		 * \param[in] t_start start position within any spline section
		 * \param[in] t_end end position within any spline section
		 * \param[in] neurite current neurite
		 * \param[in] startSec first spline section of neurite
		 */
		void calculate_segment_axial_positions_at_support_nodes
		(
			std::vector<number>& segAxPosOut,
			number t_start,
			number t_end,
			const NeuriteProjector::Neurite& neurite,
			size_t startSec
		);

		/*!
		 * \brief Force evaluations positions to be spaced equally in arclength
		 *
		 * The methods forces the evaluation positions along the spline curve
		 * to be spaced equally wrt to arclength of the spline curve itself for
		 * a given neurite
		 *
		 * \param[out] segAxPosOut evaluation points
		 * \param[in] t_start start position within any spline section
		 * \param[in] t_end end position within any spline section
		 * \param[in] neurite current neurite
		 * \param[in] startSec first spline section of neurite
		 * \param[in] segLength segment length
		 */
		void calculate_segment_axial_positions_constant_seg_length
	    (
	        std::vector<number>& segAxPosOut,
	        number t_start,
	        number t_end,
	        const NeuriteProjector::Neurite& neurite,
	        size_t startSec,
	        number segLength
	    );

		/*!
		 * \brief Force evaluations positions to be spaced equally in arclength
		 *
		 * The methods forces the evaluation positions along the spline curve
		 * to be spaced equally wrt to arclength of the spline curve itself for
		 * a given neurite
		 *
		 * Note this is a variant of the above methods for sole purpose of testing
		 *
		 * \param[out] segAxPosOut evaluation points
		 * \param[in] t_start start position within any spline section
		 * \param[in] t_end end position within any spline section
		 * \param[in] neurite current neurite
		 * \param[in] startSec first spline section of neurite
		 * \param[in] segLength segment length
		 */
        void calculate_segment_axial_positions_constant_seg_length_variant
        (
   	        std::vector<number>& segAxPosOut,
   	        number t_start,
   	        number t_end,
   	        const NeuriteProjector::Neurite& neurite,
   	        size_t startSec,
   	        number segLength
   	    );

		/*!
		 * \brief create neurite 1d
		 *
		 * \param[in] vNeurites
		 * \param[in] vPos
		 * \param[in] vR
		 * \param[in] nid
		 * \param[in] g
		 * \param[in] aaPos
		 * \param[in] aaSurfParams
		 * \param[in] aaDiam
		 * \param[in] connectingVrts
		 */
		 void create_neurite_1d(
			const std::vector<NeuriteProjector::Neurite>& vNeurites,
			const std::vector<std::vector<vector3> >& vPos,
			const std::vector<std::vector<number> >& vR,
			size_t nid,
			number anisotropy, Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<
			Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid::VertexAttachmentAccessor<Attachment<number> >& aaDiam,
			Vertex* connectingVrt
		);

		 /*!
		  * \brief creates the root neurite's vertices, edges and faces
		  *
		  * Note that the inner root neurite (ER) vertices and outer root
		  * neurite (PM) vertices as well the corresponding radii are stored
		  * in the variables outVerts, outVertsInner, outRads and outRadsInner.
		  * \param[in] vNeurites
		  * \param[in] vPos
		  * \param[in] vR
		  * \param[in] nid
		  * \param[in,out] g
		  * \param[in,out] sh
		  * \param[in] erScaleFactor
		  * \param[in,out] aaPos
		  * \param[out] outVerts
		  * \param[out] outVertsInner
		  * \param[out] outRads
		  * \param[out] outRadsInner
		  * \param[in] withER
		  *
		  */
		 void create_neurite_root_vertices
		 (
			const std::vector<NeuriteProjector::Neurite>& vNeurites,
			const std::vector<std::vector<vector3> >& vPos,
			const std::vector<std::vector<number> >& vR,
			size_t nid,
			Grid& g,
			SubsetHandler& sh,
			number erScaleFactor,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			std::vector<Vertex*>* outVerts,
			std::vector<Vertex*>* outVertsInner,
			std::vector<number>* outRads,
			std::vector<number>* outRadsInner,
			bool withER=true
		);

		 /*!
		  * \brief smoothes 1d line-graph structure constrained
		  *
		  * Depending on the largest radius the first vertices before or after a
		  * branch are ignored if they do not belong to the root neurite branch
		  * \param[in] vPoints points of graph
		  * \param[in] n number of iterations
		  * \param[in] h smoothing parameter
		  * \param[in] gamma smoothing parameter
		  * \param[in] minAngle min angle for root branch detection
		  * \param[in] maxRadiusRatio max radii ratio for root branch detection
		  */
		void ConstrainedSmoothingAlongRootBranch
		(
			std::vector<SWCPoint>& vPointsIn,
			size_t n,
			number h,
			number gamma,
			number minAngle,
			number maxRadiusRatio
		);

		/*!
		 * \brief regularizes branching points
		 *
		 * \param[in,out] vPoints
		 * \param[in] orthogonalize
		 *  1. Find BP.
		 *  2. Identify root branch and branching children by min angle criterion
		 *  3. Collect point before branching point (P), branching point itself
		 *     (B) and point on root branch after branching point (Q) and
		 *     collect all other points after branching point in a vector Rs
		 *  4. Save branching point B, create edge between P and Q as PQ
	     *  5. Project B onto PQ as B' and erase B afterwards
	     *  If orthogonalize:
	     *  6. Create additional point A normal to PQ and connect to B'
	     *  7. Before connecting orthogonal point A is rotated (orthogonally)
	     *  to find the best connection angle towards B'
	     *  8. Connect A to each point in vector Rs
		 */
		void RegularizeBranchingPoints
		(
			std::vector<SWCPoint>& vPoints,
			bool orthogonalize=false
		);

		/*!
		 * \brief Inserts a vertex at root branching neurites
		 *
		 * \param[in, out] vPoints list of SWC points
		 */
		void MitigateRootBranchingNeurites
		(
			std::vector<SWCPoint>& vPoints
		);

		/*!
		 * \brief Retriangulate the connecting regions
		 *
		 * \param[in, out] sh subset handler
		 * \param[in, out] grid
		 * \param[in, out] aaPos position accessor
		 * \param[in] si subset index
		 * \param[in] minAngle minimal dihedral angle
		 */
		void RetriangulateConnectingRegions
		(
			SubsetHandler& sh,
			Grid& grid,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			size_t si,
			number minAngle
		);

	} // end namespace neuro_collection
} // end namespace ug


#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_GRID_GENERATION_H
