/*!
 * neurite_grid_generation.h
 *
 *  Created on: May 6, 2019
 *      Author: stephan
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_GRID_GENERATION_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_GRID_GENERATION_H
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
#include "neurite_util.h"

namespace ug {
	namespace neuro_collection {
		/*!
		 * \brief creates neurites with ER
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
			SubsetHandler& sh,
			std::vector<Vertex*>* connectingVrts = NULL,
			std::vector<Edge*>* connectingEdges = NULL,
			std::vector<Face*>* connectingFaces = NULL,
			number initialOffset = 0.0
		);


		/*!
		 * \brief creates neurite surface
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
		 * \brief calculate length over radius
		 * \param[in] t_start
		 * \param[in] t_end
		 * \param[in] neurite
		 * \param[in] startSec
		 */
		number calculate_length_over_radius
		(
			number t_start,
			number t_end,
			const NeuriteProjector::Neurite& neurite,
			size_t startSec
		);


		/*!
		 * \brief calculates segment axial positions
		 * \param[in] segAxPosOut
		 * \param[in] t_start
		 * \param[in] t_end
		 * \param[in] neurite
		 * \param[in] startSec
		 * \param[in] segLength
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
		 * \brief create neurite 1d
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
	}
}


#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_GRID_GENERATION_H
