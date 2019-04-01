/*
 * neurite_axial_refinement_marker.h
 *
 *  Created on: 2019-03-15
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__NEURITE_AXIAL_REFINEMENT_MARKER_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__NEURITE_AXIAL_REFINEMENT_MARKER_H


#include "common/util/smart_pointer.h"  // for SmartPtr
#include "lib_disc/domain.h"  // for Domain3d
#include "lib_grid/attachments/attachment_pipe.h"  // for Attachment
#include "lib_grid/grid/grid.h"  // for AttachmentAccessors
#include "lib_grid/grid/grid_base_objects.h"  // for Edge
#include "lib_grid/refinement/projectors/neurite_projector.h"  // for SurfaceParams
#include "lib_grid/refinement/refiner_interface.h"  // for IRefiner

// configuration file for compile options
#include "nc_config.h"

#ifdef NC_WITH_PARMETIS
#include "../../Parmetis/src/unificator_interface.h"  // for IUnificator
#endif

namespace ug {
namespace neuro_collection {


class NeuriteAxialRefinementMarker
#ifdef NC_WITH_PARMETIS
	: public parmetis::IUnificator<Volume>
#endif
{
	public:
		NeuriteAxialRefinementMarker(SmartPtr<Domain3d> dom);

		virtual ~NeuriteAxialRefinementMarker();

		void mark(SmartPtr<IRefiner> refiner);

#ifdef NC_WITH_PARMETIS
		typedef Volume::side side_t;
		typedef Attachment<int> AElemIndex;
		typedef Attachment<std::vector<int> > AElemIndices;

		// inherited from IUnificator
		virtual void unify
		(
			MultiGrid* mg,
			int lvl,
			int localOffset,
			const Grid::AttachmentAccessor<Volume, AElemIndex>& aaElemInd, // local indices!
			const Grid::AttachmentAccessor<side_t, AElemIndices>& aaSideElemInd, // global indices!
			std::vector<std::pair<int, int> >& unificationPairs // global indices!
		) const;
#endif

	protected:
		void find_bp_volumes();

		inline bool is_bp_volume(Volume* vol) const;

	private:
		void mark_bp_volumes(MultiGrid* mg, int lvl) const;
		bool is_central_bp_vol(Volume*) const;

	protected:
		Attachment<bool> m_aBP;
		Grid::VolumeAttachmentAccessor<Attachment<bool> > m_aaBP;

		typedef NeuriteProjector::SurfaceParams NPSP;
		Grid::VertexAttachmentAccessor<Attachment<NPSP> > m_aaSurfParams;

		SmartPtr<Domain3d> m_spDom;

};


}  // namespace neuro_collection
}  // namespace ug

#endif  // UG__PLUGINS__NEURO_COLLECTION__UTIL__NEURITE_AXIAL_REFINEMENT_MARKER_H
