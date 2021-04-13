/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-03-15
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


/**
 * @brief This class is used to refine neurites anisotropically in axial direction.
 *
 * It is aimed at the refinement of coarse grids created using the neurites_from_swc
 * functions, which contain anisotropic hexahedra in the neurites, but isotropic hexahedra
 * in the branching points.
 * The goal is to automatically detect whether an element belongs to a neurite or a branching
 * point and refine anisotropically (only axially) in the neurites and by copying (no refinement)
 * in the branching points.
 *
 * The class tries to identify the branching point volumes and does that in a fashion that strongly
 * depends on the way the coarse grids are created in the neurites_from_swc routines.
 *
 * This identification method also requires radial neighbors to stick together during redistribution,
 * which is why the class inherits from parmetis::IUnificator<Volume>, so that it can be used as such
 * by a ClusteredDualGraphManager during Parmetis redistribution.
 *
 * @note The class only works on neurites_from_swc-created geometries with hexahedral elements.
 *       It is not perfect. Sometimes, the refinement will not be properly anisotropic.
 * @note The class should work both on ER-containing and ER-less geometries.
 */
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

		mutable SmartPtr<Domain3d> m_spDom;

};


}  // namespace neuro_collection
}  // namespace ug

#endif  // UG__PLUGINS__NEURO_COLLECTION__UTIL__NEURITE_AXIAL_REFINEMENT_MARKER_H
