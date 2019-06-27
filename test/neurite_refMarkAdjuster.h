/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-01-20
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_REFMARKADJUSTER_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_REFMARKADJUSTER_H

#include "lib_grid/refinement/ref_mark_adjuster_interface.h"
#include "lib_grid/refinement/hanging_node_refiner_multi_grid.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"


namespace ug {
namespace neuro_collection {


/**
 * @brief Refinement mark adjuster for anisotropic refinement of neurites
 * This refinement adjuster, unless disabled, will make refinement anisotropic along neurites.
 * Given a neurite projector containing spline data describing the neurites in question,
 * the refinement mark adjuster will mark all faces along the neurites RM_ANISOTROPIC,
 * as well as all edges tangential to the axial direction of the neurite. Only edges in axial
 * direction will be marked RM_REFINE. At branching points, it will use RM_COPY as well as at somata.
 *
 * At some point, the original anisotropy of the neurite mesh will have been "refined out"
 * and this refinement mark adjuster can be switched off using disable().
 *
 * @todo: Have the refinement mark adjuster decide on when it switches itself off.
 */
class NeuriteRefMarkAdjuster
: public IRefMarkAdjuster
{
	public:
        NeuriteRefMarkAdjuster
        (
            ConstSmartPtr<NeuriteProjector> np,
            ConstSmartPtr<ISubsetHandler> ssh,
            const Grid::VertexAttachmentAccessor<Attachment<vector3> >& aaPos
        ) : IRefMarkAdjuster(), m_spNP(np), m_ssh(ssh), m_aaPos(aaPos) {}

		virtual ~NeuriteRefMarkAdjuster() {}

        void disable()
        {enable(false);}

		virtual void ref_marks_changed
		(
			IRefiner& ref,
			const std::vector<Vertex*>& vrts,
			const std::vector<Edge*>& edges,
			const std::vector<Face*>& faces,
			const std::vector<Volume*>& vols
		);

	protected:
        void change_edge_mark(IRefiner& ref, Edge* e, const vector3* dir = NULL) const;
        void change_face_mark(IRefiner& ref, Grid& grid, Face* f, const vector3* dir = NULL) const;
        void change_vol_mark(IRefiner& ref, Grid& grid, Volume* v) const;

	private:
		ConstSmartPtr<NeuriteProjector> m_spNP;
		ConstSmartPtr<ISubsetHandler> m_ssh;
        const Grid::VertexAttachmentAccessor<Attachment<vector3> > m_aaPos;
};

void add_neurite_ref_mark_adjuster(IRefiner* ref, SmartPtr<NeuriteRefMarkAdjuster> nrma);


} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_REFMARKADJUSTER_H
