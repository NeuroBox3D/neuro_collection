/*
 * neurite_refMarkAdjuster.h
 *
 *  Created on: 20.01.2017
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__NEURITE_REFMARKADJUSTER_H_
#define UG__PLUGINS__NEURO_COLLECTION__NEURITE_REFMARKADJUSTER_H_

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
 * TODO: Have the refinement mark adjuster decide on when it switches itself off.
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

#endif // UG__PLUGINS__NEURO_COLLECTION__NEURITE_REFMARKADJUSTER_H_
