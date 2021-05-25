/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2018-01-30
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

#include "misc_util.h"

#include "lib_disc/domain_traits.h"                                 // for domain_traits
#include "lib_disc/function_spaces/approximation_space.h"           // for ApproximationSpace
#include "lib_grid/algorithms/geom_obj_util/anisotropy_util.h"      // for is_anisotropic
#include "lib_grid/algorithms/debug_util.h"                         // for ElementDebugInfo
#include "lib_grid/file_io/file_io.h"                               // for SaveGridToFile
#include "lib_grid/global_attachments.h"                            // for GlobalAttachments
#include "lib_grid/grid/neighborhood_util.h"                        // for GetConnectedNeighbor
#include "lib_grid/refinement/projectors/neurite_projector.h"       // for NeuriteProjector
#include "lib_grid/tools/grid_level.h"                              // for GridLevel
#include "lib_grid/tools/surface_view.h"                            // for SurfaceView
#include "lib_grid/refinement/projectors/projection_handler.h"      // for ProjectionHandler


namespace ug {
namespace neuro_collection {



template <typename TDomain>
void adjust_attachments
(
	SmartPtr<TDomain> domain
)
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename SurfaceView::traits<elem_type>::const_iterator const_iterator;

	SmartPtr<MultiGrid> grid = domain->grid();
	typename TDomain::position_accessor_type aaPos = domain->position_accessor();

	// access to neurite projector params
	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
			"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(*grid.get(), aSP);

	SurfaceView sv(domain->subset_handler());
	const_iterator iter = sv.begin<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_iterator iterEnd = sv.end<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);

	for (; iter != iterEnd; ++iter)
	{
		elem_type* elem = *iter;

		const size_t numVertices = elem->num_vertices();
		for (size_t i = 0; i < numVertices; i++) {
			GridObject* go = grid.get()->get_parent(elem->vertex(i));
			if (go != NULL) {
				/// VERTEX parent
				if (go->base_object_id() == VERTEX) {
					const Vertex* const v = dynamic_cast<Vertex*>(go);
					aaSurfParams[elem->vertex(i)] = aaSurfParams[v];
				}

				/// EDGE parent
				if (go->base_object_id() == EDGE) {
					const Edge* const e = dynamic_cast<Edge*>(go);
					const NPSP* const npsp1 = &aaSurfParams[e->vertex(0)];
					const NPSP* const npsp2 = &aaSurfParams[e->vertex(1)];
					aaSurfParams[elem->vertex(i)].axial =  (npsp1->axial+npsp2->axial)*0.5;
					aaSurfParams[elem->vertex(i)].radial = (npsp1->radial+npsp2->radial)*0.5;
				}

				/// FACE parent
				if (go->base_object_id() == FACE) {
					const Face* const f = dynamic_cast<Face*>(go);
					number axial = 0;
					number radial = 0;
					for (size_t i = 0; i < f->num_vertices(); i++) {
						axial += aaSurfParams[f->vertex(i)].axial;
						radial += aaSurfParams[f->vertex(i)].radial;
					}
					axial = axial / f->num_vertices();
					radial = radial / f->num_vertices();
				}

				/// VOLUME parent
				if (go->base_object_id() == VOLUME) {
					const Volume* const v = dynamic_cast<Volume*>(go);
					number axial = 0;
					number radial = 0;
					for (size_t i = 0; i < v->num_vertices(); i++) {
						axial += aaSurfParams[v->vertex(i)].axial;
						radial += aaSurfParams[v->vertex(i)].radial;
					}
					axial = axial / v->num_vertices();
					radial = radial / v->num_vertices();
				}
			}
		}
	}
}

template <typename TDomain>
void mark_anisotropic_in_local_neurite_direction
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	number thresholdRatio
)
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename SurfaceView::traits<elem_type>::const_iterator const_iterator;

	Grid& grid = *refiner->grid();
	typename TDomain::position_accessor_type aaPos = domain->position_accessor();

	// access to neurite projector params
	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
			"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(grid, aSP);

	// get surface view and prepare loop over all surface elements
	SurfaceView sv(domain->subset_handler());
	const_iterator iter = sv.begin<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_iterator iterEnd = sv.end<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);

	// loop elements for marking
	std::vector<Edge*> longEdges;
	for (; iter != iterEnd; ++iter)
	{
		longEdges.clear();
		AnisotropyState state = long_edges_of_anisotropic_elem(*iter, grid, aaPos, thresholdRatio, longEdges);
		if (state == ISOTROPIC)
			continue;

		UG_COND_THROW(!longEdges.size(), "Element is anisotropic, but no long edges present.");

		// check whether edges point in local neurite direction
		const Edge* const longEdge = longEdges[0];
		number axialDistance = fabs(aaSurfParams[longEdge->vertex(1)].axial-aaSurfParams[longEdge->vertex(0)].axial);

		if (axialDistance > 1e-6)
		{
			// check this is not a BP, ignore BPs
			uint32_t nid = aaSurfParams[(*iter)->vertex(0)].neuriteID & ((1 << 20) - 1);
			uint32_t thisNID = aaSurfParams[(*iter)->vertex(0)].neuriteID & ((1 << 20) - 1);
			for (size_t i = 1; i < (*iter)->num_vertices(); ++i)
				thisNID = std::max(thisNID, aaSurfParams[(*iter)->vertex(i)].neuriteID & ((1 << 20) - 1));
			if (thisNID != nid)
			{
				continue;
			}

			// mark elem
			refiner->mark(*iter, RM_CLOSURE);

			// mark long edges
			const size_t nEdges = longEdges.size();
			for (size_t e = 0; e < nEdges; ++e)
				refiner->mark(longEdges[e], RM_FULL);

			// mark all sides
			typename Grid::traits<side_type>::secure_container sl;
			grid.associated_elements(sl, *iter);
			const size_t slSz = sl.size();
			for (size_t s = 0; s < slSz; ++s)
				if (refiner->get_mark(sl[s]) != RM_FULL)
					refiner->mark(sl[s], RM_CLOSURE);
		}
	}
}


void MarkNeuriteForAxialRefinement(SmartPtr<IRefiner> refiner, SmartPtr<Domain3d> domain)
{
	typedef SurfaceView::traits<Volume>::const_iterator const_vol_it;

	Grid& grid = *refiner->grid();
	Domain3d::position_accessor_type aaPos = domain->position_accessor();

#ifdef UG_PARALLEL
	const DistributedGridManager& dgm = *domain->distributed_grid_manager();
#endif

	// access to neurite projector params
	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
			"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	if (grid.num_volumes())
	{
		UG_COND_THROW(!grid.has_vertex_attachment(aSP),
			"Grid does not have the required attachment 'npSurfParams'.");
		aaSurfParams.access(grid, aSP);
	}

	// volume attachment to work with here
	Attachment<int> aState;
	grid.attach_to_volumes(aState);
	Grid::VolumeAttachmentAccessor<Attachment<int> > aaState;
	aaState.access(grid, aState);


	// get surface view and prepare loop over all surface elements
	SurfaceView sv(domain->subset_handler());
	const_vol_it iter = sv.begin<Volume>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_vol_it iterEnd = sv.end<Volume>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);

	std::vector<Volume*> volumes;
	std::vector<Edge*> edges;

	grid.begin_marking();
	for (; iter != iterEnd; ++iter)
	{
		if (aaState[*iter])  // already undecided
			continue;

		// start a queue of elements all at the same axial position
		// to decide whether we are on a neurite segment or a BP
		volumes.clear();
		edges.clear();

		const size_t nInitialVrt = (*iter)->num_vertices();
		UG_COND_THROW(nInitialVrt != 8, "Volume element with " << nInitialVrt << " vertices found in grid.\n"
			"But grid must only contain hexahedra for this function to work.\n"
				<< ElementDebugInfo(grid, *iter));

		uint32_t nid = aaSurfParams[(*iter)->vertex(0)].neuriteID & ((1 << 20) - 1);
		for (size_t i = 1; i < nInitialVrt; ++i)
			nid = std::max(nid, aaSurfParams[(*iter)->vertex(i)].neuriteID & ((1 << 20) - 1));

//UG_LOGN(ElementDebugInfo(grid, *iter) << " with NID " << nid << ":");

		bool isBP = false;
		std::queue<Volume*> q;
		q.push(*iter);
		while (!q.empty())
		{
			Volume* elem = q.front();
			q.pop();

			int& state = aaState[elem];
			if (state > 1)  // BP or neurite
			{
//UG_LOGN("  BP because of neighbor " << ElementDebugInfo(grid, elem) << " which is already assigned.");
				isBP = true;
				break;
			}

			if (state == 1) // visited, but undecided
				continue;

			state = 1;
			volumes.push_back(elem);

			const size_t nVrt = elem->num_vertices();
			UG_COND_THROW(nVrt != 8, "Volume element with " << nVrt << " vertices found in grid.\n"
				"But grid must only contain hexahedra for this function to work.\n"
				<< ElementDebugInfo(grid, elem));

			std::vector<number> vrtAxPos(8);

			// check this is not a BP
			uint32_t thisNID = aaSurfParams[elem->vertex(0)].neuriteID & ((1 << 20) - 1);
			for (size_t i = 1; i < nVrt; ++i)
				thisNID = std::max(thisNID, aaSurfParams[elem->vertex(i)].neuriteID & ((1 << 20) - 1));
			if (thisNID != nid)
			{
				isBP = true;
//UG_LOGN("  BP because of " << ElementDebugInfo(grid, *iter) << " with NID " << thisNID);
				break;
			}

			// axial positions need to be clearly discriminable
			for (size_t i = 0; i < nVrt; ++i)
			{
				if ((aaSurfParams[elem->vertex(i)].neuriteID & ((1 << 20) - 1)) < thisNID)
					vrtAxPos[i] = 0.0;
				else
					vrtAxPos[i] = aaSurfParams[elem->vertex(i)].axial;
			}

			std::sort(vrtAxPos.begin(), vrtAxPos.end());
			if (vrtAxPos[0] == 0.0 && vrtAxPos[3] > 0.0)
				for (size_t i = 0; i < 3; ++i)
					vrtAxPos[i] =  vrtAxPos[3];

			number axLength = vrtAxPos[4] - vrtAxPos[3];
			if (vrtAxPos[3] - vrtAxPos[0] >= axLength || vrtAxPos[7] - vrtAxPos[4] >= axLength)
			{
				isBP = true;
/*
UG_LOGN("  BP because of " << ElementDebugInfo(grid, *iter) << " which is malformed:");
UG_LOG("    axial positions: ");
for (size_t j = 0; j < nVrt; ++j)
	UG_LOG(vrtAxPos[j] << " ");
UG_LOGN("");
*/
				break;
			}

			// identify lateral sides and edges
			Grid::traits<Face>::secure_container sl;
			grid.associated_elements(sl, elem);
			const size_t slSz = sl.size();
			for (size_t s = 0; s < slSz; ++s)
			{
				Face* side = sl[s];

				if (grid.is_marked(side))
				{
//UG_LOGN("  side " << s << " already marked");
					continue;
				}

				grid.mark(side);

				// if a side separates two partitions, then this must be a BP
				// (neurites should not be cut axially)
#ifdef UG_PARALLEL
				if (dgm.is_in_horizontal_interface(side))
				{
					isBP = true;
					break;
				}
#endif
				// check orientation
				size_t nAxialEdges = 0;
				Grid::traits<Edge>::secure_container el;
				grid.associated_elements(el, side);
				const size_t elSz = el.size();
				UG_COND_THROW(elSz != 4, "Side element with " << elSz << " edges found in grid.\n"
					"But grid must only contain quadrilateral side elements for this function to work.");
				for (size_t e = 0; e < elSz; ++e)
				{
					Edge* ed = el[e];
					if ((((aaSurfParams[ed->vertex(0)].neuriteID & ((1 << 20) - 1)) < thisNID || aaSurfParams[ed->vertex(0)].axial <= vrtAxPos[3])
							&& ((aaSurfParams[ed->vertex(1)].neuriteID & ((1 << 20) - 1)) == thisNID && aaSurfParams[ed->vertex(1)].axial >= vrtAxPos[4]))
						|| (((aaSurfParams[ed->vertex(1)].neuriteID & ((1 << 20) - 1)) < thisNID || aaSurfParams[ed->vertex(1)].axial <= vrtAxPos[3])
							&& ((aaSurfParams[ed->vertex(0)].neuriteID & ((1 << 20) - 1)) == thisNID && aaSurfParams[ed->vertex(0)].axial >= vrtAxPos[4])))
					{
						if (!grid.is_marked(ed))
						{
							edges.push_back(ed);
							grid.mark(ed);
						}
						++nAxialEdges;
					}
				}
//UG_LOGN("  nAxialEdges = " << nAxialEdges);
				if (nAxialEdges == 2)
				{
					Volume* opp = GetConnectedNeighbor(grid, side, elem);
					if (opp && aaState[opp] != 1)
						q.push(opp);
				}
			}
		}

		if (isBP)
		{
			const size_t nVol = volumes.size();
			for (size_t i = 0; i < nVol; ++i)
			{
				aaState[volumes[i]] = 3;  // BP
				refiner->mark(volumes[i], RM_CLOSURE);
			}
		}
		else
		{
			const size_t nVol = volumes.size();
			for (size_t i = 0; i < nVol; ++i)
			{
				aaState[volumes[i]] = 2;  // regular neurite
				refiner->mark(volumes[i], RM_CLOSURE);
//UG_LOGN("  " << ElementDebugInfo(grid, volumes[i]));
			}
			const size_t nEdge = edges.size();
			for (size_t i = 0; i < nEdge; ++i)
				refiner->mark(edges[i], RM_REFINE);
		}
	}
	grid.end_marking();

	grid.detach_from_volumes(aState);


	typedef SurfaceView::traits<Face>::const_iterator const_face_it;
	const_face_it itFace = sv.begin<Face>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_face_it itFaceEnd = sv.end<Face>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; itFace != itFaceEnd; ++itFace)
		refiner->mark(*itFace, RM_CLOSURE);

	typedef SurfaceView::traits<Edge>::const_iterator const_edge_it;
	const_edge_it itEdge = sv.begin<Edge>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_edge_it itEdgeEnd = sv.end<Edge>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; itEdge != itEdgeEnd; ++itEdge)
	{
		if (refiner->get_mark(*itEdge) != RM_REFINE)
			refiner->mark(*itEdge, RM_CLOSURE);
	}
}


template <typename TDomain>
void RemoveAllNonDefaultRefinementProjectors(SmartPtr<TDomain> dom)
{
	SmartPtr<RefinementProjector> rp = dom->refinement_projector();
	ProjectionHandler* ph = dynamic_cast<ProjectionHandler*>(rp.get());
	UG_COND_THROW(!ph, "Domain refinement projector is not a projection handler.");

	const size_t nProj = ph->num_projectors();
	for (size_t i = 0; i < nProj; ++i)
		ph->set_projector(i, ph->default_projector());
}


////////////////////////////////////////////////////////////////////////
/// GetCoordinatesFromVertexByIndex
////////////////////////////////////////////////////////////////////////
const vector3* GetCoordinatesFromVertexByIndex(Grid& grid, const int index)
{
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	if(!grid.has_vertex_attachment(aPosition)) { grid.attach_to_vertices(aPosition); }

	ConstVertexIterator vit = grid.vertices_begin();
	int i = 0;
	while (vit++ != grid.vertices_end()) { if (index == i++) { break; } }

	return new vector3(aaPos[*vit]);
}

// explicit template specializations
#ifdef UG_DIM_1
	template void mark_anisotropic_in_local_neurite_direction<Domain1d>(SmartPtr<IRefiner>, SmartPtr<Domain1d>, number);
	template void RemoveAllNonDefaultRefinementProjectors(SmartPtr<Domain1d>);
	template void adjust_attachments(SmartPtr<Domain1d>);
#endif
#ifdef UG_DIM_2
	template void mark_anisotropic_in_local_neurite_direction<Domain2d>(SmartPtr<IRefiner>, SmartPtr<Domain2d>, number);
	template void RemoveAllNonDefaultRefinementProjectors(SmartPtr<Domain2d>);
	template void adjust_attachments(SmartPtr<Domain2d>);
#endif
#ifdef UG_DIM_3
	template void mark_anisotropic_in_local_neurite_direction<Domain3d>(SmartPtr<IRefiner>, SmartPtr<Domain3d>, number);
	template void RemoveAllNonDefaultRefinementProjectors(SmartPtr<Domain3d>);
	template void adjust_attachments(SmartPtr<Domain3d>);
#endif


}  // namespace neuro_collection
}  // namespace ug
