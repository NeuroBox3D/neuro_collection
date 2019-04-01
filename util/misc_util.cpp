/*
 * misc_util.cpp
 *
 *  Created on: 2018-01-30
 *      Author: mbreit
 */

#include "misc_util.h"

#include "lib_disc/domain_traits.h"                                 // for domain_traits
#include "lib_disc/function_spaces/approximation_space.h"           // for ApproximationSpace
#include "lib_grid/algorithms/geom_obj_util/anisotropy_util.h"      // for is_anisotropic
#include "lib_grid/algorithms/debug_util.h"                         // for ElementDebugInfo
#include "lib_grid/global_attachments.h"                            // for GlobalAttachments
#include "lib_grid/grid/neighborhood_util.h"                        // for GetConnectedNeighbor
#include "lib_grid/refinement/projectors/neurite_projector.h"       // for NeuriteProjector
#include "lib_grid/tools/grid_level.h"                              // for GridLevel
#include "lib_grid/tools/surface_view.h"                            // for SurfaceView


namespace ug {
namespace neuro_collection {


template <typename TDomain>
void mark_global(SmartPtr<IRefiner> refiner, SmartPtr<TDomain> domain)
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename SurfaceView::traits<elem_type>::const_iterator const_iterator;

	// get surface view
	SurfaceView sv(domain->subset_handler());

	// loop elements for marking
	const_iterator iter = sv.begin<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_iterator iterEnd = sv.end<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; iter != iterEnd; ++iter)
		refiner->mark(*iter, RM_FULL);
}



template <typename TDomain>
void mark_anisotropic
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
		if (state != ISOTROPIC)
		{
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


template <typename TDomain>
void mark_anisotropic_onlyX
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

		// check whether edges point in x direction
		Edge* longEdge = longEdges[0];
		MathVector<TDomain::dim> dir;
		VecSubtract(dir, aaPos[longEdge->vertex(1)], aaPos[longEdge->vertex(0)]);
		VecNormalize(dir, dir);
		if (fabs(dir[0]) > 0.9)
		{
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




// explicit template specializations
#ifdef UG_DIM_1
	template void mark_global<Domain1d>(SmartPtr<IRefiner>, SmartPtr<Domain1d>);
	template void mark_anisotropic<Domain1d>(SmartPtr<IRefiner>, SmartPtr<Domain1d>, number);
	template void mark_anisotropic_onlyX<Domain1d>(SmartPtr<IRefiner>, SmartPtr<Domain1d>, number);
#endif
#ifdef UG_DIM_2
	template void mark_global<Domain2d>(SmartPtr<IRefiner>, SmartPtr<Domain2d>);
	template void mark_anisotropic<Domain2d>(SmartPtr<IRefiner>, SmartPtr<Domain2d>, number);
	template void mark_anisotropic_onlyX<Domain2d>(SmartPtr<IRefiner>, SmartPtr<Domain2d>, number);
#endif
#ifdef UG_DIM_3
	template void mark_global<Domain3d>(SmartPtr<IRefiner>, SmartPtr<Domain3d>);
	template void mark_anisotropic<Domain3d>(SmartPtr<IRefiner>, SmartPtr<Domain3d>, number);
	template void mark_anisotropic_onlyX<Domain3d>(SmartPtr<IRefiner>, SmartPtr<Domain3d>, number);
#endif


}  // namespace neuro_collection
}  // namespace ug
