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

#include "neurite_axial_refinement_marker.h"

#include "common/error.h"  // for UG_THROW
#include "lib_grid/global_attachments.h"  // for GlobalAttachments
#include "lib_grid/grid/grid_util.h"  // for CompareVertices
#include "lib_grid/tools/surface_view.h"  // for SurfaceView
#ifdef NC_WITH_PARMETIS
#include "lib_grid/grid/grid_base_object_traits.h"  // for geometry_traits
#include "lib_grid/grid/neighborhood_util.h"  // for GetConnectedNeighbor
#endif

#include <set>

namespace ug {
namespace neuro_collection {


NeuriteAxialRefinementMarker::NeuriteAxialRefinementMarker(SmartPtr<Domain3d> dom)
: m_spDom(dom)
{
	SmartPtr<MultiGrid> mg = m_spDom->grid();

	// attach and access helper attachment
	if (!mg->has_volume_attachment(m_aBP))
		mg->attach_to_volumes_dv(m_aBP, false);
	else
		UG_THROW("Grid already has helper attachment attached.");

	m_aaBP.access(*mg, m_aBP);

	// access neurite projector local coordinates attachment
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
		"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	UG_COND_THROW(mg->num_vertices() && !mg->has_vertex_attachment(aSP),  // make sure the attachment has been read from file
		"npSurfParams attachment not attached.");
	if (!mg->has_vertex_attachment(aSP))  // only the file reading proc MUST already have it attached
		mg->attach_to_vertices(aSP);
	m_aaSurfParams.access(*mg, aSP);
}


NeuriteAxialRefinementMarker::~NeuriteAxialRefinementMarker()
{
	SmartPtr<MultiGrid> mg = m_spDom->grid();
	mg->detach_from_volumes(m_aBP);
}


void NeuriteAxialRefinementMarker::mark(SmartPtr<IRefiner> refiner)
{
	// identify BP edges
	find_bp_volumes();

	// get surface view and grid
	SurfaceView sv(m_spDom->subset_handler());
	SmartPtr<MultiGrid> mg = m_spDom->grid();

	// mark all volumes and faces RM_CLOSURE
	typedef SurfaceView::traits<Volume>::const_iterator const_vol_it;
	const_vol_it itVol = sv.begin<Volume>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_vol_it itVolEnd = sv.end<Volume>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; itVol != itVolEnd; ++itVol)
		refiner->mark(*itVol, RM_CLOSURE);

	typedef SurfaceView::traits<Face>::const_iterator const_face_it;
	const_face_it itFace = sv.begin<Face>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_face_it itFaceEnd = sv.end<Face>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; itFace != itFaceEnd; ++itFace)
		refiner->mark(*itFace, RM_CLOSURE);

	// mark all axial edges RM_REFINE
	std::vector<EdgeDescriptor> edgeDescs1(4);
	std::vector<EdgeDescriptor> edgeDescs2(4);
	std::vector<EdgeDescriptor> edgeDescs3(4);
	float axPos[8];

	Grid::traits<Edge>::secure_container el;

	itVol = sv.begin<Volume>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; itVol != itVolEnd; ++itVol)
	{
		Volume* vol = *itVol;

		if (is_bp_volume(vol))
			continue;

		// find the axial edges
		Hexahedron* hex = dynamic_cast<Hexahedron*>(vol);
		UG_COND_THROW(!hex, "Found volume that is not a hexahedron.\n"
			"This implementation can only handle hexahedron grids.");

		uint32_t nid = 0;
		for (size_t i = 0; i < 8; ++i)
			nid = std::max(nid, m_aaSurfParams[hex->vertex(i)].neuriteID & ((1 << 20) - 1));

		// axial pos is only valid of neurite ID is that of the parent for BP vertices
		// otherwise take 0.0
		for (size_t i = 0; i < 8; ++i)
		{
			if ((m_aaSurfParams[hex->vertex(i)].neuriteID & ((1 << 20) - 1)) == nid)
				axPos[i] = m_aaSurfParams[hex->vertex(i)].axial;
			else
				axPos[i] = 0.0;
		}

		number length1 = 0.0;
		edgeDescs1[0].set_vertex(0, hex->vertex(0));
		edgeDescs1[0].set_vertex(1, hex->vertex(1));
		length1 += fabs(axPos[0] - axPos[1]);
		edgeDescs1[1].set_vertex(0, hex->vertex(2));
		edgeDescs1[1].set_vertex(1, hex->vertex(3));
		length1 += fabs(axPos[2] - axPos[3]);
		edgeDescs1[2].set_vertex(0, hex->vertex(4));
		edgeDescs1[2].set_vertex(1, hex->vertex(5));
		length1 += fabs(axPos[4] - axPos[5]);
		edgeDescs1[3].set_vertex(0, hex->vertex(6));
		edgeDescs1[3].set_vertex(1, hex->vertex(7));
		length1 += fabs(axPos[6] - axPos[7]);

		number length2 = 0.0;
		edgeDescs2[0].set_vertex(0, hex->vertex(0));
		edgeDescs2[0].set_vertex(1, hex->vertex(3));
		length2 += fabs(axPos[0] - axPos[3]);
		edgeDescs2[1].set_vertex(0, hex->vertex(1));
		edgeDescs2[1].set_vertex(1, hex->vertex(2));
		length2 += fabs(axPos[1] - axPos[2]);
		edgeDescs2[2].set_vertex(0, hex->vertex(4));
		edgeDescs2[2].set_vertex(1, hex->vertex(7));
		length2 += fabs(axPos[4] - axPos[7]);
		edgeDescs2[3].set_vertex(0, hex->vertex(5));
		edgeDescs2[3].set_vertex(1, hex->vertex(6));
		length2 += fabs(axPos[5] - axPos[6]);

		number length3 = 0.0;
		edgeDescs3[0].set_vertex(0, hex->vertex(0));
		edgeDescs3[0].set_vertex(1, hex->vertex(4));
		length3 += fabs(axPos[0] - axPos[4]);
		edgeDescs3[1].set_vertex(0, hex->vertex(1));
		edgeDescs3[1].set_vertex(1, hex->vertex(5));
		length3 += fabs(axPos[1] - axPos[5]);
		edgeDescs3[2].set_vertex(0, hex->vertex(2));
		edgeDescs3[2].set_vertex(1, hex->vertex(6));
		length3 += fabs(axPos[2] - axPos[6]);
		edgeDescs3[3].set_vertex(0, hex->vertex(3));
		edgeDescs3[3].set_vertex(1, hex->vertex(7));
		length3 += fabs(axPos[3] - axPos[7]);

		std::vector<EdgeDescriptor>* edgeDescs = &edgeDescs3;
		if (length1 > length2)
		{
			if (length1 > length3)
				edgeDescs = &edgeDescs1;
		}
		else if (length2 > length3)
			edgeDescs = &edgeDescs2;

		mg->associated_elements(el, hex);
		const size_t elSz = el.size();
		for (size_t i = 0; i < 4; ++i)
		{
			for (size_t e = 0; e < elSz; ++e)
			{
				Edge* edge = el[e];
				if (CompareVertices(edge, &(*edgeDescs)[i]))
				{
					refiner->mark(edge, RM_REFINE);
					break;
				}
			}
		}
	}
}


#ifdef NC_WITH_PARMETIS
void NeuriteAxialRefinementMarker::unify
(
	MultiGrid* mg,
	int lvl,
	int localOffset,
	const Grid::AttachmentAccessor<Volume, AElemIndex>& aaElemInd, // local indices!
	const Grid::AttachmentAccessor<side_t, AElemIndices>& aaSideElemInd, // global indices!
	std::vector<std::pair<int, int> >& unificationPairs // global indices!
) const
{
	Grid::traits<Volume>::secure_container vl;
	Grid::traits<side_t>::secure_container sl;
	Grid::traits<Edge>::secure_container el;

	float axPos[8];
	std::vector<EdgeDescriptor> edgeDescs1(4);
	std::vector<EdgeDescriptor> edgeDescs2(4);
	std::vector<EdgeDescriptor> edgeDescs3(4);

	// mark BP volumes
	mg->begin_marking();
	mark_bp_volumes(mg, lvl);

	typedef geometry_traits<Volume>::const_iterator const_vol_it;
	const_vol_it itVol = mg->begin<Volume>(lvl);
	const_vol_it itVolEnd = mg->end<Volume>(lvl);
	for (; itVol != itVolEnd; ++itVol)
	{
		Volume* vol = *itVol;

		// exclude ghosts
		if (aaElemInd[vol] == -1)
			continue;

		Hexahedron* hex = dynamic_cast<Hexahedron*>(vol);
		UG_COND_THROW(!hex, "Found volume that is not a hexahedron.\n"
			"This implementation can only handle hexahedron grids.");

		std::set<Volume*> volsToUnify;

		// case 1: BP volume
		if (mg->is_marked(vol))
		{
			// Is this a central BP volume?
			if (!is_central_bp_vol(vol))
				continue;

			// Is this a neurite without ER?
			const bool withER = m_aaSurfParams[vol->vertex(0)].radial < 1.0;
			if (!withER)
				break;  // no need to unify anything anywhere

			// If it is, we have to find the surrounding volumes
			// and find all pairs connected in the dual graph.
			for (size_t i = 0; i < 8; ++i)
			{
				mg->associated_elements(vl, vol->vertex(i));
				const size_t vlSz = vl.size();
				for (size_t v = 0; v < vlSz; ++v)
					volsToUnify.insert(vl[v]);
			}
		}

		// case 2: neurite volume
		else
		{
			// we only look for the central volume
			mg->associated_elements(sl, vol);
			const size_t slSz = sl.size();
			if (slSz != 6)
				continue;

			// now identify the four lateral faces //

			// axial pos is only valid of neurite ID is that of the parent for BP vertices
			// otherwise take 0.0

			uint32_t nid = 0;
			for (size_t i = 0; i < 8; ++i)
				nid = std::max(nid, m_aaSurfParams[hex->vertex(i)].neuriteID & ((1 << 20) - 1));

			for (size_t i = 0; i < 8; ++i)
			{
				if ((m_aaSurfParams[vol->vertex(i)].neuriteID & ((1 << 20) - 1)) == nid)
					axPos[i] = m_aaSurfParams[hex->vertex(i)].axial;
				else
					axPos[i] = 0.0;
			}

			number length1 = 0.0;
			edgeDescs1[0].set_vertex(0, hex->vertex(0));
			edgeDescs1[0].set_vertex(1, hex->vertex(1));
			length1 += fabs(axPos[0] - axPos[1]);
			edgeDescs1[1].set_vertex(0, hex->vertex(2));
			edgeDescs1[1].set_vertex(1, hex->vertex(3));
			length1 += fabs(axPos[2] - axPos[3]);
			edgeDescs1[2].set_vertex(0, hex->vertex(4));
			edgeDescs1[2].set_vertex(1, hex->vertex(5));
			length1 += fabs(axPos[4] - axPos[5]);
			edgeDescs1[3].set_vertex(0, hex->vertex(6));
			edgeDescs1[3].set_vertex(1, hex->vertex(7));
			length1 += fabs(axPos[6] - axPos[7]);

			number length2 = 0.0;
			edgeDescs2[0].set_vertex(0, hex->vertex(0));
			edgeDescs2[0].set_vertex(1, hex->vertex(3));
			length2 += fabs(axPos[0] - axPos[3]);
			edgeDescs2[1].set_vertex(0, hex->vertex(1));
			edgeDescs2[1].set_vertex(1, hex->vertex(2));
			length2 += fabs(axPos[1] - axPos[2]);
			edgeDescs2[2].set_vertex(0, hex->vertex(4));
			edgeDescs2[2].set_vertex(1, hex->vertex(7));
			length2 += fabs(axPos[4] - axPos[7]);
			edgeDescs2[3].set_vertex(0, hex->vertex(5));
			edgeDescs2[3].set_vertex(1, hex->vertex(6));
			length2 += fabs(axPos[5] - axPos[6]);

			number length3 = 0.0;
			edgeDescs3[0].set_vertex(0, hex->vertex(0));
			edgeDescs3[0].set_vertex(1, hex->vertex(4));
			length3 += fabs(axPos[0] - axPos[4]);
			edgeDescs3[1].set_vertex(0, hex->vertex(1));
			edgeDescs3[1].set_vertex(1, hex->vertex(5));
			length3 += fabs(axPos[1] - axPos[5]);
			edgeDescs3[2].set_vertex(0, hex->vertex(2));
			edgeDescs3[2].set_vertex(1, hex->vertex(6));
			length3 += fabs(axPos[2] - axPos[6]);
			edgeDescs3[3].set_vertex(0, hex->vertex(3));
			edgeDescs3[3].set_vertex(1, hex->vertex(7));
			length3 += fabs(axPos[3] - axPos[7]);

			std::vector<EdgeDescriptor>* edgeDescs = &edgeDescs3;
			if (length1 > length2)
			{
				if (length1 > length3)
					edgeDescs = &edgeDescs1;
			}
			else if (length2 > length3)
				edgeDescs = &edgeDescs2;

			mg->associated_elements(el, hex);
			const size_t elSz = el.size();
			for (size_t i = 0; i < 4; ++i)
			{
				for (size_t e = 0; e < elSz; ++e)
				{
					Edge* edge = el[e];
					if (CompareVertices(edge, &(*edgeDescs)[i]))
					{
						mg->associated_elements(vl, edge);
						const size_t vlSz = vl.size();
						for (size_t v = 0; v < vlSz; ++v)
							volsToUnify.insert(vl[v]);
					}
				}
			}
		}


		// loop all volumes to unify
		std::set<Volume*>::const_iterator it = volsToUnify.begin();
		std::set<Volume*>::const_iterator itEnd = volsToUnify.end();
		for (; it != itEnd; ++it)
		{
			Volume* uVol = *it;
			int ind = aaElemInd[uVol] + localOffset;

			// get sides
			mg->associated_elements(sl, uVol);
			const size_t slSz = sl.size();
			for (size_t s = 0; s < slSz; ++s)
			{
				side_t* side = sl[s];

				// check whether connected volume is also BP and add to unification pairs if so
				Volume* connVol = GetConnectedNeighbor(*mg, side, uVol);
				if (connVol && volsToUnify.find(connVol) != volsToUnify.end())
				{
					int connInd = aaElemInd[connVol] + localOffset;
					unificationPairs.push_back(std::make_pair(ind, connInd));
				}
			}
		}
	}

	mg->end_marking();
}

#endif


void NeuriteAxialRefinementMarker::find_bp_volumes()
{
	// get surface view and grid
	SurfaceView sv(m_spDom->subset_handler());
	SmartPtr<MultiGrid> mg = m_spDom->grid();

	Grid::traits<Volume>::secure_container vl;

	// mark bp edges
	typedef SurfaceView::traits<Volume>::const_iterator const_vol_it;
	const_vol_it itVol = sv.begin<Volume>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_vol_it itVolEnd = sv.end<Volume>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; itVol != itVolEnd; ++itVol)
	{
		Volume* vol = *itVol;

		// Is this a central BP volume?
		if (!is_central_bp_vol(vol))
			continue;

		// mark volume
		m_aaBP[vol] = true;

		// Is this a neurite without ER?
		const bool withER = m_aaSurfParams[vol->vertex(0)].radial < 1.0;
		if (!withER)
			continue;

		// If it is, we have to find the surrounding volumes and mark them as well.
		const size_t nVrt = vol->num_vertices();
		for (size_t i = 0; i < nVrt; ++i)
		{
			mg->associated_elements(vl, vol->vertex(i));
			const size_t vlSz = vl.size();
			for (size_t v = 0; v < vlSz; ++v)
				m_aaBP[vl[v]] = true;
		}
	}
}


void NeuriteAxialRefinementMarker::mark_bp_volumes(MultiGrid* mg, int lvl) const
{
	Grid::traits<Volume>::secure_container vl;

	// mark bp edges
	typedef geometry_traits<Volume>::iterator const_vol_it;
	const_vol_it itVol = mg->begin<Volume>(lvl);
	const_vol_it itVolEnd = mg->end<Volume>(lvl);
	for (; itVol != itVolEnd; ++itVol)
	{
		Volume* vol = *itVol;

		// Is this a central BP volume?
		if (!is_central_bp_vol(vol))
			continue;

		// mark volume
		mg->mark(vol);

		// Is this a neurite without ER?
		const bool withER = m_aaSurfParams[vol->vertex(0)].radial < 1.0;
		if (!withER)
			continue;

		// If it is, we have to find the surrounding volumes and mark them as well.
		const size_t nVrt = vol->num_vertices();
		for (size_t i = 0; i < nVrt; ++i)
		{
			mg->associated_elements(vl, vol->vertex(i));
			const size_t vlSz = vl.size();
			for (size_t v = 0; v < vlSz; ++v)
				mg->mark(vl[v]);
		}
	}
}




bool NeuriteAxialRefinementMarker::is_bp_volume(Volume* vol) const
{
	return m_aaBP[vol];
}


bool NeuriteAxialRefinementMarker::is_central_bp_vol(Volume* vol) const
{
	const size_t nVrt = vol->num_vertices();
	UG_COND_THROW(nVrt != 8, "Found volume with " << nVrt << " vertices.\n"
		"This implementation can only handle hexahedron grids.");

	// Is this a central BP volume?
	// This is the case if and only if four vertices belong to one neurite (the parent)
	// and the other four belong to the same parent, but also a child.
	uint32_t nid[2] = {0, 0};
	size_t cnt[2] = {0, 0};
	for (size_t i = 0; i < nVrt; ++i)
	{
		uint32_t thisNid = m_aaSurfParams[vol->vertex(i)].neuriteID;
		if (!cnt[0])
		{
			nid[0] = thisNid;
			++cnt[0];
		}
		else if (thisNid == nid[0])
		{
			++cnt[0];
		}
		else if (!cnt[1])
		{
			nid[1] = thisNid;
			++cnt[1];
		}
		else if (thisNid == nid[1])
		{
			++cnt[1];
		}
	}

	if (cnt[0] != 4 || cnt[1] != 4)
		return false;
	if ((nid[0] & ((1 << 20) - 1)) != (nid[1] & ((1 << 20) - 1)))
		return false;
	if ((nid[0] < (1 << 20) && nid[1] < (1 << 20)) || (nid[0] >= (1 << 20) && nid[1] >= (1 << 20)))
		return false;

	return true;
}


}  // namespace neuro_collection
}  // namespace ug
