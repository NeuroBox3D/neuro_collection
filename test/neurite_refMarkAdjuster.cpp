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

#include "neurite_refMarkAdjuster.h"


namespace ug {
namespace neuro_collection {


static void mark_face_copy(IRefiner& ref, Grid& grid, Face* f)
{
    ref.mark(f, RM_COPY);

    MultiGrid::traits<Edge>::secure_container el;
    grid.associated_elements(el, f);

    size_t el_sz = el.size();
    for (size_t k = 0; k < el_sz; ++k)
       ref.mark(el[k], RM_COPY);
}

static void mark_vol_copy(IRefiner& ref, Grid& grid, Volume* v)
{
    ref.mark(v, RM_COPY);

    MultiGrid::traits<Face>::secure_container fl;
    grid.associated_elements(fl, v);

    size_t fl_sz = fl.size();
    for (size_t k = 0; k < fl_sz; ++k)
       ref.mark(fl[k], RM_COPY);
}

void NeuriteRefMarkAdjuster::change_edge_mark(IRefiner& ref, Edge* e, const vector3* dir) const
{
    vector3 neuriteDir;
    if (!dir)
    {
        try {m_spNP->direction_at_grid_object(neuriteDir, e);}
        catch (UGError& err)
        {
            // no direction could be calculated since edge axial param is out of bounds,
            // i.e. <0 (soma) or >1 (tip)
            // we do not refine here
            ref.mark(e, RM_COPY);
            return;
        }
        VecNormalize(neuriteDir, neuriteDir);
    }
    else neuriteDir = *dir;

    vector3 edgeDir;
    VecSubtract(edgeDir, m_aaPos[e->vertex(0)], m_aaPos[e->vertex(1)]);
    VecNormalize(edgeDir, edgeDir);

    // mark for refinement if angle between edge and neurite is less than 45 degrees
    // else mark anisotropic
    number scp = VecProd(edgeDir, neuriteDir);
    if (fabs(scp) > 0.707)
        ref.mark(e, RM_REFINE);
    else
        ref.mark(e, RM_ANISOTROPIC);
}


void NeuriteRefMarkAdjuster::change_face_mark(IRefiner& ref, Grid& grid, Face* f, const vector3* dir) const
{
    MultiGrid::traits<Edge>::secure_container el;
    grid.associated_elements(el, f);
    size_t el_sz = el.size();

    vector3 neuriteDir;
    if (!dir)
    {
        try {m_spNP->direction_at_grid_object(neuriteDir, f);}
        catch (UGError& err)
        {
            // no direction could be calculated since edge axial param is out of bounds,
            // i.e. <0 (soma) or >1 (tip)
            // we do not refine here
            ref.mark(f, RM_COPY);

            for (size_t k = 0; k < el_sz; ++k)
                ref.mark(el[k], RM_COPY);

            return;
        }
        VecNormalize(neuriteDir, neuriteDir);
    }
    else neuriteDir = *dir;

    // mark the face for anisotropic refinement
    ref.mark(f, RM_ANISOTROPIC);

    // mark the face edges appropriately
    for (size_t k = 0; k < el_sz; ++k)
        change_edge_mark(ref, el[k], &neuriteDir);
}


void NeuriteRefMarkAdjuster::change_vol_mark(IRefiner& ref, Grid& grid, Volume* v) const
{
    // mark the face edges appropriately
    MultiGrid::traits<Face>::secure_container fl;
    grid.associated_elements(fl, v);
    size_t fl_sz = fl.size();

    // get neurite direction
    vector3 neuriteDir;
    try {m_spNP->direction_at_grid_object(neuriteDir, v);}
    catch (UGError& err)
    {
        // no direction could be calculated since volume axial param is out of bounds,
        // i.e. <0 (soma) or >1 (tip)
        // we do not refine here
        ref.mark(v, RM_COPY);
        for (size_t k = 0; k < fl_sz; ++k)
            mark_face_copy(ref, grid, fl[k]);

        return;
    }

    // mark the face for anisotropic refinement
    ref.mark(v, RM_ANISOTROPIC);

    for (size_t k = 0; k < fl_sz; ++k)
        change_face_mark(ref, grid, fl[k], &neuriteDir);
}


void NeuriteRefMarkAdjuster::ref_marks_changed
(
	IRefiner& ref,
	const std::vector<Vertex*>& vrts,
	const std::vector<Edge*>& edges,
	const std::vector<Face*>& faces,
	const std::vector<Volume*>& vols
)
{
	UG_COND_THROW(!m_ssh.valid(), "No SubsetHandler set in InterfaceRefMarkAdjuster.");

	if (!ref.grid()) return;
	Grid& grid = *ref.grid();

	Grid::edge_traits::secure_container		assEdges;
	Grid::face_traits::secure_container		assFaces;
	Grid::volume_traits::secure_container 	assVols;

	// TODO: treat somata

#if 0
	// EDGES //
	size_t sz = edges.size();
	for (size_t i = 0; i < sz; ++i)
        change_edge_mark(ref, edges[i]);
#endif

	// FACES //
	size_t sz = faces.size();
	for (size_t i = 0; i < sz; ++i)
	{
	    // only refine quads (if anisotropic)
	    Quadrilateral* q = dynamic_cast<Quadrilateral*>(faces[i]);
	    if (!q)
	    {
	        mark_face_copy(ref, grid, faces[i]);
	        continue;
	    }

	    // calculate anisotropy of face
	    vector3 sidevec;
	    VecSubtract(sidevec, m_aaPos[q->vertex(1)], m_aaPos[q->vertex(0)]);
        number a = VecLength(sidevec);
        VecSubtract(sidevec, m_aaPos[q->vertex(2)], m_aaPos[q->vertex(1)]);
        number aniso = a / VecLength(sidevec);
        if (aniso < 1.0) aniso = 1.0 / aniso;

        // isotropic face: only copy, also mark edges copy
        if (aniso < 1.4142)
            mark_face_copy(ref, grid, q);
        // anisotropic face: mark anisotropic and process edges
        else
            change_face_mark(ref, grid, q);
	}

	// VOLUMES //
	sz = vols.size();
	for (size_t i = 0; i < sz; ++i)
	{
        // only refine hexahedra (if anisotropic)
        Hexahedron* h = dynamic_cast<Hexahedron*>(vols[i]);
        if (!h)
        {
            mark_vol_copy(ref, grid, vols[i]);
            continue;
        }

		// calculate anisotropy of volume
        vector3 sidevec;
        VecSubtract(sidevec, m_aaPos[h->vertex(1)], m_aaPos[h->vertex(0)]);
        number a = VecLength(sidevec);
        VecSubtract(sidevec, m_aaPos[h->vertex(2)], m_aaPos[h->vertex(1)]);
        number b = VecLength(sidevec);
        VecSubtract(sidevec, m_aaPos[h->vertex(4)], m_aaPos[h->vertex(0)]);
        number c = VecLength(sidevec);
        number aniso = std::max(a, std::max(b,c)) / std::min(a, std::min(b,c));

        // isotropic volume: only copy, also mark faces and edges copy
        if (aniso < 1.4142)
            mark_vol_copy(ref, grid, h);
        // anisotropic volume: mark volume and volume faces appropriately
        else
            change_vol_mark(ref, grid, h);
	}
}



void add_neurite_ref_mark_adjuster(IRefiner* ref, SmartPtr<NeuriteRefMarkAdjuster> nrma)
{
	HangingNodeRefiner_MultiGrid* href = dynamic_cast<HangingNodeRefiner_MultiGrid*>(ref);
	UG_COND_THROW(!href, "An interface refinement mark adjuster can only be added to an instance of HangingNodeRefiner_MultiGrid.");

	href->add_ref_mark_adjuster(nrma);
}




} // namespace neuro_collection
} // namespace ug
