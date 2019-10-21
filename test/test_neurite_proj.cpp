/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-12-27
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

/// plugin and configuration
#include "test_neurite_proj.h"
#include "neurite_refMarkAdjuster.h"
#include "neurite_util.h"
#include "neurite_grid_generation.h"
#include "../util/misc_util.h"

/// ug
#include "lib_grid/refinement/projectors/projection_handler.h" // ProjectionHandler
#include "lib_grid/refinement/projectors/neurite_projector.h" // NeuriteProjector
#include "lib_grid/file_io/file_io_ugx.h"  // GridWriterUGX
#include "lib_grid/file_io/file_io.h"  // SaveGridHierarchyTransformed
#include "lib_grid/grid/geometry.h" // MakeGeometry3d
#include "lib_grid/refinement/global_multi_grid_refiner.h" // GlobalMultigridRefiner
#include "lib_grid/global_attachments.h" // GlobalAttachments
#include "lib_grid/algorithms/extrusion/extrusion.h" // Extrude
#include "lib_grid/refinement/regular_refinement.h"  // Refine
#include "lib_grid/algorithms/geom_obj_util/face_util.h" // CalculateNormal
#include "lib_grid/algorithms/element_side_util.h" // GetOpposingSide
#include "lib_grid/algorithms/grid_generation/icosahedron.h" // Icosahedron
#include "lib_grid/algorithms/smoothing/manifold_smoothing.h" // TangentialSmoothing
#include "lib_grid/algorithms/remeshing/resolve_intersections.h" // ResolveTriangleIntersection
#include "lib_disc/function_spaces/error_elem_marking_strategy.h" // GlobalMarking
#include "lib_disc/domain_util.h"   // LoadDomain
#include "lib_disc/quadrature/gauss_legendre/gauss_legendre.h" // Gauss-Legendre
#include "lib_algebra/small_algebra/small_algebra.h" // Invert
#include "common/util/string_util.h"  // TrimString
#include "common/util/file_util.h"  // FindFileInStandardPaths

/// other
#include <boost/lexical_cast.hpp>

/// debug id
ug::DebugID NC_TNP("NC_DID.TNP");

namespace ug {
	namespace neuro_collection {

	////////////////////////////////////////////////////////////////////////
	/// import_swc
	////////////////////////////////////////////////////////////////////////
	void import_swc
	(
		const std::string& fileName,
		std::vector<SWCPoint>& vPointsOut,
		number scale=1.0
	)
	{
    vPointsOut.clear();

    std::string inFileName = FindFileInStandardPaths(fileName.c_str());

    std::ifstream inFile(inFileName.c_str());
    UG_COND_THROW(!inFile, "SWC input file '" << fileName << "' could not be opened for reading.");

    // read line by line
    std::string line;
    size_t lineCnt = 0;
    size_t curInd = 0;
    std::map<int, size_t> indexMap;
    while (std::getline(inFile, line))
    {
        ++lineCnt;

        // trim whitespace
        line = TrimString(line);

        // ignore anything from possible '#' onwards
        size_t nChar = line.size();
        for (size_t i = 0; i < nChar; ++i)
        {
            if (line.at(i) == '#')
            {
                line = line.substr(0, i);
                break;
            }
        }

        // empty lines can be ignored
        if (line.empty()) continue;

        // split the line into tokens
        std::istringstream buf(line);
        std::istream_iterator<std::string> beg(buf), end;
        std::vector<std::string> strs(beg, end);

        // assert number of tokens is correct
        UG_COND_THROW(strs.size() != 7, "Error reading SWC file '" << inFileName
            << "': Line " << lineCnt << " does not contain exactly 7 values.");

        // collect data for SWC point
        vPointsOut.resize(vPointsOut.size()+1);
        SWCPoint& pt = vPointsOut.back();

        // get index from file and map to our index
        indexMap[boost::lexical_cast<int>(strs[0])] = curInd;

        // type
        int type = boost::lexical_cast<int>(strs[1]);
        switch (type)
        {
        	case 0: pt.type = SWC_UNDF; break;
            case 1: pt.type = SWC_SOMA; break;
            case 2: pt.type = SWC_AXON; break;
            case 3: pt.type = SWC_DEND; break;
            case 4: pt.type = SWC_APIC; break;
            case 5: pt.type = SWC_FORK; break;
            case 6: pt.type = SWC_END; break;
            default: pt.type = SWC_CUSTOM;
        }

        // coordinates
        pt.coords.x() = boost::lexical_cast<number>(strs[2]) * scale;
        pt.coords.y() = boost::lexical_cast<number>(strs[3]) * scale;
        pt.coords.z() = boost::lexical_cast<number>(strs[4]) * scale;

        // radius
        pt.radius = boost::lexical_cast<number>(strs[5]) * scale;

        // connections
        int conn = boost::lexical_cast<int>(strs[6]);
        if (conn >= 0)
        {
            std::map<int, size_t>::const_iterator it = indexMap.find(conn);
            UG_COND_THROW(it == indexMap.end(), "Error reading SWC file '" << inFileName
                << "': Line " << lineCnt << " refers to unknown parent index " << conn << ".");

            size_t parentID = indexMap[conn];
            pt.conns.push_back(parentID);
            vPointsOut[parentID].conns.push_back(curInd);
        }

        // increase current point index
        ++curInd;
    	}
	}

	////////////////////////////////////////////////////////////////////////
	/// import_swc
	////////////////////////////////////////////////////////////////////////
	void import_swc_old
	(
    const std::string& fileName,
    std::vector<SWCPoint>& vPointsOut,
    bool correct,
    number scale)
	{
    vPointsOut.clear();

    std::ifstream inFile(fileName.c_str());
    UG_COND_THROW(!inFile, "SWC input file '" << fileName << "' could not be opened for reading.");

    // read line by line
    std::string line;
    size_t lineCnt = 0;
    size_t curInd = 0;
    std::map<int, size_t> indexMap;
    while (std::getline(inFile, line))
    {
        ++lineCnt;

        // trim whitespace
        line = TrimString(line);

        // ignore anything from possible '#' onwards
        size_t nChar = line.size();
        for (size_t i = 0; i < nChar; ++i)
        {
            if (line.at(i) == '#')
            {
                line = line.substr(0, i);
                break;
            }
        }

        // empty lines can be ignored
        if (line.empty()) continue;

        // split the line into tokens
        std::istringstream buf(line);
        std::istream_iterator<std::string> beg(buf), end;
        std::vector<std::string> strs(beg, end);

        // assert number of tokens is correct
        UG_COND_THROW(strs.size() != 7, "Error reading SWC file '" << fileName
            << "': Line " << lineCnt << " does not contain exactly 7 values.");

        // collect data for SWC point
        vPointsOut.resize(vPointsOut.size()+1);
        SWCPoint& pt = vPointsOut.back();

        // get index from file and map to our index
        indexMap[boost::lexical_cast<int>(strs[0])] = curInd;

        // type
        int type = boost::lexical_cast<int>(strs[1]);
        switch (type)
        {
        	case 0: pt.type = SWC_UNDF; break;
            case 1: pt.type = SWC_SOMA; break;
            case 2: pt.type = SWC_AXON; break;
            case 3: pt.type = SWC_DEND; break;
            case 4: pt.type = SWC_APIC; break;
            case 5: pt.type = SWC_FORK; break;
            case 6: pt.type = SWC_END; break;
            default: pt.type = SWC_CUSTOM;
        }

        // coordinates
        pt.coords.x() = boost::lexical_cast<number>(strs[2]);
        pt.coords.y() = boost::lexical_cast<number>(strs[3]);
        pt.coords.z() = boost::lexical_cast<number>(strs[4]);

        // radius
        pt.radius = boost::lexical_cast<number>(strs[5]) * scale;

        // connections
        int conn = boost::lexical_cast<int>(strs[6]);
        if (conn >= 0)
        {
            std::map<int, size_t>::const_iterator it = indexMap.find(conn);
            UG_COND_THROW(it == indexMap.end(), "Error reading SWC file '" << fileName
                << "': Line " << lineCnt << " refers to unknown parent index " << conn << ".");

            size_t parentID = indexMap[conn];
            pt.conns.push_back(parentID);
            vPointsOut[parentID].conns.push_back(curInd);
        }

        // increase current point index
        ++curInd;
    }

    if (correct) {
    	for (size_t conn = 0; conn < curInd; conn++) {
      		size_t parentId = indexMap[conn];
      		swc_type type = vPointsOut[parentId].type;
       		if ((type != SWC_SOMA) && (type != SWC_UNDF)) {
       			if (vPointsOut[parentId].conns.size() == 3) { /// branch
       				std::cout << "Correcting branch no: " << conn << std::endl;
       				ug::vector3& p1 = vPointsOut[vPointsOut[parentId].conns[0]].coords;
       				ug::vector3& p2 = vPointsOut[vPointsOut[parentId].conns[1]].coords;
       				vPointsOut[parentId].coords = ug::vector3(p1[0]/2+p2[0]/2, p2[1]/2+p1[1]/2, p1[2]/2+p2[2]/2);
       			} else {
       				UG_THROW("More than two branches detected. Current implementation does not support this.")
       			}
       		}
       	}
    }
}

////////////////////////////////////////////////////////////////////////
/// smoothing
////////////////////////////////////////////////////////////////////////
void smoothing(std::vector<SWCPoint>& vPointsInOut, size_t n, number h, number gamma)
{
	// find neurite root vertices
	const size_t nP = vPointsInOut.size();
	std::vector<size_t> rootVrts;
	std::vector<bool> treated(nP, false);
	for (size_t i = 0; i < nP; ++i)
	{
		if (treated[i]) continue;
		treated[i] = true;
		if (vPointsInOut[i].type != SWC_SOMA) continue;

		// here, we have a soma point;
		// find first non-soma point in all directions
		std::queue<size_t> q;
		q.push(i);

		while (!q.empty())
		{
			const size_t ind = q.front();
			const SWCPoint& pt = vPointsInOut[ind];
			q.pop();

			if (pt.type == SWC_SOMA)
			{
				const size_t nConn = pt.conns.size();
				for (size_t j = 0; j < nConn; ++j)
					if (!treated[pt.conns[j]])
						q.push(pt.conns[j]);
			}
			else
				rootVrts.push_back(ind);

			treated[ind] = true;
		}
	}

	// starting at root vertices, smooth the entire tree(s),
	// but leave out soma vertices as well as branching points
	std::vector<vector3> newPos(nP);
	for (size_t i = 0; i < n; ++i)
	{
		treated.clear();
		treated.resize(nP, false);

		std::stack<size_t> stack;
		for (size_t rv = 0; rv < rootVrts.size(); ++rv)
			stack.push(rootVrts[rv]);

		while (!stack.empty())
		{
			size_t ind = stack.top();
			stack.pop();
			const SWCPoint& pt = vPointsInOut[ind];
			const vector3& x = pt.coords;
			vector3& x_new = newPos[ind];

			UG_COND_THROW(treated[ind], "Circle detected in supposedly tree-shaped neuron!\n"
				"Position: " << vPointsInOut[ind].coords);
			treated[ind] = true;

			// somata are not smoothed and not iterated over
			if (pt.type == SWC_SOMA)
			{
				x_new = x;
				continue;
			}

			// branching points are not smoothed, but iterated over
			size_t connSz = pt.conns.size();
			for (size_t c = 0; c < connSz; ++c)
				if (!treated[pt.conns[c]])
					stack.push(pt.conns[c]);

			if (connSz != 2)
			{
				x_new = x;
				continue;
			}

			// here we have a non-branching, non-end, non-soma point: smooth
			const vector3& x1 = vPointsInOut[pt.conns[0]].coords;
			const vector3& x2 = vPointsInOut[pt.conns[1]].coords;

			number d1 = VecDistanceSq(x1, x);
			number d2 = VecDistanceSq(x2, x);
			number w1 = std::exp(-d1/(h*h));
			number w2 = std::exp(-d2/(h*h));

			// only really smooth if both adjacent edges are short
			number w = std::min(w1, w2);

			// correction
			vector3 corr;
			VecScaleAdd(corr, w, x1, -2*w, x, w, x2);
			VecScale(corr, corr, 1.0 / (1.0 + 2*w));

			// take only the part orthogonal to x1 - x2,
			// we do not want to shift x towards the nearer neighbor
			VecSubtract(x_new, x1, x2); // using x_new as intermediate variable
			number normSq = VecNormSquared(x_new);
			VecScaleAdd(corr, 1.0, corr, - VecProd(corr, x_new) / normSq, x_new);
			VecScaleAdd(x_new, 1.0, x, gamma, corr);
		}

		// assign new positions
		for (size_t p = 0; p < nP; ++p)
			if (treated[p]) // soma points may not have been treated
				vPointsInOut[p].coords = newPos[p];
	}

}

////////////////////////////////////////////////////////////////////////
/// collapse_short_edges
////////////////////////////////////////////////////////////////////////
void collapse_short_edges(Grid& g, SubsetHandler& sh)
{
	// get access to positions
	UG_COND_THROW(!g.has_vertex_attachment(aPosition), "Position attachment not attached to grid.")
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	// get access to diameter attachment
	ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter");
	UG_COND_THROW(!g.has_vertex_attachment(aDiam), "No diameter attachment attached to grid.");
	Grid::AttachmentAccessor<Vertex, ANumber> aaDiam(g, aDiam);

	// a short edge is one that is shorter than its diameter
	// sort all short edges in a priority queue
	std::priority_queue<std::pair<Edge*, number>, std::vector<std::pair<Edge*, number> >, EdgeLengthCompare> pq;
	EdgeIterator eit = g.begin<Edge>();
	EdgeIterator edge_end = g.end<Edge>();
	for (; eit != edge_end; ++eit)
	{
		Edge* e = *eit;
		number length = EdgeLengthSq(e, aaPos);
		number diam = std::max(aaDiam[e->vertex(0)], aaDiam[e->vertex(1)]);
		diam = diam*diam;
		if (length < diam)
			pq.push(std::make_pair(e, length));
	}

	while (!pq.empty())
	{
		std::pair<Edge*, number> elp = pq.top();
		pq.pop();

		// edge length might not be up to date; if so: re-insert with correct length (if necessary)
		Edge* curEdge = elp.first;
		number curLen = EdgeLengthSq(curEdge, aaPos);
		if (curLen != elp.second)
		{
			number curDiam = std::max(aaDiam[curEdge->vertex(0)], aaDiam[curEdge->vertex(1)]);
			curDiam = curDiam*curDiam;
			if (curLen < curDiam)
			{
				elp.second = curLen;
				pq.push(elp);
			}
			continue;
		}

		// do not consider short edges connecting two branching points
		Vertex* v1 = curEdge->vertex(0);
		Vertex* v2 = curEdge->vertex(1);
		Grid::AssociatedEdgeIterator it = g.associated_edges_begin(v1);
		Grid::AssociatedEdgeIterator it_end = g.associated_edges_end(v1);
		size_t nAssV1 = 0;
		for (; it != it_end; ++it)
			++nAssV1;

		it = g.associated_edges_begin(v2);
		it_end = g.associated_edges_end(v2);
		size_t nAssV2 = 0;
		for (; it != it_end; ++it)
			++nAssV2;

		if (nAssV1 > 2 && nAssV2 > 2)
			continue;

		// otherwise, collapse edge

		// (a) calculate position (and radius) for new vertex --
		// if the old edge was parallel to one of the adjacent edges,
		// then add the complete edge length to that adjacent edge
		number newDiam;
		vector3 newPos;
		vector3 x1 = aaPos[v1];
		vector3 x2 = aaPos[v2];
		vector3 d0, d1, d2;
		VecSubtract(d0, x2, x1);
		VecNormalize(d0, d0);

		// never move branching points
		if (nAssV1 > 2)
		{
			newDiam = aaDiam[v1];
			newPos = x1;
		}
		else if (nAssV2 > 2)
		{
			newDiam = aaDiam[v2];
			newPos = x2;
		}

		// if the edge is a terminal edge, set new vertex at previous terminal vertex
		else if (nAssV1 == 1)
		{
			newDiam = aaDiam[v1];
			newPos = x1;
		}
		else if (nAssV2 == 1)
		{
			newDiam = aaDiam[v2];
			newPos = x2;
		}
		else
		{
			// calculate directions of adjacent edges
			it = g.associated_edges_begin(v1);
			if (*it != curEdge)
				VecSubtract(d1, x1, aaPos[GetOpposingSide(g, *it, v1)]);
			else
				VecSubtract(d1, x1, aaPos[GetOpposingSide(g, *(++it), v1)]);

			it = g.associated_edges_begin(v2);
			if (*it != curEdge)
				VecSubtract(d2, aaPos[GetOpposingSide(g, *it, v2)], x2);
			else
				VecSubtract(d2, aaPos[GetOpposingSide(g, *(++it), v2)], x2);

			VecNormalize(d1, d1);
			VecNormalize(d2, d2);
			number w1 = 1.0 - fabs(VecProd(d0, d1));
			number w2 = 1.0 - fabs(VecProd(d0, d2));

			// if all three directions are practically co-linear, choose middle
			if (w1 < 0.05 && w2 < 0.05)	// corresponds to a deviation of about 18 degrees
			{
				newDiam = 0.5 * (aaDiam[v1] + aaDiam[v2]);
				VecScaleAdd(newPos, 0.5, x1, 0.5, x2);
			}
			// otherwise, weighted sum
			else
			{
				newDiam = w1 * aaDiam[v1] + w2 * aaDiam[v2];
				newDiam /= (w1 + w2);
				VecScaleAdd(newPos, w1, x1, w2, x2);
				VecScale(newPos, newPos, 1.0 / (w1 + w2));
			}
		}

		// (b) actual collapse
		Vertex* newVrt = *g.create<RegularVertex>();
		sh.assign_subset(newVrt, sh.get_subset_index(curEdge));
		CollapseEdge(g, curEdge, newVrt);

		// (c) assign the new vertex its position and diameter
		aaPos[newVrt] = newPos;
		aaDiam[newVrt] = newDiam;
	}
}


////////////////////////////////////////////////////////////////////////
/// convert_pointlist_to_neuritelist
////////////////////////////////////////////////////////////////////////
void convert_pointlist_to_neuritelist
(
    const std::vector<SWCPoint>& vPoints,
    std::vector<SWCPoint>& vSomaPoints,
    std::vector<std::vector<vector3> >& vPosOut,
    std::vector<std::vector<number> >& vRadOut,
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > >& vBPInfoOut,
    std::vector<size_t>& vRootNeuriteIndsOut
)
{
    // clear out vectors
    vPosOut.clear();
    vRadOut.clear();
    vBPInfoOut.clear();
    vRootNeuriteIndsOut.clear();

	size_t nPts = vPoints.size();
	std::vector<bool> ptProcessed(nPts, false);
	size_t nProcessed = 0;
	size_t curNeuriteInd = 0;

	while (nProcessed != nPts)
	{
		// find first soma's root point in geometry and save its index as i
		size_t i = 0;
		for (; i < nPts; ++i)
		{
			if (vPoints[i].type == SWC_SOMA && !ptProcessed[i])
				break;
		}
		UG_COND_THROW(i == nPts, "No soma contained in (non-empty) list of unprocessed SWC points, \n"
				"i.e., there is at least one SWC point not connected to any soma.");
		vSomaPoints.push_back(vPoints[i]);

		// collect neurite root points
		std::vector<std::pair<size_t, size_t> > rootPts;
		std::queue<std::pair<size_t, size_t> > soma_queue;
		soma_queue.push(std::make_pair((size_t)-1, i));
		while (!soma_queue.empty())
		{
			size_t pind = soma_queue.front().first;
			size_t ind = soma_queue.front().second;
			soma_queue.pop();

			const SWCPoint& pt = vPoints[ind];

			if (pt.type == SWC_SOMA)
			{
				ptProcessed[ind] = true;
				++nProcessed;

				size_t nConn = pt.conns.size();
				for (size_t j = 0; j < nConn; ++j)
					if (pt.conns[j] != pind)
						soma_queue.push(std::make_pair(ind, pt.conns[j]));
			}
			else
				rootPts.push_back(std::make_pair(pind, ind));
		}

		vPosOut.resize(vPosOut.size() + rootPts.size());
		vRadOut.resize(vRadOut.size() + rootPts.size());
		vBPInfoOut.resize(vBPInfoOut.size() + rootPts.size());

		std::stack<std::pair<size_t, size_t> > processing_stack;
		for (size_t i = 0; i < rootPts.size(); ++i)
			processing_stack.push(rootPts[i]);

		vRootNeuriteIndsOut.push_back(curNeuriteInd);

		// helper map to be used to correctly save BPs:
		// maps branch root parent ID to neurite ID and BP ID
		std::map<size_t, std::pair<size_t, size_t> > helperMap;

		while (!processing_stack.empty())
		{
			size_t pind = processing_stack.top().first;
			size_t ind = processing_stack.top().second;
			processing_stack.pop();

			ptProcessed[ind] = true;
			++nProcessed;

			const SWCPoint& pt = vPoints[ind];

			UG_COND_THROW(pt.type == SWC_SOMA, "Detected neuron with more than one soma.");

			// push back coords and radius information to proper neurite
			vPosOut[curNeuriteInd].push_back(pt.coords);
			vRadOut[curNeuriteInd].push_back(pt.radius);

			size_t nConn = pt.conns.size();

			// branching point
			if (nConn > 2)
			{
				// branch with minimal angle will continue current branch
				vector3 parentDir;
				VecSubtract(parentDir, pt.coords, vPoints[pind].coords);
				VecNormalize(parentDir, parentDir);

				size_t parentToBeDiscarded = 0;
				size_t minAngleInd = 0;
				number minAngle = std::numeric_limits<number>::infinity();

				for (size_t i = 0; i < nConn; ++i)
				{
					if (pt.conns[i] == pind)
					{
						parentToBeDiscarded = i;
						continue;
					}

					vector3 dir;
					VecSubtract(dir, vPoints[pt.conns[i]].coords, pt.coords);
					VecNormalize(dir, dir);

					number angle = acos(VecProd(dir, parentDir));
					if (angle < minAngle)
					{
						minAngle = angle;
						minAngleInd = i;
					}
				}

				// branching point info
				std::pair<size_t, std::vector<size_t> > bp;
				bp.first = vPosOut[curNeuriteInd].size()-1; // BP is located at this index of the neurite

				// resize out vectors to accommodate new neurites starting here
				size_t newSize = vPosOut.size() + nConn - 2;
				vPosOut.resize(newSize);
				vRadOut.resize(newSize);
				vBPInfoOut.resize(newSize);

				for (size_t i = 0; i < nConn; ++i)
				{
					if (i == parentToBeDiscarded || i == minAngleInd)
						continue;

					// push new neurite starting point index to stack
					processing_stack.push(std::make_pair(ind, pt.conns[i]));

					// save current point ID with current neurite ID and BP ID
					// for later assignment of branching neurite ID to BP
					// NOT: bp.second.push_back(++futureNeuriteInd);
					helperMap[ind] = std::make_pair(curNeuriteInd, vBPInfoOut[curNeuriteInd].size());
				}

				// push next index of the current neurite to stack
				processing_stack.push(std::make_pair(ind, pt.conns[minAngleInd]));

				// add BP to BP vector
				vBPInfoOut[curNeuriteInd].push_back(bp);
			}

			// end point
			else if (nConn == 1)
			{
				// if the stack is not empty, the next ID on it will start a new neurite
				if (!processing_stack.empty())
				{
					++curNeuriteInd;

					size_t nextParentID = processing_stack.top().first;

					// is this a root neurite? (this is the case if helper map does not contain next parent)
					std::map<size_t, std::pair<size_t, size_t> >::const_iterator it = helperMap.find(nextParentID);
					if (it != helperMap.end())
					{
						// push back parent position and radius to new neurite
						vPosOut[curNeuriteInd].push_back(vPoints[nextParentID].coords);
						vRadOut[curNeuriteInd].push_back(vPoints[nextParentID].radius);

						// push back new neurite ID to BP at parent
						size_t nextParentNeuriteID = helperMap[nextParentID].first;
						size_t nextParentBPID = helperMap[nextParentID].second;
						vBPInfoOut[nextParentNeuriteID][nextParentBPID].second.push_back(curNeuriteInd);
					}
					// else: the next point is the root point of a root neurite
					else
					{
						vRootNeuriteIndsOut.push_back(curNeuriteInd);
						// Identify soma/neurite BPs
						/*
						std::pair<size_t, std::vector<size_t> > bp;
						bp.first = 0; // first point is at soma
						vBPInfoOut[curNeuriteInd].push_back(bp);
						*/
					}
				}
			}

			// normal point
			else
			{
				for (size_t i = 0; i < nConn; ++i)
				{
					if (pt.conns[i] != pind)
						processing_stack.push(std::make_pair(ind, pt.conns[i]));
				}
			}
		}

		// last neurite of each neuron does not increase counter
		++curNeuriteInd;
	}

	#if 0
    	size_t numSomaPoints = vSomaPoints.size();
    	UG_DLOGN(NC_TNP, 0, "Number of soma points: " << numSomaPoints);
    	for (size_t i = 0; i < numSomaPoints; i++) {
    		UG_DLOGN(NC_TNP, 0, "Coordinates for soma: " << vSomaPoints[i].coords);
    	}
	#endif
}



////////////////////////////////////////////////////////////////////////
/// create_spline_data_for_neurites
////////////////////////////////////////////////////////////////////////
void create_spline_data_for_neurites
(
    std::vector<NeuriteProjector::Neurite>& vNeuritesOut,
    const std::vector<std::vector<vector3> >& vPos,
    const std::vector<std::vector<number> >& vR,
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > >* vBPInfo = NULL
)
{
    size_t nNeurites = vPos.size();
    vNeuritesOut.resize(nNeurites);

    // first: reserve memory for branching region vectors
    // (we will point to their elements in BranchingPoints and do not want the vectors to reallocate!)
    if (vBPInfo)
        for (size_t n = 0; n < nNeurites; ++n)
            vNeuritesOut[n].vBR.reserve((*vBPInfo)[n].size()+1);

    for (size_t n = 0; n < nNeurites; ++n)
    {
        NeuriteProjector::Neurite& neuriteOut = vNeuritesOut[n];
        const std::vector<vector3>& pos = vPos[n];
        std::vector<number> r = vR[n];
        std::vector<std::pair<size_t, std::vector<size_t> > >* bpInfo = NULL;

        if (vBPInfo)
             bpInfo = &((*vBPInfo)[n]);

        // parameterize to achieve constant velocity on piece-wise linear geom
        size_t nVrt = pos.size();
        std::vector<number> tSuppPos(nVrt);
        std::vector<number> dt(nVrt);
        number totalLength = 0.0;
        for (size_t i = 0; i < nVrt-1; ++i)
        {
            tSuppPos[i] = totalLength;
            totalLength += VecDistance(pos[i], pos[i+1]);
        }
        for (size_t i = 0; i < nVrt-1; ++i)
            tSuppPos[i] /= totalLength;
        tSuppPos[nVrt-1] = 1.0;

        for (size_t i = 0; i < nVrt-1; ++i)
            dt[i+1] = tSuppPos[i+1] - tSuppPos[i];


        // now calculate cubic splines in all dimensions and for radius
        // using moments (as in script: Wittum, Numerik 0)
        DenseMatrix<VariableArray2<number> > mat;
        DenseVector<VariableArray1<number> > x0, x1, x2, xr, rhs;
        mat.resize(nVrt, nVrt);
        x1.resize(nVrt);
        rhs.resize(nVrt);

        for (size_t i = 0; i < nVrt; ++i)
            mat(i,i) = 2.0;
        for (size_t i = 1; i < nVrt-1; ++i)
        {
            number h2 = tSuppPos[i+1] - tSuppPos[i-1];
            mat(i,i+1) = dt[i+1] / h2;
            mat(i,i-1) = dt[i] / h2;
        }
        UG_COND_THROW(!Invert(mat), "Failed to invert moment matrix for spline calculation.")

        for (size_t i = 1; i < nVrt-1; ++i)
            rhs[i] = 6.0 / (tSuppPos[i+1] - tSuppPos[i-1]) *
                     ((pos[i+1][0] - pos[i][0]) / dt[i+1]
                     - (pos[i][0] - pos[i-1][0]) / dt[i]) ;
        x0 = mat*rhs;

        for (size_t i = 1; i < nVrt-1; ++i)
            rhs[i] = 6.0 / (tSuppPos[i+1] - tSuppPos[i-1]) *
                     ((pos[i+1][1] - pos[i][1]) / dt[i+1]
                     - (pos[i][1] - pos[i-1][1]) / dt[i]) ;
        x1 = mat*rhs;

        for (size_t i = 1; i < nVrt-1; ++i)
            rhs[i] = 6.0 / (tSuppPos[i+1] - tSuppPos[i-1]) *
                     ((pos[i+1][2] - pos[i][2]) / dt[i+1]
                     - (pos[i][2] - pos[i-1][2]) / dt[i]) ;
        x2 = mat*rhs;

        for (size_t i = 1; i < nVrt-1; ++i)
            rhs[i] = 6.0 / (tSuppPos[i+1] - tSuppPos[i-1]) *
                     ((r[i+1] - r[i]) / dt[i+1]
                     - (r[i] - r[i-1]) / dt[i]) ;
        xr = mat*rhs;

        // FIXME: find suitable permissible render vector
        vector3 neuriteDir;
        VecSubtract(neuriteDir, pos[nVrt-1], pos[0]);
        VecNormalize(neuriteDir, neuriteDir);
        if (fabs(neuriteDir[0]) < fabs(neuriteDir[1]))
        {
        	if (fabs(neuriteDir[0]) < fabs(neuriteDir[2]))
            	neuriteOut.refDir = vector3(1,0,0);
        	else
            	neuriteOut.refDir = vector3(0,0,1);
        }
        else
        {
			if (fabs(neuriteDir[1]) < fabs(neuriteDir[2]))
				neuriteOut.refDir = vector3(0,1,0);
			else
				neuriteOut.refDir = vector3(0,0,1);
		}

        //neuriteOut.refDir = vector3(0,1/sqrt(2),1/sqrt(2));
		neuriteOut.vSec.reserve(nVrt-1);

        // this will be 0 for root branches and 1 otherwise
        size_t brInd = neuriteOut.vBR.size();
        std::vector<std::pair<size_t, std::vector<size_t> > >::const_iterator brIt;
        std::vector<std::pair<size_t, std::vector<size_t> > >::const_iterator brIt_end;
        if (bpInfo)
        {
            // vBR may already contain an initial BR; now resize to accommodate all others
            neuriteOut.vBR.resize(brInd + bpInfo->size());
            brIt = bpInfo->begin();
            brIt_end = bpInfo->end();
        }

        for (size_t i = 0; i < nVrt-1; ++i)
        {
            NeuriteProjector::Section sec(tSuppPos[i+1]);
            number* param = &sec.splineParamsX[0];
            param[0] = (x0[i]-x0[i+1]) / (6.0 * dt[i+1]);
            param[1] = 0.5 * x0[i+1];
            param[2] = -(dt[i+1]/6.0 * (x0[i] + 2.0*x0[i+1]) + (pos[i+1][0] - pos[i][0]) / dt[i+1]);
            param[3] = pos[i+1][0];
            param = &sec.splineParamsY[0];
            param[0] = (x1[i]-x1[i+1]) / (6.0 * dt[i+1]);
            param[1] = 0.5 * x1[i+1];
            param[2] = -(dt[i+1]/6.0 * (x1[i] + 2.0*x1[i+1]) + (pos[i+1][1] - pos[i][1]) / dt[i+1]);
            param[3] = pos[i+1][1];
            param = &sec.splineParamsZ[0];
            param[0] = (x2[i]-x2[i+1]) / (6.0 * dt[i+1]);
            param[1] = 0.5 * x2[i+1];
            param[2] = -(dt[i+1]/6.0 * (x2[i] + 2.0*x2[i+1]) + (pos[i+1][2] - pos[i][2]) / dt[i+1]);
            param[3] = pos[i+1][2];
            param = &sec.splineParamsR[0];
            param[0] = (xr[i]-xr[i+1]) / (6.0 * dt[i+1]);
            param[1] = 0.5 * xr[i+1];
            param[2] = -(dt[i+1]/6.0 * (xr[i] + 2.0*xr[i+1]) + (r[i+1] - r[i]) / dt[i+1]);
            param[3] = r[i+1];

            // branching points?
            if (bpInfo && brIt != brIt_end && brIt->first == i+1)
            {
                NeuriteProjector::BranchingRegion& br = neuriteOut.vBR[brInd];
                std::vector<size_t>::const_iterator itBranch = brIt->second.begin();
                std::vector<size_t>::const_iterator itBranch_end = brIt->second.end();

                // create BP for parent neurite's BR
                br.bp = make_sp(new NeuriteProjector::BranchingPoint());

                // register parent BR at BP
                br.bp->vNid.push_back(n);
                br.bp->vRegions.push_back(&br);

                br.t = tSuppPos[i+1];

                // loop all child neurites starting at this BP
                for (; itBranch != itBranch_end; ++itBranch)
                {
                    size_t childID = *itBranch;

                    // save pointer to BP at child neurite's BR
                    NeuriteProjector::BranchingRegion newChildBR;
                    newChildBR.bp = br.bp;
                    vNeuritesOut[childID].vBR.push_back(newChildBR);
                    NeuriteProjector::BranchingRegion& childBR = vNeuritesOut[childID].vBR[0];

                    // register child BR at BP
                    br.bp->vNid.push_back(childID);
                    br.bp->vRegions.push_back(&childBR);

                    childBR.t = 0;
                }
                ++brInd;
                ++brIt;
            }

            neuriteOut.vSec.push_back(sec);
        	}
    	}
	}

	////////////////////////////////////////////////////////////////////////
	/// create_neurite_old
	////////////////////////////////////////////////////////////////////////
	void create_neurite_old
	(
		const std::vector<NeuriteProjector::Neurite>& vNeurites,
		const std::vector<std::vector<vector3> >& vPos,
		const std::vector<std::vector<number> >& vR,
		size_t nid,
		Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		std::vector<Vertex*>* connectingVrts = NULL,
		std::vector<Edge*>* connectingEdges = NULL,
		std::vector<Vertex*>* outVerts = NULL,
		std::vector<number>* outRads = NULL,
		bool bWithER = false
	)
	{
    const NeuriteProjector::Neurite& neurite = vNeurites[nid];
    const std::vector<vector3>& pos = vPos[nid];
    std::vector<number> r = vR[nid];

    number neurite_length = 0.0;
    for (size_t i = 1; i < pos.size(); ++i)
    	neurite_length += VecDistance(pos[i], pos[i-1]);

    size_t nSec = neurite.vSec.size();

    const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
    std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit = vBR.begin();
    std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end = vBR.end();

    std::vector<Vertex*> vVrt;
    std::vector<Edge*> vEdge;
    vVrt.resize(4);
    vEdge.resize(4);

    vector3 vel;
    UG_COND_THROW(nSec == 0, "Number of sections > 0 required. FIX: Don't collapse root edges of neurites.");
    UG_DLOGN(NC_TNP, 0, "nSec: " << nSec);
    const NeuriteProjector::Section& sec = neurite.vSec[0];
    number h = sec.endParam;
    vel[0] = -3.0*sec.splineParamsX[0]*h*h - 2.0*sec.splineParamsX[1]*h - sec.splineParamsX[2];
    vel[1] = -3.0*sec.splineParamsY[0]*h*h - 2.0*sec.splineParamsY[1]*h - sec.splineParamsY[2];
    vel[2] = -3.0*sec.splineParamsZ[0]*h*h - 2.0*sec.splineParamsZ[1]*h - sec.splineParamsZ[2];

    vector3 projRefDir;
    VecNormalize(vel, vel);
    number fac = VecProd(neurite.refDir, vel);
    VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
    VecNormalize(projRefDir, projRefDir);

    vector3 thirdDir;
    VecCross(thirdDir, vel, projRefDir);

    number angleOffset = 0.0;

    if (connectingVrts && connectingEdges)
    {
        vVrt = *connectingVrts;
        vEdge = *connectingEdges;

        // calculate start angle offset
        vector3 center(0.0);
        for (size_t i = 0; i < 4; ++i)
            VecAdd(center, center, aaPos[(*connectingVrts)[i]]);
        center /= 4;

        vector3 centerToFirst;
        VecSubtract(centerToFirst, aaPos[(*connectingVrts)[0]], center);

        vector2 relCoord;
        VecScaleAdd(centerToFirst, 1.0, centerToFirst, -VecProd(centerToFirst, vel), vel);
        relCoord[0] = VecProd(centerToFirst, projRefDir);
        VecScaleAdd(centerToFirst, 1.0, centerToFirst, -relCoord[0], projRefDir);
        relCoord[1] = VecProd(centerToFirst, thirdDir);
        VecNormalize(relCoord, relCoord);

        if (fabs(relCoord[0]) < 1e-8)
            angleOffset = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
        else
            angleOffset = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);
        if (angleOffset < 0) angleOffset += 2.0*PI;

        // ignore first branching region (the connecting region)
        ++brit;
    }
    else
    {
        // create first layer of vertices/edges for membrane
        for (size_t i = 0; i < 4; ++i)
        {
            Vertex* v = *g.create<RegularVertex>();
            vVrt[i] = v;
            number angle = 0.5*PI*i;
            VecScaleAdd(aaPos[v], 1.0, pos[0], r[0]*cos(angle), projRefDir, r[0]*sin(angle), thirdDir);

            aaSurfParams[v].neuriteID = nid;
            aaSurfParams[v].axial = 0.0;
            aaSurfParams[v].angular = angle;
            outVerts->push_back(v);
            UG_DLOGN(NC_TNP, 0, "aaPos[v]: " << aaPos[v]);
        }
        outRads->push_back(r[0]);

        for (size_t i = 0; i < 4; ++i) {
            vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));
        }
    }

    // Now create dendrite to the next branching point and iterate this process.
    // We want to create each of the segments with approx. the same aspect ratio.
    // To that end, we first calculate the length of the section to be created (in units of radius)
    // and then divide this number by 2^n where n is the number of anisotropic refinements to
    // be performed to make all segments (more or less) isotropic. The result is the number
    // of segments to be used for the section.
    number t_start = 0.0;
    number t_end = 0.0;
    vector3 lastPos = pos[0];
    size_t curSec = 0;

    while (true)
    {
    	t_start = t_end;

    	// if there is another BP, store its extensions here
		number bp_start = 1.0;
		number bp_end = 0.0;

    	// last section: create until tip
    	if (brit == brit_end)
    		t_end = 1.0;

    	// otherwise: section goes to next branching point
    	else
    	{
    		// calculate the exact position of the branching point,
    		// i.e., the axial position of the intersection of the branching neurite's
    		// spline with the surface of the current neurite
    		// this is necessary esp. when the branching angle is very small
			const std::vector<uint32_t>& vBranchInd = brit->bp->vNid;
			size_t nBranches = vBranchInd.size();

    		for (size_t br = 1; br < nBranches; ++br)
    		{
    			uint32_t brInd = vBranchInd[br];

    			// get position and radius of first point of branch
    			const number brRadSeg1 = vR[brInd][0];

    			// get position and radius of branching point
    			const number bpTPos = brit->t;
    			size_t brSec = curSec;
    			for (; brSec < nSec; ++brSec)
    			{
    				const NeuriteProjector::Section& sec = neurite.vSec[brSec];
    				if (bpTPos - sec.endParam < 1e-6*bpTPos)
    					break;
    			}
    			UG_COND_THROW(brSec == nSec, "Could not find section containing branching point "
    				"at t = " << bpTPos << ".");
    			const number bpRad = vR[nid][brSec+1];

    			// calculate branch and neurite directions
				vector3 branchDir;
				const NeuriteProjector::Section& childSec = vNeurites[brInd].vSec[0];
				number te = childSec.endParam;

				const number* s = &childSec.splineParamsX[0];
				number& v0 = branchDir[0];
				v0 = -3.0*s[0]*te - 2.0*s[1];
				v0 = v0*te - s[2];

				s = &childSec.splineParamsY[0];
				number& v1 = branchDir[1];
				v1 = -3.0*s[0]*te - 2.0*s[1];
				v1 = v1*te - s[2];

				s = &childSec.splineParamsZ[0];
				number& v2 = branchDir[2];
				v2 = -3.0*s[0]*te - 2.0*s[1];
				v2 = v2*te - s[2];

				VecNormalize(branchDir, branchDir);

				vector3 neuriteDir;
    			const NeuriteProjector::Section& sec = neurite.vSec.at(brSec);
    			vel[0] = -sec.splineParamsX[2];
    			vel[1] = -sec.splineParamsY[2];
    			vel[2] = -sec.splineParamsZ[2];
    			number velNorm = sqrt(VecNormSquared(vel));
    			VecScale(neuriteDir, vel, 1.0/velNorm);

				// calculate offset of true branch position, which is r1*cot(alpha)
				number brScProd = VecProd(neuriteDir, branchDir);
				number surfBPoffset = bpRad * brScProd / sqrt(1.0 - brScProd*brScProd);

				// calculate true length of BP, which is 2*r2*sin(alpha)
				number surfBPhalfLength = brRadSeg1 * sqrt(1.0 - brScProd*brScProd);

				// finally set bp start and end
				bp_start = std::min(bp_start, bpTPos + (surfBPoffset - surfBPhalfLength) / neurite_length);
				bp_end = std::max(bp_end, bpTPos + (surfBPoffset + surfBPhalfLength) / neurite_length);
    		}

    		t_end = bp_start;
    	}

    	// calculate total length in units of radius
    	// = integral from t_start to t_end over: ||v(t)|| / r(t) dt
    	number lengthOverRadius = calculate_length_over_radius(t_start, t_end, neurite, curSec);
    	size_t nSeg = (size_t) floor(lengthOverRadius / 8);
    	// at least one segment is required to create a neurite
    	if (nSeg == 0) { nSeg = 1; }
    	UG_COND_THROW(nSeg == 0, "Number of segments > 0 required.");
    	number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
    	UG_DLOGN(NC_TNP, 0, "segLength: " << segLength);
    	UG_DLOGN(NC_TNP, 0, "nSeg: " << nSeg);
    	std::vector<number> vSegAxPos(nSeg);
    	calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite, curSec, segLength);

    	// add the branching point to segment list (if present)
    	if (brit != brit_end)
    	{
    		vSegAxPos.resize(nSeg+1);
    		vSegAxPos[nSeg] = bp_end;
    		++nSeg;
    	}

    	// create mesh for segments
    	Selector sel(g);
    	for (size_t s = 0; s < nSeg; ++s)
    	{
    		// get exact position, velocity and radius of segment end
    		number segAxPos = vSegAxPos[s];
    		for (; curSec < nSec; ++curSec)
    		{
				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				if (sec.endParam >= segAxPos)
					break;
    		}

    		const NeuriteProjector::Section& sec = neurite.vSec[curSec];
    		vector3 curPos;
    		number monom = sec.endParam - segAxPos;
    		const number* sp = &sec.splineParamsX[0];
			number& p0 = curPos[0];
			number& v0 = vel[0];
			p0 = sp[0]*monom + sp[1];
			p0 = p0*monom + sp[2];
			p0 = p0*monom + sp[3];
			v0 = -3.0*sp[0]*monom -2.0*sp[1];
			v0 = v0*monom - sp[2];

    		sp = &sec.splineParamsY[0];
			number& p1 = curPos[1];
			number& v1 = vel[1];
			p1 = sp[0]*monom + sp[1];
			p1 = p1*monom + sp[2];
			p1 = p1*monom + sp[3];
			v1 = -3.0*sp[0]*monom -2.0*sp[1];
			v1 = v1*monom - sp[2];

    		sp = &sec.splineParamsZ[0];
			number& p2 = curPos[2];
			number& v2 = vel[2];
			p2 = sp[0]*monom + sp[1];
			p2 = p2*monom + sp[2];
			p2 = p2*monom + sp[3];
			v2 = -3.0*sp[0]*monom -2.0*sp[1];
			v2 = v2*monom - sp[2];

    		sp = &sec.splineParamsR[0];
			number radius;
			radius = sp[0]*monom + sp[1];
			radius = radius*monom + sp[2];
			radius = radius*monom + sp[3];

			VecNormalize(vel, vel);

			// calculate reference dir projected to normal plane of velocity
			number fac = VecProd(neurite.refDir, vel);
			VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
			VecNormalize(projRefDir, projRefDir);

			VecCross(thirdDir, vel, projRefDir);

			// extrude from last pos to new pos
			if (s == nSeg-1 && brit != brit_end)
				sel.enable_autoselection(true); // for last segment (BP), select new elems

			vector3 extrudeDir;
			VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, lastPos);
			Extrude(g, &vVrt, &vEdge, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

			sel.enable_autoselection(false);

			// set new positions and param attachments; also ensure correct face orientation
			for (size_t j = 0; j < 4; ++j)
			{
				number angle = 0.5*PI*j + angleOffset;
				if (angle > 2*PI) angle -= 2*PI;
				Vertex* v = vVrt[j];
				vector3 radialVec;
				VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
				VecAdd(aaPos[v], curPos, radialVec);

				aaSurfParams[v].neuriteID = nid;
				aaSurfParams[v].axial = segAxPos;
				aaSurfParams[v].angular = angle;

				Grid::traits<Face>::secure_container faceCont;
				g.associated_elements(faceCont, vEdge[j]);  // faceCont must contain exactly one face
				vector3 normal;
				CalculateNormal(normal, faceCont[0], aaPos);
				if (VecProd(normal, radialVec) < 0)
					g.flip_orientation(faceCont[0]);
			}
			lastPos = curPos;
	    }

    	// connect branching neurites if present
    	if (brit != brit_end)
    	{
    		// find branching child neurite
			SmartPtr<NeuriteProjector::BranchingPoint> bp = brit->bp;
			UG_COND_THROW(bp->vNid.size() > 2,
				"This implementation can only handle branching points with one branching child.");

			size_t child_nid;
			if (bp->vNid[0] != nid)
				child_nid = bp->vNid[0];
			else
				child_nid = bp->vNid[1];

			// find out branching child neurite initial direction
			vector3 childDir;
			const NeuriteProjector::Section& childSec = vNeurites[child_nid].vSec[0];
			number te = childSec.endParam;

			const number* s = &childSec.splineParamsX[0];
			number& v0 = childDir[0];
			v0 = -3.0*s[0]*te - 2.0*s[1];
			v0 = v0*te - s[2];

			s = &childSec.splineParamsY[0];
			number& v1 = childDir[1];
			v1 = -3.0*s[0]*te - 2.0*s[1];
			v1 = v1*te - s[2];

			s = &childSec.splineParamsZ[0];
			number& v2 = childDir[2];
			v2 = -3.0*s[0]*te - 2.0*s[1];
			v2 = v2*te - s[2];

			// now choose best side of hexahedron to connect to
			vector3 normal;
			Face* best = NULL;
			number bestProd = 0.0;
			Selector::traits<Face>::iterator fit = sel.faces_begin();
			Selector::traits<Face>::iterator fit_end = sel.faces_end();
			for (; fit != fit_end; ++fit)
			{
				CalculateNormal(normal, *fit, aaPos);
				number prod = VecProd(normal, childDir);
				if (prod > bestProd)
				{
					best = *fit;
					bestProd = prod;
				}
			}
			UG_COND_THROW(!best, "None of the branching point faces pointed in a suitable direction.")
			sel.deselect(sel.faces_begin(), sel.faces_end());

			// remove face and call creation of child neurite (recursion)
			std::vector<Vertex*> vrts(4);
			UG_COND_THROW(best->num_vertices() != 4, "Hexaeder face does not have 4 vertices!");
			for (size_t j = 0; j < 4; ++j)
				vrts[j] = best->vertex(j);
			std::vector<Edge*> edges(4);

			Grid::traits<Edge>::secure_container edgeCont;
			g.associated_elements(edgeCont, best);
			size_t esz = edgeCont.size();

			for (size_t j = 0; j < 4; ++j)
			{
				Vertex* first = vrts[j];
				Vertex* second = vrts[(j+1) % 4];

				size_t k = 0;
				for (; k < esz; ++k)
				{
					if ((edgeCont[k]->vertex(0) == first && edgeCont[k]->vertex(1) == second)
						|| (edgeCont[k]->vertex(0) == second && edgeCont[k]->vertex(1) == first))
					{
						edges[j] = edgeCont[k];
						break;
					}
				}
				UG_COND_THROW(k == esz, "Connecting edges for child neurite could not be determined.");
			}

			g.erase(best);

			UG_DLOGN(NC_TNP, 0, "Creating child")
			/// Note:  create prism to connect to in case the branching angle is small or big
			create_neurite_old(vNeurites, vPos, vR, child_nid, g, aaPos, aaSurfParams, &vrts, &edges, NULL, NULL, bWithER);
    	}

    	// update t_end and curSec
    	if (brit != brit_end)
    		t_end = bp_end;

    	for (; curSec < nSec; ++curSec)
    	{
			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			if (sec.endParam >= t_end)
				break;
    	}

    	// check whether tip has been reached
    	if (brit == brit_end)
    		break;
    	else
			++brit;
    }

    // close the tip of the neurite
    const NeuriteProjector::Section& lastSec = neurite.vSec[nSec-1];
    vel = vector3(-lastSec.splineParamsX[2], -lastSec.splineParamsY[2], -lastSec.splineParamsZ[2]);
    number radius = lastSec.splineParamsR[3];
    VecScale(vel, vel, radius/sqrt(VecProd(vel, vel)));
    Extrude(g, &vVrt, &vEdge, NULL, vel, aaPos, EO_CREATE_FACES, NULL);
    vector3 center = CalculateBarycenter(vVrt.begin(), vVrt.end(), aaPos);
    MergeMultipleVertices(g, vVrt.begin(), vVrt.end());

    Vertex* v = *vVrt.begin();
    aaPos[v] = center;

    aaSurfParams[v].neuriteID = nid;
    aaSurfParams[v].axial = 2.0;
    aaSurfParams[v].angular = 0.0;
	}


	////////////////////////////////////////////////////////////////////////
	/// test_split_geom
	////////////////////////////////////////////////////////////////////////
	void test_split_geom(number percentage) {
		Grid g;
		SubsetHandler sh(g);
		sh.set_default_subset_index(0);
		g.attach_to_vertices(aPosition);
		Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
		Selector sel(g);

		std::vector<ug::vector3> vCoords;
		vCoords.push_back(ug::vector3(0,0,0));
		vCoords.push_back(ug::vector3(0,1,0));
		vCoords.push_back(ug::vector3(1,1,0));
		vCoords.push_back(ug::vector3(1,0,0));
		std::vector<Vertex*> vVrt;
		std::vector<Edge*> vEdge;
		vVrt.resize(4);
		vEdge.resize(4);
		for (size_t i = 0; i < 4; ++i)
		{
		   Vertex* v = *g.create<RegularVertex>();
		   vVrt[i] = v;
		   aaPos[v] = vCoords[i];
		   sel.select(v);
		}

		 for (size_t i = 0; i < 4; ++i) {
		     vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));
		 }

		 IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "test_shrunk_geom2_before.ugx");
		 ug::vector3 diffVec;
		 VecSubtract(diffVec, aaPos[vVrt[0]], aaPos[vVrt[(1)%4]]);
		 std::vector<ug::Vertex*> vertices;
		 std::vector<ug::Edge*> edges;
		 split_quadrilateral_along_edges(vVrt, g, aaPos, percentage, diffVec, vertices, edges);
		 IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "test_shrunk_geom2_after.ugx");
	}

	////////////////////////////////////////////////////////////////////////
	/// test_shrink_geom_copy
	////////////////////////////////////////////////////////////////////////
	void test_shrink_geom_copy(number length=0.1) {
		Grid g;
		SubsetHandler sh(g);
	    sh.set_default_subset_index(0);
	    g.attach_to_vertices(aPosition);
	    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	    Selector sel(g);

	    std::vector<ug::vector3> vCoords;
	    vCoords.push_back(ug::vector3(0,0,0));
	    vCoords.push_back(ug::vector3(0,1,0));
	    vCoords.push_back(ug::vector3(1,1,0));
	    vCoords.push_back(ug::vector3(1,0,0));
	    std::vector<Vertex*> vVrt;
	    std::vector<Edge*> vEdge;
	    vVrt.resize(4);
	    vEdge.resize(4);
	    for (size_t i = 0; i < 4; ++i)
	    {
	         Vertex* v = *g.create<RegularVertex>();
	         vVrt[i] = v;
	         aaPos[v] = vCoords[i];
	         sel.select(v);
	    }

	    for (size_t i = 0; i < 4; ++i) {
	        vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));
	    }

	    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "test_shrunk_geom_copy_before.ugx");
	    std::vector<ug::Vertex*> vVrtOut;
	    std::vector<ug::Edge*> vEdgeOut;
	    shrink_quadrilateral_copy(vVrt, vVrtOut, vVrtOut, vEdgeOut, g, aaPos, length);
	    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "test_shrunk_geom_copy_after.ugx");
	}

	////////////////////////////////////////////////////////////////////////
	/// test_shrink_geom
	////////////////////////////////////////////////////////////////////////
	void test_shrink_geom
	(
		number length=0.01
	)
	{
	Grid g;
	SubsetHandler sh(g);
    sh.set_default_subset_index(0);
    g.attach_to_vertices(aPosition);
    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
    Selector sel(g);

    std::vector<ug::vector3> vCoords;
    vCoords.push_back(ug::vector3(0,0,0));
    vCoords.push_back(ug::vector3(0,1,0));
    vCoords.push_back(ug::vector3(1,1,0));
    vCoords.push_back(ug::vector3(1,0,0));
    std::vector<Vertex*> vVrt;
    std::vector<Edge*> vEdge;
    vVrt.resize(4);
    vEdge.resize(4);
    for (size_t i = 0; i < 4; ++i)
    {
         Vertex* v = *g.create<RegularVertex>();
         vVrt[i] = v;
         aaPos[v] = vCoords[i];
         sel.select(v);
    }

    for (size_t i = 0; i < 4; ++i) {
        vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));
    }

    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "test_shrunk_geom_before.ugx");

    ug::vector3 center;
    center = CalculateBarycenter(sel.vertices_begin(), sel.vertices_end(), aaPos);
    sel.clear();
    for (size_t i = 0; i < 4; ++i)
    {
       ug::vector3 dir;
       VecSubtract(dir, aaPos[vVrt[i]], center);
       UG_DLOGN(NC_TNP, 0, "dir:" << dir)
       VecScaleAdd(aaPos[vVrt[i]], 1.0, aaPos[vVrt[i]], length, dir);
       if (VecLength(dir) > VecDistance(center, aaPos[vVrt[i]])) {
    	   UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
       }
    }

    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "test_shrunk_geom_after.ugx");
	}

	////////////////////////////////////////////////////////////////////////
	/// create_neurite
	////////////////////////////////////////////////////////////////////////
	static void create_neurite
	(
		const std::vector<NeuriteProjector::Neurite>& vNeurites,
		const std::vector<std::vector<vector3> >& vPos,
		const std::vector<std::vector<number> >& vR,
		size_t nid,
		number anisotropy,
		Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		std::vector<Vertex*>* connectingVrts = NULL,
		std::vector<Edge*>* connectingEdges = NULL,
		std::vector<Face*>* connectingFaces = NULL,
		number initialOffset = 0.0
	)
	{
	const NeuriteProjector::Neurite& neurite = vNeurites[nid];
	const std::vector<vector3>& pos = vPos[nid];
	const std::vector<number>& r = vR[nid];

	number neurite_length = 0.0;
	for (size_t i = 1; i < pos.size(); ++i)
		neurite_length += VecDistance(pos[i], pos[i-1]);

	size_t nSec = neurite.vSec.size();

	const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit = vBR.begin();
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end = vBR.end();

	std::vector<Vertex*> vVrt;
	std::vector<Edge*> vEdge;
	std::vector<Face*> vFace;
	vVrt.resize(4);
	vEdge.resize(4);
	vFace.resize(1);

	vector3 vel;
	const NeuriteProjector::Section& sec = neurite.vSec[0];
	number h = sec.endParam;
	vel[0] = -3.0*sec.splineParamsX[0]*h*h - 2.0*sec.splineParamsX[1]*h - sec.splineParamsX[2];
	vel[1] = -3.0*sec.splineParamsY[0]*h*h - 2.0*sec.splineParamsY[1]*h - sec.splineParamsY[2];
	vel[2] = -3.0*sec.splineParamsZ[0]*h*h - 2.0*sec.splineParamsZ[1]*h - sec.splineParamsZ[2];

	vector3 projRefDir;
	VecNormalize(vel, vel);
	number fac = VecProd(neurite.refDir, vel);
	VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
	VecNormalize(projRefDir, projRefDir);

	vector3 thirdDir;
	VecCross(thirdDir, vel, projRefDir);

	number angleOffset = 0.0;

	number t_start = 0.0;
	number t_end = 0.0;

	if (connectingVrts && connectingEdges && connectingFaces)
	{
		vVrt = *connectingVrts;
		vEdge = *connectingEdges;
		vFace = *connectingFaces;

		// calculate start angle offset
		vector3 center(0.0);
		for (size_t i = 0; i < 4; ++i)
			VecAdd(center, center, aaPos[(*connectingVrts)[i]]);
		center /= 4;

		vector3 centerToFirst;
		VecSubtract(centerToFirst, aaPos[(*connectingVrts)[0]], center);

		vector2 relCoord;
		VecScaleAdd(centerToFirst, 1.0, centerToFirst, -VecProd(centerToFirst, vel), vel);
		relCoord[0] = VecProd(centerToFirst, projRefDir);
		VecScaleAdd(centerToFirst, 1.0, centerToFirst, -relCoord[0], projRefDir);
		relCoord[1] = VecProd(centerToFirst, thirdDir);
		VecNormalize(relCoord, relCoord);

		if (fabs(relCoord[0]) < 1e-8)
			angleOffset = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
		else
			angleOffset = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);
		if (angleOffset < 0) angleOffset += 2.0*PI;

		// ignore first branching region (the connecting region)
		++brit;

		// apply initial offset (to prevent first segment being shorter than the others)
		t_end = initialOffset / neurite_length;
	}
	else
	{
		// create first layer of vertices/edges
		// surface vertices first
		for (size_t i = 0; i < 4; ++i)
		{
			Vertex* v = *g.create<RegularVertex>();
			vVrt[i] = v;
			number angle = 0.5*PI*i;
			VecScaleAdd(aaPos[v], 1.0, pos[0], r[0]*cos(angle), projRefDir, r[0]*sin(angle), thirdDir);

			aaSurfParams[v].neuriteID = nid;
			aaSurfParams[v].axial = 0.0;
			aaSurfParams[v].angular = angle;
			aaSurfParams[v].radial = 1.0;
		}

		// edges
		for (size_t i = 0; i < 4; ++i)
			vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));

		// faces
		vFace[0] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3]));
	}

	// Now create dendrite to the next branching point and iterate this process.
	// We want to create each of the segments with approx. the same aspect ratio.
	// To that end, we first calculate the length of the section to be created (in units of radius)
	// and then divide this number by 2^n where n is the number of anisotropic refinements to
	// be performed to make all segments (more or less) isotropic. The result is the number
	// of segments to be used for the section.
	vector3 lastPos = pos[0];
	size_t curSec = 0;

	while (true)
	{
		t_start = t_end;

		// if there is another BP, store its extensions here
		number bp_start = 1.0;
		number bp_end = 0.0;

		// initial branch offsets (to prevent the first segment being shorter than the following)
		std::vector<number> branchOffset;
		number surfBPoffset;

		// last section: create until tip
		if (brit == brit_end)
			t_end = 1.0;

		// otherwise: section goes to next branching point
		else
		{
			// calculate the exact position of the branching point,
			// i.e., the axial position of the intersection of the branching neurite's
			// spline with the surface of the current neurite
			// this is necessary esp. when the branching angle is very small
			const std::vector<uint32_t>& vBranchInd = brit->bp->vNid;
			size_t nBranches = vBranchInd.size();

			branchOffset.resize(brit->bp->vNid.size(), 0.0);
			for (size_t br = 1; br < nBranches; ++br)
			{
				uint32_t brInd = vBranchInd[br];

				// get position and radius of first point of branch
				const number brRadSeg1 = vR[brInd][0];

				// get position and radius of branching point
				const number bpTPos = brit->t;
				size_t brSec = curSec;
				for (; brSec < nSec; ++brSec)
				{
					const NeuriteProjector::Section& sec = neurite.vSec[brSec];
					if (bpTPos - sec.endParam < 1e-6*bpTPos)
						break;
				}
				UG_COND_THROW(brSec == nSec, "Could not find section containing branching point "
					"at t = " << bpTPos << ".");
				const number bpRad = vR[nid][brSec+1];

				// calculate branch and neurite directions
				vector3 branchDir;
				const NeuriteProjector::Section& childSec = vNeurites[brInd].vSec[0];
				number te = childSec.endParam;

				const number* s = &childSec.splineParamsX[0];
				number& v0 = branchDir[0];
				v0 = -3.0*s[0]*te - 2.0*s[1];
				v0 = v0*te - s[2];

				s = &childSec.splineParamsY[0];
				number& v1 = branchDir[1];
				v1 = -3.0*s[0]*te - 2.0*s[1];
				v1 = v1*te - s[2];

				s = &childSec.splineParamsZ[0];
				number& v2 = branchDir[2];
				v2 = -3.0*s[0]*te - 2.0*s[1];
				v2 = v2*te - s[2];

				VecNormalize(branchDir, branchDir);

				vector3 neuriteDir;
				const NeuriteProjector::Section& sec = neurite.vSec.at(brSec);
				vel[0] = -sec.splineParamsX[2];
				vel[1] = -sec.splineParamsY[2];
				vel[2] = -sec.splineParamsZ[2];
				number velNorm = sqrt(VecNormSquared(vel));
				VecScale(neuriteDir, vel, 1.0/velNorm);

				// calculate offset of true branch position, which is r1/sqrt(2)*cot(alpha)
				const number brScProd = VecProd(neuriteDir, branchDir);
				const number sinAlphaInv = 1.0 / sqrt(1.0 - brScProd*brScProd);
				surfBPoffset = 0.5*sqrt(2.0) * bpRad * brScProd * sinAlphaInv;

				// calculate offset of new branch, which is r1/sqrt(2)/sin(alpha)
				branchOffset[br] = 0.5*sqrt(2.0) * bpRad * sinAlphaInv;

				// calculate true half-length of BP, which is r2/sin(alpha)
				number surfBPhalfLength = brRadSeg1 * sinAlphaInv;

				// finally set bp start and end
				bp_start = std::min(bp_start, bpTPos + (surfBPoffset - surfBPhalfLength) / neurite_length);
				bp_end = std::max(bp_end, bpTPos + (surfBPoffset + surfBPhalfLength) / neurite_length);
			}

			t_end = bp_start;
		}

		// calculate total length in units of radius
		// = integral from t_start to t_end over: ||v(t)|| / r(t) dt
		number lengthOverRadius = calculate_length_over_radius(t_start, t_end, neurite, curSec);

		// to reach the desired anisotropy on the surface in the refinement limit,
		// it has to be multiplied by pi/2 h
		size_t nSeg = (size_t) floor(lengthOverRadius / (anisotropy*0.5*PI));
		if (!nSeg)
			nSeg = 1;
		number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
		std::vector<number> vSegAxPos(nSeg);
		calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite, curSec, segLength);

		// add the branching point to segment list (if present)
		if (brit != brit_end)
		{
			vSegAxPos.resize(nSeg+1);
			vSegAxPos[nSeg] = bp_end;
			++nSeg;
		}


		// in case we construct to a BP, find out the branching angle
		// in order to adjust this neurites offset bit by bit
		number addOffset = 0.0;
		size_t child_nid;
		size_t connFaceInd = 0;
		if (brit != brit_end)
		{
			// find branching child neurite
			SmartPtr<NeuriteProjector::BranchingPoint> bp = brit->bp;
			UG_COND_THROW(bp->vNid.size() > 2,
				"This implementation can only handle branching points with one branching child.");

			if (bp->vNid[0] != nid)
				child_nid = bp->vNid[0];
			else
				child_nid = bp->vNid[1];

			// find out branching child neurite initial direction
			vector3 childDir;
			const NeuriteProjector::Section& childSec = vNeurites[child_nid].vSec[0];
			number te = childSec.endParam;

			const number* sp = &childSec.splineParamsX[0];
			number& vc0 = childDir[0];
			vc0 = -3.0*sp[0]*te - 2.0*sp[1];
			vc0 = vc0*te - sp[2];

			sp = &childSec.splineParamsY[0];
			number& vc1 = childDir[1];
			vc1 = -3.0*sp[0]*te - 2.0*sp[1];
			vc1 = vc1*te - sp[2];

			sp = &childSec.splineParamsZ[0];
			number& vc2 = childDir[2];
			vc2 = -3.0*sp[0]*te - 2.0*sp[1];
			vc2 = vc2*te - sp[2];


			// find out neurite direction in next BP
			number bpAxPos = vSegAxPos[nSeg-1];
			size_t tmpSec = curSec;
			for (; tmpSec < nSec; ++tmpSec)
			{
				const NeuriteProjector::Section& sec = neurite.vSec[tmpSec];
				if (sec.endParam >= bpAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[tmpSec];
			number monom = sec.endParam - bpAxPos;
			sp = &sec.splineParamsX[0];
			number& v0 = vel[0];
			v0 = -3.0*sp[0]*monom -2.0*sp[1];
			v0 = v0*monom - sp[2];

			sp = &sec.splineParamsY[0];
			number& v1 = vel[1];
			v1 = -3.0*sp[0]*monom -2.0*sp[1];
			v1 = v1*monom - sp[2];

			sp = &sec.splineParamsZ[0];
			number& v2 = vel[2];
			v2 = -3.0*sp[0]*monom -2.0*sp[1];
			v2 = v2*monom - sp[2];

			VecNormalize(vel, vel);

			// calculate offset
			number fac = VecProd(neurite.refDir, vel);
			VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
			VecNormalize(projRefDir, projRefDir);
			VecCross(thirdDir, vel, projRefDir);

			/*
			vector2 relCoord;
			VecScaleAppend(childDir, -VecProd(childDir, vel), vel);
			relCoord[0] = VecProd(childDir, projRefDir);
			VecScaleAppend(childDir, -relCoord[0], projRefDir);
			relCoord[1] = VecProd(childDir, thirdDir);
			VecNormalize(relCoord, relCoord);
			*/
			vector2 relCoord(VecProd(childDir, projRefDir), VecProd(childDir, thirdDir));
			VecNormalize(relCoord, relCoord);

			number branchOffset = 0.0;
			if (fabs(relCoord[0]) < 1e-8)
				branchOffset = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
			else
				branchOffset = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);

			addOffset = branchOffset - angleOffset;
			connFaceInd = floor(std::fmod(addOffset+4*PI, 2*PI) / (PI/2));
			addOffset = std::fmod(addOffset - (connFaceInd*PI/2 + PI/4) + 4*PI, 2*PI);
			if (addOffset > PI)
				addOffset -= 2*PI;
			addOffset /= nSeg - 1;
		}

		// create mesh for segments
		Selector sel(g);
		for (size_t s = 0; s < nSeg; ++s)
		{
			// get exact position, velocity and radius of segment end
			number segAxPos = vSegAxPos[s];
			for (; curSec < nSec; ++curSec)
			{
				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				if (sec.endParam >= segAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			vector3 curPos;
			number monom = sec.endParam - segAxPos;
			const number* sp = &sec.splineParamsX[0];
			number& p0 = curPos[0];
			number& v0 = vel[0];
			p0 = sp[0]*monom + sp[1];
			p0 = p0*monom + sp[2];
			p0 = p0*monom + sp[3];
			v0 = -3.0*sp[0]*monom -2.0*sp[1];
			v0 = v0*monom - sp[2];

			sp = &sec.splineParamsY[0];
			number& p1 = curPos[1];
			number& v1 = vel[1];
			p1 = sp[0]*monom + sp[1];
			p1 = p1*monom + sp[2];
			p1 = p1*monom + sp[3];
			v1 = -3.0*sp[0]*monom -2.0*sp[1];
			v1 = v1*monom - sp[2];

			sp = &sec.splineParamsZ[0];
			number& p2 = curPos[2];
			number& v2 = vel[2];
			p2 = sp[0]*monom + sp[1];
			p2 = p2*monom + sp[2];
			p2 = p2*monom + sp[3];
			v2 = -3.0*sp[0]*monom -2.0*sp[1];
			v2 = v2*monom - sp[2];

			sp = &sec.splineParamsR[0];
			number radius;
			radius = sp[0]*monom + sp[1];
			radius = radius*monom + sp[2];
			radius = radius*monom + sp[3];

			VecNormalize(vel, vel);

			// calculate reference dir projected to normal plane of velocity
			number fac = VecProd(neurite.refDir, vel);
			VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
			VecNormalize(projRefDir, projRefDir);

			VecCross(thirdDir, vel, projRefDir);

			// usual segment: extrude
			if (s != nSeg - 1 || brit == brit_end)
			{
				// apply additional offset
				angleOffset = std::fmod(angleOffset + addOffset + 2*PI, 2*PI);

				// extrude from last pos to new pos
				vector3 extrudeDir;
				VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, lastPos);
				std::vector<Volume*> vVol;
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos, EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j)
				{
					number angle = 0.5*PI*j + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// ensure correct volume orientation
				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);
			}

			// BP segment: create BP with tetrahedra/pyramids and create whole branch
			else
			{
				std::vector<Vertex*> vNewVrt(4, NULL);

				// create all needed vertices (except branch mid point)
				for (size_t j = 0; j < 4; ++j)
				{
					Vertex* v = vNewVrt[j] = *g.create<RegularVertex>();

					number angle = 0.5*PI*j + angleOffset;
					if (angle > 2*PI)
						angle -= 2*PI;
					vector3 radialVec;
					VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// correct offsets of non-connecting vertices (opposite direction than branching ones)
				VecScaleAppend(aaPos[vVrt[(connFaceInd+2)%4]], -2.0*surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd+3)%4]], -2.0*surfBPoffset, vel);
				VecScaleAppend(aaPos[vNewVrt[(connFaceInd+2)%4]], -2.0*surfBPoffset, vel);
				VecScaleAppend(aaPos[vNewVrt[(connFaceInd+3)%4]], -2.0*surfBPoffset, vel);

				aaSurfParams[vVrt[(connFaceInd+2)%4]].axial -= 2.0*surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd+3)%4]].axial -= 2.0*surfBPoffset / neurite_length;
				aaSurfParams[vNewVrt[(connFaceInd+2)%4]].axial -= 2.0*surfBPoffset / neurite_length;
				aaSurfParams[vNewVrt[(connFaceInd+3)%4]].axial -= 2.0*surfBPoffset / neurite_length;

				// prepare connectingVrts, -Edges, -Faces for branch
				std::vector<Vertex*> vBranchVrts(4);
				vBranchVrts[0] = vVrt[(connFaceInd+1)%4];
				vBranchVrts[1] = vNewVrt[(connFaceInd+1)%4];
				vBranchVrts[2] = vNewVrt[connFaceInd];
				vBranchVrts[3] = vVrt[connFaceInd];

				std::vector<Edge*> vBranchEdges(4);
				std::vector<Face*> vBranchFaces(1);
				for (size_t j = 0; j < 4; ++j)
				{
					if (j != 3)
						vBranchEdges[j] = *g.create<RegularEdge>(EdgeDescriptor(vBranchVrts[j], vBranchVrts[(j+1)%4]));
					else
						vBranchEdges[j] = vEdge[connFaceInd];
				}
				vBranchFaces[0] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vBranchVrts[0], vBranchVrts[1], vBranchVrts[2], vBranchVrts[3]));

				// add branch neurite ID to its initial vertices
				for (size_t j = 0; j < 4; ++j)
				{
					Vertex* brVrt = vBranchVrts[j];
					aaSurfParams[brVrt].neuriteID += (brit - vBR.begin()) << 20;  // add branching region index
					aaSurfParams[brVrt].neuriteID += 1 << 28;  // add child ID (always 0, since there can only be one child here)
				}

				// recursively build branch
				create_neurite(vNeurites, vPos, vR, child_nid, anisotropy,
					g, aaPos, aaSurfParams, &vBranchVrts, &vBranchEdges, &vBranchFaces, branchOffset[1]);


				// prepare connectingVrts, -Edges, -Faces for remainder of own neurite
				for (size_t j = 0; j < 4; ++j)
				{
					if (j != connFaceInd)
						vEdge[j] = *g.create<RegularEdge>(EdgeDescriptor(vNewVrt[j], vNewVrt[(j+1)%4]));
					else
						vEdge[j] = vBranchEdges[1];
				}

				vFace[0] = *g.create<Quadrilateral>(QuadrilateralDescriptor(vNewVrt[0], vNewVrt[1], vNewVrt[2], vNewVrt[3]));


				// create all inner BP elements
				g.create<Hexahedron>(HexahedronDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3],
					vNewVrt[0], vNewVrt[1], vNewVrt[2], vNewVrt[3]));

				/*
				Vertex* innerVrt = *g.create<RegularVertex>();

				VecScaleAdd(aaPos[innerVrt], 0.5, curPos, 0.5, lastPos);
				aaSurfParams[innerVrt].neuriteID = nid;
				aaSurfParams[innerVrt].neuriteID += (brit - vBR.begin()) << 20;  // add branching region index
				aaSurfParams[innerVrt].neuriteID += 1 << 28;  // add child ID (always 0, since there can only be one child here)
				aaSurfParams[innerVrt].axial = 0.5*vSegAxPos[s] + 0.5*vSegAxPos[s-1];
				aaSurfParams[innerVrt].angular = 0.0;
				aaSurfParams[innerVrt].radial = 0.0;

				g.create<Pyramid>(PyramidDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3], innerVrt));
				g.create<Pyramid>(PyramidDescriptor(vNewVrt[3], vNewVrt[2], vNewVrt[1], vNewVrt[0], innerVrt));
				for (size_t j = 0; j < 4; ++j)
					g.create<Pyramid>(PyramidDescriptor(vVrt[j], vNewVrt[j], vNewVrt[(j+1)%4], vVrt[(j+1)%4], innerVrt));
				*/

				vVrt.swap(vNewVrt);
			}

			lastPos = curPos;
		}

		// update t_end and curSec
		if (brit != brit_end)
			t_end = bp_end;

		for (; curSec < nSec; ++curSec)
		{
			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			if (sec.endParam >= t_end)
				break;
		}

		// check whether tip has been reached
		if (brit == brit_end)
			break;
		else
			++brit;
	}

	/*
	// close the tip of the neurite
	const NeuriteProjector::Section& lastSec = neurite.vSec[nSec-1];
	vel = vector3(-lastSec.splineParamsX[2], -lastSec.splineParamsY[2], -lastSec.splineParamsZ[2]);
	number radius = lastSec.splineParamsR[3];
	VecScale(vel, vel, radius/sqrt(VecProd(vel, vel)));
	Vertex* tip = *g.create<RegularVertex>();
	for (size_t i = 0; i < 4; ++i)
		VecScaleAppend(aaPos[tip], 0.25, aaPos[vVrt[i]]);
	VecAdd(aaPos[tip], aaPos[tip], vel);
	g.create<Pyramid>(PyramidDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3], tip));

	aaSurfParams[tip].neuriteID = nid;
	aaSurfParams[tip].axial = 2.0;
	aaSurfParams[tip].angular = 0.0;
	aaSurfParams[tip].radial = 1.0;
	*/
}

	////////////////////////////////////////////////////////////////////////
	/// export_to_ugx
	///////////////////////////////////////////////////////////////////////
	void export_to_ugx
	(
		Grid& g,
		SubsetHandler& sh,
		const std::string& fileName
	)
	{
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	if (!ugxWriter.write_to_file(fileName.c_str()))
		UG_THROW("Grid could not be written to file '" << fileName << "'.");
	}

	////////////////////////////////////////////////////////////////////////
	/// export_to_swc
	///////////////////////////////////////////////////////////////////////
	void export_to_swc
	(
		Grid& g,
		SubsetHandler& sh,
		const std::string& fileName
	)
	{
	// get access to positions
	UG_COND_THROW(!g.has_vertex_attachment(aPosition), "Position attachment not attached to grid.")
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	// get access to diameter attachment
	ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter");
	UG_COND_THROW(!g.has_vertex_attachment(aDiam), "No diameter attachment attached to grid.");
	Grid::AttachmentAccessor<Vertex, ANumber> aaDiam(g, aDiam);

    // analyze subset names to find out corresponding swc-types
    size_t nss = sh.num_subsets();
    std::vector<size_t> vType(nss);
    bool soma_subset_present = false;
    for (size_t i = 0; i < nss; ++i)
    {
        std::string name(sh.get_subset_name(i));
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (name.find("SOMA") != std::string::npos)
        {
            soma_subset_present = true;
            vType[i] = 1;
        }
        else if (name.find("AXON") != std::string::npos)
            vType[i] = 2;
        else if (name.find("APIC") != std::string::npos)
            vType[i] = 4;
        else if (name.find("DEND") != std::string::npos)
            vType[i] = 3;
        else vType[i] = 0;
    }

    if (!soma_subset_present)
        UG_DLOGN(NC_TNP, 0, "Warning: No somatic subset could be identified.")

	if (g.begin<Vertex>() == g.end<Vertex>())
	{
		UG_DLOGN(NC_TNP, 0, "Warning: No vertices contained in grid.")
		return;
	}

    // find soma vertex (if identifiable)
    Vertex* start = *g.begin<Vertex>();
    if (soma_subset_present)
    {
        g.begin_marking();
        std::queue<Vertex*> q; // corresponds to breadth-first
        q.push(start);
        while (!q.empty())
        {
            Vertex* v = q.front();
            if (vType[sh.get_subset_index(v)] == 1) break;
            g.mark(v);
            q.pop();

            // push neighboring elems to queue
            Grid::traits<Edge>::secure_container edges;
            g.associated_elements(edges, v);

            size_t sz = edges.size();
            for (size_t e = 0; e < sz; ++e)
            {
                Vertex* otherEnd = GetOpposingSide(g, edges[e], v);
                if (!g.is_marked(otherEnd))
                    q.push(otherEnd);
            }
        }
        g.end_marking();

        if (q.empty())
        {
            UG_DLOGN(NC_TNP, 0, "Warning: No soma vertex could be found in the requested neuron.");
        }
        else
            start = q.front();
    }

    // write the neuron to file
    std::ofstream outFile(fileName.c_str(), std::ios::out);
    UG_COND_THROW(!outFile.is_open(), "Could not open output file '" << fileName << "'.");

    outFile << "# This file has been generated by UG4." << std::endl;

    std::stack<std::pair<Vertex*, int> > stack; // corresponds to depth-first
    stack.push(std::make_pair(start, -1));

    g.begin_marking();
    int ind = 0;   // by convention, swc starts with index 1
    bool all_types_identified = true;
    while (!stack.empty())
    {
        // get all infos regarding vertex
        std::pair<Vertex*, int>& info = stack.top();
        Vertex* v = info.first;
        int conn = info.second;
        stack.pop();

        // mark curr vrt
        g.mark(v);

        size_t type = vType[sh.get_subset_index(v)];
        if (!type) all_types_identified = false;

        const Domain3d::position_type& coord = aaPos[v];

        number radius = 0.5*aaDiam[v];

        // write line to file
        outFile << ++ind << " " << type << " "
            << coord[0] << " " << coord[1] << " " << coord[2] << " "
            << radius << " " << conn << std::endl;

        // push neighboring elems to queue
        Grid::traits<Edge>::secure_container edges;
        g.associated_elements(edges, v);

        size_t sz = edges.size();
        for (size_t e = 0; e < sz; ++e)
        {
            Vertex* otherEnd = GetOpposingSide(g, edges[e], v);
            if (!g.is_marked(otherEnd))
                stack.push(std::make_pair(otherEnd, ind));
        }
    }
    g.end_marking();

    if (!all_types_identified)
        UG_DLOGN(NC_TNP, 0, "WARNING: Some vertex type(s) - soma, dendrite, axon, etc. -\n"
            "could not be identified using the subset names.\n"
            << "To ensure correct types in the resulting swc file, the ugx subset names\n"
            "need to contain one of the strings \"SOMA\", \"AXON\", \"DEND\", \"APIC\",\n"
            "upper/lower case can be ignored.");

    outFile.close();
	}


	////////////////////////////////////////////////////////////////////////
	/// swc_points_to_grid
	///////////////////////////////////////////////////////////////////////
	void swc_points_to_grid
	(
		const std::vector<SWCPoint>& vPts,
		Grid& g,
		SubsetHandler& sh,
		number scale_length = 1.0
	)
	{
	if (!g.has_vertex_attachment(aPosition))
		g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	ANumber aDiam = GlobalAttachments::attachment<ANumber>("diameter");
	if (!g.has_vertex_attachment(aDiam))
		g.attach_to_vertices(aDiam);
	Grid::AttachmentAccessor<Vertex, ANumber> aaDiam(g, aDiam);

	// create grid
	const size_t nP = vPts.size();
	std::vector<Vertex*> vrts(nP, NULL);
	for (size_t i = 0; i < nP; ++i)
	{
		const SWCPoint& pt = vPts[i];

		// create vertex and save
		Vertex* v = vrts[i] = *g.create<RegularVertex>();
		VecScale(aaPos[v], pt.coords, scale_length);
		sh.assign_subset(v, pt.type - 1);
		aaDiam[v] = 2 * pt.radius * scale_length;

		// create edge connections to already created vertices
		for (size_t j = 0; j < pt.conns.size(); ++j)
		{
			if (pt.conns[j] < i)
			{
				Edge* e = *g.create<RegularEdge>(EdgeDescriptor(vrts[pt.conns[j]], v));
				sh.assign_subset(e, vPts[pt.conns[j]].type - 1);
			}
		}
	}

	// final subset managment
	AssignSubsetColors(sh);
	sh.set_subset_name("soma", 0);
	sh.set_subset_name("axon", 1);
	sh.set_subset_name("dend", 2);
	sh.set_subset_name("apic", 3);
	sh.set_subset_name("fork", 4);
	sh.set_subset_name("end", 5);
	sh.set_subset_name("custom", 6);
	EraseEmptySubsets(sh);
	}


	////////////////////////////////////////////////////////////////////////
	/// read_swc
	///////////////////////////////////////////////////////////////////////
	void test_convert_swc_to_ugx(
		const std::string& fileName
	)
	{
		std::vector<SWCPoint> vPoints;
		import_swc(fileName, vPoints);

		// export original cell to ugx
		Grid g;
		SubsetHandler sh(g);
		swc_points_to_grid(vPoints, g, sh);
		std::string fn_noext = FilenameWithoutExtension(fileName);
		std::string fn = fn_noext + "_orig.ugx";
		export_to_ugx(g, sh, fn);
	}

	////////////////////////////////////////////////////////////////////////
	/// test_smoothing
	///////////////////////////////////////////////////////////////////////
	void test_smoothing(
		const std::string& fileName,
		size_t n,
		number h,
		number gamma
	)
	{
	std::vector<SWCPoint> vPoints;
	import_swc(fileName, vPoints);

	// export original cell to ugx
	Grid g;
	SubsetHandler sh(g);
	swc_points_to_grid(vPoints, g, sh);

	std::string fn_noext = FilenameWithoutExtension(fileName);
	std::string fn = fn_noext + "_orig.ugx";
	export_to_ugx(g, sh, fn);

	// smoothing
	smoothing(vPoints, n, h, gamma);

	// export smoothed cell to ugx
	g.clear_geometry();
	swc_points_to_grid(vPoints, g, sh);
	fn = fn_noext + "_smooth.ugx";
	export_to_ugx(g, sh, fn);

	// collapse short edges
	collapse_short_edges(g, sh);

	// export collapsed edges cell to ugx and swc
	fn = fn_noext + "_collapse.ugx";
	export_to_ugx(g, sh, fn);
	fn = fn_noext + "_precond.swc";
	export_to_swc(g, sh, fn);
	}

	////////////////////////////////////////////////////////////////////////
	/// test_import_swc
	///////////////////////////////////////////////////////////////////////
	void test_import_swc
	(
		const std::string& fileName,
		number anisotropy,
		size_t numRefs
	)
	{
	// preconditioning
	// test_smoothing(fileName, 5, 1.0, 1.0);

	// read in file to intermediate structure
	std::vector<SWCPoint> vPoints;
	std::vector<SWCPoint> vSomaPoints;
	import_swc(fileName, vPoints);

	// convert intermediate structure to neurite data
	std::vector<std::vector<vector3> > vPos;
	std::vector<std::vector<number> > vRad;
	std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
	std::vector<size_t> vRootNeuriteIndsOut;

	convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

	// prepare grid and projector
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	Selector sel(g);


	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
			"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	if (!g.has_vertex_attachment(aSP))
		g.attach_to_vertices(aSP);

	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(g, aSP);

	ProjectionHandler projHandler(&sh);
	SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
	projHandler.set_geometry(geom3d);

	SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
	projHandler.set_projector(0, neuriteProj);

	// create spline data
	std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
	create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

	// create coarse grid
	for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i)
		create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
			anisotropy, g, aaPos, aaSurfParams, NULL, NULL);

	// at branching points, we have not computed the correct positions yet,
	// so project the complete geometry using the projector
	VertexIterator vit = g.begin<Vertex>();
	VertexIterator vit_end = g.end<Vertex>();
	for (; vit != vit_end; ++vit)
		neuriteProj->project(*vit);

	// assign subset
	AssignSubsetColors(sh);
	sh.set_subset_name("neurites", 0);

	// output
	std::string outFileName = FilenameWithoutPath(std::string("testNeuriteProjector.ugx"));
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	ugxWriter.add_projection_handler(projHandler, "defPH", 0);
	if (!ugxWriter.write_to_file(outFileName.c_str()))
		UG_THROW("Grid could not be written to file '" << outFileName << "'.");


	// refinement
	Domain3d dom;
	try {LoadDomain(dom, outFileName.c_str());}
	UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");

	std::string curFileName("testNeuriteProjector.ugx");
	number offset = 5.0;

	curFileName = outFileName.substr(0, outFileName.size()-4) + "_refined_0.ugx";
	try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
	UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

	GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
	for (uint i = 0; i < numRefs; ++i)
	{
		ref.refine();

		std::ostringstream oss;
		oss << "_refined_" << i+1 << ".ugx";
		curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
		try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
		UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
	}
	}

	////////////////////////////////////////////////////////////////////////
	/// test_import_swc_with_er
	///////////////////////////////////////////////////////////////////////
	void test_import_swc_with_er
	(
		const std::string& fileNameIn,
		const std::string& fileNameOut,
		number erScaleFactor,
		number anisotropy,
		size_t numRefs
	)
	{
	// preconditioning
	//test_smoothing(fileNameIn, 5, 1.0, 1.0);


	// read in file to intermediate structure
	std::vector<SWCPoint> vPoints;
	std::vector<SWCPoint> vSomaPoints;
	//std::string fn_noext = FilenameWithoutExtension(fileNameIn);
	//std::string fn_precond = fn_noext + "_precond.swc";
	//import_swc(fn_precond, vPoints);
	import_swc(fileNameIn, vPoints);

	// convert intermediate structure to neurite data
	std::vector<std::vector<vector3> > vPos;
	std::vector<std::vector<number> > vRad;
	std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
	std::vector<size_t> vRootNeuriteIndsOut;

	convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

	// prepare grid and projector
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);


	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
			"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	if (!g.has_vertex_attachment(aSP))
		g.attach_to_vertices(aSP);

	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(g, aSP);

	SubsetHandler psh(g);
	psh.set_default_subset_index(0);

	ProjectionHandler projHandler(&psh);
	SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
	projHandler.set_geometry(geom3d);

	SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
	projHandler.set_projector(0, neuriteProj);

	// create spline data
	std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
	create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

	// create coarse grid
	for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i)
		create_neurite_with_er(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
			erScaleFactor, anisotropy, g, aaPos, aaSurfParams, sh);

	// assign subset
	AssignSubsetColors(sh);
	sh.set_subset_name("cyt", 0);
	sh.set_subset_name("er", 1);
	sh.set_subset_name("pm", 2);
	sh.set_subset_name("erm", 3);

	// output before projection
	std::string outFileNameBase = FilenameAndPathWithoutExtension(fileNameOut);
	IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, fileNameOut);

	// at branching points, we have not computed the correct positions yet,
	// so project the complete geometry using the projector
	VertexIterator vit = g.begin<Vertex>();
	VertexIterator vit_end = g.end<Vertex>();
	for (; vit != vit_end; ++vit)
		neuriteProj->project(*vit);

	// output after projection (overwrites the outFileName if projection successful)
	std::string outFileName = outFileNameBase + ".ugx";
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	ugxWriter.add_subset_handler(psh, "projSH", 0);
	ugxWriter.add_projection_handler(projHandler, "defPH", 0);
	if (!ugxWriter.write_to_file(outFileName.c_str()))
		UG_THROW("Grid could not be written to file '" << outFileName << "'.");


	// refinement
	Domain3d dom;
	dom.create_additional_subset_handler("projSH");
	try {LoadDomain(dom, outFileName.c_str());}
	UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");

	number offset = 5.0;
	std::string curFileName = outFileNameBase + "_refined_0.ugx";

	try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
	UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

	if (numRefs == 0)
		return;

	GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
	for (uint i = 0; i < numRefs; ++i)
	{
		ref.refine();
		std::ostringstream oss;
		oss << "_refined_" << i+1 << ".ugx";
		curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
		try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
		UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
	}
	}

	////////////////////////////////////////////////////////////////////////
	/// test_import_swc_surf
	///////////////////////////////////////////////////////////////////////
	void test_import_swc_surf(
		const std::string& fileName
	)
	{
	// preconditioning
	test_smoothing(fileName, 5, 1.0, 1.0);

	// read in file to intermediate structure
	std::vector<SWCPoint> vPoints;
	std::vector<SWCPoint> vSomaPoints;
	std::string fn_noext = FilenameWithoutExtension(fileName);
	std::string fn_precond = fn_noext + "_precond.swc";
	import_swc(fn_precond, vPoints);
	//import_swc(fileName, vPoints);

	// convert intermediate structure to neurite data
	std::vector<std::vector<vector3> > vPos;
	std::vector<std::vector<number> > vRad;
	std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
	std::vector<size_t> vRootNeuriteIndsOut;
	convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

	// prepare grid and projector
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	Selector sel(g);

	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
			"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	if (!g.has_vertex_attachment(aSP))
		g.attach_to_vertices(aSP);

	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(g, aSP);

	ProjectionHandler projHandler(&sh);
	SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
	projHandler.set_geometry(geom3d);

	SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
	projHandler.set_projector(0, neuriteProj);

	// create spline data
	std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
	create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

	// create coarse grid
	const number anisotropy = 2.0;
	for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i)
		create_neurite_surf(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
			anisotropy, g, aaPos, aaSurfParams, NULL, NULL);

	// at branching points, we have not computed the correct positions yet,
	// so project the complete geometry using the projector
	VertexIterator vit = g.begin<Vertex>();
	VertexIterator vit_end = g.end<Vertex>();
	for (; vit != vit_end; ++vit)
		neuriteProj->project(*vit);

	// create soma
	sel.clear();
	sh.set_default_subset_index(1);
	create_soma(vSomaPoints, g, aaPos);
	sh.set_default_subset_index(0);

	// refinement
	AssignSubsetColors(sh);
	sh.set_subset_name("neurites", 0);
	sh.set_subset_name("soma", 1);

	std::string outFileName = FilenameWithoutPath(fn_noext + "_surf.ugx");
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	ugxWriter.add_projection_handler(projHandler, "defPH", 0);
	if (!ugxWriter.write_to_file(outFileName.c_str()))
		UG_THROW("Grid could not be written to file '" << outFileName << "'.");

	Domain3d dom;
	try {LoadDomain(dom, outFileName.c_str());}
	UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");

	std::string curFileName("testNeuriteProjector.ugx");
	number offset = 2;
	try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
	UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

	GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
	for (size_t i = 0; i < 2; ++i)
	{
		ref.refine();
		std::ostringstream oss;
		oss << "_refined_" << i+1 << ".ugx";
		curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
		try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
		UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
	}
	}

	////////////////////////////////////////////////////////////////////////
	/// test_import_swc_1d
	///////////////////////////////////////////////////////////////////////
	void test_import_swc_1d(
		const std::string& fileName,
		number anisotropy,
		size_t numRefs,
		number scale
	)
	{
	// preconditioning
	//test_smoothing(fileName, 5, 1.0, 1.0);

	// read in file to intermediate structure
	std::vector<SWCPoint> vPoints;
	std::vector<SWCPoint> vSomaPoints;
	//std::string fn_noext = FilenameWithoutExtension(fileName);
	//std::string fn_precond = fn_noext + "_precond.swc";
	//import_swc(fn_precond, vPoints, false);
	import_swc(fileName, vPoints, scale);

	// convert intermediate structure to neurite data
	std::vector<std::vector<vector3> > vPos;
	std::vector<std::vector<number> > vRad;
	std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
	std::vector<size_t> vRootNeuriteIndsOut;

	convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

	// prepare grid and projector
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	Selector sel(g);


	typedef NeuriteProjector::SurfaceParams NPSP;
	UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
		"GlobalAttachment 'npSurfParams' not declared.");
	Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	if (!g.has_vertex_attachment(aSP))
		g.attach_to_vertices(aSP);

	Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	aaSurfParams.access(g, aSP);

	UG_COND_THROW(!GlobalAttachments::is_declared("diameter"),
		"GlobalAttachment 'diameter' not declared.");
	Attachment<number> aDiam = GlobalAttachments::attachment<Attachment<number> >("diameter");
	if (!g.has_vertex_attachment(aDiam))
		g.attach_to_vertices(aDiam);

	Grid::VertexAttachmentAccessor<Attachment<number> > aaDiam;
	aaDiam.access(g, aDiam);
	ProjectionHandler projHandler(&sh);
	SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
	projHandler.set_geometry(geom3d);

	SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
	projHandler.set_projector(0, neuriteProj);

	// create spline data
	std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
	create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

	// create coarse grid
	for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i)
		create_neurite_1d(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
			anisotropy, g, aaPos, aaSurfParams, aaDiam, NULL);


	// subsets
	AssignSubsetColors(sh);
	sh.set_subset_name("neurites", 0);

	// export geometry
	std::string outFileName = FilenameWithoutPath(std::string("testNeuriteProjector.ugx"));
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	ugxWriter.add_projection_handler(projHandler, "defPH", 0);
	if (!ugxWriter.write_to_file(outFileName.c_str()))
		UG_THROW("Grid could not be written to file '" << outFileName << "'.");


	// refinements
	Domain3d dom;
	try {LoadDomain(dom, outFileName.c_str());}
	UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");

	std::string curFileName("testNeuriteProjector.ugx");
	number offset = 2.0;
	curFileName = outFileName.substr(0, outFileName.size()-4) + "_refined_0.ugx";
	try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
	UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

	GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
	for (size_t i = 0; i < numRefs; ++i)
	{
		ref.refine();
		std::ostringstream oss;
		oss << "_refined_" << i+1 << ".ugx";
		curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
		try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
		UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
	}
	}

	////////////////////////////////////////////////////////////////////////
	/// test neurite projector with four section tube
	////////////////////////////////////////////////////////////////////////
	void test_neurite_projector_with_four_section_tube()
	{
    // we set up a neuron going
    // (0,0,0) -> (1,0,0) -> (3,1,0) -> (5,1,1) -> (7,0,0)
    // with radii
    //    0.05    ->    0.1  ->    0.2  ->    0.15 ->    0.05

    // first step: calculate corresponding spline data
    std::vector<std::vector<vector3> > vPos(1);
    std::vector<vector3>& pos = vPos[0];
    pos.resize(5);
    pos[0] = vector3(0,0,0);
    pos[1] = vector3(1,0,0);
    pos[2] = vector3(3,1,0);
    pos[3] = vector3(5,1,1);
    pos[4] = vector3(7,0,0);

    std::vector<std::vector<number> > vR(1);
    std::vector<number>& r = vR[0];
    r.resize(5);
    r[0] = 0.05;
    r[1] = 0.1;
    r[2] = 0.2;
    r[3] = 0.15;
    r[4] = 0.05;


    // second step: create coarse grid
    Grid g;
    SubsetHandler sh(g);
    sh.set_default_subset_index(0);
    g.attach_to_vertices(aPosition);
    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
    Selector sel(g);


    typedef NeuriteProjector::SurfaceParams NPSP;
    UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
            "GlobalAttachment 'npSurfParams' not declared.");
    Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
    if (!g.has_vertex_attachment(aSP))
        g.attach_to_vertices(aSP);

    Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
    aaSurfParams.access(g, aSP);


    ProjectionHandler projHandler(&sh);
    SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
    projHandler.set_geometry(geom3d);

    SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
    projHandler.set_projector(0, neuriteProj);

    std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
    create_spline_data_for_neurites(vNeurites, vPos, vR);

    const number anisotropy = 8.0;
    create_neurite_surf(vNeurites, vPos, vR, 0, anisotropy, g, aaPos, aaSurfParams);


    // third step: refinement
    AssignSubsetColors(sh);
    sh.set_subset_name("surf", 0);

    std::string fileName = FilenameWithoutPath(std::string("testNeuriteProjector.ugx"));
    GridWriterUGX ugxWriter;
    ugxWriter.add_grid(g, "defGrid", aPosition);
    ugxWriter.add_subset_handler(sh, "defSH", 0);
    ugxWriter.add_projection_handler(projHandler, "defPH", 0);
    if (!ugxWriter.write_to_file(fileName.c_str()))
        UG_THROW("Grid could not be written to file '" << fileName << "'.");


    Domain3d dom;
    try {LoadDomain(dom, fileName.c_str());}
    UG_CATCH_THROW("Failed loading domain from '" << fileName << "'.");

    GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
    for (size_t i = 0; i < 4; ++i)
    {
        ref.refine();

        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        std::string curFileName = fileName.substr(0, fileName.size()-4) + oss.str();
        number offset = 1;
        try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
        UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
    }
	}

	////////////////////////////////////////////////////////////////////////
	/// test neurite projector with four section tube and branch point
	////////////////////////////////////////////////////////////////////////
	void test_neurite_projector_with_four_section_tube_and_branch_point()
	{
    // we set up a neurite going
    // (0,0,0) -> (1,0,0) -> (3,1,0) -> (5,1,1) -> (7,0,0)
    // with radii
    //    0.05 ->    0.1  ->    0.2  ->    0.15 ->    0.05
    // and branching from that at (3,1,0) another neurite going
    // (3,1,0) -> (3,3,-1) -> (2,3,-2)
    // with radii
    //    0.1  ->    0.15  ->    0.05

    // first step: calculate corresponding spline data
    std::vector<std::vector<vector3> > vPos(2);
    std::vector<vector3>& pos0 = vPos[0];
    pos0.resize(5);
    pos0[0] = vector3(0,0,0);
    pos0[1] = vector3(1,0,0);
    pos0[2] = vector3(3,1,0);
    pos0[3] = vector3(5,1,1);
    pos0[4] = vector3(7,0,0);

    std::vector<vector3>& pos1 = vPos[1];
    pos1.resize(3);
    pos1[0] = vector3(3,1,0);
    pos1[1] = vector3(3,3,-1);
    pos1[2] = vector3(2,3,-2);

    std::vector<std::vector<number> > vR(2);
    std::vector<number>& r0 = vR[0];
    r0.resize(5);
    r0[0] = 0.05;
    r0[1] = 0.1;
    r0[2] = 0.2;
    r0[3] = 0.15;
    r0[4] = 0.05;

    std::vector<number>& r1 = vR[1];
    r1.resize(3);
    r1[0] = 0.1;
    r1[1] = 0.15;
    r1[2] = 0.05;

    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > bpInfo(2);
    bpInfo[0].resize(1);
    bpInfo[0].begin()->first = 2;   // at this vertex in neurite 0 starts another neurite (or more)
    bpInfo[0].begin()->second = std::vector<size_t>(1,1);  // these neurites start here

    // second step: create coarse grid
    Grid g;
    SubsetHandler sh(g);
    sh.set_default_subset_index(0);
    g.attach_to_vertices(aPosition);
    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
    Selector sel(g);

    typedef NeuriteProjector::SurfaceParams NPSP;
    UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
            "GlobalAttachment 'npSurfParams' not declared.");
    Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
    if (!g.has_vertex_attachment(aSP))
        g.attach_to_vertices(aSP);

    Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
    aaSurfParams.access(g, aSP);

    ProjectionHandler projHandler(&sh);
    SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
    projHandler.set_geometry(geom3d);

    SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
    projHandler.set_projector(0, neuriteProj);

    std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
    create_spline_data_for_neurites(vNeurites, vPos, vR, &bpInfo);

    const number anisotropy = 8.0;
    create_neurite_surf(vNeurites, vPos, vR, 0, anisotropy, g, aaPos, aaSurfParams, NULL, NULL);

    // at branching points, we have not computed the correct positions yet,
    // so project the complete geometry using the projector
    VertexIterator vit = g.begin<Vertex>();
    VertexIterator vit_end = g.end<Vertex>();
    for (; vit != vit_end; ++vit)
		neuriteProj->project(*vit);

    // third step: refinement
    AssignSubsetColors(sh);
    sh.set_subset_name("surf", 0);

    std::string fileName = FilenameWithoutPath(std::string("testNeuriteProjector.ugx"));
    GridWriterUGX ugxWriter;
    ugxWriter.add_grid(g, "defGrid", aPosition);
    ugxWriter.add_subset_handler(sh, "defSH", 0);
    ugxWriter.add_projection_handler(projHandler, "defPH", 0);
    if (!ugxWriter.write_to_file(fileName.c_str()))
        UG_THROW("Grid could not be written to file '" << fileName << "'.");


    Domain3d dom;
    try {LoadDomain(dom, fileName.c_str());}
    UG_CATCH_THROW("Failed loading domain from '" << fileName << "'.");

    SmartPtr<MultiGrid> mg = dom.grid();

    ProjectionHandler* ph = dynamic_cast<ProjectionHandler*>(dom.refinement_projector().get());
    UG_COND_THROW(!ph, "Refinement projector in domain is not a ProjectionHandler.");
    SmartPtr<RefinementProjector> spRP = ph->projector(0);
    NeuriteProjector* np = dynamic_cast<NeuriteProjector*>(spRP.get());
    UG_COND_THROW(!np, "Refinement projector in projection handler is not a NeuriteProjector.");
    SmartPtr<NeuriteProjector> spNP = SPNULL;
    spNP = SmartPtr<NeuriteProjector>(np, spRP.refcount_ptr());
    spNP->set_geometry(dom.geometry3d());
    HangingNodeRefiner_MultiGrid ref(*dom.grid(), dom.refinement_projector());
    SmartPtr<NeuriteRefMarkAdjuster> nrma(new NeuriteRefMarkAdjuster(spNP, dom.subset_handler(),
                                                                     dom.position_accessor()));
    add_neurite_ref_mark_adjuster(&ref, nrma);

    for (size_t i = 0; i < 6; ++i)
    {
        // only 3 anisotropic refinements
        if (i == 3) nrma->disable();

        int topLv = mg->num_levels() - 1;
        FaceIterator it = dom.grid()->begin<Face>(topLv);
        FaceIterator it_end = dom.grid()->end<Face>(topLv);
        for (; it!= it_end; ++it)
            ref.mark(*it, RM_REFINE);

        UG_DLOGN(NC_TNP, 0, "refinement step " << i);
        ref.refine();
        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        std::string curFileName = fileName.substr(0, fileName.size()-4) + oss.str();
        number offset = 1;
        try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
        UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
    }
	}

	////////////////////////////////////////////////////////////////////////////
	/// test_import_swc_general
	////////////////////////////////////////////////////////////////////////////
	void test_import_swc_general(
		const std::string& fileName,
		bool correct,
		number erScaleFactor,
		bool withER,
		number anisotropy,
		size_t numRefs,
		bool regularize) {
	{
		// read in file to intermediate structure
		std::vector<SWCPoint> vPoints;
		std::vector<SWCPoint> vSomaPoints;
		std::string fn_noext = FilenameWithoutExtension(fileName);
		// presmooth
		test_smoothing(fileName, 5, 1.0, 1.0);
		std::string fn_precond = fn_noext + "_precond.swc";

		///std::string fn_precond = fn_noext + ".swc";
	    std::string fn_precond_with_soma = fn_noext + "_precond_with_soma.swc";
		import_swc_old(fn_precond, vPoints, correct, 1.0);

		std::vector<ug::vector3> vSurfacePoints;
		std::vector<ug::vector3> vPosSomaClosest;
	    size_t lines;
	    get_closest_points_to_soma(fn_precond, vPosSomaClosest, lines);
	    UG_DLOGN(NC_TNP, 0, "Found #" << vPosSomaClosest.size() << " soma points.");

		// convert intermediate structure to neurite data
		std::vector<std::vector<vector3> > vPos;
		std::vector<std::vector<number> > vRad;
		std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
		std::vector<size_t> vRootNeuriteIndsOut;
	    std::vector<ug::vector3> vPointSomaSurface;
		convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);
	    std::vector<SWCPoint> somaPoint = vSomaPoints;
	    std::vector<SWCPoint> savedSomaPoint = vSomaPoints;

		// prepare grid and projector
		Grid g;
		SubsetHandler sh(g);
		sh.set_default_subset_index(0);
		g.attach_to_vertices(aPosition);
		Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

		Selector sel(g);
	    typedef NeuriteProjector::SurfaceParams NPSP;
	    UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
	            "GlobalAttachment 'npSurfParams' not declared.");
	    Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	    if (!g.has_vertex_attachment(aSP))
	        g.attach_to_vertices(aSP);
	    Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	    aaSurfParams.access(g, aSP);

		/// FIXME: In new implementation radius should be possible to scale with
	    /// 1.0 not 1.05 -> might create intersections which lead to segfault
	    somaPoint[0].radius *= 1.05;
	    UG_DLOGN(NC_TNP, 0, "Creating soma...");
	    create_soma(somaPoint, g, aaPos, sh, 1);
	    UG_DLOGN(NC_TNP, 0, " done.");
	    /// get_closest_points_on_soma(vPosSomaClosest, vPointSomaSurface, g, aaPos, sh, 1);
	    std::vector<ug::Vertex*> vPointSomaSurface2;
	    get_closest_vertices_on_soma(vPosSomaClosest, vPointSomaSurface2, g, aaPos, sh, 1);
	    UG_DLOGN(NC_TNP, 0, "Got closest points on soma with size: " << vPointSomaSurface.size());
	    /// add_soma_surface_to_swc(lines, fn_precond, fn_precond_with_soma, vPointSomaSurface);
	    std::vector<ug::vector3> newVerts = find_quad_verts_on_soma(g, aaPos, vPointSomaSurface2, vRad, 1, sh, 1.0, vPos.size());
	    replace_first_root_neurite_vertex_in_swc(lines, fn_precond, fn_precond_with_soma, newVerts);
	    UG_DLOGN(NC_TNP, 0, "Replaced soma points for neurites to SWC file.")
	    g.clear_geometry();
	    import_swc_old(fn_precond_with_soma, vPoints, correct, 1.0);

	    /// TODO: Smooth before or after regularize?
    	Grid g3;
    	SubsetHandler sh3(g3);
    	swc_points_to_grid(vPoints, g3, sh3);
    	export_to_ugx(g3, sh3, "before_regularize.ugx");
    	UG_DLOGN(NC_TNP, 0, "Regularizing branching points...")
    	/// smoothing(vPoints, n, h, gamma);
    	regularize_bps(vPoints, regularize);
    	Grid g2;
    	SubsetHandler sh2(g2);
    	swc_points_to_grid(vPoints, g2, sh2);
    	export_to_ugx(g2, sh2, "after_regularize.ugx");
	   	/// constrained_smoothing(vPoints, vRootNeuriteIndsOut.size(), 0.1, 0.1, 10, 0.1);

	    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

		SubsetHandler psh(g);
		psh.set_default_subset_index(0);

		ProjectionHandler projHandler(&psh);
		SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
		projHandler.set_geometry(geom3d);

		SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
		projHandler.set_projector(0, neuriteProj);

		// create spline data
		std::vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
		create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

		// assign soma region
		std::vector<NeuriteProjector::SomaRegion> somaRegions;
		somaRegions.push_back(NeuriteProjector::SomaRegion(0));
		somaRegions.push_back(NeuriteProjector::SomaRegion(-1));
		for (size_t i = 0; i < vNeurites.size(); i++) {
			vNeurites[i].vSR = somaRegions;
		}

		std::vector<Vertex*> outVerts;
		std::vector<number> outRads;
		std::vector<Vertex*> outVertsInner;
		std::vector<number> outRadsInner;

	    UG_DLOGN(NC_TNP, 0, "Generating neurites...");
        ///for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
	    for (size_t i = 0; i < 1; ++i) {
	    	create_neurite_root_vertices(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
	    		g, sh, erScaleFactor, aaPos, &outVerts, &outVertsInner, &outRads,
	    		&outRadsInner, withER);
	    }

	    UG_DLOGN(NC_TNP, 0, " done.");
	    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites.ugx");

	    UG_LOGN("neurites added");

	    /// Outer soma
	    // Note: axisVectors (outer soma) and axisVectorsInner (inner soma) save the
	    // cylinder center, diameter and length parameters for the CylinderProjectors
	    UG_DLOG(NC_TNP, 0, "Creating soma...");
	    sh.set_default_subset_index(4); /// soma starts now at 4
	    somaPoint = vSomaPoints;
	    create_soma(somaPoint, g, aaPos, sh, 4, 3);
	    UG_DLOGN(NC_TNP, 0, " done.");
	    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_testNeuriteProjector_after_adding_neurites_and_connecting_inner_soma_to_outer_ER.ugxafter_adding_neurites_and_soma.ugx");
	    std::vector<Vertex*> outQuadsInner;
	    std::vector<std::pair<size_t, std::pair<ug::vector3, ug::vector3> > > axisVectors;
	    std::vector<std::vector<ug::Vertex*> > connectingVertices(vRootNeuriteIndsOut.size());
	   	std::vector<std::vector<ug::Vertex*> > connectingVerticesInner(vRootNeuriteIndsOut.size());
	   	std::vector<std::vector<ug::Edge*> > connectingEdges(vRootNeuriteIndsOut.size());
	    std::vector<std::vector<ug::Edge*> > connectingEdgesInner(vRootNeuriteIndsOut.size());

	    /// connect soma with neurites (actually just finds the surface quads on the outer soma)
	    connect_neurites_with_soma(g, aaPos, aaSurfParams, outVerts, outVertsInner, outRads, outQuadsInner,
	    		4, sh, fileName, 1.0, axisVectors, vNeurites, connectingVertices, connectingVerticesInner,
	    		connectingEdges, connectingEdgesInner, true, 0.01, 10, 0.00001, 0.5, 12, 1);

	    ///UG_LOGN("Connected neurites with soma");

	    /// A test for new strategy
	    /*
	   connect_neurites_with_soma_var(g, aaPos, aaSurfParams, outVerts, outVertsInner, outRads, outQuadsInner,
	    	    		4, sh, fileName, erScaleFactor, axisVectors, vNeurites, connectingVertices, connectingVerticesInner,
	    	    		connectingEdges, connectingEdgesInner, true, 0.01, 10, 0.00001, 0.5, 12, 1);
	    	    		*/



	    /// delete old vertices from incorrect neurite starts
	    g.erase(outVerts.begin(), outVerts.end());
	    g.erase(outVertsInner.begin(), outVertsInner.end());
	    outVerts.clear();
	    outVertsInner.clear();
	    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_finding_initial_edges.ugx");
	    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_finding_initial_edges.ugx");

	    sh.set_default_subset_index(0);
	    UG_LOGN("Generating neurites...")
	    UG_DLOGN(NC_TNP, 0, "Generating neurites...")
	    ///for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
	    for (size_t i = 0; i < 1; ++i) {
	    	if (withER) {
	    		create_neurite_with_er(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
	    			erScaleFactor, anisotropy, g, aaPos, aaSurfParams, sh, &outVerts,
	    			&outVertsInner, &outRads, &outRadsInner);
	    	} else {
	    		create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
	    				anisotropy, g, aaPos, aaSurfParams);
	    	}
	    }
	    UG_DLOGN(NC_TNP, 0, " done.");
	    UG_LOGN("Generating inner soma");

	    /// Inner soma
	    UG_DLOGN(NC_TNP, 0, "Creating inner soma...")
	    number scaleER = 0.5;
	    somaPoint.front().radius = somaPoint.front().radius * scaleER;
	    size_t newSomaIndex = sh.num_subsets();
	    create_soma(somaPoint, g, aaPos, sh, newSomaIndex, 3);
	    UG_DLOGN(NC_TNP, 0, " done.");
	    std::vector<Vertex*> outQuadsInner2;
	    UG_DLOGN(NC_TNP, 0, "Size of outQuadsInner: " << outQuadsInner.size())
	    std::vector<std::pair<size_t, std::pair<ug::vector3, ug::vector3> > > axisVectorsInner;
	    /// connect inner soma with neurites (actually just finds/creates the surface quads on the inner soma)
	    connect_neurites_with_soma(g, aaPos, aaSurfParams, outVerts, outVertsInner, outRadsInner, outQuadsInner2, newSomaIndex, sh, fileName, scaleER, axisVectorsInner, vNeurites, connectingVertices, connectingVerticesInner, connectingEdges, connectingEdgesInner, false);
	    UG_LOGN("Connected inner soma");
	    SaveGridToFile(g, sh, "testNeuriteProjector_after_finding_surface_quads.ugx");

	    connect_new(g, sh, aaPos, newSomaIndex, 1, aaSurfParams);
	    SaveGridToFile(g, sh, "testNeuriteProjector_after_finding_surface_quads_and_connect_new.ugx");
	    return;

		// assign subset
		AssignSubsetColors(sh);
		sh.set_subset_name("cyt", 0);
		sh.set_subset_name("er", 1);
		sh.set_subset_name("pm", 2);
		sh.set_subset_name("erm", 3);
		sh.set_subset_name("soma (outer)", 4);
		// ER which should be always subset index 5 (newSomaIndex can be replace with 5?)
		sh.set_subset_name("soma (inner)", newSomaIndex);

		/// Note: Neurite connections get their names respectively subset assignment above in the connect methods
		for (int i = newSomaIndex+1; i < sh.num_subsets(); i++) {
		  	std::stringstream ss;
		   	ss << "inner-connex #" << i;
		   	sh.set_subset_name(ss.str().c_str(), i);
		 }
		IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_renaming.ugx");
		UG_LOGN("After inner connex");
		SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_renaming.ugx");


		/// Double Vertices might occur during Qhull gen faces -> remove these here
		/// RemoveDoubles<3>(g, g.begin<Vertex>(), g.end<Vertex>(), aaPos, 0.0001);
		IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_renaming_and_double_freed.ugx");

	    EraseEmptySubsets(sh);
	    AssignSubsetColors(sh);
	    for (size_t i = newSomaIndex+vRootNeuriteIndsOut.size(); i < (size_t)sh.num_subsets(); i++) {
	    	std::stringstream ss;
	     	ss << "inter-soma-connex #" << i;
	     	sh.set_subset_name(ss.str().c_str(), i);
	    }

	    /// TODO: Debug from here for segfault during connect method...
	    /// connect now inner soma to neurites (Note: Old strategy which uses distance not angle based criterion. Change this?)
	    /// TODO: how to find surface quads to connect: AdaptSurface destroys good properties somehow -> Verdrehung
	    UG_DLOGN(NC_TNP, 0, "Soma's inner index: " << newSomaIndex);
	    UG_LOGN("somaIndex: " << newSomaIndex);
	    if (withER) {
	    	connect_inner_neurites_to_inner_soma(newSomaIndex, 1, g, aaPos, sh, aaSurfParams, erScaleFactor);
	    }
	    SavePreparedGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_inner_soma_to_outer_ER.ugx");

		UG_LOGN("After adding neurites and connecting inner soma to outer ER");

	    /// Note: Called two times for inner and outer polygon on outer soma but use the same plane defined by the outer soma's inner quad vertices
	    connect_outer_and_inner_root_neurites_to_outer_soma_variant(4, vRootNeuriteIndsOut.size(), g, aaPos, sh, outVerts, 12, true);

	    /// Extrude ER volume a little bit further into normal direction towards inner soma, like the pyramids to close outer soma, to avoid intersections
	    if (withER) {
	    	extend_ER_within(g, sh, aaPos, aaSurfParams, newSomaIndex, 1, erScaleFactor, outVertsInner);
	    	SavePreparedGridToFile(g, sh, "after_extend_ER_and_before_connect_outer.ugx");
	    	connect_outer_and_inner_root_neurites_to_outer_soma_variant(4, vRootNeuriteIndsOut.size(), g, aaPos, sh, outVertsInner, 4, true);
	    	EraseEmptySubsets(sh);
	    	AssignSubsetColors(sh);
	    }

	    /// Reassign elements for connecting parts ER and somata to erm subset
	    sel.clear();
	    SelectSubset(sel, sh, 3, true);
	    CloseSelection(sel);
	    AssignSelectionToSubset(sel, sh, 3);
	    SavePreparedGridToFile(g, sh, "before_tetrahedralize_and_after_reassigned.ugx");

		/// assign correct axial parameters for "somata"
		set_somata_axial_parameters(g, sh, aaSurfParams, 4, 5);

	    /// tetrahedralizes somata with specified and fixed indices 4 and 5
	    tetrahedralize_soma(g, sh, aaPos, aaSurfParams, 4, 5, savedSomaPoint);

	    /// reassign soma volumes to appropriate subsets
	    reassign_volumes(g, sh, 4, 5, scaleER, savedSomaPoint[0], aaPos);

	    /// assign SurfParams for all vertices from tetrahedralize call (now in subset 4 and 5)
	    fix_axial_parameters(g, sh, aaSurfParams, aaPos, 4, 5, savedSomaPoint[0], scaleER);

	    UG_LOGN("After tetrahedralize");

		SavePreparedGridToFile(g, sh, "after_tetrahedralize_soma.ugx");
		/// After merge doubles might occur, delete them. Boundary faces are retained,
		/// however triangles occur now at boundary interface and quadrilaterals, thus
		/// delete the triangles to keep the quadrilaterals from the start of neurites
		RemoveDoubles<3>(g, g.begin<Vertex>(), g.end<Vertex>(), aaPos, 0.00001);
		SavePreparedGridToFile(g, sh, "after_tetrahedralize_soma_and_removed_doubles.ugx");
		DeleteInnerEdgesFromQuadrilaterals(g, sh, 4);
		DeleteInnerEdgesFromQuadrilaterals(g, sh, 1);
		SavePreparedGridToFile(g, sh, "after_tetrahedralize_soma_and_conversion.ugx");

		for (int i = 0; i <= sh.num_subsets(); i++) {
			for (VertexIterator iter = sh.begin<Vertex>(i); iter != sh.end<Vertex>(i); iter++) {
				UG_DLOGN(NC_TNP, 0, "attachment value (aSP) for subset " << i << ": " << aaSurfParams[*iter]);
			}
		}

	    // at branching points, we have not computed the correct positions yet,
		// so project the complete geometry using the projector
		VertexIterator vit = g.begin<Vertex>();
		VertexIterator vit_end = g.end<Vertex>();
		for (; vit != vit_end; ++vit) {
			neuriteProj->project(*vit);
		}

	    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_all.ugx");
	    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_all.ugx");

		// output
		std::string outFileNameBase = FilenameAndPathWithoutExtension(fileName);
		std::string outFileName = outFileNameBase + ".ugx";
		GridWriterUGX ugxWriter;
		ugxWriter.add_grid(g, "defGrid", aPosition);
		ugxWriter.add_subset_handler(sh, "defSH", 0);
		ugxWriter.add_subset_handler(psh, "projSH", 0);
		ugxWriter.add_projection_handler(projHandler, "defPH", 0);
		if (!ugxWriter.write_to_file(outFileName.c_str()))
			UG_THROW("Grid could not be written to file '" << outFileName << "'.");

		if (numRefs == 0)
			return;

		// refinement
		Domain3d dom;
		dom.create_additional_subset_handler("projSH");
		try {LoadDomain(dom, outFileName.c_str());}
		UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");

		number offset = 5.0;
		std::string curFileName = outFileNameBase + "_refined_0.ugx";

		try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
		UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

		GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
		for (uint i = 0; i < numRefs; ++i)
		{
			ref.refine();
			std::ostringstream oss;
			oss << "_refined_" << i+1 << ".ugx";
			curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
			try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
			UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
		}
	}
	}


	////////////////////////////////////////////////////////////////////////////
	/// test_import_swc_general
	////////////////////////////////////////////////////////////////////////////
	void test_import_swc_general_var(
		const std::string& fileName,
		bool correct,
		number erScaleFactor,
		bool withER,
		number anisotropy,
		size_t numRefs,
		bool regularize
	) {
		using namespace std;
		// read in file to intermediate structure
		vector<SWCPoint> vPoints;
		vector<SWCPoint> vSomaPoints;
		string fn_noext = FilenameWithoutExtension(fileName);

		// presmooth whole structure
		test_smoothing(fileName, 5, 1.0, 1.0);
		string fn_precond = fn_noext + "_precond.swc";

		// import preconditioned swc structure
		std::string fn_precond_with_soma = fn_noext + "_precond_with_soma.swc";
		import_swc_old(fn_precond, vPoints, correct, 1.0);

		// this finds the closest points towards soma
		vector<vector3> vSurfacePoints;
		vector<vector3> vPosSomaClosest;
		size_t lines;
	    get_closest_points_to_soma(fn_precond, vPosSomaClosest, lines);
	    UG_DLOGN(NC_TNP, 0, "We found #" << vPosSomaClosest.size() << " soma " <<
	    		"points. This means there are " << vPosSomaClosest.size() <<
	    		"neurites starting from soma");

		// convert intermediate structure to neurite data
		vector<vector<vector3> > vPos;
		vector<vector<number> > vRad;
		vector<vector<pair<size_t, vector<size_t> > > > vBPInfo;
		vector<size_t> vRootNeuriteIndsOut;
		vector<vector3> vPointSomaSurface;
		convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);
	    vector<SWCPoint> somaPoint = vSomaPoints;
	    vector<SWCPoint> savedSomaPoint = vSomaPoints;

		// prepare grid and projector
		Grid g;
		SubsetHandler sh(g);
		sh.set_default_subset_index(0);
		g.attach_to_vertices(aPosition);
		Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

		Selector sel(g);
		typedef NeuriteProjector::SurfaceParams NPSP;
		UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
		              "GlobalAttachment 'npSurfParams' not declared.");
		Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
		if (!g.has_vertex_attachment(aSP)) {
			g.attach_to_vertices(aSP);
		}

		Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
		aaSurfParams.access(g, aSP);

		/// TODO: Check that this new scaling does not interfere with grid generation, then remove this noop
		somaPoint[0].radius *= 1.00;
		UG_DLOGN(NC_TNP, 0, "Creating (outer sphere) soma in subset 1");
	    create_soma(somaPoint, g, aaPos, sh, 1);

	    // Get closest _vertices_ on soma surface
		vector<Vertex*> vPointSomaSurface2;
		get_closest_vertices_on_soma(vPosSomaClosest, vPointSomaSurface2, g, aaPos, sh, 1);
	    UG_DLOGN(NC_TNP, 0, "Got " << "# " << vPointSomaSurface.size() << "closest vertices on (outer sphere) soma");

		// add_soma_surface_to_swc(lines, fn_precond, fn_precond_with_soma, vPointSomaSurface);
	    // Replace first vertex of root neurite with soma surface vertex
	    std::vector<ug::vector3> newVerts = find_quad_verts_on_soma(g, aaPos, vPointSomaSurface2, vRad, 1, sh, 1.0, vPos.size());
	    replace_first_root_neurite_vertex_in_swc(lines, fn_precond, fn_precond_with_soma, newVerts);
	    UG_DLOGN(NC_TNP, 0, "Replaced soma points for root neurites to SWC file.")
	    g.clear_geometry();

	    // re-read the now correct SWC file with soma
	    import_swc_old(fn_precond_with_soma, vPoints, correct, 1.0);

	    // Regularize, smooth possibly again, then convert to neuritelist
	    /// TODO: Smooth before or after regularize?
	    regularize_bps(vPoints, regularize);
		/// constrained_smoothing(vPoints, vRootNeuriteIndsOut.size(), 0.1, 0.1, 10, 0.1);
	    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

	    // Projection handling
		SubsetHandler psh(g);
		psh.set_default_subset_index(0);
		ProjectionHandler projHandler(&psh);
		SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
		projHandler.set_geometry(geom3d);
		SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
		projHandler.set_projector(0, neuriteProj);

		// create spline data
		vector<NeuriteProjector::Neurite>& vNeurites = neuriteProj->neurites();
		create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

		// assign soma regions for later projection refinement
		vector<NeuriteProjector::SomaRegion> somaRegions;
		somaRegions.push_back(NeuriteProjector::SomaRegion(0));
		somaRegions.push_back(NeuriteProjector::SomaRegion(-1));
		for (size_t i = 0; i < vNeurites.size(); i++) {
			vNeurites[i].vSR = somaRegions;
		}


		// helper vectors to store radii and verts for soma/neurite connection
		vector<Vertex*> outVerts;
		vector<number> outRads;
		vector<Vertex*> outVertsInner;
		vector<number> outRadsInner;

	    UG_DLOGN(NC_TNP, 0, "Generating neurites...");
	    /// for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
	    for (size_t i = 0; i < 1; ++i) {
		   	create_neurite_root_vertices(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
		   		g, sh, erScaleFactor, aaPos, &outVerts, &outVertsInner, &outRads,
		   		&outRadsInner, withER);
		 }

		 UG_DLOGN(NC_TNP, 0, " done.");
		 IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites.ugx");

	    /// (Outer sphere) Soma
		// Note: axisVectors (outer soma) and axisVectorsInner (inner soma) save the
		// cylinder center, diameter and length parameters for the CylinderProjectors
	    UG_DLOG(NC_TNP, 0, "Creating soma...");
	    sh.set_default_subset_index(4); /// soma starts now at 4
	    somaPoint = vSomaPoints;
	    create_soma(somaPoint, g, aaPos, sh, 4, 3);
	    UG_DLOGN(NC_TNP, 0, " done.");
	    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_testNeuriteProjector_after_adding_neurites_and_connecting_inner_soma_to_outer_ER.ugxafter_adding_neurites_and_soma.ugx");
	    vector<Vertex*> outQuadsInner;
	    vector<pair<size_t, pair<vector3, vector3> > > axisVectors;
	    vector<vector<Vertex*> > connectingVertices(vRootNeuriteIndsOut.size());
	   	vector<vector<Vertex*> > connectingVerticesInner(vRootNeuriteIndsOut.size());
	   	vector<vector<Edge*> > connectingEdges(vRootNeuriteIndsOut.size());
	    vector<vector<Edge*> > connectingEdgesInner(vRootNeuriteIndsOut.size());

	    /// TODO: Remove obsolete parameters
		connect_neurites_with_soma_var(g, aaPos, aaSurfParams, outVerts, outVertsInner, outRads, outQuadsInner,
		    	    		4, sh, fileName, erScaleFactor, axisVectors, vNeurites, connectingVertices, connectingVerticesInner,
		    	    		connectingEdges, connectingEdgesInner, true, 0.01, 10, 0.00001, 0.5, 12, 1);

		UG_DLOGN(NC_TNP, 0, "Found (outer sphere) soma dodecagons to connect neurite with");

	    // delete old vertices from incorrect neurite starts
	    g.erase(outVerts.begin(), outVerts.end());
	    g.erase(outVertsInner.begin(), outVertsInner.end());
		outVerts.clear();
		outVertsInner.clear();
		IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_finding_initial_edges.ugx");
		SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_finding_initial_edges.ugx");

	    sh.set_default_subset_index(0);
	    UG_LOGN("Generating neurites...")
	    UG_DLOGN(NC_TNP, 0, "Generating neurites...")
		/// for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
		for (size_t i = 0; i < 1; ++i) {
		 	if (withER) {
		   		create_neurite_with_er(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
		   			erScaleFactor, anisotropy, g, aaPos, aaSurfParams, sh, &outVerts,
		   			&outVertsInner, &outRads, &outRadsInner);
		   	} else {
		   		create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i],
		   				anisotropy, g, aaPos, aaSurfParams);
		   	}
		}
	    UG_DLOGN(NC_TNP, 0, " done.");
	    UG_LOGN("Generating inner soma");

	    /// (Inner sphere) ER
	    UG_DLOGN(NC_TNP, 0, "Creating (inner sphere) ER");
	    /// TODO: Make this a parameter again
	    number scaleER = 0.5;
	    somaPoint.front().radius = somaPoint.front().radius * scaleER;
	    size_t newSomaIndex = sh.num_subsets();
	    create_soma(somaPoint, g, aaPos, sh, newSomaIndex, 3);
	    UG_DLOGN(NC_TNP, 0, " done.");
	    vector<Vertex*> outQuadsInner2;
		UG_DLOGN(NC_TNP, 0, "Size of outQuadsInner: " << outQuadsInner.size())
		vector<pair<size_t, pair<vector3, vector3> > > axisVectorsInner;
		// Find surface quads on the (inner sphere) ER to connect with to ER
		connect_neurites_with_soma(g, aaPos, aaSurfParams, outVerts, outVertsInner, outRadsInner, outQuadsInner2, newSomaIndex, sh, fileName, scaleER, axisVectorsInner, vNeurites, connectingVertices, connectingVerticesInner, connectingEdges, connectingEdgesInner, false);
		UG_LOGN("Connected inner soma");
		SaveGridToFile(g, sh, "testNeuriteProjector_after_finding_surface_quads.ugx");

		/// Connects (inner sphere) ER with ER part of dendrite
		if (withER) {
			connect_new(g, sh, aaPos, newSomaIndex, 1, aaSurfParams);
		}
		SaveGridToFile(g, sh, "testNeuriteProjector_after_finding_surface_quads_and_connect_new.ugx");
		///return;

		// assign subset
		AssignSubsetColors(sh);
		sh.set_subset_name("cyt", 0);
		sh.set_subset_name("er", 1);
		sh.set_subset_name("pm", 2);
		sh.set_subset_name("erm", 3);
		sh.set_subset_name("soma (outer)", 4);
		// ER which should be always subset index 5 (newSomaIndex can be replace with 5?)
		sh.set_subset_name("soma (inner)", newSomaIndex);

		/// Note: Neurite connections get their names respectively subset assignment above in the connect methods
		for (int i = newSomaIndex+1; i < sh.num_subsets(); i++) {
		  	stringstream ss;
		   	ss << "inner-connex #" << i;
		   	sh.set_subset_name(ss.str().c_str(), i);
		 }
		IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_renaming.ugx");
		UG_LOGN("After inner connex");
		SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_renaming.ugx");

		/// Double Vertices might occur during Qhull gen faces -> remove these here
		/// RemoveDoubles<3>(g, g.begin<Vertex>(), g.end<Vertex>(), aaPos, 0.0001);
		IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_renaming_and_double_freed.ugx");

	    EraseEmptySubsets(sh);
	    AssignSubsetColors(sh);
	    for (size_t i = newSomaIndex+vRootNeuriteIndsOut.size(); i < (size_t)sh.num_subsets(); i++) {
	    	stringstream ss;
		   	ss << "inter-soma-connex #" << i;
		   	sh.set_subset_name(ss.str().c_str(), i);
		 }

		SavePreparedGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_inner_soma_to_outer_ER.ugx");
		UG_LOGN("After adding neurites and connecting inner soma to outer ER");

		UG_LOGN("somaIndex before connect pm with soma: " << newSomaIndex);
		std::vector<std::vector<ug::Vertex*> > outVertsClean;
		for (size_t i = 0; i < outVerts.size() / 12; i++) {
			outVertsClean.push_back(std::vector<ug::Vertex*>(outVerts.begin() + (i*12), outVerts.begin() + (i+1)*12));
		}

		// This connects the outer dodecagon with the soma surface (plasma membrane)
		// The indices start from 5 (4 is outer sphere) up to 5+numRootNeurites (5+numRootNeurites+1 is inner sphere)
	    /// connect_outer_and_inner_root_neurites_to_outer_soma_variant(4, vRootNeuriteIndsOut.size(), g, aaPos, sh, outVerts, 12, true);
		connect_pm_with_soma(newSomaIndex, g, aaPos, sh, outVertsClean, false, 0);
		SaveGridToFile(g, sh, "after_connect_pm_with_soma.ugx");


	    // This connects the inner quad with the soma surface (ER)
	    /// Extrude ER volume a little bit further into normal direction towards
		/// inner soma, like the pyramids to close outer soma, to avoid intersections
	    if (withER) {
	    	/// TODO: Check that this still works, probably it will:
	    	/// TODO: Check also that aaSurfParams are set accordingly for neurite projection later
	    	extend_ER_within(g, sh, aaPos, aaSurfParams, newSomaIndex, 1, erScaleFactor, outVertsInner);
	    	SavePreparedGridToFile(g, sh, "after_extend_ER_and_before_connect_outer.ugx");
			std::vector<std::vector<ug::Vertex*> > outVertsInnerClean;
			for (size_t i = 0; i < outVertsInner.size() / 4; i++) {
				outVertsInnerClean.push_back(std::vector<ug::Vertex*>(outVertsInner.begin() + (i*4), outVertsInner.begin() + (i+1)*4));
			}
			UG_LOGN("newSomaIndex (er with er): " << newSomaIndex);
	    	// The indices start from 5+2*numRootNeurites+1 (newSomaIndex is 5+numRootNeurites) up to 5+3*numRootNeurites
	    	/// connect_outer_and_inner_root_neurites_to_outer_soma_variant(4, vRootNeuriteIndsOut.size(), g, aaPos, sh, outVertsInner, 4, true);
			size_t numQuads = outVertsInnerClean.size();
	    	connect_er_with_er(newSomaIndex+2*numQuads, g, aaPos, sh, outVertsInnerClean, 2*numQuads-1, true);
	    	SaveGridToFile(g, sh, "after_connect_er_with_er.ugx");
	    	EraseEmptySubsets(sh);
	    	AssignSubsetColors(sh);
	    }

	    /// Reassign elements for connecting parts ER and somata to erm subset
	    sel.clear();
	    SelectSubset(sel, sh, 3, true);
	    CloseSelection(sel);
	    AssignSelectionToSubset(sel, sh, 3);
	    SavePreparedGridToFile(g, sh, "before_tetrahedralize_and_after_reassigned.ugx");

		/// TODO: tested until here
	    return;

		/// assign correct axial parameters for "somata"
		set_somata_axial_parameters(g, sh, aaSurfParams, 4, 5);

	    /// tetrahedralizes somata with specified and fixed indices 4 and 5
	    tetrahedralize_soma(g, sh, aaPos, aaSurfParams, 4, 5, savedSomaPoint);

	    /// reassign soma volumes to appropriate subsets
	    reassign_volumes(g, sh, 4, 5, scaleER, savedSomaPoint[0], aaPos);

	    /// assign SurfParams for all vertices from tetrahedralize call (now in subset 4 and 5)
	    fix_axial_parameters(g, sh, aaSurfParams, aaPos, 4, 5, savedSomaPoint[0], scaleER);

	    UG_LOGN("After tetrahedralize");

		SavePreparedGridToFile(g, sh, "after_tetrahedralize_soma.ugx");
		/// After merge doubles might occur, delete them. Boundary faces are retained,
		/// however triangles occur now at boundary interface and quadrilaterals, thus
		/// delete the triangles to keep the quadrilaterals from the start of neurites
		RemoveDoubles<3>(g, g.begin<Vertex>(), g.end<Vertex>(), aaPos, 0.00001);
		SavePreparedGridToFile(g, sh, "after_tetrahedralize_soma_and_removed_doubles.ugx");
		DeleteInnerEdgesFromQuadrilaterals(g, sh, 4);
		DeleteInnerEdgesFromQuadrilaterals(g, sh, 1);
		SavePreparedGridToFile(g, sh, "after_tetrahedralize_soma_and_conversion.ugx");

		for (int i = 0; i <= sh.num_subsets(); i++) {
			for (VertexIterator iter = sh.begin<Vertex>(i); iter != sh.end<Vertex>(i); iter++) {
				UG_DLOGN(NC_TNP, 0, "attachment value (aSP) for subset " << i << ": " << aaSurfParams[*iter]);
			}
		}

	    // at branching points, we have not computed the correct positions yet,
		// so project the complete geometry using the projector
		VertexIterator vit = g.begin<Vertex>();
		VertexIterator vit_end = g.end<Vertex>();
		for (; vit != vit_end; ++vit) {
			neuriteProj->project(*vit);
		}

	    IF_DEBUG(NC_TNP, 0) SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_all.ugx");
	    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_all.ugx");

		// output
		string outFileNameBase = FilenameAndPathWithoutExtension(fileName);
		string outFileName = outFileNameBase + ".ugx";
		GridWriterUGX ugxWriter;
		ugxWriter.add_grid(g, "defGrid", aPosition);
		ugxWriter.add_subset_handler(sh, "defSH", 0);
		ugxWriter.add_subset_handler(psh, "projSH", 0);
		ugxWriter.add_projection_handler(projHandler, "defPH", 0);
		if (!ugxWriter.write_to_file(outFileName.c_str()))
			UG_THROW("Grid could not be written to file '" << outFileName << "'.");

		if (numRefs == 0)
			return;

		// refinement
		Domain3d dom;
		dom.create_additional_subset_handler("projSH");
		try {LoadDomain(dom, outFileName.c_str());}
		UG_CATCH_THROW("Failed loading domain from '" << outFileName << "'.");

		number offset = 5.0;
		string curFileName = outFileNameBase + "_refined_0.ugx";

		try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
		UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

		GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
		for (uint i = 0; i < numRefs; ++i)
		{
			ref.refine();
			ostringstream oss;
			oss << "_refined_" << i+1 << ".ugx";
			curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
			try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
			UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
		}
	}


	////////////////////////////////////////////////////////////////////////////
	/// test_shrinkage
	////////////////////////////////////////////////////////////////////////////
	void test_shrinkage() {
		std::vector<ug::Vertex*> verts;
		Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> > aaSurfParams;
		Grid::VertexAttachmentAccessor<APosition> aaPos;
		correct_axial_offset(verts, aaSurfParams, aaPos, 0.5);
	}
	} // namespace neuro_collection
} // namespace ug
