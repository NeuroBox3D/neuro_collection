/*
 * Test run: ../bin/ugshell -call "test_import_swc_general(\"smith.swc\", false, 0.5, true)"
 * test_neurite_proj.cpp
 *
 *  Created on: 27.12.2016
 *      Author: mbreit
 */

#include "test_neurite_proj.h"
#include "lib_grid/refinement/projectors/projection_handler.h"
#include "../../ugcore/ugbase/lib_grid/refinement/projectors/cylinder_volume_projector.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"
#include "lib_grid/refinement/projectors/cylinder_projector.h"
#include "lib_grid/refinement/projectors/sphere_projector.h"
#include "lib_grid/refinement/projectors/soma_projector.h"
#include "lib_grid/file_io/file_io_ugx.h"  // GridWriterUGX
#include "lib_grid/file_io/file_io.h"  // SaveGridHierarchyTransformed
#include "lib_grid/grid/geometry.h" // MakeGeometry3d
#include "lib_grid/refinement/global_multi_grid_refiner.h" // GlobalMultigridRefiner
#include "lib_grid/global_attachments.h"
#include "lib_grid/algorithms/extrusion/extrusion.h" // Extrude
#include "lib_grid/refinement/regular_refinement.h"  // Refine
#include "lib_grid/algorithms/geom_obj_util/face_util.h" // CalculateNormal
#include "lib_disc/domain_util.h"   // LoadDomain
#include "lib_algebra/small_algebra/small_algebra.h" // Invert
#include "common/util/string_util.h"  // TrimString etc.
#include "neurite_refMarkAdjuster.h"
#include "lib_disc/function_spaces/error_elem_marking_strategy.h" // GlobalMarking
#include "lib_grid/algorithms/element_side_util.h" // GetOpposingSide
#include "lib_disc/quadrature/gauss_legendre/gauss_legendre.h"
#include "lib_grid/algorithms/grid_generation/icosahedron.h" // icosahedron
#include "common/math/ugmath_types.h"
#include "lib_grid/algorithms/subset_color_util.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h" // AdaptSurfaceGridToCylinder
#include "../../ugcore/ugbase/bridge/domain_bridges/selection_bridge.cpp"
#include "lib_grid/algorithms/smoothing/manifold_smoothing.h" // TangentialSmoothing
#include "lib_grid/algorithms/remeshing/resolve_intersections.h" // ResolveTriangleIntersection
#include <boost/geometry.hpp>
#include "../../ElementQualityStatistics/element_quality_statistics.h"

#include <boost/lexical_cast.hpp>

#include <istream>
#include <sstream>
#include <list>
#include <vector>
#include <queue>
#include <stack>
#include <vector>

// configuration file for compile options
#include "config.h"

/// use our implementation quickhull but prefer qhull.org if available
#ifdef NC_WITH_QHULL
	#include "qhull.cpp"
#else
	#include "quickhull.cpp"
#endif

namespace ug {
namespace neuro_collection {
	/**
	 * @brief imports a SWC file
	 */
	void import_swc
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
            case 1: pt.type = SWC_SOMA; break;
            case 2: pt.type = SWC_AXON; break;
            case 3: pt.type = SWC_DEND; break;
            case 4: pt.type = SWC_APIC; break;
            default: pt.type = SWC_UNDF;
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

	/**
	 * @brief smooeths out the SWC file
	 */
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

	/**
	 * @brief function object to compare edge length
	 */
	struct EdgeLengthCompare
{
	bool operator()(const std::pair<Edge*, number> e1, const std::pair<Edge*, number> e2)
	{return e1.second > e2.second;}
};

	/**
	 * @brief collapses too short edges, i.e. smaller then diameter of section
	 */
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

	// Soma (Default subset index 0) should be collapsed into one single point
   	if (sh.num_elements<ug::Edge>(0) != 0) {
   		UG_WARNING("Soma not properly collapsed into one single point. #Edges > 0.");
   	}
}

	/**
	 * @brief SWC list to NeuriteList
	 */
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
	vSomaPoints.clear();
    vPosOut.clear();
    vRadOut.clear();
    vBPInfoOut.clear();
    vRootNeuriteIndsOut.clear();

    // find first soma's root point in geometry and save its index as i
    size_t nPts = vPoints.size();
    size_t i = 0;
    int firstSoma = -1;
    size_t numSomatas = 0;
    for (; i < nPts; ++i) {
        if (vPoints[i].type == SWC_SOMA) {
        	if (firstSoma == -1) {
        		firstSoma = i;
        	}
        	numSomatas++;
        }
    }

    i = firstSoma;
    UG_COND_THROW(i == nPts, "No soma contained in swc point list.")
    UG_COND_THROW(numSomatas > 1, "Too many somatas contained in swc point list.")

    // collect neurite root points
    std::vector<std::pair<size_t, size_t> > rootPts;
    std::queue<std::pair<size_t, size_t> > soma_queue;
    soma_queue.push(std::make_pair((size_t)-1,i));
    while (!soma_queue.empty())
    {
        size_t pind = soma_queue.front().first;
        size_t ind = soma_queue.front().second;
        soma_queue.pop();

        const SWCPoint& pt = vPoints[ind];

        if (pt.type == SWC_SOMA) {
        	vSomaPoints.push_back(vPoints[i]);
            size_t nConn = pt.conns.size();
            for (size_t i = 0; i < nConn; ++i)
                if (pt.conns[i] != pind) {
                    soma_queue.push(std::make_pair(ind, pt.conns[i]));
                }
        } else {
        	rootPts.push_back(std::make_pair(pind, ind));
        }
    }

    vPosOut.resize(rootPts.size());
    vRadOut.resize(rootPts.size());
    vBPInfoOut.resize(rootPts.size());

    std::stack<std::pair<size_t, size_t> > processing_stack;
    for (size_t i = 0; i < rootPts.size(); ++i)
        processing_stack.push(rootPts[i]);

    size_t curNeuriteInd = 0;
    vRootNeuriteIndsOut.push_back(0);

    // helper map to be used to correctly save BPs:
    // maps branch root parent ID to neurite ID and BP ID
    std::map<size_t, std::pair<size_t, size_t> > helperMap;

    while (!processing_stack.empty())
    {
        size_t pind = processing_stack.top().first;
        size_t ind = processing_stack.top().second;
        processing_stack.pop();

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

    size_t numSomaPoints = vSomaPoints.size();
    UG_LOGN("Number of soma points: " << numSomaPoints);
    for (size_t i = 0; i < numSomaPoints; i++) {
    	UG_LOGN("Coordinates for soma point " << i << ": " << vSomaPoints[i].coords);
    }

}

	/**
	 * @brief create spline data operates on NeuriteList to create splines
	 */
	static void create_spline_data_for_neurites
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

        // Note: find suitable permissible render vector
        neuriteOut.refDir = vector3(0,0,1);
        neuriteOut.vSec.reserve(nVrt-1);

        // this will be 0 for root branches and 1 otherwise
        size_t brInd = neuriteOut.vBR.size();
        std::vector<std::pair<size_t, std::vector<size_t> > >::const_iterator brIt;
        std::vector<std::pair<size_t, std::vector<size_t> > >::const_iterator brIt_end;
        if (bpInfo)
        {
            // fill in missing tend info for first BR
            if (brInd)
            {
                number threeRad = 5.0 * r[0];
                number dist = 0.0;
                size_t k = 0;
                while (dist < threeRad && ++k < nVrt)
                    dist += VecDistance(pos[k-1], pos[k]);
                if (dist <= threeRad)
                    neuriteOut.vBR[0].tend = 1.0;
                else
                {
                    neuriteOut.vBR[0].tend = tSuppPos[k] - (dist-threeRad) / VecDistance(pos[k-1], pos[k])
                                                           * (tSuppPos[k] - tSuppPos[k-1]);
                }
            }

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

                // calculate t_start and t_end for parent's BR
                // we surround by approx. 3*radius in all directions
                // Note: This might need some patching up
                //       for branching points with almost parallel child and parent,
                //       where the branching point mollifying effect persists even farther away
                number threeRad = 5.0 * r[i+1];
                number dist = 0.0;
                size_t k = i+2;
                while (dist < threeRad && --k > 0)
                    dist += VecDistance(pos[k-1], pos[k]);
                if (dist <= threeRad)
                    br.tstart = 0.0;
                else
                {
                    br.tstart = tSuppPos[k-1] + (dist-threeRad) / VecDistance(pos[k-1], pos[k])
                                                * (tSuppPos[k] - tSuppPos[k-1]);
                }
                dist = 0.0;
                k = i+1;
                while (dist < threeRad && ++k < nVrt)
                    dist += VecDistance(pos[k-1], pos[k]);
                if (dist <= threeRad)
                    br.tend = 1.0;
                else
                {
                    br.tend = tSuppPos[k] - (dist-threeRad) / VecDistance(pos[k-1], pos[k])
                                                * (tSuppPos[k] - tSuppPos[k-1]);
                }

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

                    // calculate t_start and t_end for child's BR
                    childBR.tstart = 0;
                }
                ++brInd;
                ++brIt;
            }

            neuriteOut.vSec.push_back(sec);
        }
    }

/*
    // debug
    // print out branching region data
    for (size_t n = 0; n < vNeuritesOut.size(); ++n)
    {
        UG_LOGN("Neurite " << n);
        const NeuriteProjector::Neurite& neurite = vNeuritesOut[n];
        const std::vector<NeuriteProjector::BranchingRegion> vBR = neurite.vBR;

        for (size_t b = 0; b < vBR.size(); ++b)
        {
            UG_LOGN("  BR " << b << ": " << vBR[b].tstart << ".." << vBR[b].tend);
            SmartPtr<NeuriteProjector::BranchingPoint> bp = vBR[b].bp;

            UG_LOGN("    associated BP data:")
            size_t bpSz = bp->vNid.size();
            if (bpSz != bp->vRegions.size())
            {
            	UG_LOGN(      "Size mismatch: vNid " << bpSz << ", vRegions " << bp->vRegions.size());
            }
            else
            {
				for (size_t i = 0; i < bpSz; ++i)
				{
					UG_LOGN("      " << bp->vNid[i] << " (" << bp->vRegions[i]->tstart
						<< ", " << bp->vRegions[i]->tend << ")");
				}
            }
        }
    }
*/
}

	/**
	 * @brief helper method to calculate length
	 */
	number calculate_length_over_radius
(
	number t_start,
	number t_end,
	const NeuriteProjector::Neurite& neurite,
	size_t startSec
)
{
	GaussLegendre gl(5);
	size_t nPts = gl.size();

	std::vector<NeuriteProjector::Section>::const_iterator sec_it = neurite.vSec.begin() + startSec;
	std::vector<NeuriteProjector::Section>::const_iterator sec_end = neurite.vSec.end();

	// check that startSec was correct
	number sec_tstart = startSec > 0 ? (sec_it - 1)->endParam : 0.0;
	number sec_tend = sec_it->endParam;

	UG_COND_THROW(sec_tend < t_start || sec_tstart > t_start,
		"Wrong section iterator given to calc_length_over_radius().\n"
		"Section goes from " << (sec_it-1)->endParam << " to " << sec_it->endParam
		<< ", but t_start is " << t_start << ".");

	number integral = 0.0;
	while (sec_it != sec_end)
	{
		// integrate from t_start to min{t_end, sec_tend}
		const NeuriteProjector::Section& sec = *sec_it;
		sec_tend = sec.endParam;
		number dt = std::min(t_end, sec_tend) - t_start;
		number sec_integral = 0.0;
		for (size_t i = 0; i < nPts; ++i)
		{
			number t = sec_tend - (t_start + dt*gl.point(i)[0]);

			vector3 vel;
			const number* s = &sec.splineParamsX[0];
			number& v0 = vel[0];
			v0 = -3.0*s[0]*t - 2.0*s[1];
			v0 = v0*t - s[2];

			s = &sec.splineParamsY[0];
			number& v1 = vel[1];
			v1 = -3.0*s[0]*t - 2.0*s[1];
			v1 = v1*t - s[2];

			s = &sec.splineParamsZ[0];
			number& v2 = vel[2];
			v2 = -3.0*s[0]*t - 2.0*s[1];
			v2 = v2*t - s[2];

			s = &sec.splineParamsR[0];
			number r = s[0]*t + s[1];
			r = r*t + s[2];
			r = r*t + s[3];

			sec_integral += gl.weight(i) * sqrt(VecNormSquared(vel)) / r;
		}

		integral += dt * sec_integral;

		// update lower bound and iterator
		t_start = sec_tend;
		if (t_start >= t_end) break;
		++sec_it;
	}

	return integral;
}

	/**
	 * @brief helper method to calculate the axial positions
	 */
	void calculate_segment_axial_positions
(
	std::vector<number>& segAxPosOut,
	number t_start,
	number t_end,
	const NeuriteProjector::Neurite& neurite,
	size_t startSec,
	number segLength
)
{
	const size_t nSeg = segAxPosOut.size();

	GaussLegendre gl(5);
	size_t nPts = gl.size();

	std::vector<NeuriteProjector::Section>::const_iterator sec_it = neurite.vSec.begin() + startSec;
	std::vector<NeuriteProjector::Section>::const_iterator sec_end = neurite.vSec.end();

	// check that startSec was correct
	number sec_tstart = startSec > 0 ? (sec_it - 1)->endParam : 0.0;
	number sec_tend = sec_it->endParam;
	UG_COND_THROW(sec_tend < t_start || sec_tstart > t_start,
		"Wrong section iterator given to calc_length_over_radius().");

	number integral = 0.0;
	size_t seg = 0;
	while (sec_it != sec_end)
	{
		// integrate from t_start to min{t_end, sec_tend}
		const NeuriteProjector::Section& sec = *sec_it;
		sec_tend = sec.endParam;
		number dt = std::min(t_end, sec_tend) - t_start;
		number sec_integral = 0.0;
		for (size_t i = 0; i < nPts; ++i)
		{
			number t = sec_tend - (t_start + dt*gl.point(i)[0]);

			vector3 vel;
			const number* s = &sec.splineParamsX[0];
			number& v0 = vel[0];
			v0 = -3.0*s[0]*t - 2.0*s[1];
			v0 = v0*t - s[2];

			s = &sec.splineParamsY[0];
			number& v1 = vel[1];
			v1 = -3.0*s[0]*t - 2.0*s[1];
			v1 = v1*t - s[2];

			s = &sec.splineParamsZ[0];
			number& v2 = vel[2];
			v2 = -3.0*s[0]*t - 2.0*s[1];
			v2 = v2*t - s[2];

			s = &sec.splineParamsR[0];
			number r = s[0]*t + s[1];
			r = r*t + s[2];
			r = r*t + s[3];

			sec_integral += gl.weight(i) * sqrt(VecNormSquared(vel)) / r;
		}
		integral += dt * sec_integral;

		// calculate exact position by linear interpolation, whenever integral has surpassed it
		while (integral >= (seg+1)*segLength)
		{
			number lastIntegral = integral - dt * sec_integral;
			segAxPosOut[seg] = t_start + ((seg+1)*segLength - lastIntegral) / sec_integral;
			++seg;
		}

		// update lower bound and iterator
		t_start = sec_tend;
		if (t_start >= t_end) break;
		++sec_it;
	}

	// maybe write last position (should not happen, but you never know)
	segAxPosOut[nSeg-1] = t_end;
}

	/**
	 * @brief creates soma
	 */
	static void create_soma
(
		const std::vector<SWCPoint>& somaPts,
		Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		SubsetHandler& sh,
		size_t si,
		size_t numRefs = 2
)
{
	UG_COND_THROW(somaPts.size() != 1, "Currently only one soma point is allowed by this implementation");
	Selector sel(g);
	GenerateIcosphere(g, somaPts.front().coords, somaPts.front().radius, numRefs, aPosition, &sel);
	AssignSelectionToSubset(sel, sh, si);
}

	static void connect_inner_soma_with_ER_by_projection
	(
			size_t somaIndex, /// inner soma index: beginning of quads is somaIndex+1
			size_t numQuads, /// number of quads
		    Grid& g,
		    Grid::VertexAttachmentAccessor<APosition>& aaPos,
		    SubsetHandler& sh
	) {

	}


	static float calculate_angle(const vector3& pos, const vector3& origin, const vector3& point, const vector3& n) {
		vector3 v1, v2;
		VecSubtract(v1, point, pos);
		VecSubtract(v2, origin, pos);
		number dot = VecDot(v1, v2);
		vector3 cross, normal;
		VecCross(cross, v1, v2);
		VecNormalize(normal, n);
		number det = VecDot(normal, cross);
		return rad_to_deg(std::atan2(det, dot));
	}



	/**
	 * @brief calculates angle offset
	 */
	static float calculate_angle(const vector3& pos, const vector3& origin, const vector3& point) {
		vector3 v1, v2;
		VecSubtract(v1, point, pos);
		VecSubtract(v2, origin, pos);
		number dot = VecDot(v1, v2);
		number v1_len = VecLengthSq(v1);
		number v2_len = VecLengthSq(v2);

		number x = (dot / (std::sqrt(v1_len*v2_len)));
		UG_LOGN("acos(" << x << ")");

		return rad_to_deg(acos(x));
	}

	/**
	 * @brief calculates angles from [-180, 180] interval to [0, 360]
	 */
	static float deg_to_full_range(float angle) {
		return fmod(angle+360, 360);
	}

	/**
	 * @brief calculates the angles according to a center/origin vertex
	 */
	static void calculate_angles(const vector3& originCenter, const std::vector<ug::vector3>& points, std::vector<number>& angles, size_t refIndex=-1) {
		/// reference point
		ug::vector3	ref = points[(refIndex > points.size() ? 0 : refIndex)];

		/// calculate each angle
		for (std::vector<ug::vector3>::const_iterator it = points.begin(); it != points.end(); ++it) {
			angles.push_back(calculate_angle(originCenter, ref, *it));
		}
	}

	/**
	 * @brief calculates the angles according to a center/origin vertex specified
	 */
	static void calculate_angles(const vector3& originCenter, const std::vector<ug::vector3>& points,  std::vector<number>& angles, std::vector<ug::vector3>& normals, const ug::vector3& ref) {
		size_t i = 0;
		/// calculate each angle
		for (std::vector<ug::vector3>::const_iterator it = points.begin(); it != points.end(); ++it) {
			angles.push_back(calculate_angle(originCenter, ref, *it, normals[i]));
			i++;
		}
	}

	/**
	 * @brief calculates the angles according to a center/origin vertex specified implicilty
	 */
	static void calculate_angles(const vector3& originCenter, const std::vector<ug::vector3>& points, std::vector<number>& angles, std::vector<ug::vector3>& normals, size_t refIndex=-1) {
		/// reference point
		ug::vector3	ref = points[(refIndex > points.size() ? 0 : refIndex)];

		size_t i = 0;
		/// calculate each angle
		for (std::vector<ug::vector3>::const_iterator it = points.begin(); it != points.end(); ++it) {
			angles.push_back(calculate_angle(originCenter, ref, *it, normals[i]));
			i++;
		}
	}

	/**
	 * sorts a vector values and stores the original indices
	 */
	static void sortIndices(const std::vector<number>& values, std::vector<size_t>& indices) {
		std::vector<std::pair<size_t, size_t> > a;
		for (size_t i = 0 ; i < values.size() ; i++) {
		    a.push_back(make_pair(values[i], i));
		}

		std::sort(a.begin(), a.end());
		for (size_t i = 0; i < a.size(); i++) {
			indices.push_back(a[i].second);
		}
	}

	/**
	 * @brief this connects the outer soma surface with the root neurites vertices ER (smaller quad) and PM (larger quad)
	 * First the root neurite vertices are projected to the outer sphere / soma surface. Then the smallest angle pairs
	 * are calculated and then the associated unprojected vertices are connected with the small and large quad soma sphere vertices
	 */
    static void connect_outer_and_inner_root_neurites_to_outer_soma
    (
		size_t somaIndex,
		size_t numQuads,
	    Grid& g,
	    Grid::VertexAttachmentAccessor<APosition>& aaPos,
	    SubsetHandler& sh,
	    std::vector<ug::Vertex*>& rootNeurites, /// Note: is just one single array => need to iterate by stride 4
	    std::vector<ug::Vertex*>& rootNeuritesInner /// TODO: Implement the same projection and connection method for innerRootNeurites vertices
    )
    {
		Selector::traits<Edge>::iterator eit;
		Selector::traits<Edge>::iterator eit_end;
		Selector sel(g);

    	/// TODO Project root neurite vertices to outer soma surface
    	/// Find minimal angle vertices... and merge them at outer soma surface
		/// (don't need to introduce factor 1.05 for outer soma, can let at 1.00)
		std::vector<std::vector<ug::Vertex*> > projectedVertices;
		std::vector<std::vector<ug::Vertex*> > projectedVertices2;
		std::vector<std::vector<ug::vector3> > projected;
		std::vector<std::vector<ug::vector3> > allNormals;
		projected.resize(numQuads);
		projectedVertices.resize(numQuads);
		projectedVertices2.resize(numQuads);
    	for (size_t i = 1; i < numQuads+1; i++) {
    		UG_LOGN("Selecting now subset: " << somaIndex+i);
    		sel.clear();
    		/// Select outer soma inner quad
    		SelectSubsetElements<Vertex>(sel, sh, somaIndex+i, true);
 			ug::Vertex* v0 = *(sel.vertices_begin());
 			UG_LOGN("First vertex of subset: " << aaPos[v0]);
    		SelectSubsetElements<Edge>(sel, sh, somaIndex+i, true);
    		eit = sel.edges_begin();
    		eit_end = sel.edges_end();
    		std::vector<std::pair<ug::Vertex*, ug::Vertex*> > es;
    		size_t count = 0;

			UG_LOGN("Trying to find edges...");
			for (; eit != eit_end; ++eit) {
				Edge* e = *eit;
				UG_LOGN("trying e->vertex(0)");
				std::pair<ug::Vertex*, ug::Vertex*> p;
				if (e->vertex(0) == v0) {
					p.first = v0;
					p.second = e->vertex(1);
					count++;
					es.push_back(p);
				}

				UG_LOGN("trying e->vertex(1)");
				if (e->vertex(1) == v0) {
					p.first = e->vertex(0);
					p.second = v0;
					es.push_back(p);
					count++;
				}
			}
			UG_LOGN("Found edges");

			UG_COND_THROW(count != 2, "Number of edges has to be two!");

			ug::vector3 v1, v2;
			VecSubtract(v1, aaPos[es[0].first], aaPos[es[0].second]);
			UG_LOGN("Subtracted v1");
			VecSubtract(v2, aaPos[es[1].first], aaPos[es[1].second]);
			UG_LOGN("Subtracted v2");
			ug::vector3 normal;
			VecCross(normal, v1, v2);
			VecNormalize(normal, normal);
			sel.clear();

			UG_LOGN("calculated cross product");

			size_t numVerts = 4;
			ug::vector3 vProjected;
			std::vector<ug::vector3> normals;
			for (size_t l = 0; l < numVerts; l++) {
				ug::Vertex* vit = rootNeurites[(i-1)*numVerts+l];
				ug::vector3 v;
				VecSubtract(v, aaPos[vit], v1);
				number dot = VecDot(v, normal);
				vProjected = aaPos[vit];
				VecScaleAdd(vProjected, 1.0, vProjected, dot, normal);
				ug::Vertex* projVert = *g.create<ug::RegularVertex>();
				const ug::vector3& n = -normal;
				ProjectPointToPlane(vProjected, aaPos[vit], aaPos[es[0].first], n);
				aaPos[projVert] = vProjected;
				projected[i-1].push_back(vProjected);
				projectedVertices[i-1].push_back(vit); /// save original vertex from which we projected (root neurite)
				projectedVertices2[i-1].push_back(projVert); /// the actual projected vertex on the soma surface
				normals.push_back(normal);
			}
			allNormals.push_back(normals);

			UG_LOGN("First projection!");

    	}

    	/// TODO: find the corresponding vertices (See comments below how to do it)
		std::vector<std::vector<number> > allAngles;
		std::vector<std::vector<number> > allAnglesInner;

		Selector::traits<Vertex>::iterator vit;
		Selector::traits<Vertex>::iterator vit_end;
		size_t j = 1;
		for (std::vector<std::vector<ug::vector3> >::const_iterator it = projected.begin(); it != projected.end(); ++it) {
			ug::vector3 centerOut;
			std::vector<number> angles;
			CalculateCenter(centerOut, &(*it)[0], it->size());
			///calculate_angles(centerOut, *it, angles);
			calculate_angles(centerOut, *it, angles, allNormals[j-1], (*it)[0]);
			allAngles.push_back(angles);

			sel.clear();
			std::vector<number> angles2;
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+j, true);
			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			std::vector<ug::vector3> verts;
			for (; vit != vit_end; ++vit) {
				verts.push_back(aaPos[*vit]);
			}
			calculate_angles(centerOut, verts, angles2, allNormals[j-1], (*it)[0]);
			allAnglesInner.push_back(angles2);
			j++;
		}

		UG_LOGN("Found angles");

		for (size_t i = 0; i < 1; i++) {
			for (size_t j = 0; j < allAngles[i].size(); j++) {
				UG_LOGN("(old angles): " << allAngles[i][j]); /// allAngles are projected inner vertices (from outer vertices)
			}
		}

		for (size_t i = 0; i < 1; i++) {
			for (size_t j = 0; j < allAnglesInner[i].size(); j++) {
				UG_LOGN("(old angles inner): " << allAnglesInner[i][j]); /// allAnglesInner are unprojected inner vertices
			}
		}



		/// TODO: Maybe have to center the projected vertices around the soma sphere's large respectively small surface quad center
		for (size_t i = 0; i < allAngles.size(); i++) {
			for (size_t j = 0; j < allAngles[i].size(); j++) {
				number a = allAngles[i][j];
				number b = allAnglesInner[i][j];
				/// Convert from -180,180 to 0,360 -> then sort array based on this
				if (a > 180) {
					a -= 360;
				}
				if (a < -180) {
					a += 360;
				}

				if (a < 0) {
					a = a + 360;
				}

				if (b > 180) {
					b -= 360;
				}
				if (b < -180) {
					b += 360;
				}

				if (b < 0) {
					b = b + 360;
				}

				allAngles[i][j] = a;
				allAnglesInner[i][j] = b;
			}
		}

		for (size_t i = 0; i < 1; i++) {
			for (size_t j = 0; j < allAngles[i].size(); j++) {
				UG_LOGN("(new angles): " << allAngles[i][j]); /// projected inner vertices (from outer vertices of root neurite)
			}
		}

		for (size_t i = 0; i < 1; i++) {
			for (size_t j = 0; j < allAnglesInner[i].size(); j++) {
				UG_LOGN("(new angles inner): " << allAnglesInner[i][j]); /// unprojected inner vertices
			}
		}

		/*
		std::vector<std::vector<size_t> >indices;
		std::vector<std::vector<size_t> >indices2;
		indices.resize(allAngles.size());
		indices2.resize(allAnglesInner.size());
		/// Sort angles to determine the pairs: Pair the smallest angle, then the second smallest angle, and so forth...
		/// This will be used to connect neurite root vertices with the vertices on the outer soma -> same strategy
		/// as before: We remember pairs of unprojected and projected vertices and connect unprojected vertices with
		/// the vertices which correspond to the minimum angle pairs we found below by using the projected vertices
		for (size_t i = 0; i < allAngles.size(); i++) {
			std::vector<number> angleCopy = allAngles[i];
			std::vector<number> angleCopy2 = allAnglesInner[i];
			sortIndices(angleCopy, indices[i]);
			sortIndices(angleCopy2, indices2[i]);
		}
		*/


		/// get mapping of outer vertices to inner (Unprojected) vertices
		std::map<Vertex*, Vertex*> innerToOuter;
		for (size_t i = 0 ; i < projected.size(); i++) {
			/// Store the unprojected (outer vertices - root neurite start vertices)
			sel.clear();
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+i+1, true);
			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			std::vector<ug::Vertex*> vertsOuterVtx; /// each inner quad vertices
			for (; vit != vit_end; ++vit) {
				vertsOuterVtx.push_back(*vit);
			}

			/// map the unprojected vertex to the corresponding outer vertex
			for (size_t j = 0; j < 4; j++) {
				//innerToOuter[vertsOuterVtx[indices2[i][j]]] = projectedVertices[i][indices[i][j]]; // (projectedVertices contains unprojected vertices -> but needs to be reshuffled because angle vertices have been sorted)
				///innerToOuter[vertsOuterVtx[j]] = projectedVertices[i][j]; // (projectedVertices contains unprojected vertices -> but needs to be reshuffled because angle vertices have been sorted)
			}

			/// this should connect the projected vertices with unprojected root neurites

			std::vector<std::pair<ug::Vertex*, number> > anglesOfProjectedInnerVertices;
			std::vector<std::pair<ug::Vertex*, number> > anglesOfOrginalSomaInnerVertices;

			/// calculate angles for inner original soma vertices and projected vertices (projectedVertices is ug::Vertex*)
			ug::vector3 centerOut;
			std::vector<number> angles;
			CalculateCenter(centerOut, &projected[i][0], 4);
			calculate_angles(centerOut, projected[i], angles, allNormals[i], projected[i][0]);
			for (size_t n = 0; n < 4; n++) {
				anglesOfProjectedInnerVertices.push_back(make_pair(projectedVertices[i][n], angles[n]));
			}

			sel.clear();
			std::vector<number> angles2;
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+i+1, true);
			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			std::vector<ug::vector3> verts;
			std::vector<ug::Vertex*> vertsVtx;
			for (; vit != vit_end; ++vit) {
				verts.push_back(aaPos[*vit]);
				vertsVtx.push_back(*vit);
			}
			calculate_angles(centerOut, verts, angles2, allNormals[i], projected[i][0]);
			for (size_t n = 0; n < 4; n++) {
				anglesOfOrginalSomaInnerVertices.push_back(make_pair(vertsVtx[n], angles2[n]));
			}

			/// Store: Projected vertex on soma surface to unprojected vertex (root neurite vertices)
			std::map<ug::Vertex*, ug::Vertex*> projectedToUnprojected;
			for (size_t l = 0; l < 4; l++) {
				projectedToUnprojected[projectedVertices[i][l]] = rootNeurites[(i)*4+l];
			}

			/// TODO: sort both vector lists: angleOfProjectedInnerVertices and anglesOfOriginalSomaInnerVertices by value

			for (size_t l = 0; l < anglesOfOrginalSomaInnerVertices.size(); l++) {
				Vertex* p1 = anglesOfOrginalSomaInnerVertices[l].first; // original soma vertex
				Vertex* temp = anglesOfProjectedInnerVertices[l].first; // projected vertex on somasurface
				Vertex* p2 = projectedToUnprojected[temp]; // look up the unprojected vertex (root neurite vertex)

				UG_COND_THROW(!p1, "P1 not found!");
				UG_COND_THROW(!p2, "P2 not found!");

				innerToOuter[p1] = p2;
			}

			/*
			/// NOTE: if we run through the original soma quad vertices -> we don't know which order -> thus the projectedVertices2 zuordnung is always not determined...
			/// Old strategy: Find closest vertex based on angle
			for (size_t l = 0; l < 4; l++) {
				ug::Vertex* vit = rootNeurites[(i)*4+l]; /// from root neurite outer (Unprojected vertex) i's l's vertex to ...
				innerToOuter[vit] = projectedVertices2[i][l]; /// the corresponding projected vertex on soma surface's quad i with vertex l
				/// TODO: for each vit do a loop with the projectedVertices corresponding to this -> find the smallest distance
			}
			*/
		}

		/*
		UG_LOGN("Sorted indices:")
		for (size_t j = 0; j < 1; j++) {
			for (size_t i = 0; i < indices[j].size(); i++) {
			UG_LOGN("i (indices):" << indices[j][i]);
			UG_LOGN("angles[i] (sorted): " << allAngles[j][indices[j][i]]);
			}
		}

		for (size_t j = 0; j < 1; j++) {
			for (size_t i = 0; i < indices2[j].size(); i++) {
			UG_LOGN("i (indices):" << indices2[j][i]);
			UG_LOGN("angles[i] (sorted): " << allAnglesInner[j][indices2[j][i]]);
			}
		}

		UG_LOGN("sorted angles");
		*/

		for (std::map<ug::Vertex*, ug::Vertex*>::iterator it = innerToOuter.begin(); it != innerToOuter.end(); ++it) {
			ug::Vertex* p1 = it->first;
			ug::Vertex* p2 = it->second;
			ug::Edge* e = *g.create<RegularEdge>(EdgeDescriptor(p1, p2));
		}

		UG_LOGN("created edges !!! for outer done !!! might be screwed because angles not sorted.")

		/*
		for (size_t i = 0; i < 1; i++) {
			sel.clear();
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+1+i, true);
			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			std::vector<ug::Vertex*> vertsInnerVtx;

			for (; vit != vit_end; ++vit) {
				vertsInnerVtx.push_back(*vit);
			}

			for (size_t j = 0; j < 2; j++) {
				ug::Vertex* p1 = rootNeurites[(i*4)+indices2[i][j]];
				ug::Vertex* p2 = vertsInnerVtx[indices[i][j]];
				ug::Edge* e = *g.create<RegularEdge>(EdgeDescriptor(p1, p2));
			}
		}
		*/

		/*
		std::vector<std::vector<std::pair<size_t, size_t > > > pairs;
				for (size_t k = 0; k < numQuads; k++) {
					std::vector<std::pair<size_t, size_t> > pair;
					for (size_t i = 0; i < allAngles[k].size(); i++) {
						number dist = std::numeric_limits<number>::infinity();
						size_t smallest = 0;
						for (size_t j = 0 ; j < allAnglesInner[k].size(); j++) {
							number altDist = allAnglesInner[k][j] - allAngles[k][i];
							altDist += (altDist>180) ? -360 : (altDist<-180) ? 360 : 0;
							altDist = std::abs(altDist);
							if (altDist < dist) {
								dist = altDist;
								smallest = j;
							}
						}
						pair.push_back(make_pair<size_t, size_t>(indices[k][i], indices2[k][smallest])); /// TODO: use indices from above to map
					}
					pairs.push_back(pair);
				}

			UG_LOGN("found pairs");
			for (size_t i = 0; i < pairs.size(); i++) {
				sel.clear();
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+1+i, true);
				vit = sel.vertices_begin();
				vit_end = sel.vertices_end();
				std::vector<ug::Vertex*> vertsInnerVtx;

				for (; vit != vit_end; ++vit) {
					vertsInnerVtx.push_back(*vit);
				}

				const std::vector<std::pair<size_t, size_t > >& temp = pairs[i]; // first quad
				for (size_t j = 0; j < temp.size(); j++) {
					const std::pair<size_t, size_t >& pair = temp[j]; // get a pair
					ug::Vertex* p1 = rootNeurites[(i*4)+pair.first];
					ug::Vertex* p2 = vertsInnerVtx[pair.second];
					ug::Edge* e = *g.create<RegularEdge>(EdgeDescriptor(p1, p2));
				}
			}
			*/

    }


    /**
     * @brief This connects the inner neurites (ER) to the inner sphere's (ER) surface quads
     * Projected vertices to the inner sphere's quad plane will correspond to an unprojected outer sphere's quad vertex.
     * The projected vertices are used to find the vertices with smallest angle difference between the inner sphere's quad vertices and the projected vertices onto that plane
     * Then the unprojected vertices and the corresponding vertex on the inner sphere's quad are connected by edges and faces.
     * For connecting inner neurites to inner soma / ER it is also possible to compare distances...
     */
	static void connect_inner_neurites_to_inner_soma
	(
			size_t somaIndex, /// inner soma index: beginning of inner sphere's quads (ER) is somaIndex+1, outer sphere's quads (ER) is somaIndex-numQuads-1
			size_t numQuads, /// number of total surface quads or neurites to connect to
		    Grid& g,
		    Grid::VertexAttachmentAccessor<APosition>& aaPos,
		    SubsetHandler& sh
	) {
		/// For the test geometry this has to be 10
		/// UG_COND_THROW(somaIndex != 10, "Soma index needs to be 10!");
		Selector::traits<Vertex>::iterator vit;
		Selector::traits<Vertex>::iterator vit_end;
		Selector::traits<Edge>::iterator eit;
		Selector::traits<Edge>::iterator eit_end;
		Selector sel(g);

		/// the projected vertices coordinates as a vector3 and the corresponding grid vertices
		std::vector<std::vector<ug::vector3> > projected;
		std::vector<std::vector<ug::Vertex*> > projectedVertices;
		std::vector<ug::vector3> normals;
		projected.resize(numQuads);
		projectedVertices.resize(numQuads);
		SaveGridToFile(g, "before_projections_inner.ugx");

		/// find all edges for each inner sphere's surface quad - take two starting at the same vertex to get two edges for normal calculation
		for (size_t i = 1; i < numQuads+1; i++) {
			UG_LOGN("Selecting now subset: " << somaIndex+i);
			sel.clear();
			/// Select inner soma quad
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+i, true);
			ug::Vertex* v0 = *(sel.vertices_begin());
			UG_LOGN("First vertex of subset: " << aaPos[v0]);
			SelectSubsetElements<Edge>(sel, sh, somaIndex+i, true);
			eit = sel.edges_begin();
			eit_end = sel.edges_end();
			std::vector<std::pair<ug::Vertex*, ug::Vertex*> > es;
			size_t count = 0;

			UG_LOGN("Trying to find edges...");
			for (; eit != eit_end; ++eit) {
				Edge* e = *eit;
				UG_LOGN("trying e->vertex(0)");
				std::pair<ug::Vertex*, ug::Vertex*> p;
				if (e->vertex(0) == v0) {
					p.first = v0;
					p.second = e->vertex(1);
					count++;
					es.push_back(p);
				}

				UG_LOGN("trying e->vertex(1)");
				if (e->vertex(1) == v0) {
					p.first = e->vertex(0);
					p.second = v0;
					es.push_back(p);
					count++;
				}
			}

			/// Two edges found starting in same vertex
			UG_LOGN("Found edges");
			UG_COND_THROW(count != 2, "Number of edges has to be two to calculate a normal");

			/// Now calculate the normal for this plane / inner sphere's quad
			ug::vector3 v1, v2;
			VecSubtract(v1, aaPos[es[0].first], aaPos[es[0].second]);
			UG_LOGN("Subtracted v1");
			VecSubtract(v2, aaPos[es[1].first], aaPos[es[1].second]);
			UG_LOGN("Subtracted v2");
			ug::vector3 normal;
			VecCross(normal, v1, v2);
			VecNormalize(normal, normal);
			sel.clear();
			UG_LOGN("calculated cross product");

			/// Project each of the outer sphere's quad (ER) to the plane described by the normal of the inner sphere's quad
			ug::vector3 vProjected;
			SelectSubsetElements<Vertex>(sel, sh, somaIndex-numQuads+i-1, true);
			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			for (; vit != vit_end; ++vit) {
				ug::vector3 v;
				VecSubtract(v, aaPos[*vit], v1);
				number dot = VecDot(v, normal);
				vProjected = aaPos[*vit];
				VecScaleAdd(vProjected, 1.0, vProjected, dot, normal);
				ug::Vertex* projVert = *g.create<ug::RegularVertex>();
				const ug::vector3& n = -normal;
				ProjectPointToPlane(vProjected, aaPos[*vit], aaPos[es[0].first], n);
				aaPos[projVert] = vProjected;
				projected[i-1].push_back(vProjected);
				projectedVertices[i-1].push_back(*vit); /// save original vertex from which we proojected
			}
			normals.push_back(normal);
			UG_LOGN("First projection!");
		}

		/// TODO: comment code
		/// TODO: cleanup Code
		/// TODO: split up large file

		/// TODO: Could also calculate an averaged plane, e.g. calculate two plane normals for each quad, average them
		/// TODO: Should get normal not from the two points of each inner quad but define the normal to be the edge through the center of the inner sphere's quad (ER) and outer sphere's quad (ER) part

		/// Note: Strategie für innere Verbindungen
		/// 1. Berechne normale (Axis) durch den Mittelpunkt der beiden Oberflächenquads (inner soma und äußeres soma) für den ER Teil
		/// 2. Projiziere äußere soma quad vertices des ER auf inneres Soma quad ebene definiert durch normal und mittelpunkt
		/// 3. Projiziere innere soma quad vertices des ER auf inneres soma quad ebene definiert durch normal und mittelpunkt
		/// 4. Zentriere projizierte äußere soma quad vertices des er um den mittelpunkt des quads auf innerem soma
		/// 5. Berechne winkel und finde paare mit minimalem winkel offset (oder closest distance geht evt. auch)
		/// 6. Move die alten inneren soma quad vertices auf die position der projizierten vertices von äußerem soma quad
		/// 7. merke die paare von projizierten vertices, starte mit jewiels 2 paaren (2 projizierte vertices auf innerem soma , 2 bewegte vertices auf dem inneren soma)
		/// 8. builde faces für 4 vertices immer -> hexaeder. ==> Verbindungshexaeder muss dann gespeichert werden für refinement! ebenso soma radien für beide branching point anschlüsse!!!

		/// Note: Strategie für äußere Verbindungen
		/// Projiziere ebenso auf ebene, jetzt aber auf den äußeren quad vertices des ERs (äußeres Soma)
		/// Projiziere Neuritenstartknoten ebenso darauf. Finde kleinsten Winkel, dann merge diese Vertices!

		/// Calculate all angles with respect to a reference point for the projected outer sphere's quad (ER) vertices to the inner sphere's quad and the inner sphere's vertices
		std::vector<std::vector<number> > allAngles;
		std::vector<std::vector<number> > allAnglesInner;
		size_t j = 1;
		for (std::vector<std::vector<ug::vector3> >::const_iterator it = projected.begin(); it != projected.end(); ++it) {
			ug::vector3 centerOut;
			std::vector<number> angles;
			CalculateCenter(centerOut, &(*it)[0], it->size());
			///calculate_angles(centerOut, *it, angles);
			calculate_angles(centerOut, *it, angles, normals);
			allAngles.push_back(angles);

			sel.clear();
			std::vector<number> angles2;
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+j, true);
			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			std::vector<ug::vector3> verts;
			for (; vit != vit_end; ++vit) {
				verts.push_back(aaPos[*vit]);
			}
			calculate_angles(centerOut, verts, angles2, normals);
			allAnglesInner.push_back(angles2);
			j++;
		}

		for (std::vector<std::vector<number> >::const_iterator it = allAngles.begin(); it != allAngles.end(); ++it) {
			UG_LOGN("Quad angles projected from outer to inner...");
			for (std::vector<number>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				UG_LOGN(*it2);
			}
			UG_LOGN("---");
		}

		for (std::vector<std::vector<number> >::const_iterator it = allAnglesInner.begin(); it != allAnglesInner.end(); ++it) {
			UG_LOGN("Quad angles inner...");
			for (std::vector<number>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				UG_LOGN(*it2);
			}
			UG_LOGN("---");
		}

		/// Remember pairs: For each projected vertices to the inner sphere's quad plane a corresponding vertices we projected from the outer sphere's quad exist these have to be connected by edges / faces to create a hexaeder
		std::vector<std::vector<std::pair<size_t, size_t > > > pairs;
		for (size_t k = 0; k < numQuads; k++) {
			std::vector<std::pair<size_t, size_t> > pair;
			for (size_t i = 0; i < allAngles[k].size(); i++) {
				number dist = std::numeric_limits<number>::infinity();
				size_t smallest = 0;
				for (size_t j = 0 ; j < allAnglesInner[k].size(); j++) {
					number altDist = allAnglesInner[k][j] - allAngles[k][i];
					altDist += (altDist>180) ? -360 : (altDist<-180) ? 360 : 0;
					altDist = std::abs(altDist);
					if (altDist < dist) {
						dist = altDist;
						smallest = j;
					}
				}
				pair.push_back(make_pair<size_t, size_t>(i, smallest));
			}
			pairs.push_back(pair);
		}

		for (std::vector<std::vector<std::pair<size_t, size_t > > >::const_iterator it = pairs.begin(); it != pairs.end(); ++it) {
			UG_LOGN("***");
			for (std::vector<std::pair<size_t, size_t> >::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				UG_LOGN("Pair " << it2->first << " -> " << it2->second);
			}
			UG_LOGN("***");
		}
		/// TODO: The angle is not calculate correctly: why? Because two calls to calculate_angles, but the REFERENCE point differs. Has to be precisely the same to make sense.

		/// find closest vertex instead of minimum angle difference: this should be save for the inner sphere and outer sphere ER part connection
		j = 1;
		std::vector<std::vector<std::pair<ug::vector3, ug::vector3> > > myPairs;
		std::map<Vertex*, Vertex*> myPairs2;
		for (size_t i = 0; i < projected.size(); i++) { /// each outer quad projected vertices (now on inner soma)
			sel.clear();
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+j, true);
			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			std::vector<ug::vector3> vertsInner; /// each inner quad vertices
			std::vector<ug::Vertex*> vertsInnerVtx; /// each inner quad vertices

			for (; vit != vit_end; ++vit) {
				vertsInner.push_back(aaPos[*vit]);
				vertsInnerVtx.push_back(*vit);
			}

			std::vector<std::pair<ug::vector3, ug::vector3> > pair;
			std::vector<std::pair<ug::Vertex*, ug::Vertex*> > pair2;
			for (size_t k = 0; k < projected[i].size(); k++) {
				number dist = std::numeric_limits<number>::infinity();
				ug::vector3 p1;
				ug::vector3 p2;
				size_t p1i;
				size_t p2i;
				for (size_t l = 0; l < vertsInner.size(); l++) {
					number altDist = VecDistance(projected[i][k], vertsInner[l]);
					if (altDist < dist) {
						dist = altDist;
						p1 = projected[i][k];
						p2 = vertsInner[l];
						p1i = k;
						p2i = l;
					}
				}
				pair.push_back(make_pair(p1, p2));
				myPairs2[vertsInnerVtx[p2i]] = projectedVertices[i][p1i];
			}
			myPairs.push_back(pair);
			j++;
		}

		UG_LOGN("vector pairs...");
		for (std::vector<std::vector<std::pair<ug::vector3, ug::vector3 > > >::const_iterator it = myPairs.begin(); it != myPairs.end(); ++it) {
			UG_LOGN("***");
			for (std::vector<std::pair<ug::vector3, ug::vector3> >::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				UG_LOGN("Pair " << it2->first << " -> " << it2->second);
			}
			UG_LOGN("***");
		}
		/// Nachdem vertices paare gefunden sind, müssen diese gemerkt werden welcher original vertex projiziert wurde.
		/// TODO: dann findet man die edges des inneren somas jeweils und für die edges die korrepsoniderten projizierten vertices => create face!
		/// Note: Winkelstrategie müsste anders behandelt werden wie oben beschrieben

		/// TODO Man kann zusätzlich auch die inneren Vertices auf die Ebene projizieren die zuvor definiert wurde
		/// (denn nicht alle liegen in der Ebene nur die zwei punkte die genommen wurden um die Normale zu berechnen)...
		/// und die projizierten äußeren Soma ER Quad verts auf das Zentrum verschieben siehe unten

		/// Iterate over all neurite connections (numQuads) and get the vertices of the inner sphere's quad edge each and find the corresponding unprojected (outer sphere's quad vertices) and form a face
		/// It is also possible to do the same procedure with the sorted angle differences above to create these faces
		for (size_t i = 1; i < numQuads+1; i++) {
			sel.clear();
			UG_LOGN("Selecting now subset: " << somaIndex+i);
			/// Select inner soma quad
			SelectSubsetElements<Edge>(sel, sh, somaIndex+i, true);
			eit = sel.edges_begin();
			eit_end = sel.edges_end();
			for (; eit != eit_end; ++eit) {
				UG_LOGN("Edge!");
				Edge* e = *eit;
				ug::Vertex* p1 = e->vertex(0);
				ug::Vertex* p2 = e->vertex(1);
				ug::Vertex* p3 = myPairs2[e->vertex(0)];
				ug::Vertex* p4 = myPairs2[e->vertex(1)];
	    		ug::Face* f = *g.create<Quadrilateral>(QuadrilateralDescriptor(p1, p3, p4, p2));
			}
		}

		/// TODO: Vor dem verbinden, verschiebe inner quad vertices auf die
		/// projizierten positionen auf innerem Soma welche das Resultat sind
		/// von der Projektion der äußeren Quad Soma vertices: Sollte nicht nötig
		/// für gutgeartete innere Quads - benötigt falls projiziertes Quad weit
		/// entfernt von dem inneren Quad ist - insebsondere um minimales Angle
		/// Differenz oder minimale Entfernung zu berechnen um Vertices zu verbinden

		SaveGridToFile(g, "after_projections_inner.ugx");


		/// Old strategy with quickhull
		/*
			/// Get vertices attached to inner soma
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+i, true);
			Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
			Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
			/// select inner vertices connected to outer soma (start vertices of ER)

			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			UG_LOGN("selected inner vertices attached to outer soma: " << sel.num());
			for (; vit != vit_end; ++vit) {
				Selector sel2(g);
				sel2.select(*vit);
				ExtendSelection(sel2, 1, true);
				if (sel2.num() == numNeighborHood) {
					temp.push_back(aaPos[*vit]);
					foo2.push_back(aaPos[*vit]);
				}
			}
			sel.clear();

			UG_COND_THROW(temp.size() != 8, "Need 8 vertices for calculating all faces.");
			#ifdef NC_WITH_QHULL
				using ug::neuro_collection::convexhull::gen_face;
				using ug::neuro_collection::convexhull::erase_face;
				gen_face(temp, g, sh, si+i, aaPos);
				erase_face(g, sh, si+i, aaPos, foo);
				erase_face(g, sh, si+i, aaPos, foo2);
			#else
				using ug::neuro_collection::quickhull::gen_face;
				gen_face(temp, temp2, g, sh, si+i, aaPos);
			#endif
		}
		*/
	}

	static std::vector<ug::vector3> find_quad_verts_on_soma
(
	   Grid& g,
	   Grid::VertexAttachmentAccessor<APosition>& aaPos,
	   std::vector<ug::Vertex*> verticesOld,
	   std::vector<std::vector<number> > outRads,
	   size_t si,
	   SubsetHandler& sh,
	   number rimSnapThresholdFactor,
	   size_t numQuads
) {
	size_t numVerts = 4;
	std::vector<ug::vector3> centers;
	UG_LOGN("3. AdaptSurfaceGridToCylinder")
	Selector sel(g);
	for (size_t i = 0; i < verticesOld.size(); i++) {
		sel.clear();
		ug::vector3 normal;
		CalculateVertexNormal(normal, g, verticesOld[i], aaPos);
		number radius = outRads[i][0];
		AdaptSurfaceGridToCylinder(sel, g, verticesOld[i], normal, radius, 1.0*rimSnapThresholdFactor, aPosition);
		UG_LOGN("Adaption done");

		sel.clear();
		sel.select(verticesOld[i]);
		ExtendSelection(sel, 1, true);
		CloseSelection(sel);
		sel.deselect(verticesOld[i]);
		g.erase(verticesOld[i]);

		UG_LOGN("num edges: " << sel.num<Edge>());
		size_t numEdges = sel.num<Edge>();
		size_t j = 0;
		while (numEdges > numVerts) {
			SubsetHandler::traits<Edge>::iterator eit = sel.begin<Edge>();
			SubsetHandler::traits<Edge>::iterator end = sel.end<Edge>();
			number bestLength = -1;
			Edge* eBest = NULL;
			for (; eit != end; ++eit) {
				const Edge* ee = *eit;
				Vertex* const* verts = ee->vertices();
				if (bestLength == -1) {
					bestLength = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
					eBest = *eit;
				} else {
					number length = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
					if (length < bestLength) {
						eBest = *eit;
						bestLength = length;
					}
				}
			}
			CollapseEdge(g, eBest, eBest->vertex(0));
			numEdges--;
			j++;
		}

		UG_LOGN("Collapsing done");


		std::vector<ug::vector3> vertices;
		Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
		Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
		UG_LOGN("Pushing vertices");
		for (; vit != vit_end; ++vit) {
			vertices.push_back(aaPos[*vit]);
		}

		UG_LOGN("Number of vertices: " << vertices.size());
		ug::vector3 centerOut;
		CalculateCenter(centerOut, &vertices[0], sel.num<ug::Vertex>());
		UG_LOGN("centerOut: " << centerOut);
		centers.push_back(centerOut);
	}
	return centers;
	}


	/**
	 * @brief connects neurites to soma
	 * TODO: Find suitable parameters for tangential smooth and resolve intersection
	 */
	static void connect_neurites_with_soma
(
	   Grid& g,
	   Grid::VertexAttachmentAccessor<APosition>& aaPos,
	   Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
	   std::vector<Vertex*> outVerts,
	   std::vector<Vertex*> outVertsInner,
	   std::vector<number> outRads,
	   std::vector<Vertex*>& smallerQuadVerts,
	   size_t si,
	   SubsetHandler& sh,
	   const std::string& fileName,
	   number rimSnapThresholdFactor,
	   std::vector<std::pair<size_t, std::pair<ug::vector3, ug::vector3> > >& axisVectors,
	   std::vector<NeuriteProjector::Neurite>& vNeurites,
	   std::vector<std::vector<ug::Vertex*> >& connectingVertices,
	   std::vector<std::vector<ug::Vertex*> >& connectingVerticesInner,
	   std::vector<std::vector<ug::Edge*> >& connectingEdges,
       std::vector<std::vector<ug::Edge*> >& connectingEdgesInner,
	   bool createInner=true,
	   number alpha=0.01,
	   int numIterations=10,
	   number resolveThreshold=0.00001,
	   number scale=0.5
) {
	UG_LOGN("1. Find the vertices representing dendrite connection to soma.");
	/// 1. Finde die 4 Vertices die den Dendritenanschluss darstellen zum Soma
	std::vector<std::vector<ug::vector3> > quads;
	std::vector<number> quadsRadii;
	size_t numVerts = 4;
	size_t numQuads = outVerts.size()/numVerts;

	for (size_t i = 0; i < numQuads; i++) {
		std::vector<ug::vector3> temp;
		for (size_t j = 0; j < numVerts; j++) {
			temp.push_back(aaPos[outVerts[(i*4)+j]]);
		}
		UG_LOGN("push a quad!");
		quads.push_back(temp);
	}

	UG_LOGN("2. Calculate center of each quad, find next surface vertex on soma.")
	/// 2. Berechne den Schwerpunkt jedes Quads und finde den nächstgelegenen
	///    Vertex des Oberflächengitters vom Soma
	std::vector<ug::vector3> centerOuts;
	std::vector<ug::vector3> centerOuts2;
	std::vector<Vertex*> bestVertices;
	for (size_t i = 0; i < numQuads; i++) {
		const ug::vector3* pointSet = &(quads[i][0]);
		ug::vector3 centerOut;
		CalculateCenter(centerOut, pointSet, numVerts);
		centerOuts.push_back(centerOut);
		Selector sel(g);
		SelectSubsetElements<Vertex>(sel, sh, si, true);
		Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
		Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
		number best = -1;
		ug::Vertex* best_vertex = NULL;
		for (; vit != vit_end; ++vit) {
			number dist = VecDistance(aaPos[*vit], centerOut);
			if (best == -1) {
				best = dist;
				best_vertex = *vit;
			} else if (dist < best) {
				best = dist;
				best_vertex = *vit;
			}
		}
		UG_COND_THROW(!best_vertex, "No best vertex found for quad >>" << i << "<<.");
		bestVertices.push_back(best_vertex);
	}

	UG_LOGN("3. AdaptSurfaceGridToCylinder")
	/// 3. Für jeden Vertex v führe AdaptSurfaceGridToCylinder mit Radius entsprechend
	///    dem anzuschließenden Dendritenende aus. Dadurch entsteht auf der Icosphere
	///    um jedes v ein trianguliertes 6- bzw. 5-Eck.
	Selector sel(g);
	for (size_t i = 0; i < bestVertices.size(); i++) {
		sel.clear();
		ug::vector3 normal;
		CalculateVertexNormal(normal, g, bestVertices[i], aaPos);
		number radius = outRads[i];
		AdaptSurfaceGridToCylinder(sel, g, bestVertices[i], normal, radius, 1.0*rimSnapThresholdFactor, aPosition);
	}

	AssignSubsetColors(sh);
	std::stringstream ss;
	ss << fileName << "_before_deleting_center_vertices.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	UG_LOGN("5. MergeVertices")
	/// 5. Wandle die stückweise linearen Ringe um die Anschlusslöcher per
	///    MergeVertices zu Vierecken um.
	sel.clear();
	for (std::vector<Vertex*>::iterator it = bestVertices.begin(); it != bestVertices.end(); ++it) {
		sel.select(*it);
		ExtendSelection(sel, 1, true);
		CloseSelection(sel);
		AssignSelectionToSubset(sel, sh, sh.num_subsets()+1);
		sel.clear();
	}

	AssignSubsetColors(sh);
	ss << fileName << "_before_getting_neighborhoods.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	UG_LOGN("4. Remove each vertex. Creates holes in soma")
	/// 4. Lösche jedes v, sodass im Soma Anschlusslöcher für die Dendriten entstehen.
	sel.clear();
	for (std::vector<Vertex*>::iterator it = bestVertices.begin(); it != bestVertices.end(); ++it) {
		sel.select(*it);
	}

	size_t numSubsets = sh.num_subsets();
	AssignSelectionToSubset(sel, sh, numSubsets);
	EraseElements<Vertex>(g, sh.begin<Vertex>(numSubsets), sh.end<Vertex>(numSubsets));
	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);
	ss << fileName << "_after_deleting_center_vertices.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	/// Collapse now edges and take smallest edges first
	size_t beginningOfQuads = si+1; // subset index where quads are stored in
	for (size_t i = 0; i < numQuads; i++) {
		size_t si = beginningOfQuads+i;
		size_t numEdges = sh.num<Edge>(si);
		size_t j = 0;
		while (numEdges > numVerts) {
			SubsetHandler::traits<Edge>::iterator eit = sh.begin<Edge>(si);
			SubsetHandler::traits<Edge>::iterator end = sh.end<Edge>(si);
			number bestLength = -1;
			Edge* eBest = NULL;
			for (; eit != end; ++eit) {
				const Edge* ee = *eit;
				Vertex* const* verts = ee->vertices();
				if (bestLength == -1) {
					bestLength = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
					eBest = *eit;
				} else {
					number length = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
					if (length < bestLength) {
						eBest = *eit;
						bestLength = length;
					}
				}
			}
			CollapseEdge(g, eBest, eBest->vertex(0));
			numEdges--;
			j++;
			std::stringstream ss;
			ss << fileName << "_after_collapse_number_" << j << "_for_quad_" << i << ".ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
		}

		SubsetHandler::traits<Vertex>::iterator vit = sh.begin<Vertex>(si);
		SubsetHandler::traits<Vertex>::iterator vend = sh.end<Vertex>(si);

		/// Outer Soma is assumed to start at -5 * outRads[i] of corresponding neurite and inner soma is assumed to start at -1
		for (; vit != vend; ++vit) {
			aaSurfParams[*vit].soma = true;
			if (createInner)
				aaSurfParams[*vit].axial = -5 * outRads[i];
			else
				aaSurfParams[*vit].axial = -1;
		}
	}

	if (createInner) {
		/// Shrink each quad on the outer soma surface
		for (size_t i = 0; i < numQuads; i++) {
			sel.clear();
			size_t si = beginningOfQuads+i;

			SelectSubsetElements<Vertex>(sel, sh, si, true);
			std::vector<Vertex*> vrts;
			sel.clear();
			UG_LOGN("verts size: " << vrts.size());

			std::vector<Edge*> edges;
			SelectSubsetElements<Edge>(sel, sh, si, true);
			edges.assign(sel.edges_begin(), sel.edges_end());
			sel.clear();

			UG_LOGN("edges size: " << edges.size());

			vrts.push_back(edges[0]->vertex(0));
			vrts.push_back(edges[0]->vertex(1));

			Vertex* prevVertex = edges[0]->vertex(1);
			size_t numIterations = edges.size()-1;
			std::vector<size_t> indices;
			edges.erase(edges.begin());
			UG_LOGN("number of edges: " << edges.size());
			while (!edges.empty()) {
				UG_LOGN("Still running: edges.size(): " << edges.size());
				for (size_t i = 0; i < edges.size(); i++) {
					Edge* nextEdge = edges[i];
					if (nextEdge->vertex(0) == prevVertex) {
						UG_LOGN("push first if");
						vrts.push_back(nextEdge->vertex(1));
						prevVertex = nextEdge->vertex(1);
						edges.erase(edges.begin()+i);
						break;
					}
					UG_LOGN("in between");
					if (nextEdge->vertex(1) == prevVertex) {
						UG_LOGN("push second if")
		            		vrts.push_back(nextEdge->vertex(0));
						prevVertex = nextEdge->vertex(0);
						edges.erase(edges.begin()+i);
						break;
					}
				}
			}

			vrts.erase(vrts.end()-1);

			std::vector<ug::Edge*> vEdgeOut;
			std::vector<ug::Vertex*> vVrtOut;

			Selector selToAssign(g);
			UG_COND_THROW(vrts.size() != 4, "Non-quadrilateral encountered. Cannot shrink a non-quadrilateral!");
			shrink_quadrilateral_copy(vrts, vVrtOut, vVrtOut, vEdgeOut, g, aaPos, -scale, false, &selToAssign, NULL);
			for (std::vector<Vertex*>::const_iterator it = vVrtOut.begin(); it != vVrtOut.end(); ++it) {
				smallerQuadVerts.push_back(*it);
			}
			AssignSelectionToSubset(selToAssign, sh, si+numQuads);
			sel.clear();
			selToAssign.clear();

			SubsetHandler::traits<Vertex>::iterator vit = sh.begin<Vertex>(si);
			SubsetHandler::traits<Vertex>::iterator vend = sh.end<Vertex>(si);

			/// Inner soma is assumed to start at -5 * outRads[i] of corresponding neurite too
			for (; vit != vend; ++vit) {
				aaSurfParams[*vit].soma = true;
				aaSurfParams[*vit].axial =  -5 * outRads[i];
				vNeurites[i].somaStart = 5 * outRads[i];
			}
		}
	}

	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);

	ss << fileName << "_after_merging_cylinder_vertices.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	UG_LOGN("8. TangentialSmooth");
	/// Note: TangentialSmooth -> alpha has to be corrected for different geometries.
	/// TangentialSmooth(g, g.vertices_begin(), g.vertices_end(), aaPos, alpha, numIterations);

	ss << fileName << "_after_merging_cylinder_vertices_and_tangential_smooth.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	UG_LOGN("6. Extrude rings along normal")
	/// 6. Extrudiere die Ringe entlang ihrer Normalen mit Höhe 0 (Extrude mit
	///    aktivierter create faces Option).
	sel.clear();
	std::vector<std::vector<Vertex*> > somaVerts;
	std::vector<std::vector<Vertex*> > allVerts;

	std::vector<std::vector<Vertex*> > somaVertsInner;
	std::vector<std::vector<Vertex*> > allVertsInner;


	for (size_t i = 0; i < numQuads; i++) {
		size_t si = beginningOfQuads+i;
		ug::vector3 normal;
		CalculateVertexNormal(normal, g, *sh.begin<Vertex>(si), aaPos);
		UG_LOGN("normal (outer): " << normal);
		ug::vector3 axisVector;
		CalculateCenter(sh.begin<Vertex>(si), sh.end<Vertex>(si), aaPos);
		/// indicate soma posiiton
		/*for (SubsetHandler::traits<Vertex>::iterator it = sh.begin<Vertex>(si); it != sh.end<Vertex>(si); ++it) {
			UG_LOGN("setting axial to -1!");
			aaSurfParams[*it].axial = -1;
		}*/

		VecAdd(axisVector, axisVector, normal);
		std::vector<Edge*> edges;
		std::vector<Vertex*> vertices;

		for (SubsetHandler::traits<Edge>::const_iterator it = sh.begin<Edge>(si); it != sh.end<Edge>(si); ++it) {
			edges.push_back(*it);
		}

		for (SubsetHandler::traits<Vertex>::const_iterator it = sh.begin<Vertex>(si); it != sh.end<Vertex>(si); ++it) {
			vertices.push_back(*it);
		}

		if (createInner) {
			for (SubsetHandler::traits<Edge>::const_iterator it = sh.begin<Edge>(si+numQuads); it != sh.end<Edge>(si+numQuads); ++it) {
				edges.push_back(*it);
			}
			for (SubsetHandler::traits<Vertex>::const_iterator it = sh.begin<Vertex>(si+numQuads); it != sh.end<Vertex>(si+numQuads); ++it) {
				vertices.push_back(*it);
			}
		}

		SelectSubsetElements<Vertex>(sel, sh, si, true);
		std::vector<Vertex*> temp;
		temp.assign(sel.vertices_begin(), sel.vertices_end());
		somaVerts.push_back(temp);
		sel.clear();
		temp.clear();

		if (createInner) {
			SelectSubsetElements<Vertex>(sel, sh, si+numQuads, true);
			temp.assign(sel.vertices_begin(), sel.vertices_end());
			somaVertsInner.push_back(temp);
			sel.clear();
			temp.clear();
		}

		SelectSubsetElements<Vertex>(sel, sh, si, true);

		if (createInner) {
			SelectSubsetElements<Vertex>(sel, sh, si+numQuads, true);
		}

		///Extrude(g, &vertices, &edges, NULL, normal, aaPos, EO_CREATE_FACES, NULL);

		/// store found soma vertices as connectingvertices for initial neurite vertices (used in create_neurite_general)
		for (size_t j= 0; j < vertices.size(); j++) {
			if (sh.get_subset_index(vertices[j]) == si) {
				connectingVertices[i].push_back(vertices[j]);
			}
			if (sh.get_subset_index(vertices[j]) == si+numQuads) {
				connectingVerticesInner[i].push_back(vertices[j]);
			}
		}

		/// store found soma edges as connectingedges for initial neurite vertices (used in create_neurite_general)
		for (size_t j= 0; j < edges.size(); j++) {
			if (sh.get_subset_index(edges[j]) == si) {
				connectingEdges[i].push_back(edges[j]);
			}
			if (sh.get_subset_index(edges[j]) == si+numQuads) {
				connectingEdgesInner[i].push_back(edges[j]);
			}
		}

		ug::vector3 centerOut2 = CalculateCenter(vertices.begin(), vertices.end(), aaPos);
		centerOuts2.push_back(centerOut2);
		sel.clear();

		/// indicate start of neurite with axial 0 and soma false explicitly: neurite start is 0 for outer soma and for inner soma it is -1 + radOut[i] * 5;
		for (std::vector<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
			aaSurfParams[*it].axial = -1 + 5 * outRads[i];
			vNeurites[i].somaStart = 5 * outRads[i];
			aaSurfParams[*it].soma = false;
		}

		SelectSubsetElements<Vertex>(sel, sh, si, true);
		std::vector<Vertex*> temp2;
		temp2.assign(sel.vertices_begin(), sel.vertices_end());
		allVerts.push_back(temp2);
		sel.clear();
		temp2.clear();

		if (createInner) {
			SelectSubsetElements<Vertex>(sel, sh, si+numQuads, true);
			temp2.assign(sel.vertices_begin(), sel.vertices_end());
			allVertsInner.push_back(temp2);
			sel.clear();
			temp2.clear();
		}

		ug::vector3 cylinderCenter = CalculateCenter(sh.begin<Vertex>(si), sh.end<Vertex>(si), aaPos);
		axisVectors.push_back(make_pair(si, make_pair(axisVector, cylinderCenter)));
	}

	ss << fileName << "_after_extruding_cylinders_before_removing_common_vertices.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	/// delete common vertices, thus keep only newly extruded vertices (used in next step 7.)
	for (size_t i = 0; i < numQuads; i++) {
		size_t numSomaVerts = somaVerts[i].size();
		for (size_t j = 0; j < numSomaVerts; j++) {
			allVerts[i].erase(std::remove(allVerts[i].begin(), allVerts[i].end(), somaVerts[i][j]), allVerts[i].end());
		}
	}

	if (createInner) {
		/// delete common vertices, thus keep only newly extruded vertices (used in next step 7.)
		for (size_t i = 0; i < numQuads; i++) {
			size_t numSomaVerts = somaVertsInner[i].size();
			for (size_t j = 0; j < numSomaVerts; j++) {
				allVertsInner[i].erase(std::remove(allVertsInner[i].begin(), allVertsInner[i].end(), somaVertsInner[i][j]), allVertsInner[i].end());
			}
		}
	}

	ss << fileName << "_after_extruding_cylinders.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();


	UG_LOGN("7. Calculate convex hull and connect")
	/// 7. Vereine per MergeVertices die Vertices der in 6. extrudierten Ringe jeweils
	///    mit den zu ihnen nächstgelegenen Vertices des entsprechenden Dendritenendes.
	si = beginningOfQuads;
	sel.clear();
	/*
	for (size_t i = 0; i < numQuads; i++) {
		std::vector<ug::vector3> temp;
		std::vector<ug::vector3> foo;
		std::vector<ug::vector3> foo2;
		std::vector<Vertex*> temp2;
		for (size_t j = 0; j < numVerts; j++) {
			temp.push_back(aaPos[outVerts[i*4+j]]);
			foo.push_back(aaPos[outVerts[i*4+j]]);
		}
		for (size_t j = 0; j < numVerts; j++) {
			temp.push_back(aaPos[allVerts[i][j]]);
			foo2.push_back(aaPos[allVerts[i][j]]);
		}

		UG_COND_THROW(temp.size() != 8, "Need 8 vertices for calculating all faces.");
		#ifdef NC_WITH_QHULL
			using ug::neuro_collection::convexhull::gen_face;
			using ug::neuro_collection::convexhull::erase_face;
			gen_face(temp, g, sh, si+i, aaPos);
			erase_face(g, sh, si+i, aaPos, foo);
			erase_face(g, sh, si+i, aaPos, foo2);
		#else
			using ug::neuro_collection::quickhull::gen_face;
			gen_face(temp, temp2, g, sh, si+i, aaPos);
		#endif
		ug::vector3 center;
		center = CalculateCenter(sh.begin<Vertex>(si+i+1000), sh.end<Vertex>(si+i+1000), aaPos);
		ug::vector3 axis;
		VecSubtract(axis, centerOuts2[i], centerOuts[i]);
		axisVectors.push_back(make_pair(si+i+numQuads*2, make_pair(axis, center)));
		/// numQuads*2 is required: n-inner quads and n-outer quads -> these quads here are the connecting quads
		/// i.e. first come all outer quads, then all inner quads, then the connecting outer quads, then the connecting inner quads
	}
	*/

	/*
	if (!createInner) {
		si = beginningOfQuads+numQuads;
		sel.clear();
		for (size_t i = 0; i < numQuads; i++) {
			UG_LOGN("First quad to connect...: " << i);
			std::vector<ug::vector3> temp;
			std::vector<Vertex*> temp2;
			std::vector<ug::vector3> foo;
			std::vector<ug::vector3> foo2;
			UG_LOGN("Accessing outVertsInner...; " << i);
			for (size_t j = 0; j < numVerts; j++) {
				temp.push_back(aaPos[outVertsInner[i*4+j]]);
				foo.push_back(aaPos[outVertsInner[i*4+j]]);
			}
			UG_LOGN("Accessing allVertsInner...; " << i);
			for (size_t j = 0; j < numVerts; j++) {
				temp.push_back(aaPos[allVertsInner[i][j]]);
				foo2.push_back(aaPos[allVertsInner[i][j]]);
			}

			UG_LOGN("Checking consistency of temp...");
			UG_COND_THROW(temp.size() != 8, "Need 8 vertices for calculating all faces.");
			#ifdef NC_WITH_QHULL
				UG_LOGN("Trying to use convexhull...");
				using ug::neuro_collection::convexhull::gen_face;
				using ug::neuro_collection::convexhull::erase_face;
				gen_face(temp, g, sh, si+i, aaPos);
				erase_face(g, sh, si+i, aaPos, foo);
				erase_face(g, sh, si+i, aaPos, foo2);
			#else
				UG_LOGN("Trying to use quickhull... ")
				using ug::neuro_collection::quickhull::gen_face;
				gen_face(temp, temp2, g, sh, si+i, aaPos);
			#endif
				UG_LOGN("Done with quickhull/convexhull...");
		}
	}
	*/

	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);
	ss << fileName << "_after_extruding_cylinders_and_merging.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	UG_LOGN("9. Resolve potentially generated intersection(s)")
	ResolveTriangleIntersections(g, g.begin<ug::Triangle>(), g.end<ug::Triangle>(), resolveThreshold, aPosition);

	ss << fileName << "_final.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());

}

	/**
	 * @brief the first create_neurite method
	 */
	static void create_neurite
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
    UG_LOGN("nSec: " << nSec);
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
            UG_LOGN("aaPos[v]: " << aaPos[v]);
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
    		const std::vector<size_t>& vBranchInd = brit->bp->vNid;
    		size_t nBranches = vBranchInd.size();

    		for (size_t br = 1; br < nBranches; ++br)
    		{
    			size_t brInd = vBranchInd[br];

    			// get position and radius of first point of branch
    			const number brRadSeg1 = vR[brInd][0];

    			// get position and radius of branching point
    			const number bpTPos = 0.5 * (brit->tend + brit->tstart);
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
    	UG_LOGN("segLength: " << segLength);
    	UG_LOGN("nSeg: " << nSeg);
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

			UG_LOGN("Creating child")
			/// Note:  create prism to connect to in case the branching angle is small or big
			create_neurite(vNeurites, vPos, vR, child_nid, g, aaPos, aaSurfParams, &vrts, &edges, NULL, NULL, bWithER);
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

	/**
	 * @brief splits quadrilateral along edge (parallel)
	 */
	void split_quadrilateral_along_edges
	(
		std::vector<Vertex*> vVrt,
		Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		number percentage,
		ug::vector3 vecDir,
		std::vector<ug::Vertex*>& vertices,
		std::vector<ug::Edge*>& edges,
		bool conservative = true
	)
	{
		/// "middle edges"
		std::vector<ug::Vertex*> from;
		std::vector<ug::Vertex*> to;
		Selector sel(g);
		vVrt.resize(4);
		size_t numPar = 0;
		for (size_t i = 0; i < 4; ++i) {
			ug::vector3 diffVec;
			VecSubtract(diffVec, aaPos[vVrt[i]], aaPos[vVrt[(i+1)%4]]);
			VecNormalize(diffVec, diffVec);
			VecNormalize(vecDir, vecDir);
			UG_LOGN("Parallel? " << VecDot(vecDir, diffVec));
			if (abs(VecDot(vecDir, diffVec)) > (1-0.1)) {
				numPar++;
				UG_LOGN("Parallel:" << VecDot(vecDir, diffVec));
				Edge* e = g.get_edge(vVrt[i], vVrt[(i+1)%4]);
				ug::RegularVertex* newVertex = SplitEdge<ug::RegularVertex>(g, e, conservative);
				ug::vector3 dir;
				VecSubtract(dir, aaPos[vVrt[i]], aaPos[vVrt[(i+1)%4]]);
				VecScaleAdd(aaPos[newVertex], 1.0, aaPos[vVrt[i]], percentage, dir);
				e = g.get_edge(newVertex, vVrt[(i+1)%4]);
				ug::RegularVertex* newVertex2 = SplitEdge<ug::RegularVertex>(g, e, conservative);
				VecScaleAdd(aaPos[newVertex2], 1.0, aaPos[vVrt[(i+1)%4]], -percentage, dir);
				from.push_back(newVertex);
				to.push_back(newVertex2);
			}
		 }


		/// Note: Verify this is correct - seems okay.
		edges.push_back(g.get_edge(to[0], from[0]));
		ug::RegularEdge* e1 = *g.create<RegularEdge>(EdgeDescriptor(to[0], from[1]));
		edges.push_back(e1);

		edges.push_back(g.get_edge(to[1], from[1]));
		ug::RegularEdge* e2 = *g.create<RegularEdge>(EdgeDescriptor(to[1], from[0]));
		edges.push_back(e2);

		vertices.push_back(from[0]);
		vertices.push_back(to[0]);
		vertices.push_back(from[1]);
		vertices.push_back(to[1]);

		UG_COND_THROW(numPar != 2, "Shrinking of connecting quadrilateral failed!");
	}

	/**
	 * @brief callback method to test split quadrilateral with ugshell
	 */
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

		    SaveGridToFile(g, sh, "test_shrunk_geom2_before.ugx");
		    ug::vector3 diffVec;
		    VecSubtract(diffVec, aaPos[vVrt[0]], aaPos[vVrt[(1)%4]]);
		    std::vector<ug::Vertex*> vertices;
		    std::vector<ug::Edge*> edges;
		    split_quadrilateral_along_edges(vVrt, g, aaPos, percentage, diffVec, vertices, edges);
		    SaveGridToFile(g, sh, "test_shrunk_geom2_after.ugx");
	}

	/**
	 * @brief shrinks quadrilateral towards center
	 */
	void shrink_quadrilateral_center
		(
				std::vector<Vertex*>& vVrt,
				Grid& g,
				Grid::VertexAttachmentAccessor<APosition>& aaPos,
				number percentage,
				ug::vector3& center
		)
		{
		    for (size_t i = 0; i < 4; ++i)
		    {
		       ug::vector3 dir;
		       VecSubtract(dir, aaPos[vVrt[i]], center);

		       UG_LOGN("dir:" << dir)
		       VecScaleAdd(aaPos[vVrt[i]], 1.0, aaPos[vVrt[i]], percentage, dir);

		       if (percentage > 1) {
		          UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
		       }
		    }
		}

	/*!
	 * @brief reorders the vertices accordingly
	 */
	void reorder_connecting_elements
	(
	   	std::vector<ug::Vertex*>& v,
		std::vector<ug::Edge*> e
	) {
		std::vector<ug::Vertex*> sorted;
		sorted.push_back(v[0]);

		for (size_t j = 1; j < v.size(); ++j) {
			Vertex* next = v[j];
			for (size_t i = 0; i < e.size(); i++) {
				if ( (e[i]->vertex(0) == next && e[i]->vertex(1) == sorted[j-1]) ||
					 (e[i]->vertex(1) == next && e[i]->vertex(0) == sorted[j-1])) {
					sorted.push_back(next);
					break;
				}
			}
		}
		UG_COND_THROW(sorted.size() != 4, "Did not find vertices to sort...");
		v = sorted;
	}

	/**
	 * @brief shrinks a quadrilateral and creates a copy of the smaller
	 */
	void shrink_quadrilateral_copy
	(
			const std::vector<Vertex*>& vVrt,
			std::vector<Vertex*>& outvVrt,
			const std::vector<Vertex*>& oldVertices,
			std::vector<Edge*>& outvEdge,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			bool createFaces=true,
			ISelector* outSel = NULL,
			ug::vector3* currentDir = NULL
	)
	{
		Selector sel(g);
	    for (size_t i = 0; i < 4; ++i) {
	    	sel.select(vVrt[i]);
		}


	    ug::vector3 center;
	    center = CalculateBarycenter(sel.vertices_begin(), sel.vertices_end(), aaPos);
	    sel.clear();
	    for (size_t i = 0; i < 4; ++i)
	    {
	       ug::Vertex* v = *g.create<RegularVertex>();
	       aaPos[v] = aaPos[vVrt[i]];
	       ug::vector3 dir;
	       VecSubtract(dir, aaPos[vVrt[i]], center);

	       /// TODO: Add the correct shrinkage into a prescribed direction
	       if (currentDir) {
	    	   /// 1. get Edge e starting from i % 4 to i+1 % 4
	    	   /// 2. Check if e is parallel or anti-parallel to currentDir
	    	   /// 3. If true then calculate dir1, dir2 from vertex i, i+1 to center
	    	   /// 4. Project dir1 to edge e if e was parallel to currentDir
	    	   ///    otherwise project dir1 to edge -e if e was antiparallel to currentDir
	    	   /// 5. Project dir2 to edge -e if e was parallel to currentDir
	    	   ///    otherwise project dir2 to edge e if e was antiparallel to currentDir
	       }

	       UG_LOGN("dir:" << dir)
	       VecScaleAdd(aaPos[v], 1.0, aaPos[v], percentage, dir);

	       if (percentage > 1) {
	          UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
	       }
	       outvVrt.push_back(v);
	       if (outSel) outSel->select<ug::Vertex>(v);
	    }

	    /// create new edges in new (small) quad
	    for (size_t i = 0; i < 4; ++i) {
	       ug::Edge* e = *g.create<RegularEdge>(EdgeDescriptor(outvVrt[i], outvVrt[(i+1)%4]));
	       outvEdge.push_back(e);
	       if (outSel) outSel->select<ug::Edge>(e);
	    }

	    if (createFaces) {
	    	/// create new faces
	    	for (size_t i = 0; i < 4; ++i) {
	    		ug::Face* f = *g.create<Quadrilateral>(QuadrilateralDescriptor(outvVrt[i], outvVrt[(i+1)%4], oldVertices[(i+1)%4], oldVertices[i]));
	    		/// Note: Do we need to flip faces here to wrt radial vector?
	    	}
	    }
	}

	/**
	 * @brief shrinks quadrilateral and overwrites old quadrilateral's vertices
	 */
	void shrink_quadrilateral
(
	std::vector<Vertex*> vVrt,
    Grid& g,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	number percentage
)
{
	Selector sel(g);
	vVrt.resize(4);
    for (size_t i = 0; i < 4; ++i) {
    	sel.select(vVrt[i]);
	}

    ug::vector3 center;
    /// Note: Using barycenter, could also use vertex or area centroid for this:
    /// https://en.wikipedia.org/wiki/Quadrilateral#Remarkable_points_and_lines_in_a_convex_quadrilateral
    center = CalculateBarycenter(sel.vertices_begin(), sel.vertices_end(), aaPos);
    sel.clear();

    for (size_t i = 0; i < 4; ++i)
    {
       ug::vector3 dir;
       VecSubtract(dir, aaPos[vVrt[i]], center);
       UG_LOGN("dir:" << dir)
       VecScaleAdd(aaPos[vVrt[i]], 1.0, aaPos[vVrt[i]], percentage, dir);
       if (percentage > 1) {
          UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
       }
    }
}

	/**
	 * @brief callback method to test shrink geometry copy with ugshell
	 */
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

	    SaveGridToFile(g, sh, "test_shrunk_geom_copy_before.ugx");
	    std::vector<ug::Vertex*> vVrtOut;
	    std::vector<ug::Edge*> vEdgeOut;
	    shrink_quadrilateral_copy(vVrt, vVrtOut, vVrtOut, vEdgeOut, g, aaPos, length);
	    SaveGridToFile(g, sh, "test_shrunk_geom_copy_after.ugx");
	}

	/**
	 * @brief callback method to test shrink geometry with ugshell
	 */
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

    SaveGridToFile(g, sh, "test_shrunk_geom_before.ugx");

    ug::vector3 center;
    center = CalculateBarycenter(sel.vertices_begin(), sel.vertices_end(), aaPos);
    sel.clear();
    for (size_t i = 0; i < 4; ++i)
    {
       ug::vector3 dir;
       VecSubtract(dir, aaPos[vVrt[i]], center);
       UG_LOGN("dir:" << dir)
       VecScaleAdd(aaPos[vVrt[i]], 1.0, aaPos[vVrt[i]], length, dir);
       if (VecLength(dir) > VecDistance(center, aaPos[vVrt[i]])) {
    	   UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
       }
    }

    SaveGridToFile(g, sh, "test_shrunk_geom_after.ugx");

}

	/**
	 * @brief generic comparator for SurfaceParams
	 */
	typedef float (NeuriteProjector::SurfaceParams::*membervar);
	template< membervar m > struct CompareBy {
		Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> > m_aaSurfParams;
	   bool operator()( const ug::Vertex* a, const ug::Vertex* b ) const {
	      return m_aaSurfParams[a].*m < m_aaSurfParams[b].*m ;
	   }
	   CompareBy(const Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams) {
		   m_aaSurfParams = aaSurfParams;
	   }
	};

	/**
     * @brief comparator for elements in vector
	 */
	template <typename TElem>
	struct ExistsInVector
	{
		ExistsInVector(const std::vector<TElem>& vec) : m_vec(vec) {
		}

		bool operator() (TElem elem) {
			return (std::find(m_vec.begin(), m_vec.end(), elem) != m_vec.end());
		}
	private:
		const std::vector<TElem>& m_vec;
	};

	/**
	 * @brief Correcting inner branching points of neurites
	 * Note: In case of very small shrinkage factor might result in intersections
	 */
	static void correct_edges
	(
		std::vector<ug::Vertex*>& verts,
		std::vector<ug::Edge*>& edges,
		std::vector<ug::Vertex*>& oldVertsSorted,
		Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		number scale
	)
	{
		 sort(verts.begin(), verts.end(), CompareBy< &NeuriteProjector::SurfaceParams::axial >(aaSurfParams) );
		 oldVertsSorted = verts;

		 Edge* e1 = g.get_edge(verts[0], verts[2]);
		 if (!e1) e1 = g.get_edge(verts[0], verts[3]);
		 Edge* e2 = g.get_edge(verts[1], verts[2]);
		 if (!e2) e2 = g.get_edge(verts[1], verts[3]);

		 /// "bottom vertices of connecting inner face": e1->vertex(0) - newVertex1 - newVertex2 - e1->vertex(1)
		 vector3 dir;
		 VecSubtract(dir, aaPos[e1->vertex(1)], aaPos[e1->vertex(0)]);
		 ug::RegularVertex* newVertex1 = *g.create<ug::RegularVertex>();
		 ug::RegularVertex* newVertex2 = *g.create<ug::RegularVertex>();
		 aaPos[newVertex1] = aaPos[e1->vertex(0)];
		 aaPos[newVertex2] = aaPos[e1->vertex(1)];
		 VecScaleAdd(aaPos[newVertex1], 1.0, aaPos[newVertex1], scale/2.0, dir);
		 aaSurfParams[newVertex1] = aaSurfParams[e1->vertex(0)];
		 aaSurfParams[newVertex1].axial = aaSurfParams[e1->vertex(0)].axial + scale/2*(aaSurfParams[e1->vertex(1)].axial - aaSurfParams[e1->vertex(0)].axial);
		 aaSurfParams[newVertex1].neuriteID = aaSurfParams[e1->vertex(0)].neuriteID;
		 aaSurfParams[newVertex1].scale = aaSurfParams[e1->vertex(0)].scale;
		 VecScaleAdd(aaPos[newVertex2], 1.0, aaPos[newVertex2], -scale/2.0, dir);
		 aaSurfParams[newVertex2] = aaSurfParams[e1->vertex(1)];
		 aaSurfParams[newVertex2].axial = aaSurfParams[e1->vertex(1)].axial - scale/2*(aaSurfParams[e1->vertex(1)].axial - aaSurfParams[e1->vertex(0)].axial);
		 aaSurfParams[newVertex2].neuriteID = aaSurfParams[e1->vertex(1)].neuriteID;
		 aaSurfParams[newVertex2].scale = aaSurfParams[e1->vertex(1)].scale;

		 /// "top vertices of connecting inner face": e2->vertex(0) - newVertex3 - newVertex4 - e2->vertex(1)
		 vector3 dir2;
		 VecSubtract(dir, aaPos[e2->vertex(1)], aaPos[e2->vertex(0)]);
		 VecSubtract(dir2, aaPos[e2->vertex(1)], aaPos[e2->vertex(0)]);
		 ug::RegularVertex* newVertex3 = *g.create<ug::RegularVertex>();
		 ug::RegularVertex* newVertex4 = *g.create<ug::RegularVertex>();
		 aaPos[newVertex3] = aaPos[e2->vertex(0)];
		 aaPos[newVertex4] = aaPos[e2->vertex(1)];
		 VecScaleAdd(aaPos[newVertex3], 1.0, aaPos[newVertex3], scale/2.0, dir);
		 aaSurfParams[newVertex3] = aaSurfParams[e2->vertex(0)];
		 aaSurfParams[newVertex3].axial =  aaSurfParams[e2->vertex(0)].axial + scale/2*(aaSurfParams[e2->vertex(1)].axial - aaSurfParams[e2->vertex(0)].axial);
		 aaSurfParams[newVertex3].neuriteID = aaSurfParams[e2->vertex(0)].neuriteID;
		 aaSurfParams[newVertex3].scale = aaSurfParams[e2->vertex(0)].scale;
		 VecScaleAdd(aaPos[newVertex4], 1.0, aaPos[newVertex4], -scale/2.0, dir);
		 aaSurfParams[newVertex4] = aaSurfParams[e2->vertex(1)];
		 aaSurfParams[newVertex4].axial = aaSurfParams[e2->vertex(1)].axial - scale/2*(aaSurfParams[e2->vertex(1)].axial - aaSurfParams[e2->vertex(0)].axial);
		 aaSurfParams[newVertex4].neuriteID = aaSurfParams[e2->vertex(1)].neuriteID;
		 aaSurfParams[newVertex4].scale = aaSurfParams[e2->vertex(1)].scale;

		 ug::RegularEdge* e31 = *g.create<RegularEdge>(EdgeDescriptor(newVertex1, newVertex3));
		 ug::Quadrilateral* q1 = *g.create<Quadrilateral>(QuadrilateralDescriptor(e1->vertex(0), newVertex1, newVertex3, e2->vertex(0)));
		 ug::RegularEdge* e24 = *g.create<RegularEdge>(EdgeDescriptor(newVertex4, newVertex2));
		 ug::Quadrilateral* q2 = *g.create<Quadrilateral>(QuadrilateralDescriptor(e1->vertex(1), newVertex2, newVertex4, e2->vertex(1)));
		 ug::RegularEdge* e12 =  *g.create<RegularEdge>(EdgeDescriptor(newVertex2, newVertex1));
		 ug::RegularEdge* e43 =  *g.create<RegularEdge>(EdgeDescriptor(newVertex3, newVertex4));

		 /// Verify edges are quasi parallel (should never happen but you never know)
		 VecNormalize(dir, dir);
		 VecNormalize(dir2, dir2);
		 number dotProd = VecDot(dir, dir2) / (VecLength(dir) * VecLength(dir2));
		 UG_COND_THROW( !( fabs(dotProd-1) < SMALL), "Edges need to be quasi parallel during splitting a hexaeder: " << dotProd);

		 /// erase old edges
		 g.erase(e1);
		 g.erase(e2);

		 /// set new face vertices for connection
		 verts.clear();
		 verts.push_back(newVertex1);
		 verts.push_back(newVertex3);
		 verts.push_back(newVertex4);
		 verts.push_back(newVertex2);

		 /// set new edge vertices for connection
		 edges.clear();
		 edges.push_back(e31);
		 edges.push_back(e43);
		 edges.push_back(e24);
		 edges.push_back(e12);
	}

	/**
	 * @brief helper method to correct one side of the quadrilateral
	 */
	static void correct_edges_all
	(
		std::vector<ug::Vertex*>& verts,
		std::vector<ug::Vertex*>& vertsOpp,
		std::vector<ug::Edge*>& edges,
		std::vector<ug::Edge*>& edgesOpp,
		Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		number scale
	)
	{
		/// the connecting vertices are needed later
		UG_LOGN("correcting edges connecting...")
		std::vector<ug::Vertex*> oldVertsSorted;
		correct_edges(verts, edges, oldVertsSorted, aaSurfParams, g, aaPos, scale);
		UG_LOGN("correcting edges opposing...")

		/// backside not needed
		std::vector<ug::Vertex*> oldVertsSortedOpp;
		correct_edges(vertsOpp, edgesOpp, oldVertsSortedOpp, aaSurfParams, g, aaPos, scale);
		g.create<Quadrilateral>(QuadrilateralDescriptor(vertsOpp[0], vertsOpp[1], vertsOpp[2], vertsOpp[3]));

		// connect the sides of splitted edges and fill top center and bottom center holes
		g.create<RegularEdge>(EdgeDescriptor(verts[0], vertsOpp[1]));
		g.create<RegularEdge>(EdgeDescriptor(verts[1], vertsOpp[0]));
		g.create<Quadrilateral>(QuadrilateralDescriptor(verts[0], verts[3], vertsOpp[2], vertsOpp[1]));
		g.create<RegularEdge>(EdgeDescriptor(verts[2], vertsOpp[3]));
		g.create<RegularEdge>(EdgeDescriptor(verts[3], vertsOpp[2]));
		g.create<Quadrilateral>(QuadrilateralDescriptor(verts[1], verts[2], vertsOpp[3], vertsOpp[0]));

		/// fill 4 more holes to the left and right on top and left and right on bottom
		g.create<Quadrilateral>(QuadrilateralDescriptor(verts[0], vertsOpp[1], oldVertsSortedOpp[1], oldVertsSorted[0]));
		g.create<Quadrilateral>(QuadrilateralDescriptor(verts[1], vertsOpp[0], oldVertsSortedOpp[0], oldVertsSorted[1]));
		g.create<Quadrilateral>(QuadrilateralDescriptor(verts[2], vertsOpp[3], oldVertsSortedOpp[3], oldVertsSorted[2]));
		g.create<Quadrilateral>(QuadrilateralDescriptor(verts[3], vertsOpp[2], oldVertsSortedOpp[2], oldVertsSorted[3]));
	}

	/**
	 * @brief corrects the axial offset at the inner branching points
	 * This means, we move the points with smaller axial value further down
	 * the current neurite and the larger axial values further back
	 */
	static void correct_axial_offset
	(
		std::vector<ug::Vertex*>& verts,
		Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		number scale
	)
	{
		// check for consistency
		UG_COND_THROW(verts.size() != 4, "Exactly 4 vertices are necessary on coarse grid level.");
		// sort to find min and max axial values
		sort(verts.begin(), verts.end(), CompareBy< &NeuriteProjector::SurfaceParams::axial >(aaSurfParams) );
		number length = aaSurfParams[verts[2]].axial - aaSurfParams[verts[0]].axial;
		UG_LOGN("length TIMES scale/2: " << length*scale/2)
		// update surface parameters
		aaSurfParams[verts[0]].axial = aaSurfParams[verts[0]].axial + length*scale/2;
		aaSurfParams[verts[1]].axial = aaSurfParams[verts[1]].axial + length*scale/2;
		aaSurfParams[verts[2]].axial = aaSurfParams[verts[2]].axial - length*scale/2;
		aaSurfParams[verts[3]].axial = aaSurfParams[verts[3]].axial - length*scale/2;
	}





	 /**
	 * @brief creates neurites with one inner layer
	 * Note: Could make this more general: Provide std::vector<Layer> instead of
	 * hard-coded one additional layer -> then can add provide additional nestings
	 */
	static void create_neurite_general
(
    const std::vector<NeuriteProjector::Neurite>& vNeurites,
    const std::vector<std::vector<vector3> >& vPos,
    const std::vector<std::vector<number> >& vR,
    size_t nid,
    Grid& g,
    Grid::VertexAttachmentAccessor<APosition>& aaPos,
    Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
    bool firstLayerOnly,
    std::vector<Vertex*>* connectingVrts = NULL,
    std::vector<Edge*>* connectingEdges = NULL,
    std::vector<Vertex*>* connectingVrtsInner = NULL,
    std::vector<Edge*>* connectingEdgesInner = NULL,
    std::vector<Vertex*>* outVerts = NULL,
    std::vector<Vertex*>* outVertsInner = NULL,
    std::vector<number>* outRads = NULL,
    std::vector<number>* outRadsInner = NULL,
    bool forcePositions = false
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

    /// outer layer
    std::vector<Vertex*> vVrt;
    std::vector<Edge*> vEdge;
    vVrt.resize(4);
    vEdge.resize(4);

    /// inner layer
    std::vector<Vertex*> vVrtInner;
    std::vector<Edge*> vEdgeInner;
    vVrtInner.resize(4);
    vEdgeInner.resize(4);

    vector3 vel;
    UG_COND_THROW(nSec == 0, "Number of sections > 0 required. FIX: Don't collapse root edges of neurites.");
    UG_LOGN("nSec: " << nSec);
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
    number angleOffsetInner = 0.0;

    if (connectingVrts && connectingEdges && connectingVrtsInner && connectingEdgesInner && !forcePositions)
    {
        vVrt = *connectingVrts;
        vEdge = *connectingEdges;
        vVrtInner = *connectingVrtsInner;
        vEdgeInner = *connectingEdgesInner;

        // calculate start angle offset
        vector3 center(0.0);
        vector3 center2(0.0);
        for (size_t i = 0; i < 4; ++i)
            VecAdd(center, center, aaPos[(*connectingVrts)[i]]);
        center /= 4;
        for (size_t i = 0; i < 4; ++i)
            VecAdd(center2, center2, aaPos[(*connectingVrtsInner)[i]]);
        center2 /= 4;

        vector3 centerToFirst;
        vector3 centerToFirst2;
        VecSubtract(centerToFirst, aaPos[(*connectingVrts)[0]], center);
        VecSubtract(centerToFirst2, aaPos[(*connectingVrtsInner)[0]], center2);

        vector2 relCoord;
        VecScaleAdd(centerToFirst, 1.0, centerToFirst, -VecProd(centerToFirst, vel), vel);
        relCoord[0] = VecProd(centerToFirst, projRefDir);
        VecScaleAdd(centerToFirst, 1.0, centerToFirst, -relCoord[0], projRefDir);
        relCoord[1] = VecProd(centerToFirst, thirdDir);
        VecNormalize(relCoord, relCoord);

        vector2 relCoord2;
        VecScaleAdd(centerToFirst2, 1.0, centerToFirst2, -VecProd(centerToFirst2, vel), vel);
        relCoord2[0] = VecProd(centerToFirst2, projRefDir);
        VecScaleAdd(centerToFirst2, 1.0, centerToFirst2, -relCoord2[0], projRefDir);
        relCoord2[1] = VecProd(centerToFirst2, thirdDir);
        VecNormalize(relCoord2, relCoord2);

        if (fabs(relCoord[0]) < 1e-8)
            angleOffset = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
        else
            angleOffset = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);
        if (angleOffset < 0) angleOffset += 2.0*PI;

        if (fabs(relCoord2[0]) < 1e-8)
            angleOffsetInner = relCoord2[1] < 0 ? 1.5*PI : 0.5*PI;
        else
            angleOffsetInner = relCoord2[0] < 0 ? PI - atan(-relCoord2[1]/relCoord2[0]) : atan(relCoord2[1]/relCoord2[0]);
        if (angleOffsetInner < 0) angleOffsetInner += 2.0*PI;

        /// note we use the angle offset from the outer, which might be reasonable, since we have the same direction
        angleOffsetInner = angleOffset;

        // ignore first branching region (the connecting region)
        ++brit;
    }
    else
    {
    		// create first layer of vertices and edges for outer layer
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
            	UG_LOGN("aaPos[v]: " << aaPos[v]);
            	if (forcePositions) {
            		aaPos[v] = aaPos[(*connectingVrts)[i]];
            		v = (*connectingVrts)[i];
            		aaSurfParams[v].angular = calculate_angle(pos[0], aaPos[vVrt[0]], aaPos[vVrt[i]]);
            		UG_LOGN("angle: " << aaSurfParams[v].angular);
            	}
        	}
        	outRads->push_back(r[0]);

        	for (size_t i = 0; i < 4; ++i) {
        		vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));
        	}

        	// create first layer of vertices and edges for inner layer
        	if (neurite.bHasER)
        	{
        		for (size_t i = 0; i < 4; ++i)
        		{
        			Vertex* v = *g.create<RegularVertex>();
        			vVrtInner[i] = v;
        			number angle = 0.5*PI*i;
        			VecScaleAdd(aaPos[v], 1.0, pos[0], r[0]*neurite.scaleER*cos(angle), projRefDir, r[0]*neurite.scaleER*sin(angle), thirdDir);
        			/// VecScale(aaPos[v], aaPos[v], neurite.scaleER);
        			UG_LOGN("scale with: " << neurite.scaleER)
        			aaSurfParams[v].neuriteID = nid;
        			aaSurfParams[v].axial = 0.0;
        			aaSurfParams[v].angular = angle;
        			aaSurfParams[v].scale = neurite.scaleER;
        			outVertsInner->push_back(v);
        			UG_LOGN("aaPos[v]: " << aaPos[v]);
        			if (forcePositions) {
        				aaPos[v] = aaPos[(*connectingVrtsInner)[i]];
        				v = (*connectingVrtsInner)[i];
        				aaSurfParams[v].angular = calculate_angle(pos[0], aaPos[vVrt[0]], aaPos[vVrt[i]]);
        			}
        		}
        		if (!forcePositions) shrink_quadrilateral(vVrtInner, g, aaPos, neurite.scaleER);
        		outRadsInner->push_back(r[0]*neurite.scaleER);

        		for (size_t i = 0; i < 4; ++i) {
        			vEdgeInner[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrtInner[i], vVrtInner[(i+1)%4]));
        		}
        	}
    	}

    if (firstLayerOnly) {
    	return;
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

    /// TODO: change angle here too if we forced positions
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
    		const std::vector<size_t>& vBranchInd = brit->bp->vNid;
    		size_t nBranches = vBranchInd.size();

    		for (size_t br = 1; br < nBranches; ++br)
    		{
    			size_t brInd = vBranchInd[br];

    			// get position and radius of first point of branch
    			const number brRadSeg1 = vR[brInd][0];

    			// get position and radius of branching point
    			const number bpTPos = 0.5 * (brit->tend + brit->tstart);
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
    	UG_LOGN("segLength: " << segLength);
    	UG_LOGN("nSeg: " << nSeg);
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
    	Selector sel2(g);
		vector3 extrudeDir;
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
			if (s == nSeg-1 && brit != brit_end) {
				sel.enable_autoselection(true); // for last segment (BP), select new elems
			}

			VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, lastPos);
			Extrude(g, &vVrt, &vEdge, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);
			sel.enable_autoselection(false);

			if (s == nSeg-1 && brit != brit_end) {
				sel2.enable_autoselection(true); // for last segment (BP), select new elems
			}


			Extrude(g, &vVrtInner, &vEdgeInner, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

			sel2.enable_autoselection(false);

			// set new positions and param attachments; also ensure correct face orientation
			for (size_t j = 0; j < 4; ++j)
			{
				number angle = 0.5*PI*j + angleOffset;
				if (angle > 2*PI) angle -= 2*PI;
				Vertex* v = vVrt[j];
				vector3 radialVec;
				VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
				if (!forcePositions) {
				VecAdd(aaPos[v], curPos, radialVec);
				}
				UG_LOGN("aaPos[v] (after extrude (outer)): " << aaPos[v])

				aaSurfParams[v].neuriteID = nid;
				aaSurfParams[v].axial = segAxPos;
				aaSurfParams[v].angular = angle;

				Grid::traits<Face>::secure_container faceCont;
				g.associated_elements(faceCont, vEdge[j]);  // faceCont must contain exactly one face
				vector3 normal;
				CalculateNormal(normal, faceCont[0], aaPos);
				if (VecProd(normal, radialVec) < 0)
					g.flip_orientation(faceCont[0]);

				/// correct angle offset
				if (forcePositions) {
					//aaPos[v] = aaPos[vVrt[j]]; /// TODO: This is wrong, aaPos[v] is wrong in the end because wrong angle used above...
					number angle = calculate_angle(curPos, aaPos[vVrt[0]], aaPos[vVrt[j]]) + angleOffset;
					if (angle > 2*PI) angle -=2*PI;
				}
			}

			// set new positions and param attachments; also ensure correct face orientation
			for (size_t j = 0; j < 4; ++j)
			{
				number angle = 0.5*PI*j + angleOffsetInner;
				if (angle > 2*PI) angle -= 2*PI;
					Vertex* v = vVrtInner[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius*neurite.scaleER*cos(angle), projRefDir, radius*neurite.scaleER*sin(angle), thirdDir);
					if (!forcePositions) {
					VecAdd(aaPos[v], curPos, radialVec);
					}
					/// VecScale(aaPos[v], aaPos[v], neurite.scaleER);
					UG_LOGN("aaPos[v] (after extrude (inner)): " << aaPos[v])

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].scale = neurite.scaleER;

					Grid::traits<Face>::secure_container faceCont;
					g.associated_elements(faceCont, vEdgeInner[j]);  // faceCont must contain exactly one face
					vector3 normal;
					CalculateNormal(normal, faceCont[0], aaPos);
					if (VecProd(normal, radialVec) < 0)
						g.flip_orientation(faceCont[0]);

					/*
					/// correct angle offsets
					if (forcePositions) {
        				aaPos[v] = aaPos[vVrtInner[j]];
        				aaSurfParams[v].angular = calculate_angle(curPos, aaPos[vVrtInner[0]], aaPos[vVrtInner[j]]);
        			}
        			*/
			}

			if (!forcePositions) shrink_quadrilateral(vVrtInner, g, aaPos, neurite.scaleER);
			lastPos = curPos;
    	}

    	UG_LOGN("After extruding...")
    	SaveGridToFile(g, "shit.ugx");




		vector3 currentDir;
    	/// connect branching neurites if present
    	if (brit != brit_end)
		/// Note: create prism to connect to in case the branching angle is small or big
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

			size_t current_nid;
			if (bp->vNid[0] != nid)
				current_nid = bp->vNid[1];
			else
				current_nid = bp->vNid[0];

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

			// now choose best side of hexahedron to connect to (outer)
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
			UG_COND_THROW(!best, "None of the branching point faces pointed in a suitable direction (outer).")
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

			// now choose best side of hexahedron to connect to (inner)
			best = NULL;
			bestProd = 0.0;
			fit = sel2.faces_begin();
			fit_end = sel2.faces_end();
			ug::vector3 hexCenter = CalculateCenter(fit, fit_end, aaPos);

			std::vector<Vertex*> vrtsOpposing;
			for (; fit != fit_end; ++fit)
			{
				CalculateNormal(normal, *fit, aaPos);
				number prod = VecProd(normal, childDir);
				if (prod > bestProd)
				{
					best = *fit;
					bestProd = prod;
				}

				/*
				Face* f = *fit;
				for (size_t i = 0; i < 4; ++i) {
					vrtsOpposing.push_back(f->operator[](i));
				}
				*/
			}
			UG_LOGN("Number of hexaeder verts: " << vrtsOpposing.size());
			UG_COND_THROW(!best, "None of the branching point faces pointed in a suitable direction (inner).")

			// remove face and call creation of child neurite (recursion)
			std::vector<Vertex*> vrtsInner(4);
			UG_COND_THROW(best->num_vertices() != 4, "Hexaeder face does not have 4 vertices!");
			for (size_t j = 0; j < 4; ++j) {
				vrtsInner[j] = best->operator[](j);
			}

			/// find opposing vertices
			fit = sel2.faces_begin();
			fit_end = sel2.faces_end();

			for (; fit != fit_end; ++fit) {
				Face* f = *fit;
				bool opposingFace = true;
				for (size_t i = 0; i < 4; ++i) {
					if(std::find(vrtsInner.begin(), vrtsInner.end(), f->vertex(i)) != vrtsInner.end()) {
						opposingFace = false;
					}
				}
				if (opposingFace) {
					for (size_t i = 0; i < 4; ++i) {
						vrtsOpposing.push_back(f->vertex(i));
					}
				}
			}

			sel2.deselect(sel2.faces_begin(), sel2.faces_end());

/*			vrtsOpposing.erase(std::remove_if(vrtsOpposing.begin(), vrtsOpposing.end(),
					ExistsInVector<Vertex*>(vrtsInner)), vrtsOpposing.end());
			vrtsOpposing.erase(std::remove_if(vrtsOpposing.begin(), vrtsOpposing.end(),
					ExistsInVector<Vertex*>(vVrtInner)), vrtsOpposing.end());
			*/
			UG_COND_THROW(vrtsOpposing.size() != 4, "Hexaeder has to have 4 vertices, but got: " << vrtsOpposing.size());

			std::vector<Edge*> edgesInner(4);
			edgeCont.clear();
			g.associated_elements(edgeCont, best);
			esz = edgeCont.size();

			for (size_t j = 0; j < 4; ++j)
			{
				Vertex* first = vrtsInner[j];
				Vertex* second = vrtsInner[(j+1) % 4];

				size_t k = 0;
				for (; k < esz; ++k)
				{
					if ((edgeCont[k]->vertex(0) == first && edgeCont[k]->vertex(1) == second)
						|| (edgeCont[k]->vertex(0) == second && edgeCont[k]->vertex(1) == first))
					{
						edgesInner[j] = edgeCont[k];
						break;
					}
				}
				UG_COND_THROW(k == esz, "Connecting edges for child neurite could not be determined.");
			}
			g.erase(best);

			/// split hexaeder and correct inner branching points
			std::vector<ug::Edge*> edgesOut;
			std::vector<ug::Edge*> edgesOutOpposing;
			correct_edges_all(vrtsInner, vrtsOpposing, edgesOut, edgesOutOpposing, aaSurfParams, g, aaPos, neurite.scaleER);
			edgesInner = edgesOut;

			UG_LOGN("Creating child(s) for inner and outer...")
			create_neurite_general(vNeurites, vPos, vR, child_nid, g, aaPos, aaSurfParams, false, &vrts, &edges, &vrtsInner, &edgesInner, NULL, NULL, NULL, NULL);

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

    // close the tip of the inner neurite
    vel = vector3(-lastSec.splineParamsX[2], -lastSec.splineParamsY[2], -lastSec.splineParamsZ[2]);
    radius = lastSec.splineParamsR[3];
    VecScale(vel, vel, radius/sqrt(VecProd(vel, vel)));
    VecScale(vel, vel, neurite.scaleER);
    Extrude(g, &vVrtInner, &vEdgeInner, NULL, vel, aaPos, EO_CREATE_FACES, NULL);
    center = CalculateBarycenter(vVrtInner.begin(), vVrtInner.end(), aaPos);
    MergeMultipleVertices(g, vVrtInner.begin(), vVrtInner.end());

    Vertex* vInner = *vVrtInner.begin();
    aaPos[vInner] = center;

    aaSurfParams[vInner].neuriteID = nid;
    aaSurfParams[vInner].axial = 2.0;
    aaSurfParams[vInner].angular = 0.0;
    aaSurfParams[vInner].scale = neurite.scaleER;
}

	/**
	 * @brief exports grid to ugx
	 */
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

	/**
	 * @brief exports grid to swc
	 */
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
        UG_LOGN("Warning: No somatic subset could be identified.")

	if (g.begin<Vertex>() == g.end<Vertex>())
	{
		UG_LOGN("Warning: No vertices contained in grid.")
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
            UG_LOGN("Warning: No soma vertex could be found in the requested neuron.")
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
        UG_LOGN("WARNING: Some vertex type(s) - soma, dendrite, axon, etc. -\n"
            "could not be identified using the subset names.\n"
            << "To ensure correct types in the resulting swc file, the ugx subset names\n"
            "need to contain one of the strings \"SOMA\", \"AXON\", \"DEND\", \"APIC\",\n"
            "upper/lower case can be ignored.");

    outFile.close();
}

	/**
	 * @brief convert swc points to a ug grid
	 */
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
	EraseEmptySubsets(sh);
}

	/**
	 * @brief test function to test smoothing
	 */
	void test_smoothing(const std::string& fileName, size_t n, number h, number gamma, number scale=1.0)
{
	std::vector<SWCPoint> vPoints;
	import_swc(fileName, vPoints, false, scale);


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

	/**
	 * @brief the grid generation test method without scaling and only correcting for small angles at branches (main test method)
	 */
	void test_import_swc(const std::string& fileName, bool correct)
{
	// preconditioning
    test_smoothing(fileName, 5, 1.0, 1.0);

	// read in file to intermediate structure
    std::vector<SWCPoint> vPoints;
    std::vector<SWCPoint> vSomaPoints;
    std::string fn_noext = FilenameWithoutExtension(fileName);
    std::string fn_precond = fn_noext + "_precond.swc";
    import_swc(fn_precond, vPoints, correct, 1.0);

    // convert intermediate structure to neurite data
    std::vector<std::vector<vector3> > vPos;
    std::vector<std::vector<number> > vRad;
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
    std::vector<size_t> vRootNeuriteIndsOut;

    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);
    std::vector<Vertex*> outVerts;
    std::vector<number> outRads;

    // create spline data
    std::vector<NeuriteProjector::Neurite> vNeurites;
    create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

    // create coarse grid
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


    UG_LOGN("do projection handling and generate geom3d")
    ProjectionHandler projHandler(&sh);
    SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
    projHandler.set_geometry(geom3d);
    UG_LOGN("done!")

    SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
    projHandler.set_projector(0, neuriteProj);

    // Note:  This has to be improved: When neurites are copied,
    //        pointers inside still point to our vNeurites array.
    //        If we destroy it, we're in for some pretty EXC_BAD_ACCESSes.
    UG_LOGN("add neurites")
    for (size_t i = 0; i < vNeurites.size(); ++i)
        neuriteProj->add_neurite(vNeurites[i]);
    UG_LOGN("done");

    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
    	create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, NULL, NULL, &outVerts, &outRads, false);
    }

    // at branching points, we have not computed the correct positions yet,
    // so project the complete geometry using the projector
    // Note: little bit dirty; provide proper method in NeuriteProjector to do this
    VertexIterator vit = g.begin<Vertex>();
    VertexIterator vit_end = g.end<Vertex>();
    for (; vit != vit_end; ++vit)
    {
        Edge* tmp = *g.create<RegularEdge>(EdgeDescriptor(*vit,*vit));
        neuriteProj->new_vertex(*vit, tmp);
        g.erase(tmp);
    }

    // create soma
    /*
    sel.clear();
    UG_LOGN("Creating soma!")
    sh.set_default_subset_index(1);
    std::vector<SWCPoint> somaPoint;
    somaPoint.push_back(vSomaPoints[0]);
    create_soma(somaPoint, g, aaPos, sh, 1);
    sh.set_default_subset_index(0);
    UG_LOGN("Done with soma!");
    std::vector<Vertex*> outQuads;
    connect_neurites_with_soma(g, aaPos, outVerts, outVerts, outRads, outQuads, 1, sh, fileName, 1.0);
    UG_LOGN("Done with connecting neurites!");
    */

    // refinement
    AssignSubsetColors(sh);
    sh.set_subset_name("neurites", 0);
    sh.set_subset_name("soma", 1);

    std::string outFileName = FilenameWithoutPath(std::string("testNeuriteProjector.ugx"));
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
    number offset=10;
    try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
    UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

    GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
    for (size_t i = 0; i < 4; ++i)
    {
        ref.refine();
        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
        try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
        UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
    }
}

	/**
	 * @brief get's the first points closest to soma
	 */
	void get_closest_points_to_soma
	(
		const std::string& fn_precond,
		std::vector<ug::vector3>& vPos,
		size_t& lines
	)
	{
		std::ifstream inFile(fn_precond.c_str());
	    UG_COND_THROW(!inFile, "SWC input file '" << fn_precond << "' could not be opened for reading.");

	    size_t lineCnt = 0;
	    std::string line;
	    int somaIndex;
	    std::vector<SWCPoint> swcPoints;
	    while (std::getline(inFile, line)) {
	    	lineCnt++;
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
		    UG_COND_THROW(strs.size() != 7, "Error reading SWC file '" << fn_precond
		           << "': Line " << lineCnt << " does not contain exactly 7 values.");

		   // type
		   if (boost::lexical_cast<int>(strs[6]) == -1) {
			   somaIndex = boost::lexical_cast<int>(strs[0]);
		   } else {
		        if (boost::lexical_cast<int>(strs[6]) == somaIndex) {
		        	number x = boost::lexical_cast<number>(strs[2]);
		        	number y = boost::lexical_cast<number>(strs[3]);
		        	number z = boost::lexical_cast<number>(strs[4]);
		        	vPos.push_back(ug::vector3(x, y, z));
		        }
		   }
	    }
	    lines = lineCnt;
	}

	void get_closest_vertices_on_soma
	(
		const std::vector<ug::vector3>& vPos,
		std::vector<ug::Vertex*>& vPointsSomaSurface, Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		SubsetHandler& sh,
		size_t si
	) {
		UG_LOGN("finding now: " << vPos.size());
				for (size_t i = 0; i < vPos.size(); i++) {
					const ug::vector3* pointSet = &vPos[i];
					ug::vector3 centerOut;
					CalculateCenter(centerOut, pointSet, 1);
					Selector sel(g);
					SelectSubsetElements<Vertex>(sel, sh, si, true);
					UG_LOGN("selected vertices: " << sel.num<Vertex>());
					Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
					Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
					number best = -1;
					ug::Vertex* best_vertex = NULL;
					for (; vit != vit_end; ++vit) {
						number dist = VecDistance(aaPos[*vit], centerOut);
						if (best == -1) {
							best = dist;
							best_vertex = *vit;
						} else if (dist < best) {
							best = dist;
							best_vertex = *vit;
						}
					}
					UG_COND_THROW(!best_vertex, "No best vertex found for root neurite >>" << i << "<<.");
					vPointsSomaSurface.push_back(best_vertex);
				}

	}

	/**
	 * @brief get the closest surface point on soma surface based on triangulation
	 */
	void get_closest_points_on_soma
	(
		const std::vector<ug::vector3>& vPos,
		std::vector<ug::vector3>& vPointsSomaSurface, Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		SubsetHandler& sh,
		size_t si
	)
	{
		UG_LOGN("finding now: " << vPos.size());
		for (size_t i = 0; i < vPos.size(); i++) {
			const ug::vector3* pointSet = &vPos[i];
			ug::vector3 centerOut;
			CalculateCenter(centerOut, pointSet, 1);
			Selector sel(g);
			SelectSubsetElements<Vertex>(sel, sh, si, true);
			UG_LOGN("selected vertices: " << sel.num<Vertex>());
			Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
			Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
			number best = -1;
			ug::Vertex* best_vertex = NULL;
			for (; vit != vit_end; ++vit) {
				number dist = VecDistance(aaPos[*vit], centerOut);
				if (best == -1) {
					best = dist;
					best_vertex = *vit;
				} else if (dist < best) {
					best = dist;
					best_vertex = *vit;
				}
			}
			UG_COND_THROW(!best_vertex, "No best vertex found for root neurite >>" << i << "<<.");
			vPointsSomaSurface.push_back(aaPos[best_vertex]);
		}
	}

	void replace_first_root_neurite_vertex_in_swc
	(
		const size_t& lines,
		const std::string& fn_precond,
		const std::string& fn_precond_with_soma,
		const std::vector<ug::vector3>& vPointsSomaSurface
	) {

		std::ifstream inFile(fn_precond.c_str());
	    UG_COND_THROW(!inFile, "SWC input file '" << fn_precond << "' could not be opened for reading.");
		std::ofstream outFile(fn_precond_with_soma.c_str());
	    UG_COND_THROW(!outFile, "SWC output file '" << fn_precond_with_soma << "' could not be opened for reading.");

	    size_t lineCnt = 1;
	    std::string line;
	    int somaIndex;
	    std::vector<SWCPoint> swcPoints;
        std::vector<number> rads;
        size_t j = 0;
	    while (std::getline(inFile, line)) {
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
			 UG_COND_THROW(strs.size() != 7, "Error reading SWC file '" << fn_precond
			          << "': Line " << lineCnt << " does not contain exactly 7 values.");

			 // type
			 if (boost::lexical_cast<int>(strs[6]) == -1) {
				   somaIndex = boost::lexical_cast<int>(strs[0]);
				    outFile << strs[0] << " " << strs[1] << " " << strs[2] << " "
				        		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
				        		<< strs[6] << std::endl;
			 } else {
				 int index = boost::lexical_cast<int>(strs[6]);
			     if (index == somaIndex) {
			        	number rad = boost::lexical_cast<number>(strs[5]);
			        	outFile << lineCnt << " 3 " << vPointsSomaSurface[j].x() << " "
			        			<< vPointsSomaSurface[j].y() << " " << vPointsSomaSurface[j].z()
			        			<< " " << rad << " " << somaIndex << std::endl;
			        	///lineCnt++;
			        	j++;
			     } else {
			    	 int newIndex = boost::lexical_cast<int>(strs[6]);
			    	 outFile << lineCnt << " " << strs[1] << " " << strs[2] << " "
 	 					   		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
   	 					   		<< newIndex << std::endl;
			     }
			   }
			 lineCnt++;
	    }

	}

	/**
	 * @brief add the new soma surface points to the precondioned swc file.
	 * Note that we assume every dendrite is connected to the ROOT soma point,
	 * and this might, depending on the reconstruction of the SWC file not be true.
	 */
	void add_soma_surface_to_swc
	(
		const size_t& lines,
		const std::string& fn_precond,
		const std::string& fn_precond_with_soma,
		const std::vector<ug::vector3>& vPointsSomaSurface
	)
	{
		std::ifstream inFile(fn_precond.c_str());
	    UG_COND_THROW(!inFile, "SWC input file '" << fn_precond << "' could not be opened for reading.");
		std::ofstream outFile(fn_precond_with_soma.c_str());
	    UG_COND_THROW(!outFile, "SWC output file '" << fn_precond_with_soma << "' could not be opened for reading.");

	    size_t lineCnt = 1;
	    std::string line;
	    int somaIndex;
	    std::vector<SWCPoint> swcPoints;
        std::vector<number> rads;
        size_t j = 0;
	    while (std::getline(inFile, line)) {
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
			 UG_COND_THROW(strs.size() != 7, "Error reading SWC file '" << fn_precond
			          << "': Line " << lineCnt << " does not contain exactly 7 values.");

			 // type
			 if (boost::lexical_cast<int>(strs[6]) == -1) {
				   somaIndex = boost::lexical_cast<int>(strs[0]);
				    outFile << strs[0] << " " << strs[1] << " " << strs[2] << " "
				        		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
				        		<< strs[6] << std::endl;
			 } else {
				 int index = boost::lexical_cast<int>(strs[6]);
			     if (index == somaIndex) {
			        	number rad = boost::lexical_cast<number>(strs[5]);
			        	outFile << lineCnt << " 3 " << vPointsSomaSurface[j].x() << " "
			        			<< vPointsSomaSurface[j].y() << " " << vPointsSomaSurface[j].z()
			        			<< " " << rad << " " << somaIndex << std::endl;
			        	outFile << lineCnt+1 << " " << strs[1] << " " << strs[2] << " "
					   		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
					   		<< lineCnt << std::endl;
			        	lineCnt++;
			        	j++;
			     } else {
			    	 int newIndex = boost::lexical_cast<int>(strs[6])+j;
			    	 outFile << lineCnt << " " << strs[1] << " " << strs[2] << " "
 	 					   		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
   	 					   		<< newIndex << std::endl;
			     }
			   }
			 lineCnt++;
	    }
	}

	/**
	 * @brief the grid generation test method with scaling and ER generation as well as correcting angle
	 */
	void test_import_swc_general(const std::string& fileName, bool correct, number scaleER, bool withER)
{
	UG_LOGN("scaling ER (inner layer) to: " << scaleER);
	UG_COND_THROW(scaleER == 1.0, "scaling to the same size as outer layer is NOT allowed.");
	// preconditioning
    test_smoothing(fileName, 5, 1.0, 1.0);

	// read in file to intermediate structure
    std::vector<SWCPoint> vPoints;
    std::vector<SWCPoint> vSomaPoints;
    std::string fn_noext = FilenameWithoutExtension(fileName);
    std::string fn_precond = fn_noext + "_precond.swc";
    std::string fn_precond_with_soma = fn_noext + "_precond_with_soma.swc";
    std::vector<ug::vector3> vSurfacePoints;
    std::vector<ug::vector3> vPosSomaClosest;
    size_t lines;
    get_closest_points_to_soma(fn_precond, vPosSomaClosest, lines);
    UG_LOGN("got closest points: " << vPosSomaClosest.size());

    // create coarse grid
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

    // convert intermediate structure to neurite data
    std::vector<std::vector<vector3> > vPos;
    std::vector<std::vector<number> > vRad;
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
    std::vector<size_t> vRootNeuriteIndsOut;

    import_swc(fn_precond, vPoints, correct, 1.0);
    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);
    UG_LOGN("converted to neuritelist 1!")

    std::vector<ug::vector3> vPointSomaSurface;
    std::vector<SWCPoint> somaPoint = vSomaPoints;
    /// In new implementation radius can be scaled with 1.0, not 1.05, since vertices are merged
    somaPoint[0].radius *= 1.05;
    create_soma(somaPoint, g, aaPos, sh, 1);
    UG_LOGN("created soma!")
    ///get_closest_points_on_soma(vPosSomaClosest, vPointSomaSurface, g, aaPos, sh, 1);
    std::vector<ug::Vertex*> vPointSomaSurface2;
    get_closest_vertices_on_soma(vPosSomaClosest, vPointSomaSurface2, g, aaPos, sh, 1);
    UG_LOGN("got closest points on soma: " << vPointSomaSurface.size());
    ///add_soma_surface_to_swc(lines, fn_precond, fn_precond_with_soma, vPointSomaSurface);
    std::vector<ug::vector3> newVerts = find_quad_verts_on_soma(g, aaPos, vPointSomaSurface2, vRad, 1, sh, 1.0, vPos.size());
    replace_first_root_neurite_vertex_in_swc(lines, fn_precond, fn_precond_with_soma, newVerts);
    UG_LOGN("added soma points to swc")
    g.clear_geometry();
    import_swc(fn_precond_with_soma, vPoints, correct, 1.0);

    UG_LOGN("converted to neuritelist 2!")
    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

    std::vector<Vertex*> outVerts;
    std::vector<number> outRads;
    std::vector<Vertex*> outVertsInner;
    std::vector<number> outRadsInner;

    // create spline data
    std::vector<NeuriteProjector::Neurite> vNeurites;
    create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

    UG_LOGN("do projection handling and generate geom3d")
    ProjectionHandler projHandler(&sh);
    SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
    projHandler.set_geometry(geom3d);
    UG_LOGN("done!")

    SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
    projHandler.set_projector(0, neuriteProj);

    /// indicate scale and if ER is present
    for (std::vector<NeuriteProjector::Neurite>::iterator it = vNeurites.begin(); it != vNeurites.end(); ++it) {
    	it->bHasER = true;
    	it->scaleER = scaleER;
    }

    // Note:  This has to be improved: When neurites are copied,
    //        pointers inside still point to our vNeurites array.
    //        If we destroy it, we're in for some pretty EXC_BAD_ACCESSes.
    UG_LOGN("add neurites")
    for (size_t i = 0; i < vNeurites.size(); ++i)
        neuriteProj->add_neurite(vNeurites[i]);
    UG_LOGN("done");

    UG_LOGN("generating neurites")
    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
    	create_neurite_general(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, true, NULL, NULL, NULL, NULL, &outVerts, &outVertsInner, &outRads, &outRadsInner);
    }

    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites.ugx");
    sel.clear();


    /// Outer soma
    /// Note: axisVectors (outer soma) and axisVectorsInner (inner soma) save the cylinder center, diameter and length parameters for the CylinderProjectors
    UG_LOGN("Creating soma!")
    sh.set_default_subset_index(1);
    somaPoint = vSomaPoints;
    create_soma(somaPoint, g, aaPos, sh, 1);
    UG_LOGN("Done with soma!");
    std::vector<Vertex*> outQuadsInner;
    std::vector<std::pair<size_t, std::pair<ug::vector3, ug::vector3> > > axisVectors;
    std::vector<std::vector<ug::Vertex*> > connectingVertices(vRootNeuriteIndsOut.size());
   	std::vector<std::vector<ug::Vertex*> > connectingVerticesInner(vRootNeuriteIndsOut.size());
   	std::vector<std::vector<ug::Edge*> > connectingEdges(vRootNeuriteIndsOut.size());
    std::vector<std::vector<ug::Edge*> > connectingEdgesInner(vRootNeuriteIndsOut.size());
    connect_neurites_with_soma(g, aaPos, aaSurfParams, outVerts, outVertsInner, outRads, outQuadsInner, 1, sh, fileName, 1.0, axisVectors, vNeurites, connectingVertices, connectingVerticesInner, connectingEdges, connectingEdgesInner,true);

    /*for (size_t i = 0; i < 2; i++) {
    	reorder_connecting_elements(connectingVertices[i], connectingEdges[i]);
    	reorder_connecting_elements(connectingVerticesInner[i], connectingEdgesInner[i]);
    }
    */

    /// delete old vertices from incorrect neurite starts
    g.erase(outVerts.begin(), outVerts.end());
    g.erase(outVertsInner.begin(), outVertsInner.end());
    outVerts.clear();
    outVertsInner.clear();
    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_finding_initial_edges.ugx");

    sh.set_default_subset_index(0);
    UG_LOGN("generating neurites")

    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
    	///create_neurite_general(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, false, &connectingVertices[i], &connectingEdges[i], &connectingVerticesInner[i], &connectingEdgesInner[i], &outVerts, &outVertsInner, &outRads, &outRadsInner, false);
    	create_neurite_general(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, false, NULL, NULL, NULL, NULL, &outVerts, &outVertsInner, &outRads, &outRadsInner, false);
    }

    /// connect_neurites_with_soma(g, aaPos, aaSurfParams, outVerts, outVertsInner, outRads, outQuadsInner, 1, sh, fileName, 1.0, axisVectors, vNeurites, connectingVertices, connectingVerticesInner, connectingEdges, connectingEdgesInner,true);

    /// Inner soma
    UG_LOGN("Done with connecting neurites!");
    UG_LOGN("Creating soma inner!")
    somaPoint.front().radius = somaPoint.front().radius * scaleER;
    size_t newSomaIndex = sh.num_subsets();
    create_soma(somaPoint, g, aaPos, sh, newSomaIndex, 2);
    UG_LOGN("Done with soma inner!");
    std::vector<Vertex*> outQuadsInner2;
    UG_LOGN("Size of outQuadsInner: " << outQuadsInner.size())
    std::vector<std::pair<size_t, std::pair<ug::vector3, ug::vector3> > > axisVectorsInner;
    connect_neurites_with_soma(g, aaPos, aaSurfParams, outVerts, outVertsInner, outRadsInner, outQuadsInner2, newSomaIndex, sh, fileName, scaleER, axisVectorsInner, vNeurites, connectingVertices, connectingVerticesInner, connectingEdges, connectingEdgesInner, false);

    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
    	vNeurites[i].somaRadius = somaPoint.front().radius;
    	vNeurites[i].somaPt = somaPoint.front().coords;
    }

    // save after connecting and assign subsets
    EraseEmptySubsets(sh);
    AssignSubsetColors(sh);
    sh.set_subset_name("neurites (all)", 0);
    sh.set_subset_name("soma (outer)", 1);
    for (size_t i = 2; i < newSomaIndex; i++) {
    	std::stringstream ss;
     	ss << "outer-connex #" << i;
     	sh.set_subset_name(ss.str().c_str(), i);
    }


    sh.set_subset_name("soma (inner)", newSomaIndex);
    for (size_t i = newSomaIndex+1; i < sh.num_subsets(); i++) {
    	std::stringstream ss;
     	ss << "inner-connex #" << i;
     	sh.set_subset_name(ss.str().c_str(), i);
    }
    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_not_forced.ugx");

    /// Double Vertices might occur (during Qhull gen faces) -> remove these here to be sure
    RemoveDoubles<3>(g, g.begin<Vertex>(), g.end<Vertex>(), aaPos, 0.0001);

    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_not_forced_without_doubles.ugx");

    EraseEmptySubsets(sh);
    AssignSubsetColors(sh);
    for (size_t i = newSomaIndex+vRootNeuriteIndsOut.size(); i < sh.num_subsets(); i++) {
    	std::stringstream ss;
     	ss << "inter-soma-connex #" << i;
     	sh.set_subset_name(ss.str().c_str(), i);
    }

    /// connect now inner soma to neurites (This is the old strategy: Neew strategy is to project)
    connect_inner_neurites_to_inner_soma(newSomaIndex, vRootNeuriteIndsOut.size(), g, aaPos, sh);

    /// connect now outer soma to root neurites
    connect_outer_and_inner_root_neurites_to_outer_soma(1, vRootNeuriteIndsOut.size(), g, aaPos, sh, outVerts, outVertsInner);
    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites_and_connecting_all.ugx");
    return;

    // Note: At branching points, we have not computed the correct positions yet.
    // By adding new edges (point, point) which collapse into a vertex we force a
    // projection refinement on the initial geometry. This is a little hack.
    // TODO: Provide proper method in NeuriteProjector to do this not here
    VertexIterator vit = sh.begin<Vertex>(0);
    VertexIterator vit_end = sh.end<Vertex>(0);
    for (; vit != vit_end; ++vit)
    {
        Edge* tmp = *g.create<RegularEdge>(EdgeDescriptor(*vit,*vit));
        neuriteProj->new_vertex(*vit, tmp);
        g.erase(tmp);
    }

    // refinement
    std::string outFileName = FilenameWithoutPath(std::string("testNeuriteProjector.ugx"));
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
    number offset=somaPoint.front().radius * 2 * 10; // 10 times the diameter of the soma
    try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
    UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

    GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
    size_t numRefs = 1;
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

	/**
	 * @brief the grid generation test method with scaling and ER generation as well as correcting angle
	 */
	void test_import_swc_general_smooth(const std::string& fileName, bool correct, number scaleER, bool withER)
{
	UG_LOGN("scaling ER (inner layer) to: " << scaleER);
	UG_COND_THROW(scaleER == 1.0, "scaling to the same size as outer layer is NOT allowed.");
	// preconditioning
    test_smoothing(fileName, 5, 1.0, 1.0);

	// read in file to intermediate structure
    std::vector<SWCPoint> vPoints;
    std::vector<SWCPoint> vSomaPoints;
    std::string fn_noext = FilenameWithoutExtension(fileName);
    std::string fn_precond = fn_noext + "_precond.swc";
    import_swc(fn_precond, vPoints, correct, 1.0);

    // convert intermediate structure to neurite data
    std::vector<std::vector<vector3> > vPos;
    std::vector<std::vector<number> > vRad;
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
    std::vector<size_t> vRootNeuriteIndsOut;

    std::vector<Vertex*> outVerts;
    std::vector<number> outRads;
    std::vector<Vertex*> outVertsInner;
    std::vector<number> outRadsInner;


    // create coarse grid
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

    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);
    UG_LOGN("Creating soma!")
    sh.set_default_subset_index(1);
    std::vector<SWCPoint> somaPoint = vSomaPoints;
    create_soma(somaPoint, g, aaPos, sh, 1);
    UG_LOGN("Done with soma!");


    // create spline data
    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);
    std::vector<NeuriteProjector::Neurite> vNeurites;
    create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

    UG_LOGN("do projection handling and generate geom3d")
    ProjectionHandler projHandler(&sh);
    SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
    projHandler.set_geometry(geom3d);
    UG_LOGN("done!")

    SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
    projHandler.set_projector(0, neuriteProj);

    /// indicate scale and if ER is present
    for (std::vector<NeuriteProjector::Neurite>::iterator it = vNeurites.begin(); it != vNeurites.end(); ++it) {
    	it->bHasER = true;
    	it->scaleER = scaleER;
    }

    // Note:  This has to be improved: When neurites are copied,
    //        pointers inside still point to our vNeurites array.
    //        If we destroy it, we're in for some pretty EXC_BAD_ACCESSes.
    UG_LOGN("add neurites")
    for (size_t i = 0; i < vNeurites.size(); ++i)
        neuriteProj->add_neurite(vNeurites[i]);
    UG_LOGN("done");

    UG_LOGN("generating neurites")
    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
    	create_neurite_general(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, false, NULL, NULL, NULL, NULL, &outVerts, &outVertsInner, &outRads, &outRadsInner);
    }

    SaveGridToFile(g, sh, "testNeuriteProjector_after_adding_neurites.ugx");

    // at branching points, we have not computed the correct positions yet,
    // so project the complete geometry using the projector -> for inner neurite this fails at some points -> needs to be addressed.
    // Note: little bit dirty; provide proper method in NeuriteProjector to do this
    VertexIterator vit = sh.begin<Vertex>(0);
    VertexIterator vit_end = sh.end<Vertex>(0);
    for (; vit != vit_end; ++vit)
    {
        Edge* tmp = *g.create<RegularEdge>(EdgeDescriptor(*vit,*vit));
        neuriteProj->new_vertex(*vit, tmp);
        g.erase(tmp);
    }

    // refinement
    std::string outFileName = FilenameWithoutPath(std::string("testNeuriteProjector.ugx"));
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
    number offset=10;
    try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
    UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

    GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
    for (size_t i = 0; i < 4; ++i)
    {
        ref.refine();
        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
        try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
        UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
    }
}

	/**
	 * @brief the grid generation test method with scaling only and correcting angles
	 */
	void test_import_swc_scale(const std::string& fileName, bool correct, number scale)
{
	// preconditioning
    test_smoothing(fileName, 5, 1.0, 1.0, 1.0);

	// read in file to intermediate structure
    std::vector<SWCPoint> vPoints;
    std::vector<SWCPoint> vSomaPoints;
    std::string fn_noext = FilenameWithoutExtension(fileName);
    std::string fn_precond = fn_noext + "_precond.swc";
    import_swc(fn_precond, vPoints, correct, 1.0);

    // convert intermediate structure to neurite data
    std::vector<std::vector<vector3> > vPos;
    std::vector<std::vector<number> > vRad;
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
    std::vector<size_t> vRootNeuriteIndsOut;

    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);
    std::vector<Vertex*> outVerts;
    std::vector<number> outRads;

    // create spline data (and scale radii before)
    std::vector<NeuriteProjector::Neurite> vNeurites;
    for (std::vector<std::vector<number> >::iterator it = vRad.begin(); it != vRad.end(); ++it) {
         	for (std::vector<number>::iterator itRad = it->begin(); itRad != it->end(); ++itRad) {
          		*itRad = *itRad * scale;
         	}
       }
    create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

    // create coarse grid
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

    UG_LOGN("do projection handling and generate geom3d")
    ProjectionHandler projHandler(&sh);
    SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
    projHandler.set_geometry(geom3d);
    UG_LOGN("done!")

    SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
    projHandler.set_projector(0, neuriteProj);

    UG_LOGN("add neurites")
    for (size_t i = 0; i < vNeurites.size(); ++i)
        neuriteProj->add_neurite(vNeurites[i]);
    UG_LOGN("done");

    UG_LOGN("generating neurites")
    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
    	create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, NULL, NULL, &outVerts, &outRads, false);
    }

    VertexIterator vit = g.begin<Vertex>();
    VertexIterator vit_end = g.end<Vertex>();
    for (; vit != vit_end; ++vit)
    {
        Edge* tmp = *g.create<RegularEdge>(EdgeDescriptor(*vit,*vit));
        neuriteProj->new_vertex(*vit, tmp);
        g.erase(tmp);
    }

    // create soma
    sel.clear();
    UG_LOGN("Creating soma!")
    sh.set_default_subset_index(1);
    create_soma(vSomaPoints, g, aaPos, sh, 1);
    sh.set_default_subset_index(0);
    UG_LOGN("Done with soma!");

    // connect soma with neurites
    std::vector<Vertex*> outQuads;
    std::vector<std::pair<size_t, std::pair<ug::vector3, ug::vector3> > > axisVectors;
    std::vector<std::vector<ug::Vertex*> > connectingVertices;
   	std::vector<std::vector<ug::Vertex*> > connectingVerticesInner;
   	std::vector<std::vector<ug::Edge*> > connectingEdges;
    std::vector<std::vector<ug::Edge*> > connectingEdgesInner;
    connect_neurites_with_soma(g, aaPos, aaSurfParams, outVerts, outVerts, outRads, outQuads, 1, sh, fileName, 1.0, axisVectors, vNeurites, connectingVertices, connectingVerticesInner, connectingEdges, connectingEdgesInner, false);
    UG_LOGN("Done with connecting neurites!");

    // refinement
    AssignSubsetColors(sh);
    sh.set_subset_name("neurites", 0);
    sh.set_subset_name("soma", 1);

    std::string outFileName = FilenameWithoutPath(std::string("testNeuriteProjector.ugx"));
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
    number offset=10;
    try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
    UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

    GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
    for (size_t i = 0; i < 4; ++i)
    {
        ref.refine();
        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
        try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
        UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");
    }
}

	/**
	 * @brief neurite projector test method for four section tube
	 */
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

    std::vector<NeuriteProjector::Neurite> vNeurites;
    create_spline_data_for_neurites(vNeurites, vPos, vR);
    NeuriteProjector::Neurite& neurite = vNeurites[0];


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
    neuriteProj->add_neurite(neurite);
    projHandler.set_projector(0, neuriteProj);

    std::vector<Vertex*> outVerts;
    create_neurite(vNeurites, vPos, vR, 0, g, aaPos, aaSurfParams);


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
        //sel.select(g.begin<Face>(), g.end<Face>());
        //Refine(g, sel, &projHandler);

        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        std::string curFileName = fileName.substr(0, fileName.size()-4) + oss.str();
        number offset = 1;
        try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
        UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

        //GridWriterUGX ugxWriter;
        //ugxWriter.add_grid(g, "defGrid", aPosition);
        //ugxWriter.add_subset_handler(sh, "defSH", 0);
        //ugxWriter.add_projection_handler(projHandler, "defPH", 0);
        //if (!ugxWriter.write_to_file(fileName.c_str()))
        //    UG_THROW("Grid could not be written to file '" << fileName << "'.");
    }
}

	/**
	 * @brief neurite projector test method for sections with tubes and branches
	 */
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

    std::vector<NeuriteProjector::Neurite> vNeurites(2);
    create_spline_data_for_neurites(vNeurites, vPos, vR, &bpInfo);

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

    for (size_t i = 0; i < vNeurites.size(); ++i)
        neuriteProj->add_neurite(vNeurites[i]);

    create_neurite(vNeurites, vPos, vR, 0, g, aaPos, aaSurfParams, NULL, NULL);
    //create_neurite_general(vNeurites, vPos, vR, 0, g, aaPos, aaSurfParams, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    // at branching points, we have not computed the correct positions yet,
    // so project the complete geometry using the projector
    // Note: little bit dirty; provide proper method in NeuriteProjector to do this
    VertexIterator vit = g.begin<Vertex>();
    VertexIterator vit_end = g.end<Vertex>();
    for (; vit != vit_end; ++vit)
    {
        Edge* tmp = *g.create<RegularEdge>(EdgeDescriptor(*vit,*vit));
        neuriteProj->new_vertex(*vit, tmp);
        g.erase(tmp);
    }

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

        UG_LOGN("refinement step " << i);
        ref.refine();
        //Refine(g, sel, &projHandler);

        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        std::string curFileName = fileName.substr(0, fileName.size()-4) + oss.str();
        number offset = 1;
        try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
        UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

        //GridWriterUGX ugxWriter;
        //ugxWriter.add_grid(g, "defGrid", aPosition);
        //ugxWriter.add_subset_handler(sh, "defSH", 0);
        //ugxWriter.add_projection_handler(projHandler, "defPH", 0);
        //if (!ugxWriter.write_to_file(fileName.c_str()))
        //    UG_THROW("Grid could not be written to file '" << fileName << "'.");
    }
}

	/**
	  * @brief top level vertex repositioning function for neurite projection
	  */
	void apply_neurite_projector(MultiGrid& mg, SmartPtr<NeuriteProjector> neuriteProj)
	{
	// define attachment accessors
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);

	// move manifold vertices to their projected positions
	VertexIterator vrtIter = mg.begin<Vertex>(mg.top_level());
	VertexIterator vrtEnd = mg.end<Vertex>(mg.top_level());
	for (; vrtIter != vrtEnd; ++vrtIter)
	{
		GridObject* go = mg.get_parent(*vrtIter);

		if (go->base_object_id() == EDGE)
		{
			Edge* par = dynamic_cast<Edge*>(go);
			UG_ASSERT(par, "Object with base object id EDGE is not an edge.");
			neuriteProj->new_vertex(*vrtIter, par);
		}
		else if (go->base_object_id() == FACE)
		{
			Face* par = dynamic_cast<Face*>(go);
			UG_ASSERT(par, "Object with base object id FACE is not a face.");
			neuriteProj->new_vertex(*vrtIter, par);
		}
		/// TODO: treat the volume case!
		}
	}

	/**
	 * @brief test method to test the experimental cylinder volume projector
	 */
	void test_cylinder_volume_projector()
{
    // grid preparation
    Grid g;
    SubsetHandler sh(g);
    sh.set_default_subset_index(0);
    g.attach_to_vertices(aPosition);
    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
    Selector sel(g);

    // projection handling
    ProjectionHandler projHandler(&sh);
    SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
    projHandler.set_geometry(geom3d);
    SmartPtr<CylinderVolumeProjector> cylVolProj(new CylinderVolumeProjector(geom3d, vector3(0,0,0), vector3(0,0,1)));
    projHandler.set_projector(0, cylVolProj);

    // create quadrilateral
    Vertex* v0 = *g.create<RegularVertex>();
    Vertex* v1 = *g.create<RegularVertex>();
    Vertex* v2 = *g.create<RegularVertex>();
    Vertex* v3 = *g.create<RegularVertex>();

    aaPos[v0] = vector3(1,1,0);
    aaPos[v1] = vector3(-1,1,0);
    aaPos[v2] = vector3(-1,-1,0);
    aaPos[v3] = vector3(1,-1,0);

    *g.create<Quadrilateral>(QuadrilateralDescriptor(v0, v1, v2, v3));

    // refinement
    AssignSubsetColors(sh);

    std::string fileName = FilenameWithoutPath(std::string("testCylVolProj.ugx"));
    GridWriterUGX ugxWriter;
    ugxWriter.add_grid(g, "defGrid", aPosition);
    ugxWriter.add_subset_handler(sh, "defSH", 0);
    ugxWriter.add_projection_handler(projHandler, "defPH", 0);
    if (!ugxWriter.write_to_file(fileName.c_str()))
        UG_THROW("Grid could not be written to file '" << fileName << "'.");

    Domain3d dom;
    try {LoadDomain(dom, fileName.c_str());}
    UG_CATCH_THROW("Failed loading domain from '" << fileName << "'.");
    HangingNodeRefiner_MultiGrid ref(*dom.grid(), dom.refinement_projector());
    SmartPtr<MultiGrid> mg = dom.grid();
    for (size_t i = 0; i < 6; ++i)
    {
        int topLv = mg->num_levels() - 1;
        FaceIterator it = dom.grid()->begin<Face>(topLv);
        FaceIterator it_end = dom.grid()->end<Face>(topLv);
        for (; it!= it_end; ++it)
            ref.mark(*it, RM_REFINE);

        UG_LOGN("refinement step " << i);
        ref.refine();
        //Refine(g, sel, &projHandler);

        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        std::string curFileName = fileName.substr(0, fileName.size()-4) + oss.str();
        number offset = 1;
        try {SaveGridHierarchyTransformed(*dom.grid(), *dom.subset_handler(), curFileName.c_str(), offset);}
        UG_CATCH_THROW("Grid could not be written to file '" << curFileName << "'.");

        //ugxWriter.add_grid(g, "defGrid", aPosition);
        //ugxWriter.add_subset_handler(sh, "defSH", 0);
        //ugxWriter.add_projection_handler(projHandler, "defPH", 0);
        //if (!ugxWriter.write_to_file(fileName.c_str()))
        //    UG_THROW("Grid could not be written to file '" << fileName << "'.");
    	}
	}

	/**
	 * @brief test axis shrinkage
	 */
	void test_shrinkage() {
		std::vector<ug::Vertex*> verts;
		Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> > aaSurfParams;
		Grid::VertexAttachmentAccessor<APosition> aaPos;
		correct_axial_offset(verts, aaSurfParams, aaPos, 0.5);
	}

	} // namespace neuro_collection
} // namespace ug
