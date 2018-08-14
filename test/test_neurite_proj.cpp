/*
 *
 * test_neurite_proj.cpp
 *
 *  Created on: 27.12.2016
 *      Author: mbreit
 */

#include "test_neurite_proj.h"
#include "lib_grid/refinement/projectors/projection_handler.h"
#include "../../ugcore/ugbase/lib_grid/refinement/projectors/cylinder_volume_projector.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"
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

void import_swc
(
    const std::string& fileName,
    std::vector<SWCPoint>& vPointsOut,
    bool correct)
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
        pt.radius = boost::lexical_cast<number>(strs[5]);

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



struct EdgeLengthCompare
{
	bool operator()(const std::pair<Edge*, number> e1, const std::pair<Edge*, number> e2)
	{return e1.second > e2.second;}
};


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

        // FIXME: find suitable permissible render vector
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
                // TODO: This might need some patching up
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
		if (integral >= (seg+1)*segLength)
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

static void create_soma
(
		const std::vector<SWCPoint>& somaPts,
		Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		SubsetHandler& sh
)
{
	UG_COND_THROW(somaPts.size() != 1, "Currently only one soma point is allowed by this implementation");
	GenerateIcosphere(g, somaPts.front().coords, somaPts.front().radius, 2, aPosition);
}

static void connect_neurites_with_soma
(
	   Grid& g,
	   Grid::VertexAttachmentAccessor<APosition>& aaPos,
	   std::vector<Vertex*> outVerts,
	   std::vector<number> outRads,
	   size_t si,
	   SubsetHandler& sh,
	   const std::string& fileName
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
		quads.push_back(temp);
	}

	UG_LOGN("2. Calculate center of each quad, find next surface vertex on soma.")
	/// 2. Berechne den Schwerpunkt jedes Quads und finde den nächstgelegenen
	///    Vertex des Oberflächengitters vom Soma
	std::vector<Vertex*> bestVertices;
	for (size_t i = 0; i < numQuads; i++) {
		const ug::vector3* pointSet = &(quads[i][0]);
		ug::vector3 centerOut;
		CalculateCenter(centerOut, pointSet, numVerts);
		Selector sel(g);
		SelectSubsetElements<Vertex>(sel, sh, 1, true);
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
		AdaptSurfaceGridToCylinder(sel, g, bestVertices[i], normal, radius, 1, aPosition);
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
	size_t beginningOfQuads = 2; // subset index where quads are stored in
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
	}
	ss << fileName << "_after_merging_cylinder_vertices.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	/// TODO: parameterize
	UG_LOGN("8. TangentialSmooth");
	TangentialSmooth(g, g.vertices_begin(), g.vertices_end(), aaPos, 0.01, 10);

	UG_LOGN("6. Extrude rings along normal")
	/// 6. Extrudiere die Ringe entlang ihrer Normalen mit Höhe 0 (Extrude mit
	///    aktivierter create faces Option).
	sel.clear();
	std::vector<std::vector<Vertex*> > somaVerts;
	std::vector<std::vector<Vertex*> > allVerts;
	for (size_t i = 0; i < numQuads; i++) {
		size_t si = beginningOfQuads+i;
		ug::vector3 normal;
		CalculateVertexNormal(normal, g, *sh.begin<Vertex>(si), aaPos);
		std::vector<Edge*> edges;
		edges.assign(sh.begin<Edge>(si), sh.end<Edge>(si));
		SelectSubsetElements<Vertex>(sel, sh, si, true);
		std::vector<Vertex*> temp;
		temp.assign(sel.vertices_begin(), sel.vertices_end());
		somaVerts.push_back(temp);
		Extrude(g, NULL, &edges, NULL, normal, aaPos, EO_CREATE_FACES, NULL);
		sel.clear();
		SelectSubsetElements<Vertex>(sel, sh, si, true);
		std::vector<Vertex*> temp2;
		temp2.assign(sel.vertices_begin(), sel.vertices_end());
		allVerts.push_back(temp2);
		sel.clear();
	}

	/// delete common vertices, thus keep only newly extruded vertices (useed in next step 7.)
	for (size_t i = 0; i < numQuads; i++) {
		size_t numSomaVerts = somaVerts[i].size();
		for (size_t j = 0; j < numSomaVerts; j++) {
			allVerts[i].erase(std::remove(allVerts[i].begin(), allVerts[i].end(), somaVerts[i][j]), allVerts[i].end());
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
	for (size_t i = 0; i < numQuads; i++) {
		std::vector<ug::vector3> temp;
		std::vector<Vertex*> temp2;
		for (size_t j = 0; j < numVerts; j++) {
			temp.push_back(aaPos[outVerts[i*4+j]]);
		}
		for (size_t j = 0; j < numVerts; j++) {
			temp.push_back(aaPos[allVerts[i][j]]);
		}

		UG_COND_THROW(temp.size() != 8, "Need 8 vertices for calculating all faces.");
		#ifdef NC_WITH_QHULL
			using ug::neuro_collection::convexhull::gen_face;
			gen_face(temp, g, sh, si+i, aaPos);
		#else
			using ug::neuro_collection::quickhull::gen_face;
			gen_face(temp, temp2, g, sh, si+i, aaPos);
		#endif
	}
	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);
	ss << fileName << "_after_extruding_cylinders_and_merging.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
	ss.str(""); ss.clear();

	/// TODO: parameterize
	UG_LOGN("9. Resolve intersections")
	ResolveTriangleIntersections(g, g.begin<ug::Triangle>(), g.end<ug::Triangle>(), 0.00001, aPosition);

	ss << fileName << "_final.ugx";
	SaveGridToFile(g, sh, ss.str().c_str());
}




/// Could just copy these method 2 times:
/// TODO: save segment lengths (calculate_length_over_radius) with unique index (or pointer to neurite)
/// then use create_neurite 2 times, and in the ER call (cell-within-cell) rely on the segment length
/// which we can get from the index map with the unique index or pointer of the neurite

/// create two times the spline data -> then call this (create_neurite) method 2 times:
/// once regularly, secondly, with the the vCellwithinCell and vNeurites, but we will
/// use the segment length from vNeurites not from vCellWithinCells! (this should work) -> same number of segments!
/// make sure that vCellWithinCell and vNeurites have the same number of sections!, if not, this idea won't work
/// before creating second spline data we have to scale the radii..

/// if it does not work, we can try to, create spline data once, and modify the radius in the create_neurite method!
/// and essentially duplicate the create_neurite code to call 2 times -> could refactor in another method and call from create_neurite,
/// e.g. a create_neurite_helper method.

/// if this does not work, use own strategy: iterate over cuboids and scale them

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
    std::vector<number>* outRads = NULL
)
{
	/// TODO: use vNeurites, vNeuritesER... (scale before vR) -> then assure length_over_radius is the same
	/// TODO: split up in two parts: neurite + er + etc -> refactor this into a general method which get's called 3 times, for neurite, er and etc
	/// If this is not working then use own strategy -> scale coarse grid of membrane to ER
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
    vVrt.resize(4);
    vEdge.resize(4);

    vector3 vel;
    UG_COND_THROW(nSec == 0, "Number of sections > 0 required. FIX: Don't collapse root edges of neurites.");
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
    	number lengthOverRadius;
    /*	if (vNeuritesInner != NULL) {
    		const NeuriteProjector::Neurite& neuriteInner= (*vNeuritesInner)[nid];
    		calculate_length_over_radius(t_start, t_end, neuriteInner, curSec);
    	} else {
    		calculate_length_over_radius(t_start, t_end, neurite, curSec);
    	}*/
   		lengthOverRadius = calculate_length_over_radius(t_start, t_end, neurite, curSec);

    	size_t nSeg = (size_t) floor(lengthOverRadius / 8);
    	// at least one segment is required to create a neurite
    	if (nSeg == 0) { nSeg = 1; }
    	UG_COND_THROW(nSeg == 0, "Number of segments > 0 required.");
    	number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
    	UG_LOGN("segLength: " << segLength);
    	UG_LOGN("nSeg: " << nSeg);
    	std::vector<number> vSegAxPos(nSeg);
    	/*
    	if (vNeuritesInner != NULL) {
    		const NeuriteProjector::Neurite& neuriteInner= (*vNeuritesInner)[nid];
    		calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neuriteInner, curSec, segLength);
    	} else {
    		calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite, curSec, segLength);
    	}*/
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
				aaSurfParams[v].axial = segAxPos * 2;
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
			// TODO: create prism to connect to in case the branching angle is small or big
			create_neurite(vNeurites, vPos, vR, child_nid, g, aaPos, aaSurfParams, &vrts, &edges);
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


void test_smoothing(const std::string& fileName, size_t n, number h, number gamma)
{
	std::vector<SWCPoint> vPoints;
	import_swc(fileName, vPoints, false);


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


void test_import_swc(const std::string& fileName, bool correct, number scaleER)
{
	// preconditioning
    test_smoothing(fileName, 5, 1.0, 1.0);

	// read in file to intermediate structure
    std::vector<SWCPoint> vPoints;
    std::vector<SWCPoint> vSomaPoints;
    std::string fn_noext = FilenameWithoutExtension(fileName);
    std::string fn_precond = fn_noext + "_precond.swc";
    import_swc(fn_precond, vPoints, correct);

    // convert intermediate structure to neurite data
    std::vector<std::vector<vector3> > vPos;
    std::vector<std::vector<number> > vRad;
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
    std::vector<size_t> vRootNeuriteIndsOut;

    convert_pointlist_to_neuritelist(vPoints, vSomaPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);
    std::vector<Vertex*> outVerts;
    std::vector<number> outRads;

/* debug
    std::cout << "BPInfo:" << std::endl;
    for (size_t i = 0; i < vBPInfo.size(); ++i)
    {
        std::cout << "nid " << i << ":  " << "#vrts = " << vPos[i].size() << std::endl;
        for (size_t j = 0; j < vBPInfo[i].size(); ++j)
        {
            std::cout << "  " << vBPInfo[i][j].first << ": ";
            for (size_t k = 0; k < vBPInfo[i][j].second.size(); ++k)
            	std::cout << vBPInfo[i][j].second[k] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
*/

    // create spline data
    std::vector<NeuriteProjector::Neurite> vNeurites;
    create_spline_data_for_neurites(vNeurites, vPos, vRad, &vBPInfo);

    std::vector<NeuriteProjector::Neurite> vNeuritesWithin;
    std::vector<std::vector<number> > vRadInner = vRad;
    /// TODO: scale vRad! -> and fix previously introduced bug...
    /// create_spline_data_for_neurites(vNeuritesWithin, vPos, vRadInner, &vBPInfo);

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

    /// indicate scale and if ER is present
    for (std::vector<NeuriteProjector::Neurite>::iterator it = vNeurites.begin(); it != vNeurites.end(); ++it) {
    	it->bHasER = true;
    	it->scaleER = scaleER;
    }

    // FIXME: This has to be improved: When neurites are copied,
    //        pointers inside still point to our vNeurites array.
    //        If we destroy it, we're in for some pretty EXC_BAD_ACCESSes.
    UG_LOGN("add neurites")
    for (size_t i = 0; i < vNeurites.size(); ++i)
        neuriteProj->add_neurite(vNeurites[i]);
    UG_LOGN("done");

    UG_LOGN("generating neurites")
    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
    	create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, NULL, NULL, &outVerts, &outRads);
    }

    UG_LOGN("generating ER structures")
    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i) {
    	///create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, NULL, NULL, &outVerts, &outRads, &vNeuritesWithin);
    }
    UG_LOGN("done!")


    // at branching points, we have not computed the correct positions yet,
    // so project the complete geometry using the projector
    // TODO: little bit dirty; provide proper method in NeuriteProjector to do this
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
    create_soma(vSomaPoints, g, aaPos, sh);
    sh.set_default_subset_index(0);
    UG_LOGN("Done with soma!");

    // connect soma with neurites
    connect_neurites_with_soma(g, aaPos, outVerts, outRads, 1, sh, fileName);
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

    // at branching points, we have not computed the correct positions yet,
    // so project the complete geometry using the projector
    // TODO: little bit dirty; provide proper method in NeuriteProjector to do this
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



/// top level vertex repositioning function for neurite projection
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

		// TODO: treat the volume case!
	}
}





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

        //GridWriterUGX ugxWriter;
        //ugxWriter.add_grid(g, "defGrid", aPosition);
        //ugxWriter.add_subset_handler(sh, "defSH", 0);
        //ugxWriter.add_projection_handler(projHandler, "defPH", 0);
        //if (!ugxWriter.write_to_file(fileName.c_str()))
        //    UG_THROW("Grid could not be written to file '" << fileName << "'.");
    }
}


} // namespace neuro_collection
} // namespace ug
