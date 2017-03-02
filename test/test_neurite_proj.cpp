/*
 * test_neurite_proj.cpp
 *
 *  Created on: 27.12.2016
 *      Author: mbreit
 */

#include "test_neurite_proj.h"
#include "lib_grid/refinement/projectors/cylinder_volume_projector.h"
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
#include "common/util/string_util.h"  // TrimString
#include "neurite_refMarkAdjuster.h"
#include "lib_disc/function_spaces/error_elem_marking_strategy.h" // GlobalMarking

#include <boost/lexical_cast.hpp>

#include <istream>
#include <sstream>
#include <list>
#include <vector>

namespace ug {
namespace neuro_collection {



void import_swc
(
    const std::string& fileName,
    std::vector<SWCPoint>& vPointsOut
)
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
}



void convert_pointlist_to_neuritelist
(
    const std::vector<SWCPoint>& vPoints,
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

    // TODO: we suppose we only have one cell here
    // find soma
    size_t nPts = vPoints.size();
    size_t i = 0;
    for (; i < nPts; ++i)
    {
        if (vPoints[i].type == SWC_SOMA)
            break;
    }
    UG_COND_THROW(i == nPts, "No soma contained in swc point list.");

    // TODO: process soma somehow; here, they are simply ignored
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
        if (pt.type == SWC_SOMA)
        {
            size_t nConn = pt.conns.size();
            for (size_t i = 0; i < nConn; ++i)
                if (pt.conns[i] != pind)
                    soma_queue.push(std::make_pair(ind, pt.conns[i]));
        }
        else
            rootPts.push_back(std::make_pair(pind, ind));
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
        neuriteOut.refDir = vector3(n,1-n,1-n);
        neuriteOut.vSec.reserve(nVrt-1);

        // this will be 0 for root branches and 1 otherwise
        size_t brInd = neuriteOut.vBR.size();
        std::vector<std::pair<size_t, std::vector<size_t> > >::const_iterator brIt;
        std::vector<std::pair<size_t, std::vector<size_t> > >::const_iterator brIt_end;
        if (bpInfo)
        {
            // fill in missing tend info for first BR
            if (neuriteOut.vBR.size())
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
            neuriteOut.vBR.resize(neuriteOut.vBR.size() + bpInfo->size());
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
            UG_LOGN("  BR " << b << ": " << vBR[b].tstart << ".." << vBR[b].tend);
    }
    */
}



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
    std::vector<Edge*>* connectingEdges = NULL
)
{
    const NeuriteProjector::Neurite& neurite = vNeurites[nid];
    const std::vector<vector3>& pos = vPos[nid];
    const std::vector<number>& r = vR[nid];

    size_t nSec = neurite.vSec.size();

    const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
    std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit = vBR.begin();
    std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end = vBR.end();

    std::vector<Vertex*> vVrt;
    std::vector<Edge*> vEdge;
    vVrt.resize(4);
    vEdge.resize(4);

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
        // create first layer of vertices/edges
        for (size_t i = 0; i < 4; ++i)
        {
            Vertex* v = *g.create<RegularVertex>();
            vVrt[i] = v;
            number angle = 0.5*PI*i;
            VecScaleAdd(aaPos[v], 1.0, pos[0], r[0]*cos(angle), projRefDir, r[0]*sin(angle), thirdDir);

            aaSurfParams[v].neuriteID = nid;
            aaSurfParams[v].axial = 0.0;
            aaSurfParams[v].angular = angle;
        }
        for (size_t i = 0; i < 4; ++i)
            vEdge[i] = *g.create<RegularEdge>(EdgeDescriptor(vVrt[i], vVrt[(i+1)%4]));
    }


    vector3 lastPos = pos[0];
    for (size_t i = 0; i < nSec; ++i)
    {
        const NeuriteProjector::Section& sec = neurite.vSec[i];

        // get velocity
        vel[0] = -sec.splineParamsX[2];
        vel[1] = -sec.splineParamsY[2];
        vel[2] = -sec.splineParamsZ[2];
        number velNorm = sqrt(VecNormSquared(vel));
        vel /= velNorm;

        // calculate reference dir projected to normal plane of velocity
        number fac = VecProd(neurite.refDir, vel);
        VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
        VecNormalize(projRefDir, projRefDir);

        VecCross(thirdDir, vel, projRefDir);

        // possible branching region:
        if (brit != brit_end && brit->tstart < sec.endParam)
        {
            // around the new position, create a hexahedron
            // extrude from last pos to start of hexahedron
            number radius = r[i+1];
            vector3 extrudeDir;
            VecScaleAdd(extrudeDir, 1.0, pos[i+1], -1.0, lastPos, -sqrt(0.5)*radius, vel);
            Extrude(g, &vVrt, &vEdge, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

            // set new positions and param attachments; also ensure correct face orientation
            for (size_t j = 0; j < 4; ++j)
            {
                number angle = 0.5*PI*j + angleOffset;
                if (angle > 2*PI) angle -= 2*PI;
                Vertex* v = vVrt[j];
                vector3 radialVec;
                VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
                VecScaleAdd(aaPos[v], 1.0, pos[i+1], -sqrt(0.5)*radius, vel, 1.0, radialVec);

                aaSurfParams[v].neuriteID = nid;
                aaSurfParams[v].axial = sec.endParam - sqrt(0.5)*radius/velNorm;
                aaSurfParams[v].angular = angle;

                Grid::traits<Face>::secure_container faceCont;
                g.associated_elements(faceCont, vEdge[j]);  // faceCont must contain exactly one face
                vector3 normal;
                CalculateNormal(normal, faceCont[0], aaPos);
                if (VecProd(normal, radialVec) < 0)
                    g.flip_orientation(faceCont[0]);
            }

            // extrude hex
            VecScale(extrudeDir, vel, sqrt(2.0)*radius);
            Selector sel(g);
            sel.enable_autoselection(true);
            Extrude(g, &vVrt, &vEdge, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

            // set new positions and param attachments; also ensure correct face orientation
            for (size_t j = 0; j < 4; ++j)
            {
                number angle = 0.5*PI*j + angleOffset;
                if (angle > 2*PI) angle -= 2*PI;
                Vertex* v = vVrt[j];
                vector3 radialVec;
                VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
                VecScaleAdd(aaPos[v], 1.0, pos[i+1], sqrt(0.5)*radius, vel, 1.0, radialVec);

                aaSurfParams[v].neuriteID = nid;
                aaSurfParams[v].axial = sec.endParam + sqrt(0.5)*radius/velNorm;
                aaSurfParams[v].angular = angle;

                Grid::traits<Face>::secure_container faceCont;
                g.associated_elements(faceCont, vEdge[j]);  // faceCont must contain exactly one face
                vector3 normal;
                CalculateNormal(normal, faceCont[0], aaPos);
                if (VecProd(normal, radialVec) < 0)
                    g.flip_orientation(faceCont[0]);
            }

            VecScaleAdd(lastPos, 1.0, pos[i+1], sqrt(0.5)*radius, vel);

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

            create_neurite(vNeurites, vPos, vR, child_nid, g, aaPos, aaSurfParams, &vrts, &edges);

            ++brit;
        }
        else
        {
            // extrude
            vector3 extrudeDir;
            VecSubtract(extrudeDir, pos[i+1], lastPos);
            Extrude(g, &vVrt, &vEdge, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);
            lastPos = pos[i+1];        number radius = r[i+1];

            // set new positions and param attachments; also ensure correct face orientation
            for (size_t j = 0; j < 4; ++j)
            {
                number angle = 0.5*PI*j + angleOffset;
                if (angle > 2*PI) angle -= 2*PI;
                Vertex* v = vVrt[j];
                vector3 radialVec;
                VecScaleAdd(radialVec, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
                VecAdd(aaPos[v], pos[i+1], radialVec);

                aaSurfParams[v].neuriteID = nid;
                aaSurfParams[v].axial = sec.endParam;
                aaSurfParams[v].angular = angle;

                Grid::traits<Face>::secure_container faceCont;
                g.associated_elements(faceCont, vEdge[j]);  // faceCont must contain exactly one face
                vector3 normal;
                CalculateNormal(normal, faceCont[0], aaPos);
                if (VecProd(normal, radialVec) < 0)
                    g.flip_orientation(faceCont[0]);
            }
        }
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


void test_import_swc(const std::string& fileName)
{
    // read in file to intermediate structure
    std::vector<SWCPoint> vPoints;
    import_swc(fileName, vPoints);

    // convert intermediate structure to neurite data
    std::vector<std::vector<vector3> > vPos;
    std::vector<std::vector<number> > vRad;
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > > vBPInfo;
    std::vector<size_t> vRootNeuriteIndsOut;
    convert_pointlist_to_neuritelist(vPoints, vPos, vRad, vBPInfo, vRootNeuriteIndsOut);

    std::cout << "BPInfo:" << std::endl;
    for (size_t i = 0; i < vBPInfo.size(); ++i)
    {
        std::cout << "nid " << i << ":  " << "#vrts = " << vPos[i].size() << std::endl;
        for (size_t j = 0; j < vBPInfo[i].size(); ++j)
        {
            std::cout << "  " << vBPInfo[i][j].first << ": ";
            for (size_t k = 0; k < vBPInfo[i][j].second.size(); ++k)
            {
                std::cout << vBPInfo[i][j].second[k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    // create spline data
    std::vector<NeuriteProjector::Neurite> vNeurites(2);
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


    ProjectionHandler projHandler(&sh);
    SmartPtr<IGeometry<3> > geom3d = MakeGeometry3d(g, aPosition);
    projHandler.set_geometry(geom3d);

    SmartPtr<NeuriteProjector> neuriteProj(new NeuriteProjector(geom3d));
    projHandler.set_projector(0, neuriteProj);

    for (size_t i = 0; i < vNeurites.size(); ++i)
        neuriteProj->add_neurite(vNeurites[i]);

    for (size_t i = 0; i < vRootNeuriteIndsOut.size(); ++i)
        create_neurite(vNeurites, vPos, vRad, vRootNeuriteIndsOut[i], g, aaPos, aaSurfParams, NULL, NULL);

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

    // refinement
    AssignSubsetColors(sh);
    sh.set_subset_name("surf", 0);

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

    GlobalMultiGridRefiner ref(*dom.grid(), dom.refinement_projector());
    for (size_t i = 0; i < 4; ++i)
    {
        ref.refine();

        std::ostringstream oss;
        oss << "_refined_" << i+1 << ".ugx";
        std::string curFileName = outFileName.substr(0, outFileName.size()-4) + oss.str();
        number offset = 1;
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
