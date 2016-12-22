/*
 * hybrid_neuron_communicator.cpp
 *
 *  Created on: 20.12.2016
 *      Author: mbreit
 */

#include "hybrid_neuron_communicator.h"
#include "lib_grid/algorithms/debug_util.h" // ElementDebugInfo

#include <algorithm>	// std::sort
#include <limits>       // std::numeric_limits


namespace ug {
namespace neuro_collection {


template <typename TDomain>
HybridNeuronCommunicator<TDomain>::HybridNeuronCommunicator
(
    SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
    SmartPtr<ApproximationSpace<TDomain> > spApprox1d
)
: m_spCE(SPNULL),
#ifdef UG_PARALLEL
  rcvSize(NULL), rcvFrom(NULL), rcvBuf(NULL),
  sendSize(NULL), sendTo(NULL), sendBuf(NULL),
#endif
  m_spApprox1d(spApprox1d), m_spApprox3d(spApprox3d)
{}


template <typename TDomain>
HybridNeuronCommunicator<TDomain>::~HybridNeuronCommunicator()
{
#ifdef UG_PARALLEL
    if (rcvSize) delete[] rcvSize;
    if (rcvFrom) delete[] rcvFrom;
    if (rcvBuf) delete[] (char*) rcvBuf;
    if (sendSize) delete[] sendSize;
    if (sendTo) delete[] sendTo;
    if (sendBuf)  delete[] (char*) sendBuf;
#endif
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::set_ce_object(ConstSmartPtr<cable_neuron::CableEquation<TDomain> > spCEDisc)
{
    m_spCE = spCEDisc;
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::set_potential_subsets(const std::vector<std::string>& vSubset)
{
    SubsetGroup ssGrp;
    try {ssGrp = SubsetGroup(m_spApprox3d->domain()->subset_handler(), vSubset);}
    UG_CATCH_THROW("Subset group creation failed.");

    for (size_t si = 0; si < ssGrp.size(); si++)
        m_vPotSubset3d.push_back(ssGrp[si]);

    reinit_potential_mappings();
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::reinit_potential_mappings()
{
    typedef typename DoFDistribution::traits<side_t>::const_iterator SideItType;
    typedef typename DoFDistribution::traits<Vertex>::const_iterator VrtItType;
    typedef typename TDomain::position_accessor_type posAccType;

    posAccType aaPos = m_spApprox3d->domain()->position_accessor();
    ConstSmartPtr<DoFDistribution> dd1 = m_spApprox1d->dof_distribution(GridLevel());
    ConstSmartPtr<DoFDistribution> dd3 = m_spApprox3d->dof_distribution(GridLevel());

    // find all elements of the 3d geometry boundary on which the potential is defined
    // store their center coords in the same order
    std::vector<posType> vLocPotElemPos;
    size_t numSs = m_vPotSubset3d.size();
    for (size_t s = 0; s < numSs; ++s)
    {
        int si = m_vPotSubset3d[s];
        SideItType it = dd1->template begin<side_t>(si);
        SideItType it_end = dd1->template end<side_t>(si);

        for (; it != it_end; ++it)
            vLocPotElemPos.push_back(CalculateCenter(*it, aaPos));
    }

    // find all 1d vertices and store them and their positions
    std::vector<Vertex*> vLocVrt;
    std::vector<posType> vLocVrtPos;

    VrtItType it = dd3->template begin<Vertex>();
    VrtItType it_end = dd3->template end<Vertex>();
    for (; it != it_end; ++it)
    {
       vLocVrt.push_back(*it);
       vLocVrtPos.push_back(aaPos[*it]);
    }

#ifdef UG_PARALLEL
    if (pcl::NumProcs() <= 1) goto serial_case;

    {
        // TODO: somehow work with neuron and neurite IDs and parameterization
        // along the neurites, where this is possible
        // this would come with the considerable benefit that only those
        // neurons would need to communicate which have the required neuron(s)

        // communicate search positions to every proc
        std::vector<posType> vGlobPotElemPos;
        std::vector<int> vOffsets;
        pcl::ProcessCommunicator procComm;
        procComm.allgatherv(vGlobPotElemPos, vLocPotElemPos, NULL, &vOffsets);

        // local nearest neighbor search for all elem centers
        std::vector<size_t> vNearest;
        std::vector<typename posType::value_type> vDist;
        int noNeighbors = nearest_neighbor_search(vGlobPotElemPos, vLocVrtPos, vNearest, vDist);

        // NN search return 1 if no neighbors found
        if (noNeighbors)
            vDist.resize(vGlobPotElemPos.size(), std::numeric_limits<typename posType::value_type>::max());

        // find global minimal dists
        std::vector<typename posType::value_type> minDist;
        procComm.allreduce(vDist, minDist, PCL_RO_MIN);

        // now all procs know the global min distances, however, they might not be unique
        std::vector<int> locRank, globRank;
        locRank.resize(minDist.size());
        size_t sz = locRank.size();
        int myRank = pcl::ProcRank();
        for (size_t i = 0; i < sz; ++i)
        {
            if (minDist[i] == vDist[i])
            {
                UG_COND_THROW(minDist[i] == std::numeric_limits<typename posType::value_type>::max(),
                    "No vertex can be mapped for elem at " << vGlobPotElemPos[i] << " on any proc.");

                locRank[i] = myRank;
            }
            else locRank[i] = 0;
        }

        // after this, globRank contains the rank of the proc containing the NN for all query points
        procComm.allreduce(locRank, globRank, PCL_RO_MAX);

        // fill m_vSendInfo and m_vReceiveInfo
        // first, sendInfo
        int receiver = 0;
        std::vector<Vertex*> vVrt;
        for (size_t i = 0; i < sz; ++i)
        {
            // increase receiver ID if necessary
            if (receiver < (int)(vOffsets.size())-1 && receiver >= vOffsets[receiver+1])
            {
                // push old receiver if any NN for it are present
                if (vVrt.size())
                {
                    m_mSendInfo[receiver] = vVrt;
                    vVrt.clear();
                }

                // increase receiver index as long as necessary
                ++receiver;
                while (receiver < (int)(vOffsets.size())-1 && receiver >= vOffsets[receiver+1])
                    ++receiver;
            }

            // this proc is sender if it holds the NN
            if (globRank[i] == myRank)
                vVrt.push_back(vLocVrt[vNearest[i]]);
        }

        // push last receiver if any NN for it are present
        if (vVrt.size())
        {
            m_mSendInfo[receiver] = vVrt;
            vVrt.clear();
        }

        // now receiveInfo
        size_t i = (size_t) vOffsets[(size_t) myRank];
        for (size_t s = 0; s < numSs; ++s)
        {
            int si = m_vPotSubset3d[s];
            SideItType it = dd1->template begin<side_t>(si);
            SideItType it_end = dd1->template end<side_t>(si);
            for (; it != it_end; ++it)
            {
                m_mReceiveInfo[globRank[i]].push_back(*it);
                ++i;
            }
        }


        // at last, prepare communication arrays
        // delete present ones of this is a re-initialization
        if (rcvSize) delete[] rcvSize;
        if (rcvFrom) delete[] rcvFrom;
        if (rcvBuf) delete[] (char*) rcvBuf;
        if (sendSize) delete[] sendSize;
        if (sendTo) delete[] sendTo;
        if (sendBuf)  delete[] (char*) sendBuf;

        // receiving setup
        int numRcv = (int) m_mReceiveInfo.size();
        rcvSize = new int[numRcv];
        rcvFrom = new int[numRcv];
        size_t rcvBytes = 0;

        i = 0;
        typename std::map<int, std::vector<side_t*> >::const_iterator itRec = m_mReceiveInfo.begin();
        typename std::map<int, std::vector<side_t*> >::const_iterator itRec_end = m_mReceiveInfo.end();
        for (; itRec != itRec_end; ++itRec)
        {
            rcvFrom[i] = itRec->first;
            rcvSize[i] = (int) (itRec->second.size() * sizeof(number));
            rcvBytes += rcvSize[i];
            ++i;
        }
        rcvBuf = new char[rcvBytes];

        // sending setup
        int numSend = (int) m_mSendInfo.size();
        sendSize = new int[numSend];
        sendTo = new int[numSend];
        size_t sendBytes = 0;

        i = 0;
        std::map<int, std::vector<Vertex*> >::const_iterator itSend = m_mSendInfo.begin();
        std::map<int, std::vector<Vertex*> >::const_iterator itSend_end = m_mSendInfo.end();
        for (; itSend != itSend_end; ++itSend)
        {
            sendTo[i] = itSend->first;
            sendSize[i] = (int) (itSend->second.size() * sizeof(number));
            sendBytes += sendSize[i];
            ++i;
        }
        sendBuf = new char[sendBytes];

    }

    return;

serial_case:
#endif
        // local nearest neighbor search for all elem centers
        std::vector<size_t> vNearest;
        std::vector<typename posType::value_type> vDist;
        int noNeighbors = nearest_neighbor_search(vLocPotElemPos, vLocVrtPos, vNearest, vDist);
        UG_COND_THROW(noNeighbors, "No 1d vertices are present to map 3d potential elements to.");

        // fill 3d->1d map
        size_t i = 0;
        for (size_t s = 0; s < numSs; ++s)
        {
            int si = m_vPotSubset3d[s];
            SideItType it = dd1->template begin<side_t>(si);
            SideItType it_end = dd1->template end<side_t>(si);
            for (; it != it_end; ++it)
            {
                m_mPotElemToVertex[*it] = vLocVrt[vNearest[i]];
                ++i;
            }
        }
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::coordinate_potential_values()
{
#ifdef UG_PARALLEL
    if (pcl::NumProcs() <= 1) goto serial_case;

    {
        // collect potential values of this proc in send buffer
         int numRcv = (int) m_mReceiveInfo.size();
         int numSend = (int) m_mSendInfo.size();
         number* curVal = (number*) sendBuf;

         std::map<int, std::vector<Vertex*> >::const_iterator itSend = m_mSendInfo.begin();
         std::map<int, std::vector<Vertex*> >::const_iterator itSend_end = m_mSendInfo.end();
         for (; itSend != itSend_end; ++itSend)
         {
             const std::vector<Vertex*>& vVrts = itSend->second;
             size_t sz = vVrts.size();
             for (size_t j = 0; j < sz; ++j)
             {
                 *curVal = m_spCE->vm(vVrts[j]);
                 ++curVal;
             }
         }

        // communicate
        pcl::ProcessCommunicator procComm;
        procComm.distribute_data
        (
            rcvBuf,     // receive buffer (for all data to be received)
            rcvSize,    // sizes of segments in receive buffer
            rcvFrom,    // processes from which data is received
            numRcv,     // number of procs from which data is received
            sendBuf,    // send buffer (for all data to be sent)
            sendSize,   // sizes of segments in send buffer
            sendTo,     // processes to send data to
            numSend     // number of procs to send data to
        );

        // save values in map
        curVal = (number*) rcvBuf;
        typename std::map<int, std::vector<side_t*> >::const_iterator itRec = m_mReceiveInfo.begin();
        typename std::map<int, std::vector<side_t*> >::const_iterator itRec_end = m_mReceiveInfo.end();
        for (; itRec != itRec_end; ++itRec)
        {
            const std::vector<side_t*>& vElems = itRec->second;
            size_t sz = vElems.size();
            for (size_t j = 0; j < sz; ++j)
            {
                m_mElemPot[vElems[j]] = *curVal;
                ++curVal;
            }
        }
    }

serial_case:
#endif

    typename std::map<side_t*, number>::iterator it = m_mElemPot.begin();
    typename std::map<side_t*, number>::iterator it_end = m_mElemPot.end();
    for (; it != it_end; ++it)
        it->second = m_spCE->vm(m_mPotElemToVertex[it->first]);
}


template <typename TDomain>
number HybridNeuronCommunicator<TDomain>::potential(side_t* elem) const
{
	typename std::map<side_t*, number>::const_iterator it = m_mElemPot.find(elem);
	UG_COND_THROW(it == m_mElemPot.end(), "No potential value available for "
        << ElementDebugInfo(*m_spCE->approx_space()->domain()->grid(), elem) << ".");

    return it->second;
}


template <typename TDomain>
int HybridNeuronCommunicator<TDomain>::nearest_neighbor_search
(
	const std::vector<posType>& queryPts,
	const std::vector<posType>& dataPts,
	std::vector<size_t>& vNNout,
	std::vector<typename posType::value_type>& vDistOut
) const
{
	// TODO: speed might benefit from octree or some other spatial tree structure
    size_t qSz = queryPts.size();
    size_t dSz = dataPts.size();

    if (!dSz)
    {
        if (!qSz) return 0;
        return 1; // return that no neighbors could be found at all
    }

    vNNout.resize(qSz);
    vDistOut.resize(qSz);

    for (size_t q = 0; q < qSz; ++q)
    {
        const posType& qp = queryPts[q];
        typename posType::value_type& minDist = vDistOut[q];
        size_t& minPt = vNNout[q];

        minPt = 0;
        minDist = VecDistanceSq(qp, dataPts[0]);
        for (size_t d = 1; d < dSz; ++d)
        {
            number dist = VecDistanceSq(qp, dataPts[d]);
            if (dist < minDist)
            {
                minDist = dist;
                minPt = d;
            }
        }
    }

    return 0;


	/* // probably too complicated, too little gain
	// for some speedup in the NN search, we first sort both point sets
	// with respect to the coordinate axis with biggest variance in the bigger set
	// then we will be able to drop some comparisons

    if (!queryPts.size() || !dataPts.size()) return;

	// compute variation in one pass
	const std::vector<posType>& biggerVec = queryPts.size() > dataPts.size() ? queryPts : dataPts;
	size_t sz = biggerVec.size();

	posType mean(biggerVec[0]);
	posType var(0.0);
	posType delta;

	for (size_t i = 1; i < sz; )
	{
		const posType& x = biggerVec[i];
		VecSubtract(delta, x, mean);
		VecScaleAdd(mean, 1.0, mean, 1.0/(++i), delta);
		for (size_t j = 0; j < posType::Size; ++j)
			var[j] += delta[j] * (x[j] - mean[j]);
	}

	size_t maxInd = 0;
	for (size_t j = 0; j < posType::Size; ++j)
		if (var[j] > var[maxInd])
			maxInd = j;

	// sort
	VecIndexCompare cmpQ(queryPts, maxInd);
	VecIndexCompare cmpD(dataPts, maxInd);

	size_t szQ = queryPts.size();
	size_t szD = dataPts.size();
	std::vector<size_t> qIndices(szQ);
	std::vector<size_t> dIndices(szD);
	for (size_t i = 0; i < szQ; ++i)
		qIndices[i] = i;
	for (size_t i = 0; i < szD; ++i)
		dIndices[i] = i;

	std::sort(qIndices.begin(), qIndices.end(), cmpQ);
	std::sort(dIndices.begin(), dIndices.end(), cmpD);
	*/
}



// explicit template specializations
#ifdef UG_DIM_1
    template class HybridNeuronCommunicator<Domain1d>;
#endif
#ifdef UG_DIM_2
    template class HybridNeuronCommunicator<Domain2d>;
#endif
#ifdef UG_DIM_3
    template class HybridNeuronCommunicator<Domain3d>;
#endif




} // namespace neuro_collection
} // namespace ug
