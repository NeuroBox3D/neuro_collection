/*
 * hybrid_neuron_communicator.cpp
 *
 *  Created on: 20.12.2016
 *      Author: mbreit
 */

#include "hybrid_neuron_communicator.h"
#include "lib_grid/algorithms/debug_util.h" // ElementDebugInfo

#include <algorithm>	// std::sort


namespace ug {
namespace neuro_collection {


template <typename TDomain>
HybridNeuronCommunicator<TDomain>::HybridNeuronCommunicator
(
    SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
    SmartPtr<ApproximationSpace<TDomain> > spApprox1d
)
: m_spCE(SPNULL), m_spApprox1d(spApprox1d), m_spApprox3d(spApprox3d) {}


template <typename TDomain>
HybridNeuronCommunicator<TDomain>::~HybridNeuronCommunicator()
{}


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

    reinit();
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::reinit()
{
    // find all elements of the 3d geometry boundary on which the potential is defined
	// store their coords in the same order
	typedef typename DoFDistribution::traits<side_t>::const_iterator SideItType;
	typedef typename DoFDistribution::traits<Vertex>::const_iterator VrtItType;
	typedef typename TDomain::position_accessor_type posAccType;

	posAccType aaPos = m_spApprox3d->domain()->position_accessor();
	ConstSmartPtr<DoFDistribution> dd1 = m_spApprox1d->dof_distribution(GridLevel());
	ConstSmartPtr<DoFDistribution> dd3 = m_spApprox3d->dof_distribution(GridLevel());

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

	// TODO: somehow work with neuron and neurite IDs and parameterization
	// along the neurites, where this is possible

	// find the procs on which reside the corresponding 1d vertices
	// use nearest neighbor search for that

	// first step: communicate search positions to every proc
	pcl::ProcessCommunicator procComm;

	std::vector<posType> vGlobPotElemPos;
	std::vector<int> vSizes;
	procComm.allgatherv(vGlobPotElemPos, vLocPotElemPos, &vSizes);

	// second step: local nearest neighbor search for all elem centers
	std::vector<posType> vLocVrtPos;

	VrtItType it = dd3->template begin<Vertex>();
	VrtItType it_end = dd3->template end<Vertex>();
	for (; it != it_end; ++it)
		vLocVrtPos.push_back(aaPos[*it]);

	std::vector<size_t> vNearest;
	std::vector<typename posType::value_type> vDist;
	nearest_neighbor_search(vGlobPotElemPos, vLocVrtPos, vNearest, vDist);


    // communicate both to one another to fill m_vSendInfo and m_vReceiveInfo
    // take extra care to preserve the order!
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::communicate_potential_values()
{
    // TODO: only do setup when receiver elems or sender vertices change

	// TODO: only communicate in parallel env

    // receiving setup
    int numRcv = (int) m_vReceiveInfo.size();
    int* rcvSize = new int[numRcv];
    int* rcvFrom = new int[numRcv];
    size_t rcvBytes = 0;
    for (size_t i = 0; i < (size_t) numRcv; ++i)
    {
        rcvFrom[i] = m_vReceiveInfo[i].sender;
        rcvSize[i] = (int) (m_vReceiveInfo[i].vElems.size() * sizeof(number));
        rcvBytes += rcvSize[i];
    }
    void* rcvBuf = new char[rcvBytes];

    // sending setup
    int numSend = (int) m_vSendInfo.size();
    int* sendSize = new int[numSend];
    int* sendTo = new int[numSend];
    size_t sendBytes = 0;
    for (size_t i = 0; i < (size_t) numSend; ++i)
    {
        sendTo[i] = m_vSendInfo[i].receiver;
        sendSize[i] = (int) (m_vSendInfo[i].vVrts.size() * sizeof(number));
        sendBytes += sendSize[i];
    }
    void* sendBuf = new char[sendBytes];

    // collect potential values of this proc in send buffer
    number* curVal = (number*) sendBuf;
    for (size_t i = 0; i < (size_t) numSend; ++i)
    {
        const std::vector<Vertex*>& vVrts = m_vSendInfo[i].vVrts;
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
    for (size_t i = 0; i < (size_t) numRcv; ++i)
    {
        const std::vector<side_t*>& vElems = m_vReceiveInfo[i].vElems;
        size_t sz = vElems.size();
        for (size_t j = 0; j < sz; ++j)
        {
            m_mElemPot[vElems[j]] = *curVal;
            ++curVal;
        }
    }

    // free memory
    delete[] rcvSize;
    delete[] rcvFrom;
    delete[] (char*) rcvBuf;
    delete[] sendSize;
    delete[] sendTo;
    delete[] (char*) sendBuf;
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
void HybridNeuronCommunicator<TDomain>::nearest_neighbor_search
(
	const std::vector<posType>& queryPts,
	const std::vector<posType>& dataPts,
	std::vector<size_t>& vNNout,
	std::vector<typename posType::value_type>& vDistOut
) const
{
	if (!queryPts.size() || !dataPts.size()) return;

	// for some speedup in the NN search, we first sort both point sets
	// with respect to the coordinate axis with biggest variance in the bigger set
	// then we will be able to drop some comparisons
	// TODO: speed might benefit from octree or some other spatial tree structure

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
