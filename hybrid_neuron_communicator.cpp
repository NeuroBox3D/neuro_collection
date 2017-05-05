/*
 * hybrid_neuron_communicator.cpp
 *
 *  Created on: 20.12.2016
 *      Author: mbreit, lreinhardt
 */

#include "hybrid_neuron_communicator.h"
#include "lib_grid/algorithms/debug_util.h" // ElementDebugInfo
#include "lib_algebra/cpu_algebra_types.h" // CPUAlgebra
#include "../cable_neuron/util/functions.h"	// neuron_identification
#include "common/util/vector_util.h" // GetDataPtr

#include <algorithm>	// std::sort
#include <limits>       // std::numeric_limits
#include <vector>

namespace ug {
namespace neuro_collection {


template <typename TDomain>
HybridNeuronCommunicator<TDomain>::HybridNeuronCommunicator
(
    SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
    SmartPtr<ApproximationSpace<TDomain> > spApprox1d
)
: m_spSynHandler(SPNULL),
  m_mSynapse3dVertex(std::map<synapse_id, Vertex*>()),
#ifdef UG_PARALLEL
  rcvSize(NULL), rcvFrom(NULL), rcvBuf(NULL),
  sendSize(NULL), sendTo(NULL), sendBuf(NULL),
#endif
  m_spApprox1d(spApprox1d), m_spApprox3d(spApprox3d),
  m_potFctInd(0),
  m_spGrid1d(m_spApprox1d->domain()->grid()), m_spGrid3d(m_spApprox3d->domain()->grid()),
  m_spMGSSH3d(m_spApprox3d->domain()->subset_handler()),
  m_scale_factor_from_3d_to_1d(1.0),
  m_aaPos1d(m_spApprox1d->domain()->position_accessor()),
  m_aaPos3d(m_spApprox3d->domain()->position_accessor()),
  m_aNID(GlobalAttachments::attachment<ANeuronID>("neuronID")),
  m_vNid(1,0),
  m_bInited(false)
{
	if (!m_spGrid1d->has_vertex_attachment(m_aNID))
		m_spGrid1d->attach_to_vertices(m_aNID);
	m_aaNID = Grid::VertexAttachmentAccessor<ANeuronID>(*m_spGrid1d, m_aNID);

	// calculate identifiers for each neuron
	cable_neuron::neuron_identification(*m_spGrid1d);
}


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
void HybridNeuronCommunicator<TDomain>::set_synapse_handler
(SmartPtr<cable_neuron::synapse_handler::SynapseHandler<TDomain> > spSH)
{
    m_spSynHandler = spSH;
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::set_potential_subsets(const std::vector<std::string>& vSubset)
{
    SubsetGroup ssGrp;
    try {ssGrp = SubsetGroup(m_spApprox3d->domain()->subset_handler(), vSubset);}
    UG_CATCH_THROW("Subset group creation failed.");

    for (size_t si = 0; si < ssGrp.size(); si++)
        m_vPotSubset3d.push_back(ssGrp[si]);
}

template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::set_current_subsets(const std::vector<std::string>& vSubset)
{
    SubsetGroup ssGrp;
    try {ssGrp = SubsetGroup(m_spApprox3d->domain()->subset_handler(), vSubset);}
    UG_CATCH_THROW("Subset group creation failed.");

    for (size_t si = 0; si < ssGrp.size(); si++)
    	m_vCurrentSubset3d.push_back(ssGrp[si]);
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::
set_solution_and_potential_index(ConstSmartPtr<GridFunction<TDomain, algebra_t> > u, size_t fctInd)
{
    m_spU = u;
    m_potFctInd = fctInd;
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::set_coordinate_scale_factor_3d_to_1d(number scale)
{
	m_scale_factor_from_3d_to_1d = scale;
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::set_neuron_ids(const std::vector<uint>& vNid)
{
	m_vNid = vNid;
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::reinit()
{
	if (!m_bInited)
	{
		// TODO: This is awkward.
		//       HNC is used for two purposes which require different members set.
		//       Think about better separation.
		if (!m_vPotSubset3d.empty()) reinit_potential_mapping();
		if (!m_vCurrentSubset3d.empty()) reinit_synapse_mapping();
		m_bInited = true;
	}
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::reinit_potential_mapping()
{
    typedef typename DoFDistribution::traits<side_t>::const_iterator SideItType;
    typedef typename DoFDistribution::traits<Vertex>::const_iterator VrtItType;
    typedef typename TDomain::position_accessor_type posAccType;

    posAccType aaPos1 = m_spApprox1d->domain()->position_accessor();
    posAccType aaPos3 = m_spApprox3d->domain()->position_accessor();
    ConstSmartPtr<DoFDistribution> dd1 = m_spApprox1d->dof_distribution(GridLevel());
    ConstSmartPtr<DoFDistribution> dd3 = m_spApprox3d->dof_distribution(GridLevel());

    // find all elements of the 3d geometry boundary on which the potential is defined
    // store their center coords in the same order
    std::vector<posType> vLocPotElemPos;
    std::vector<side_t*> vLocPotElems;
    size_t numSs = m_vPotSubset3d.size();
    for (size_t s = 0; s < numSs; ++s)
    {
        int si = m_vPotSubset3d[s];
        SideItType it = dd3->template begin<side_t>(si);
        SideItType it_end = dd3->template end<side_t>(si);

        for (; it != it_end; ++it)
        {
        	vLocPotElems.push_back(*it);
            vLocPotElemPos.push_back(CalculateCenter(*it, aaPos3).operator*=(m_scale_factor_from_3d_to_1d));
        }
    }

    // find all 1d vertices for neuron IDs of interest and store them and their positions
    // TODO: think about finding edges instead -> this would allow linear interpolation
    std::vector<Vertex*> vLocVrt;
    std::vector<posType> vLocVrtPos;

    size_t nNid = m_vNid.size();
    VrtItType it = dd1->template begin<Vertex>();
    VrtItType it_end = dd1->template end<Vertex>();
    for (; it != it_end; ++it)
    {
    	for (size_t n = 0; n < nNid; ++n)
    	{
    		if (m_aaNID[*it] == m_vNid[n])
    		{
    			vLocVrt.push_back(*it);
    			vLocVrtPos.push_back(aaPos1[*it]);
				break;
    		}
    	}
    }

#ifdef UG_PARALLEL
    if (pcl::NumProcs() <= 1) goto serial_case;

    {
        // TODO: somehow work with parameterization along the neurites
    	//       or do something like Voronoi tesselation of the 3d geom w.r.t. 1d vertices

    	// communicate 1d positions to every proc (as there are probably much fewer of them)
        std::vector<posType> vGlobVrtPos;
        std::vector<int> vOffsets;
        pcl::ProcessCommunicator procComm;
        procComm.allgatherv(vGlobVrtPos, vLocVrtPos, NULL, &vOffsets);

        // local nearest neighbor search for all elem centers
        std::vector<size_t> vNearest;
        std::vector<typename posType::value_type> vDist;
        int noNeighbors = nearest_neighbor_search(vLocPotElemPos, vGlobVrtPos, vNearest, vDist);

        // NN search returns 1 if no neighbors found
        UG_COND_THROW(noNeighbors, "No 1d neighbors could be found for any potential element.\n"
        	"This means there are no 1d vertices for the given neuron IDs of interest.");

        // construct (using vOffset) which proc holds the nearest 1d vertex and under which index
        // at the same time, fill recvInfo
        size_t nPE = vLocPotElemPos.size();
        size_t nProcs = vOffsets.size();
        std::map<size_t, std::vector<int> > mIndex;
        m_mReceiveInfo.clear();
        for (size_t i = 0; i < nPE; ++i)
        {
        	size_t pos = vNearest[i];

        	// perform a binary search for the largest entry in offset
        	// that is lower than or equal to pos
        	size_t low = 0;
        	size_t high = nProcs-1;
        	while (low != high)
        	{
        	    size_t mid = (low + high) / 2 + (low + high) % 2;
        	    if ((size_t) vOffsets[mid] > pos)
        	    	high = mid-1;
        	    else
        	    	low = mid;
        	}

        	mIndex[low].push_back((int) pos - vOffsets[low]);

        	// add to recvInfo
        	m_mReceiveInfo[(int) low].push_back(vLocPotElems[i]);
        }

        // fill m_vSendInfo (we need to communicate for that)
        std::vector<int> recvBuffer;
		std::vector<int> recvSizes;
		std::vector<int> senderProcs;
		int numSenderProcs = 0;
		std::vector<int> sendBuffer;
		std::vector<int> sendSizes;
		std::vector<int> recverProcs;
		int numRecverProcs = 0;

		// step 1: who has how much for whom?
        std::vector<int> vNumTo(nProcs, 0);
        for (size_t p = 0; p < nProcs; ++p)
        {
        	std::map<size_t, std::vector<int> >::const_iterator it = mIndex.find(p);
        	if (it != mIndex.end())
        	{
        		int sz = it->second.size();
        		vNumTo[p] = sz;

        		++numRecverProcs;
        		recverProcs.push_back(p);
        		sendSizes.push_back(sz * sizeof(int));
        		for (size_t i = 0; i < (size_t) sz; ++i)
        			sendBuffer.push_back(it->second[i]);
        	}
        }
        std::vector<int> vNumFrom(nProcs);
        procComm.alltoall(&vNumTo[0], 1, PCL_DT_INT, &vNumFrom[0], 1, PCL_DT_INT);

        // step 2: exchange information on who has which minDist vertices of whom
		size_t nRcv = 0;
		for (size_t p = 0; p < nProcs; ++p)
		{
			if (vNumFrom[p])
			{
				++numSenderProcs;
				senderProcs.push_back(p);
				recvSizes.push_back(vNumFrom[p] * sizeof(int));
				nRcv += vNumFrom[p];
			}
		}
		recvBuffer.resize(nRcv);

		procComm.distribute_data
		(
			GetDataPtr(recvBuffer),     // receive buffer (for all data to be received)
			GetDataPtr(recvSizes),      // sizes of segments in receive buffer (in bytes)
			GetDataPtr(senderProcs),    // processes from which data is received
			numSenderProcs,             // number of procs from which data is received
			GetDataPtr(sendBuffer),     // send buffer (for all data to be sent)
			GetDataPtr(sendSizes),      // sizes of segments in send buffer (in bytes)
			GetDataPtr(recverProcs),    // processes to send data to
			numRecverProcs              // number of procs to send data to
		);

		// step 3: fill the actual senderInfo
		sendBuffer.clear();
		sendSizes.clear();
		recverProcs.clear();

		m_mSendInfo.clear();
		size_t offset = 0;
		for (int p = 0; p < numSenderProcs; ++p)
		{
			size_t sz = recvSizes[p] / sizeof(int);
			std::vector<Vertex*>& senderVrts = m_mSendInfo[senderProcs[p]];
			senderVrts.resize(sz);
			for (size_t i = 0; i < sz; ++i)
			{
				size_t locVrtInd = recvBuffer[offset + i];
				senderVrts[i] = vLocVrt[locVrtInd];
			}

			offset += sz;
		}

		recvBuffer.clear();
		recvSizes.clear();
		senderProcs.clear();



        // at last, prepare communication arrays
        // delete present ones if this is a re-initialization
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

        size_t i = 0;
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
		SideItType it = dd3->template begin<side_t>(si);
		SideItType it_end = dd3->template end<side_t>(si);
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
	// perform 3d elem -> 1d vertex and 1d synapse -> 3d vertex mappings
	// if not yet done
	reinit();

    ConstSmartPtr<DoFDistribution> dd1 = m_spApprox1d->dof_distribution(GridLevel(), false);

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
                 // get DoFIndex for vertex
                 std::vector<DoFIndex> vIndex;
                 dd1->inner_dof_indices(vVrts[j], m_potFctInd, vIndex, false);

                 UG_COND_THROW(!vIndex.size(), "Potential function (index: "
                     << m_potFctInd << ") is not defined for "
                     << ElementDebugInfo(*this->m_spApprox1d->domain()->grid(), vVrts[j]) << ".")

                 UG_ASSERT(vIndex.size() == 1, "Apparently, shape functions different from P1 "
                     "are used in the 1d approximation space.\nThis is not supported.");

                 // save value in buffer
                 *curVal = DoFRef(*m_spU, vIndex[0]);
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

    return;

serial_case:
#endif

    typename std::map<side_t*, Vertex*>::iterator it = m_mPotElemToVertex.begin();
    typename std::map<side_t*, Vertex*>::iterator it_end = m_mPotElemToVertex.end();
    std::vector<DoFIndex> vIndex;
    for (; it != it_end; ++it)
    {
        Vertex* vrt = it->second;

        // get DoFIndex for vertex
        dd1->inner_dof_indices(vrt, m_potFctInd, vIndex, true);

        UG_COND_THROW(!vIndex.size(), "Potential function (index: "
            << m_potFctInd << ") is not defined for "
            << ElementDebugInfo(*this->m_spApprox1d->domain()->grid(), vrt) << ".")

        UG_ASSERT(vIndex.size() == 1, "Apparently, shape functions different from P1 "
            "are used in the 1d approximation space.\nThis is not supported.");

        // save value in buffer
        m_mElemPot[it->first] = DoFRef(*m_spU, vIndex[0]);
    }
}


template <typename TDomain>
number HybridNeuronCommunicator<TDomain>::potential(side_t* elem) const
{
	typename std::map<side_t*, number>::const_iterator it = m_mElemPot.find(elem);
	UG_COND_THROW(it == m_mElemPot.end(), "No potential value available for "
        << ElementDebugInfo(*m_spApprox3d->domain()->grid(), elem) << ".");

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

    vNNout.clear();
    vDistOut.clear();
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
}



template <typename TDomain>
uint HybridNeuronCommunicator<TDomain>::get_postsyn_neuron_id(synapse_id id)
{
	Edge* e;
	try {e = m_spSynHandler->postsyn_edge(id);}
	UG_CATCH_THROW("Could not determine neuron ID for synapse " << id << ".");

	return m_aaNID[e->vertex(0)];
}

template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::reinit_synapse_mapping()
{
	std::vector<MathVector<dim> > syn_coords_local;
	std::vector<synapse_id> syn_ids_local;

	std::vector<Vertex*> v3dVertices;
	std::vector<MathVector<dim> > v3dVertexPos;

	// gather plasma membrane surface vertices:
    typedef typename DoFDistribution::traits<Vertex>::const_iterator VrtItType;
	ConstSmartPtr<DoFDistribution> dd3 = m_spApprox3d->dof_distribution(GridLevel());
	size_t numSs = m_vCurrentSubset3d.size();
	for (size_t s = 0; s < numSs; ++s)
	{
		int si = m_vCurrentSubset3d[s];
		VrtItType it = dd3->template begin<Vertex>(si);
		VrtItType it_end = dd3->template end<Vertex>(si);

		for (; it != it_end; ++it)
		{
			v3dVertices.push_back(*it);
			v3dVertexPos.push_back(m_aaPos3d[*it]);
			v3dVertexPos.back().operator*=(m_scale_factor_from_3d_to_1d);
		}
	}

	// gather local synapse coords, that are interesting
	typedef cable_neuron::synapse_handler::IPostSynapse PostSynapse;
	typedef cable_neuron::synapse_handler::SynapseIter<PostSynapse> PostSynIter;
	PostSynIter synIt = m_spSynHandler->template begin<PostSynapse>();
	PostSynIter synItEnd = m_spSynHandler->template end<PostSynapse>();
	size_t nNeuron = m_vNid.size();
	for (; synIt != synItEnd; ++synIt)
	{
		PostSynapse* syn = *synIt;
		uint nid = get_postsyn_neuron_id(syn->id());
		for (size_t n = 0; n < nNeuron; ++n)
		{
			if (nid == m_vNid[n])
			{
				MathVector<dim> coords;
				get_postsyn_coordinates(syn->id(), coords);
				syn_coords_local.push_back(coords);
				syn_ids_local.push_back(syn->id());
				break;
			}
		}
	}

#ifdef UG_PARALLEL
	size_t nProcs = pcl::NumProcs();
	if (nProcs > 1)
	{
		pcl::ProcessCommunicator com;
		std::vector<MathVector<dim> > syn_coords_global; //synapses coords of neuron with given id
		std::vector<synapse_id> syn_ids_global;
		std::vector<int> vSizes, vOffsets;

		// communicate all interesting synapse coords to all procs
		com.allgatherv(syn_coords_global, syn_coords_local, &vSizes, &vOffsets);

		// compute min distances of all local 3d vertices to the global 1d synapse positions
        std::vector<size_t> vNearest;
    	std::vector<number> vDistances;
		int failure = nearest_neighbor_search(syn_coords_global, v3dVertexPos, vNearest, vDistances);

		if (failure) // this happens if no 3d vertices are present on this proc
			vDistances.resize(syn_coords_global.size(), std::numeric_limits<typename posType::value_type>::max());


		// communicate local min dists back to original 1d-synapse holders
		number* globDist = new number[nProcs * syn_coords_local.size()];
		for (size_t p = 0; p < nProcs; ++p)
		{
			if (!vSizes[p]) continue;
			com.gather((void*) &vDistances[vOffsets[p]], vSizes[p], PCL_DT_DOUBLE, (void*) globDist, vSizes[p], PCL_DT_DOUBLE, (int) p);
		}

		// compute minimizing proc for each 1d synapse
		size_t nSyn1d = syn_coords_local.size();
		std::vector<size_t> vMinProc(nSyn1d, 0);
		for (size_t s = 0; s < nSyn1d; ++s)
		{
			size_t& minProc = vMinProc[s];
			number minDist = globDist[s];
			for (size_t p = 1; p < nProcs; ++p)
			{
				if (globDist[s + p*nSyn1d] < minDist)
				{
					minDist = globDist[s + p*nSyn1d];
					minProc = p;
				}
			}
			UG_COND_THROW(minDist == std::numeric_limits<typename posType::value_type>::max(),
				"No 3d vertex in defined plasma membrane subset present on any proc.")
		}

		delete[] globDist;

		// inform minProcs about their representing a synapse
		// and also provide them with the synapse ID (for ease of use)

		// step 1: fill send buffers
		std::vector<BinaryBuffer> sendBuffers;
		std::vector<int> receiverProcs;
		std::map<size_t, size_t> mapProcToBufferIndex;
		for (size_t s = 0; s < nSyn1d; ++s)
		{
			size_t proc = vMinProc[s];
			std::map<size_t, size_t>::const_iterator mapIt = mapProcToBufferIndex.find(proc);
			if (mapIt == mapProcToBufferIndex.end())
			{
				mapProcToBufferIndex[proc] = sendBuffers.size();
				receiverProcs.push_back(proc);
				sendBuffers.resize(sendBuffers.size() + 1);
				BinaryBuffer& buf = sendBuffers.back();
				buf.write((char*) &s, sizeof(size_t));
				buf.write((char*) &syn_ids_local[s], sizeof(synapse_id));
			}
			else
			{
				BinaryBuffer& buf = sendBuffers[mapIt->second];
				buf.write((char*) &s, sizeof(size_t));
				buf.write((char*) &syn_ids_local[s], sizeof(synapse_id));
			}
		}

		// step 2: prepare (communicate who will send to whom)
		std::vector<int> vSendTo(nProcs, 0);
		int numRecProcs = receiverProcs.size();
		for (int r = 0; r < numRecProcs; ++r)
			vSendTo[receiverProcs[r]] = 1;	// set all receivers 1; the others stay 0

		std::vector<int> vNumFrom(nProcs);
		com.alltoall(&vSendTo[0], 1, PCL_DT_INT, &vNumFrom[0], 1, PCL_DT_INT);

		std::vector<int> senderProcs;
		for (size_t p = 0; p < nProcs; ++p)
			if (vNumFrom[p])
				senderProcs.push_back(p);

		// step 3: distribute data
		int numSendProcs = senderProcs.size();
		std::vector<BinaryBuffer> recBuffers(numSendProcs);

		com.distribute_data
		(
			GetDataPtr(recBuffers),
			GetDataPtr(senderProcs),
			numSendProcs,
			GetDataPtr(sendBuffers),
			GetDataPtr(receiverProcs),
			numRecProcs
		);

		// step 4: use received data to construct 1d-post-synapse -> 3d vertex map
		for (size_t s = 0; s < (size_t) numSendProcs; ++s)
		{
			BinaryBuffer& recBuf = recBuffers[s];
			int sendProc = senderProcs[s];
			int offset = vOffsets[sendProc];

			while (!recBuf.eof())
			{
				size_t ind;
				synapse_id sid;
				recBuf.read((char*) &ind, sizeof(size_t));
				recBuf.read((char*) &sid, sizeof(synapse_id));

				UG_ASSERT(m_mSynapse3dVertex.find(sid) == m_mSynapse3dVertex.end(),
					"Synapse ID " << sid << " already mapped to a vertex.");

				m_mSynapse3dVertex[sid] = v3dVertices[vNearest[offset + ind]];
			}
		}
	}
	else
	{
#endif
		std::vector<size_t> vNearest;
		std::vector<typename posType::value_type> vDistances;
		int failure = nearest_neighbor_search(syn_coords_local, v3dVertexPos, vNearest, vDistances);
		UG_COND_THROW(failure, "No 3d vertices in the defined plasma membrane subset present to map 1d synapses to.");

		for (size_t i = 0; i < syn_coords_local.size(); ++i)
			m_mSynapse3dVertex[syn_ids_local[i]] = v3dVertices[vNearest[i]];
#ifdef UG_PARALLEL
	}
#endif
}

template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::gather_synaptic_currents
(
	std::vector<Vertex*>& vActSynOut,
	std::vector<number>& vSynCurrOut,
	number time
)
{
	// reinit mappings if necessary
	reinit();

	vActSynOut.clear();
	vSynCurrOut.clear();

	// get locally active synapses
	std::vector<synapse_id> vLocActSyn;
	std::vector<number> vLocSynCurr;
	m_spSynHandler->active_postsynapses_and_currents(vLocActSyn, vLocSynCurr, m_vNid, m_aaNID, time);

	// todo: in the parallel case, do not gather all active synapse ids and currents
	//       on every proc; instead, send active synapses only to the procs that have
	//       the corresponding 3d vertex; use information computed in reinit_synapse_mapping()
#ifdef UG_PARALLEL
	if (pcl::NumProcs() == 1)
		goto serial_case;

	{
		std::vector<synapse_id> vGlobActSyn;
		std::vector<number> vGlobSynCurr;

		pcl::ProcessCommunicator com;
		com.allgatherv(vGlobActSyn, vLocActSyn, NULL, NULL);
		com.allgatherv(vGlobSynCurr, vLocSynCurr, NULL, NULL);

		for (size_t i=0; i<vGlobActSyn.size(); ++i)
		{
			synapse_id sid = vGlobActSyn[i];
			std::map<synapse_id, Vertex*>::const_iterator it;
			if ((it = m_mSynapse3dVertex.find(sid)) != m_mSynapse3dVertex.end())
			{
				vActSynOut.push_back(it->second);
				vSynCurrOut.push_back(vGlobSynCurr[i]);
			}
		}
	}

	return;
#endif

serial_case:
	for (size_t i=0; i<vLocActSyn.size(); ++i)
	{
		synapse_id sid = vLocActSyn[i];
		std::map<synapse_id, Vertex*>::const_iterator it;
		if ((it = m_mSynapse3dVertex.find(sid)) != m_mSynapse3dVertex.end())
		{
			vActSynOut.push_back(it->second);
			vSynCurrOut.push_back(vLocSynCurr[i]);
		}
	}
}



template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::get_postsyn_coordinates(synapse_id id, MathVector<dim>& vCoords)
{
	Edge* e;
	try {e = m_spSynHandler->postsyn_edge(id);}
	UG_CATCH_THROW("Could not determine edge for synapse " << id << ".");

	Vertex* v0 = (*e)[0];
	Vertex* v1 = (*e)[1];
	number localcoord = m_spSynHandler->post_synapse(id)->location();
	VecScaleAdd(vCoords, 1.0 - localcoord, m_aaPos1d[v0], localcoord, m_aaPos1d[v1]);
}


template <typename TDomain, typename TAlgebra>
HybridSynapseCurrentAssembler<TDomain, TAlgebra>::
HybridSynapseCurrentAssembler
(
	SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
	SmartPtr<ApproximationSpace<TDomain> > spApprox1d,
	SmartPtr<cable_neuron::synapse_handler::SynapseHandler<TDomain> > spSH,
	const std::vector<std::string>& PlasmaMembraneSubsetName,
	const std::string& fct,
	const std::string& fct_ip3
)
: m_fctInd(0), m_fctInd_ip3(0), m_ip3_set(true), m_F(96485.309), m_valency(2), m_current_percentage(0.1), m_spHNC(new hnc_type(spApprox3d, spApprox1d)),
  m_scaling_3d_to_1d_amount_of_substance(1e-15), m_scaling_3d_to_1d_electric_charge(1.0), m_scaling_3d_to_1d_coordinates(1e-6), m_scaling_3d_to_1d_ip3(1e-15),
  m_j_ip3_max(6e-19), m_j_ip3_decayRate(1.188), m_j_ip3_duration(3.0 / m_j_ip3_decayRate)
{
	// get function index of whatever it is that the current carries (in our case: calcium)
	FunctionGroup fctGrp(spApprox3d->function_pattern());
	try {fctGrp.add(fct);}
	UG_CATCH_THROW("Function " << fct << " could not be identified in given approximation space.");
	m_fctInd = fctGrp.unique_id(0);

	//get function index if use of ip3 is enabled
	try {fctGrp.add(fct_ip3);}
	UG_CATCH_THROW("Function " << fct_ip3 << " could not be identified in given approximation space.");
	m_fctInd_ip3 = fctGrp.unique_id(0);

	// init HNC
	m_spHNC->set_synapse_handler(spSH);
	m_spHNC->set_coordinate_scale_factor_3d_to_1d(m_scaling_3d_to_1d_coordinates);
	m_spHNC->set_current_subsets(PlasmaMembraneSubsetName);
}


template <typename TDomain, typename TAlgebra>
HybridSynapseCurrentAssembler<TDomain, TAlgebra>::
HybridSynapseCurrentAssembler
(
	SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
	SmartPtr<ApproximationSpace<TDomain> > spApprox1d,
	SmartPtr<cable_neuron::synapse_handler::SynapseHandler<TDomain> > spSH,
	const std::vector<std::string>& PlasmaMembraneSubsetName,
	const std::string& fct
)
: m_fctInd(0), m_fctInd_ip3(0), m_ip3_set(false), m_F(96485.309), m_valency(2), m_current_percentage(0.1), m_spHNC(new hnc_type(spApprox3d, spApprox1d)),
  m_scaling_3d_to_1d_amount_of_substance(1e-15), m_scaling_3d_to_1d_electric_charge(1.0), m_scaling_3d_to_1d_coordinates(1e-6), m_scaling_3d_to_1d_ip3(1e-15),
  m_j_ip3_max(6e-19), m_j_ip3_decayRate(1.188), m_j_ip3_duration(3.0 / m_j_ip3_decayRate)
{
	// get function index of whatever it is that the current carries (in our case: calcium)
	FunctionGroup fctGrp(spApprox3d->function_pattern());
	try {fctGrp.add(fct);}
	UG_CATCH_THROW("Function " << fct << " could not be identified in given approximation space.");
	m_fctInd = fctGrp.unique_id(0);

	// init HNC
	m_spHNC->set_synapse_handler(spSH);
	m_spHNC->set_coordinate_scale_factor_3d_to_1d(m_scaling_3d_to_1d_coordinates);
	m_spHNC->set_current_subsets(PlasmaMembraneSubsetName);
}


template <typename TDomain, typename TAlgebra>
number HybridSynapseCurrentAssembler<TDomain, TAlgebra>::get_ip3(Vertex* const v, number time)
{
	//check wether v is newly active
	typename std::map<Vertex*, IP3Timing>::iterator it = m_mSynapseActivationTime.find(v);
		if(it == m_mSynapseActivationTime.end()) {//true if v is not in the map yet
			struct IP3Timing t;
			t.t_start = time;
			t.t_end = time + m_j_ip3_duration;
			m_mSynapseActivationTime[v] = t;

			return m_j_ip3_max;

		} else { //v was found and it points to it
			//check wether v is still active
			number t_onset = it->second.t_start, t_end = it->second.t_end;
			if(time >= t_end) {
				//erase v
				m_mSynapseActivationTime.erase(v);
				return 0;

			} else { //v is still active
				return m_j_ip3_max * std::exp(m_j_ip3_decayRate*(t_onset - time));
			}
		}
}

template <typename TDomain, typename TAlgebra>
void HybridSynapseCurrentAssembler<TDomain, TAlgebra>::adjust_defect
(
    vector_type& d,
    const vector_type& u,
    ConstSmartPtr<DoFDistribution> dd,
    int type,
    number time,
    ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
    const std::vector<number>* vScaleMass,
    const std::vector<number>* vScaleStiff
)
{
	// we want to add inward currents to the defect
	// at all vertices representing an active synapse (or more)

	// Using the m_spHNC, get a list of all (3d) vertices on this proc
	// that are mapped to an active (1d) synapse (on any proc);
	// also get the corresponding current values at the same time.
	std::vector<Vertex*> vActiveList;
	std::vector<number> vCurrent;
	m_spHNC->gather_synaptic_currents(vActiveList, vCurrent, time);


	// calculate dt
	// if the following happens, get dt from elsewhere (to be set and updated by user)
	UG_COND_THROW(!vScaleStiff, "No stiffness scales given.");

	// sum of stiffness factors should be dt (shouldn't it!?)
	number dt = 0.0;
	size_t ntp = vScaleStiff->size();
	for (size_t tp = 0; tp < ntp; ++tp)
		dt += (*vScaleStiff)[tp];


	// loop active synapses
	size_t sz = vActiveList.size();
	for (size_t i = 0; i < sz; ++i)
	{
		Vertex* v = vActiveList[i];

		// get the DoFIndex for this vertex
		std::vector<DoFIndex> vIndex;
		dd->inner_dof_indices(v, m_fctInd, vIndex, false);

		UG_COND_THROW(!vIndex.size(), "Function given by 'set_flowing_substance_name' is not defined for "
			<< ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), v) << ".")

		UG_ASSERT(vIndex.size() == 1, "Apparently, you are using shape functions different from P1, this is not supported.");

		const DoFIndex& dofInd = vIndex[0];

		// Evaluate current (use the same current for all time points).
		// This way, we enforce explicit Euler scheme for current discretization.
		// The currents have to be scaled properly in case the units
		// of the 1d and 3d simulations do not match!
		// Add synaptic current * dt to defect.
		// Currents are outward in the synapse handler, so we add to defect.

		// if the potential rises high, synaptic currents are reversed;
		// we need to exclude calcium from this effect
		if (vCurrent[i] > 0.0) vCurrent[i] = 0.0;
		DoFRef(d, dofInd) += dt * vCurrent[i] * m_current_percentage / (m_valency*m_F) / m_scaling_3d_to_1d_amount_of_substance;

		if(!m_ip3_set) return;
		//get DoFIndex for this vertex (ip3)
		std::vector<DoFIndex> vIndex_ip3;
		dd->inner_dof_indices(v, m_fctInd_ip3, vIndex_ip3, false);

		UG_COND_THROW(!vIndex_ip3.size(), "Function for IP3 is not defined for "
			<< ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), v) << ".")

		UG_ASSERT(vIndex_ip3.size() == 1, "Apparently, you are using shape functions different from P1, this is not supported.");

		const DoFIndex& dofInd_ip3 = vIndex_ip3[0];
		DoFRef(d, dofInd_ip3) += get_ip3(v, time) * dt / m_scaling_3d_to_1d_ip3;
	}
}





// explicit template specializations
#ifdef UG_DIM_1
	template class HybridNeuronCommunicator<Domain1d>;
	template class HybridSynapseCurrentAssembler<Domain1d, CPUAlgebra>;
#endif
#ifdef UG_DIM_2
    template class HybridNeuronCommunicator<Domain2d>;
	template class HybridSynapseCurrentAssembler<Domain2d, CPUAlgebra>;
#endif
#ifdef UG_DIM_3
    template class HybridNeuronCommunicator<Domain3d>;
	template class HybridSynapseCurrentAssembler<Domain3d, CPUAlgebra>;
#endif




} // namespace neuro_collection
} // namespace ug
