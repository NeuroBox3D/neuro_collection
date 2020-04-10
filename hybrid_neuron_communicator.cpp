/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Lukas Reinhardt
 * Creation date: 2016-12-20
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

#include "hybrid_neuron_communicator.h"
#include "lib_grid/algorithms/debug_util.h" // ElementDebugInfo
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
: m_spGridDistributionCallbackID(SPNULL),
  m_spSynHandler(SPNULL),
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
  m_bPotentialMappingNeedsUpdate(true),
  m_bSynapseMappingNeedsUpdate(true)
{
	if (!m_spGrid1d->has_vertex_attachment(m_aNID))
		m_spGrid1d->attach_to_vertices(m_aNID);
	m_aaNID = Grid::VertexAttachmentAccessor<ANeuronID>(*m_spGrid1d, m_aNID);

	// calculate identifiers for each neuron
	cable_neuron::neuron_identification(*m_spGrid1d);

	// set this object as listener for distribution events
	m_spGridDistributionCallbackID = m_spGrid1d->message_hub()->register_class_callback(this,
		&HybridNeuronCommunicator<TDomain>::grid_distribution_callback);
	m_spGridDistributionCallbackID = m_spGrid3d->message_hub()->register_class_callback(this,
		&HybridNeuronCommunicator<TDomain>::grid_distribution_callback);

	m_spGridAdaptionCallbackID = m_spGrid1d->message_hub()->register_class_callback(this,
		&HybridNeuronCommunicator<TDomain>::grid_adaption_callback);
	m_spGridAdaptionCallbackID = m_spGrid3d->message_hub()->register_class_callback(this,
		&HybridNeuronCommunicator<TDomain>::grid_adaption_callback);
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

    m_spGrid1d->message_hub()->unregister_callback(m_spGridDistributionCallbackID);
    m_spGrid1d->message_hub()->unregister_callback(m_spGridAdaptionCallbackID);
    m_spGrid3d->message_hub()->unregister_callback(m_spGridDistributionCallbackID);
    m_spGrid3d->message_hub()->unregister_callback(m_spGridAdaptionCallbackID);
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
void HybridNeuronCommunicator<TDomain>::reinit_potential_mapping()
{
	if (!m_bPotentialMappingNeedsUpdate)
		return;

    typedef typename DoFDistribution::traits<vm_grid_object>::const_iterator VmElemItType;
    typedef typename DoFDistribution::traits<Vertex>::const_iterator VrtItType;
    typedef typename TDomain::position_accessor_type posAccType;

    const posAccType& aaPos1 = m_spApprox1d->domain()->position_accessor();
    const posAccType& aaPos3 = m_spApprox3d->domain()->position_accessor();
    ConstSmartPtr<DoFDistribution> dd1 = m_spApprox1d->dof_distribution(GridLevel());
    ConstSmartPtr<DoFDistribution> dd3 = m_spApprox3d->dof_distribution(GridLevel());

	// delete old potential value mappings
	m_mElemPot.clear();

    // find all elements of the 3d geometry boundary on which the potential is defined
    // store their center coords in the same order
    std::vector<posType> vLocPotElemPos;
    std::vector<vm_grid_object*> vLocPotElems;
    size_t numSs = m_vPotSubset3d.size();
    for (size_t s = 0; s < numSs; ++s)
    {
        int si = m_vPotSubset3d[s];
        VmElemItType it = dd3->template begin<vm_grid_object>(si);
        VmElemItType it_end = dd3->template end<vm_grid_object>(si);

        for (; it != it_end; ++it)
        {
        	vLocPotElems.push_back(*it);
            vLocPotElemPos.push_back(aaPos3[*it]);
            vLocPotElemPos.back() *= m_scale_factor_from_3d_to_1d;
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
        vOffsets.resize(nProcs+1, vGlobVrtPos.size());
        std::map<size_t, std::vector<int> > mIndex;
        m_mReceiveInfo.clear();
        for (size_t i = 0; i < nPE; ++i)
        {
        	size_t pos = vNearest[i];

        	// perform a binary search for the largest entry in offset
        	// that is lower than or equal to pos
        	size_t proc = std::distance(vOffsets.begin(),
        		std::upper_bound(vOffsets.begin(), vOffsets.end(), pos)) - 1;

        	mIndex[proc].push_back((int) pos - vOffsets[proc]);

        	// add to recvInfo
        	m_mReceiveInfo[(int) proc].push_back(vLocPotElems[i]);
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
        std::map<size_t, std::vector<int> >::const_iterator it = mIndex.begin();
        std::map<size_t, std::vector<int> >::const_iterator itEnd = mIndex.end();
		for (; it != itEnd; ++it)
		{
			const size_t p = it->first;
			int sz = it->second.size();
			vNumTo[p] = sz;

			++numRecverProcs;
			recverProcs.push_back(p);
			sendSizes.push_back(sz * sizeof(int));
			for (size_t i = 0; i < (size_t) sz; ++i)
				sendBuffer.push_back(it->second[i]);
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


        // at last, prepare communication arrays
        // delete present ones if this is a re-initialization
        if (rcvSize) delete[] rcvSize;
        if (rcvFrom) delete[] rcvFrom;
        if (rcvBuf) delete[] (char*) rcvBuf;
        if (sendSize) delete[] sendSize;
        if (sendTo) delete[] sendTo;
        if (sendBuf)  delete[] (char*) sendBuf;

        // receiving setup
        size_t numRcv = m_mReceiveInfo.size();
        rcvSize = new int[numRcv];
        rcvFrom = new int[numRcv];
        size_t rcvBytes = 0;

        size_t i = 0;
        typename std::map<int, std::vector<vm_grid_object*> >::const_iterator itRec = m_mReceiveInfo.begin();
        typename std::map<int, std::vector<vm_grid_object*> >::const_iterator itRec_end = m_mReceiveInfo.end();
        for (; itRec != itRec_end; ++itRec)
        {
            rcvFrom[i] = itRec->first;
            rcvSize[i] = (int) (itRec->second.size() * sizeof(number));
            rcvBytes += rcvSize[i];
            ++i;
        }
        rcvBuf = new char[rcvBytes];

        // sending setup
        size_t numSend = m_mSendInfo.size();
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

    m_bPotentialMappingNeedsUpdate = false;
    return;

serial_case:
#endif
	// local nearest neighbor search for all elem centers
	std::vector<size_t> vNearest;
	std::vector<typename posType::value_type> vDist;
	int noNeighbors = nearest_neighbor_search(vLocPotElemPos, vLocVrtPos, vNearest, vDist);
	UG_COND_THROW(noNeighbors, "No 1d vertices are present to map 3d potential elements to.");

	// delete old mapping
	m_mPotElemToVertex.clear();

	// fill 3d->1d map
	size_t i = 0;
	for (size_t s = 0; s < numSs; ++s)
	{
		int si = m_vPotSubset3d[s];
		VmElemItType it = dd3->template begin<vm_grid_object>(si);
		VmElemItType it_end = dd3->template end<vm_grid_object>(si);
		for (; it != it_end; ++it)
		{
		    m_mPotElemToVertex[*it] = vLocVrt[vNearest[i]];
			++i;
		}
	}

	m_bPotentialMappingNeedsUpdate = false;
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::coordinate_potential_values()
{
	// perform 3d elem -> 1d vertex and 1d synapse -> 3d vertex mappings
	// if not yet done
	reinit_potential_mapping();

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
        typename std::map<int, std::vector<vm_grid_object*> >::const_iterator itRec = m_mReceiveInfo.begin();
        typename std::map<int, std::vector<vm_grid_object*> >::const_iterator itRec_end = m_mReceiveInfo.end();
        for (; itRec != itRec_end; ++itRec)
        {
            const std::vector<vm_grid_object*>& vElems = itRec->second;
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

    typename std::map<vm_grid_object*, Vertex*>::iterator it = m_mPotElemToVertex.begin();
    typename std::map<vm_grid_object*, Vertex*>::iterator it_end = m_mPotElemToVertex.end();
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
number HybridNeuronCommunicator<TDomain>::potential(vm_grid_object* elem) const
{
	typename std::map<vm_grid_object*, number>::const_iterator it = m_mElemPot.find(elem);
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
	if (!m_bSynapseMappingNeedsUpdate)
		return;

	// clear previous synapse mapping
	m_mSynapse3dCoords.clear();

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
				syn_coords_local.push_back(coords /= m_scale_factor_from_3d_to_1d);
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

		const size_t nGlobSyn = syn_coords_global.size();

		// compute min distances of all local 3d vertices to the global 1d synapse positions
        std::vector<size_t> vNearest;
    	std::vector<number> vDistances;
		int failure = nearest_neighbor_search(syn_coords_global, v3dVertexPos, vNearest, vDistances);

		std::vector<number> vLocNearestPos(nGlobSyn*dim);
		for (size_t s = 0; s < nGlobSyn; ++s)
		{
			if (!failure)
				for (int d = 0; d < dim; ++d)
					vLocNearestPos[s*dim+d] = v3dVertexPos[vNearest[s]][d];
			else
				for (int d = 0; d < dim; ++d)
					vLocNearestPos[s*dim+d] = std::numeric_limits<typename posType::value_type>::quiet_NaN();
		}


		// communicate local min dists back to original 1d-synapse holders
		const size_t nSyn1d = syn_coords_local.size();
		number* globMinPos = new number[nProcs * nSyn1d * dim];
		for (size_t p = 0; p < nProcs; ++p)
		{
			if (!vSizes[p]) continue;
			com.gather((void*) &vLocNearestPos[vOffsets[p]*dim], vSizes[p]*dim, PCL_DT_DOUBLE,
				(void*) globMinPos, vSizes[p]*dim, PCL_DT_DOUBLE, (int) p);
		}

		// compute minimizing proc for each 1d synapse
		for (size_t s = 0; s < nSyn1d; ++s)
		{
			//size_t& minProc = vMinProc[s];
			MathVector<dim>& minPos = m_mSynapse3dCoords[syn_ids_local[s]];
			MathVector<dim> pos;
			number minDistSq = std::numeric_limits<number>::max();
			for (size_t p = 0; p < nProcs; ++p)
			{
				for (int d = 0; d < dim; ++d)
					pos[d] = globMinPos[p*nSyn1d*dim + s*dim + d];

				const number distSq = VecDistanceSq(pos, syn_coords_local[s]);
				if (distSq < minDistSq)
				{
					minDistSq = distSq;
					minPos = pos;
				}
			}
			UG_COND_THROW(minDistSq == std::numeric_limits<number>::max(),
				"No 3d vertex in defined plasma membrane subset present on any proc.");
		}

		delete[] globMinPos;
	}
	else
	{
#endif
		std::vector<size_t> vNearest;
		std::vector<typename posType::value_type> vDistances;
		int failure = nearest_neighbor_search(syn_coords_local, v3dVertexPos, vNearest, vDistances);
		UG_COND_THROW(failure, "No 3d vertices in the defined plasma membrane subset present to map 1d synapses to.");

		for (size_t i = 0; i < syn_coords_local.size(); ++i)
			m_mSynapse3dCoords[syn_ids_local[i]] = v3dVertexPos[vNearest[i]];
#ifdef UG_PARALLEL
	}
#endif

	m_bSynapseMappingNeedsUpdate = false;
}

template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::gather_synaptic_currents
(
	std::vector<MathVector<dim> >& vActSynPosOut,
	std::vector<number>& vSynCurrOut,
	std::vector<synapse_id>& vSynIDOut,
	number time
)
{
	// reinit mappings if necessary
	reinit_synapse_mapping();

	vActSynPosOut.clear();
	vSynCurrOut.clear();
	vSynIDOut.clear();

	// get locally active synapses
	std::vector<synapse_id> vLocActSyn;
	std::vector<number> vLocSynCurr;
	m_spSynHandler->active_postsynapses_and_currents(vLocActSyn, vLocSynCurr, m_vNid, m_aaNID, time);

	const size_t numLocActSyn = vLocActSyn.size();
	for (size_t i = 0; i < numLocActSyn; ++i)
	{
		synapse_id sid = vLocActSyn[i];
		typename std::map<synapse_id, MathVector<dim> >::const_iterator it;
		if ((it = m_mSynapse3dCoords.find(sid)) != m_mSynapse3dCoords.end())
		{
			vActSynPosOut.push_back(it->second);
			vSynCurrOut.push_back(vLocSynCurr[i]);
			vSynIDOut.push_back(sid);
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


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::grid_adaption_callback(const GridMessage_Adaption& gma)
{
	// after grid adaption, mappings need to be force-updated
	if (gma.adaption_ends())
	{
		m_bPotentialMappingNeedsUpdate = true;
		m_bSynapseMappingNeedsUpdate = true;
	}
}


template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::grid_distribution_callback(const GridMessage_Distribution& gmd)
{
	// after grid distribution, mappings need to be force-updated
	if (gmd.msg() == GMDT_DISTRIBUTION_STOPS)
	{
		m_bPotentialMappingNeedsUpdate = true;
		m_bSynapseMappingNeedsUpdate = true;
	}
}





// explicit template specializations
#ifdef UG_DIM_3
    template class HybridNeuronCommunicator<Domain3d>;
#endif




} // namespace neuro_collection
} // namespace ug
