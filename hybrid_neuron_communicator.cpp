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
    size_t numSs = m_vPotSubset3d.size();
    for (size_t s = 0; s < numSs; ++s)
    {
        int si = m_vPotSubset3d[s];
        SideItType it = dd3->template begin<side_t>(si);
        SideItType it_end = dd3->template end<side_t>(si);

        for (; it != it_end; ++it)
            vLocPotElemPos.push_back(CalculateCenter(*it, aaPos3).operator*=(m_scale_factor_from_3d_to_1d));
    }

    // find all 1d vertices for neuron IDs of interest and store them and their positions
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
    		}
    	}
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

        // NN search returns 1 if no neighbors found
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
            if (receiver < (int)(vOffsets.size())-1 && (int)i >= vOffsets[receiver+1])
            {
                // push old receiver if any NN for it are present
                if (vVrt.size())
                {
                    m_mSendInfo[receiver] = vVrt;
                    vVrt.clear();
                }

                // increase receiver index as long as necessary
                ++receiver;
                while (receiver < (int)(vOffsets.size())-1 && (int)i >= vOffsets[receiver+1])
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
            SideItType it = dd3->template begin<side_t>(si);
            SideItType it_end = dd3->template end<side_t>(si);
            for (; it != it_end; ++it)
            {
                m_mReceiveInfo[globRank[i]].push_back(*it);
                ++i;
            }
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

    // TODO: The potentials will have to be scaled properly in case the units
    // of the 1d and 3d simulations do not match!
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



template <typename TDomain>
uint HybridNeuronCommunicator<TDomain>::get_neuron_id(synapse_id id)
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
	std::vector<Vertex*> vMap;
	std::vector<number> vDistances;

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
			v3dVertices.push_back(*it);
	}

	// gather local synapse coords, that are interesting
	typedef cable_neuron::synapse_handler::IPostSynapse PostSynapse;
	typedef cable_neuron::synapse_handler::SynapseIter<PostSynapse> PostSynIter;
	PostSynIter synIt = m_spSynHandler->template begin<PostSynapse>();
	PostSynIter synItEnd = m_spSynHandler->template end<PostSynapse>();
	size_t nNeuron = m_vNid.size();
	for (; synIt != synItEnd; ++synIt)
	{
		cable_neuron::synapse_handler::IPostSynapse* syn = *synIt;
		for (size_t n = 0; n < nNeuron; ++n)
		{
			if (get_neuron_id(syn->id()) == m_vNid[n])
			{
				MathVector<dim> coords;
				get_coordinates(syn->id(), coords);
				syn_coords_local.push_back(coords);
				syn_ids_local.push_back(syn->id());
			}
		}
	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		std::vector<MathVector<dim> > syn_coords_global; //synapses coords of neuron with given id
		std::vector<synapse_id> syn_ids_global;
		std::vector<int> syn_sizes;
		std::vector<int> syn_offsets;

		//communicate all interesting synapse coords to all procs
		com.allgatherv(syn_coords_global, syn_coords_local, &syn_sizes, &syn_offsets); //syn_coords_global contains all synapses, which are interesting for a 3d mapping
		com.allgatherv(syn_ids_global, syn_ids_local, NULL, NULL);

		// compute min distances of all local 3d vertices to the global 1d synapse positions
		nearest_neighbor(syn_coords_global, v3dVertices, vMap, vDistances);

		// communicate all distances to all processes
		std::vector<number> vGloblDistances;
		std::vector<std::vector<number> > mGloblDistances;

		// TODO: for big process numbers, try to minimize communication:
		//       using gather (in a loop over 1d synapse holder procs),
		//       only send min_distances back to the holders of 1d synapses
		//       and let them decide on global minimum
		com.allgatherv(vGloblDistances, vDistances);

		//compute minimizing proc for each 1d synapse
		std::vector<int> vMinimizingProc(vDistances.size());
		for (size_t syn = 0; syn < vDistances.size(); ++syn)
		{
			number min = vGloblDistances[syn]; //init global min with proc0's min
			int min_proc = 0;					//minimizing proc
			for (size_t proc = 1; proc < com.size(); ++proc)
			{
				number& dist = vGloblDistances[proc * vDistances.size() + syn] ;
				if (vGloblDistances[proc * vDistances.size() + syn] < min)
				{
					min = dist;
					min_proc = proc;
				}
			}
			vMinimizingProc[syn] = min_proc;
		}

		/* Tried another approach below commented snippet.
		//construct index array sorted by proc (ascending order)

		std::vector<size_t> synapse_index;
		std::vector<size_t> sizes_synapse_index;
		std::vector<size_t> offsets_synapse_index;

		size_t offset=0;
		for(size_t i=0; i<vMinimizingProc.size(); ++i) {
			size_t size=0;										//number of elements for proc i
			offsets_synapse_index.push_back(offset);
			for(size_t j=0; j<vMinimizingProc.size(); ++j) {
				if( static_cast<size_t>(vMinimizingProc[j]) == i) {
					synapse_index.push_back(j);
					size++;
				}
			}
			sizes_synapse_index.push_back(size);
			offset += size;
		}


		// every proc constructs its minimizing vector of vertices
		// with offsets and sizes for the range of proc_i then follows:
		// range_i = [offsets_synapse_index, offsets_synapse_index + sizes_synapse_index)
		vMinimizing3dVertices.clear();
		vMinimizingSynapseCoords.clear();

		size_t loc_start = offsets_synapse_index[pcl::ProcRank()];
		size_t loc_end = offsets_synapse_index[pcl::ProcRank()] + sizes_synapse_index[pcl::ProcRank()];
		for(size_t i=loc_start; i<loc_end; ++i) {
			vMinimizing3dVertices.push_back(vMap[ synapse_index[i] ]);
			vMinimizingSynapseCoords.push_back( syn_coords_global[synapse_index[i]] );
		}
		//vMinimizing3dVertices should finally contain the nearest 3d Vertex to the corresponding 1d synapse at the same position in vMinimizingSynapseCoords
		*/

		//Create local HNC map out of vMap and vMinimizingProc
		for(size_t i=0; i<vMap.size(); ++i) {
			if(pcl::ProcRank() == vMinimizingProc[i] )
			{
				m_mSynapse3dVertex[ syn_ids_global[i] ] = vMap[i];
				//UG_LOGN("Synapse " << syn_ids_global[i] << " mapped to position "
				//		<< m_aaPos3d[vMap[i]] << "." << std::endl);
			}
		}

		// vMinimizing3dVertices should finally contain the nearest 3d Vertex to the corresponding 1d synapse at the same position in vMinimizingSynapseCoords
	}
	else
	{
#endif
		nearest_neighbor(syn_coords_local, v3dVertices, vMap, vDistances);

		for(size_t i=0; i<vMap.size(); ++i) {
			m_mSynapse3dVertex[ syn_ids_local[i] ] = vMap[i];
		}
#ifdef UG_PARALLEL
	}
#endif
}

template <typename TDomain>
int HybridNeuronCommunicator<TDomain>::nearest_neighbor
(
	const std::vector<MathVector<dim> >& v1dCoords,
	const std::vector<Vertex*>& v3dVertices,
	std::vector<Vertex*>& vMap,
	std::vector<number>& vDistances
)
{
	vMap.clear();
	vDistances.clear();

	// iterate over 1dSynapses
	for(size_t i=0; i<v1dCoords.size(); i++) {

		MathVector<dim> scaled3dVec;	// TODO: avoid possible segfault if v3dvertices is empty
		VecScale(scaled3dVec, m_aaPos3d[v3dVertices[0]], m_scale_factor_from_3d_to_1d);
        number min_distance = VecDistanceSq(scaled3dVec, v1dCoords[i]);
		vMap.push_back(v3dVertices[0]);

		// iterate over 3dVertices
		for(size_t j=1; j<v3dVertices.size(); ++j)
		{
			VecScale(scaled3dVec, m_aaPos3d[v3dVertices[j]], m_scale_factor_from_3d_to_1d);
			number distance = VecDistanceSq(scaled3dVec, v1dCoords[i]);
			if(distance < min_distance) {
				min_distance = distance;
				vMap[i] = v3dVertices[j];
			}
		}
		vDistances.push_back(min_distance); //after comparisons with all 3dVertices add min_distance to the output vector
	}
	return 1;
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

	// todo: only get the synapses and currents located on the neurons of interest;
	//       we need to pass the neuron_ids to the syn_handler for this
	m_spSynHandler->active_postsynapses_and_currents(vLocActSyn, vLocSynCurr, time);
	UG_LOGN("no of locally active post-synapses: " << vLocActSyn.size());

#ifdef UG_PARALLEL
	if (pcl::NumProcs() == 1)
		goto serial_case;

	{
		std::vector<synapse_id> vGlobActSyn;
		std::vector<number> vGlobSynCurr;

		// todo: only send to the procs that hold the 3d vertex to any active synapses
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
				UG_LOG_ALL_PROCS("Found an active synapse!" << std::endl
					<< "with ID " << it->first << " and current " << vGlobSynCurr[i]
					<< "." << std::endl);
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
void HybridNeuronCommunicator<TDomain>::get_coordinates(synapse_id id, MathVector<dim>& vCoords)
{
	// TODO: This implementation is unfavorable; avoid the edge iteration in every call.

	// Search for Edge
	for(EdgeIterator eIter = m_spGrid1d->begin<Edge>(); eIter != m_spGrid1d->end<Edge>(); ++eIter) {
		Edge* e = *eIter;
		std::vector<cable_neuron::synapse_handler::IBaseSynapse*> v = m_spSynHandler->get_synapses_on_edge(e);
		for(size_t j=0; j<v.size(); ++j) {
			if(v[j]->id() == id) {
				Vertex* v0 = (*e)[0];
				Vertex* v1 = (*e)[1];

				number localcoord = v[j]->location();

				VecScaleAdd(vCoords, 1.0 - localcoord, m_aaPos1d[v0], localcoord, m_aaPos1d[v1]);

				return;
			}
		}
	}
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
: m_fctInd(0), m_F(96485.309), m_valency(2), m_current_percentage(0.1), m_spHNC(new hnc_type(spApprox3d, spApprox1d)),
  m_scaling_3d_to_1d_amount_of_substance(1e-15), m_scaling_3d_to_1d_electric_charge(1.0), m_scaling_3d_to_1d_coordinates(1e-6)
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
		DoFRef(d, dofInd) += dt * vCurrent[i] * m_current_percentage / (m_valency*m_F) / m_scaling_3d_to_1d_amount_of_substance;
		//UG_LOGN("Adjusted defect by the amount " << dt * vCurrent[i] * m_current_percentage / (m_valency*m_F) / m_scaling_3d_to_1d_amount_of_substance);
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
