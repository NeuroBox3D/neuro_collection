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
#include <vector>

namespace ug {
namespace neuro_collection {


template <typename TDomain>
HybridNeuronCommunicator<TDomain>::HybridNeuronCommunicator
(
    SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
    SmartPtr<ApproximationSpace<TDomain> > spApprox1d
)
: m_spCE(SPNULL),
  m_spSynHandler(SPNULL),
  m_mSynapse3dVertex(std::map<synapse_id, Vertex*>()),
#ifdef UG_PARALLEL
  rcvSize(NULL), rcvFrom(NULL), rcvBuf(NULL),
  sendSize(NULL), sendTo(NULL), sendBuf(NULL),
#endif
  m_spApprox1d(spApprox1d), m_spApprox3d(spApprox3d),
  m_spGrid1d(SPNULL), m_spGrid3d(SPNULL), m_spMGSSH3d(SPNULL),
  m_scale_factor_from_3d_to_1d(1.0),
  m_aNID(GlobalAttachments::attachment<ANeuronID>("neuronID"))
{
	m_spGrid1d = m_spApprox1d->domain()->grid();
	m_spGrid3d = m_spApprox3d->domain()->grid();
	m_spMGSSH3d = m_spApprox3d->domain()->subset_handler();

	if (!m_spGrid1d->has_vertex_attachment(m_aNID))
		m_spGrid1d->attach_to_vertices(m_aNID);
	m_aaNID = Grid::VertexAttachmentAccessor<ANeuronID>(*m_spGrid1d, m_aNID);

	m_aaPos1d = m_spApprox1d->domain()->position_accessor();
	m_aaPos3d = m_spApprox1d->domain()->position_accessor();
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
void HybridNeuronCommunicator<TDomain>::set_ce_object(ConstSmartPtr<cable_neuron::CableEquation<TDomain> > spCEDisc)
{
    m_spCE = spCEDisc;
    m_spSynHandler = m_spCE->synapse_handler();
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
void HybridNeuronCommunicator<TDomain>::reinit_potential_mappings()
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

    // TODO: Only find those of a neuron with given ID (or from a set of IDs).
    // find all 1d vertices and store them and their positions
    std::vector<Vertex*> vLocVrt;
    std::vector<posType> vLocVrtPos;

    VrtItType it = dd1->template begin<Vertex>();
    VrtItType it_end = dd1->template end<Vertex>();
    for (; it != it_end; ++it)
    {
       vLocVrt.push_back(*it);
       vLocVrtPos.push_back(aaPos1[*it]);
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
        vIndex.clear();
        dd1->inner_dof_indices(vrt, m_potFctInd, vIndex, false);

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
void HybridNeuronCommunicator<TDomain>::neuron_identification()
{
	// FIXME: use only base level here and propagate to higher levels afterwards
	int nid = 0;
	//initialize with -1
	for(VertexIterator vIter = m_spGrid1d->begin<Vertex>(); vIter != m_spGrid1d->end<Vertex>(); ++vIter ) {
		Vertex* v = *vIter;
		m_aaNID[v] = -1;
	}

	for(VertexIterator vIter = m_spGrid1d->begin<Vertex>(); vIter != m_spGrid1d->end<Vertex>(); ++vIter ) {
		UG_COND_THROW(nid < 0, "Too many neurons");
		Vertex* v = *vIter;
		if(deep_first_search(v, nid)) {
			nid++;
		}
	}
	std::cout << "\nNIDs: " << nid + 1 << std::endl << std::endl; //dbg: prints out number of neurons (+1 because of 0 being first index)
}

template <typename TDomain>
void HybridNeuronCommunicator<TDomain>::prep_timestep(const number& t, const int& id, std::vector<number>& vCurr, std::vector<synapse_id>& vSid)
{
	m_spSynHandler->get_currents(t, id, vCurr, vSid);
}

template <typename TDomain>
int HybridNeuronCommunicator<TDomain>::deep_first_search(Vertex* v, int id)
{
	if(m_aaNID[v] >= 0) {
		return 0; //vertex v is already identified
	} else {
		m_aaNID[v] = id; //v gets id
	}

	Grid::traits<Edge>::secure_container edges; //get all edges incident to v
	m_spGrid1d->associated_elements(edges, v);

	for(size_t i=0; i<edges.size(); ++i) {
		deep_first_search( (*edges[i])[0], id ); //dfs for every connected vertex
		deep_first_search( (*edges[i])[1], id );
	}
	return 1;
}

template <typename TDomain>
int HybridNeuronCommunicator<TDomain>::get_neuron_id(synapse_id id)
{
	for(EdgeIterator eIter = m_spGrid1d->begin<Edge>(); eIter != m_spGrid1d->end<Edge>(); ++eIter) {
		Edge* e = *eIter;
		std::vector<cable_neuron::synapse_handler::IBaseSynapse*> v = m_spSynHandler->get_synapses_on_edge(e);
		for(size_t j=0; j<v.size(); ++j) {
			if(v[j]->id() == id) {
				return m_aaNID[(*e)[0]];
			}
		}
	}
	return -1;
}

template <typename TDomain>
int HybridNeuronCommunicator<TDomain>::Mapping3d(int neuron_id,
												 std::vector<Vertex*>& vMinimizing3dVertices,
												 std::vector<MathVector<dim> >& vMinimizingSynapseCoords)
{
	std::vector<MathVector<dim> > syn_coords_local;
	std::vector<synapse_id> syn_ids_local;

	std::vector<Vertex*> v3dVertices;
	std::vector<Vertex*> vMap;
	std::vector<number> vDistances;

	//gather plasma membrane surface vertices:
	for(unsigned int i=m_vCurrentSubset3d[0]; i<m_vCurrentSubset3d.size(); ++i) {
		geometry_traits<Vertex>::const_iterator vrtIt = m_spMGSSH3d->begin<Vertex>(i, 0);
		geometry_traits<Vertex>::const_iterator vrtItend = m_spMGSSH3d->end<Vertex>(i, 0);

		while(vrtIt != vrtItend) {
			v3dVertices.push_back(*vrtIt);
		}
	}


	//gather local synapse coords, that are interesting
	cable_neuron::synapse_handler::SynapseIter<cable_neuron::synapse_handler::IBaseSynapse> synIt = m_spSynHandler->template begin<cable_neuron::synapse_handler::IBaseSynapse>();
	cable_neuron::synapse_handler::SynapseIter<cable_neuron::synapse_handler::IBaseSynapse> synItEnd = m_spSynHandler->template end<cable_neuron::synapse_handler::IBaseSynapse>();
	for (; synIt != synItEnd; ++synIt)
	{
		cable_neuron::synapse_handler::IBaseSynapse* syn = *synIt;
		if( get_neuron_id(syn->id()) == neuron_id) {

			MathVector<dim> coords;
			get_coordinates(syn->id(), coords);
			syn_coords_local.push_back(coords);
			syn_ids_local.push_back(syn->id());
		}
	}

#ifdef UG_PARALLEL

	pcl::ProcessCommunicator com;
	std::vector<MathVector<dim> > syn_coords_global; //synapses coords of neuron with given id
	std::vector<synapse_id> syn_ids_global;
	std::vector<int> syn_sizes;
	std::vector<int> syn_offsets;

	//communicate all interesting synapse coords to all procs
	com.allgatherv(syn_coords_global, syn_coords_local, &syn_sizes, &syn_offsets); //syn_coords_global contains all synapses, which are interesting for a 3d mapping
	com.allgatherv(syn_ids_global, syn_ids_local, &syn_sizes, &syn_offsets);


	// compute local min distances
	nearest_neighbor(syn_coords_global, v3dVertices, vMap, vDistances);

	// communicate all distances to all processes
	std::vector<number> vGloblDistances;
	std::vector<std::vector<number> > mGloblDistances;

	// TODO: for big process numbers, try to minimize communication:
	//       using gather (in a loop over 1d synapse holder procs),
	//       only send min_distances back to the holders of 1d synapses
	//       and let them decide on global minimum
	com.allgatherv(vGloblDistances, vDistances);
	for(size_t i=0; i<com.size(); ++i) { //assemble global distances matrix
		std::vector<number> tmp;
		for(size_t j=0; j<vDistances.size(); ++j) {
			tmp.push_back(vGloblDistances[i]);
		}
		mGloblDistances.push_back(tmp);
	}

	//compute minimizing proc for each 1d synapse
	std::vector<int> vMinimizingProc;
	for(size_t j=0; j<vDistances.size(); ++j) {//iterate over columns, ie synapses
		number min = mGloblDistances[j][0]; //init global min with proc0's min
		int min_proc = 0;					//minimizing proc
		for(size_t i=0; i<com.size(); ++i) {//iterate over rows, ie local min distances
			if(mGloblDistances[j][i] < min) {
				min = mGloblDistances[j][i];
				min_proc = i;
			}
		}
		vMinimizingProc.push_back(min_proc);
	}

	/* Tried another aproach below commented snippet.
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
			m_mSynapse3dVertex[ syn_ids_global[i] ] = vMap[i];
	}

	// vMinimizing3dVertices should finally contain the nearest 3d Vertex to the corresponding 1d synapse at the same position in vMinimizingSynapseCoords





#else
	nearest_neighbor(syn_coords_local, v3dVertices, vMap, vDistances);

	for(size_t i=0; i<vMap.size(); ++i) {
		m_mSynapse3dVertex[ syn_ids_local[i] ] = vMap[i];
	}
#endif
	return 1;
}

template <typename TDomain>
int HybridNeuronCommunicator<TDomain>::nearest_neighbor(	const std::vector<MathVector<dim> >& v1dCoords,
													const std::vector<Vertex*>& v3dVertices,
													std::vector<Vertex*>& vMap,
													std::vector<number>& vDistances)
{
	vMap.clear();
	vDistances.clear();

	// iterate over 1dSynapses
	for(size_t i=0; i<v1dCoords.size(); i++) {

        number min_distance = VecDistanceSq(m_aaPos3d[v3dVertices[0]], v1dCoords[0]);
		vMap.push_back(v3dVertices[0]);

		// iterate over 3dVertices
		for(size_t j=1; j<v3dVertices.size(); ++j)
		{
			number distance = VecDistanceSq(m_aaPos3d[v3dVertices[j]], v1dCoords[i]);
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
void HybridNeuronCommunicator<TDomain>::get_coordinates(synapse_id id, MathVector<dim>& vCoords)
{

	// Search for Edge
	for(EdgeIterator eIter = m_spGrid1d->begin<Edge>(); eIter != m_spGrid1d->end<Edge>(); ++eIter) {
		Edge* e = *eIter;
		std::vector<cable_neuron::synapse_handler::IBaseSynapse*> v = m_spSynHandler->get_synapses_on_edge(e);
		for(size_t j=0; j<v.size(); ++j) {
			if(v[j]->id() == id) {
				Vertex* v0 = (*e)[0];
				Vertex* v1 = (*e)[1];

				number localcoord = v[j]->location();

				VecScaleAdd(vCoords, localcoord, m_aaPos1d[v0], localcoord, m_aaPos1d[v1]);
			}
		}
	}


}

template <typename TDomain, typename TAlgebra>
HybridSynapseCurrentAssembler<TDomain, TAlgebra>::HybridSynapseCurrentAssembler(ConstSmartPtr<HybridNeuronCommunicator<TDomain> > hnc, const std::string& fct)
:m_fctInd(0), m_F(96485.309), m_valency(2), m_current_percentage(0.1), m_spHNC(hnc)
{
	// get function index of whatever it is that the current carries (in our case: calcium)
	FunctionGroup fctGrp(this->m_spApproxSpace->function_pattern());
	try {fctGrp.add(fct);}
	UG_CATCH_THROW("Function " << fct << " could not be identified in given approximation space.");
	m_fctInd = fctGrp.unique_id(0);
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

	// TODO: Using the m_spHNC, get a list of all (3d) vertices on this proc
	//       that are mapped to an active (1d) synapse (on any proc);
	//       also get the corresponding current values at the same time.
	std::vector<Vertex*> vActiveList;
	std::vector<number> vCurrent;

	const std::vector<synapse_id> vActiveSynapsesId = m_spHNC->synapse_handler()->active_presynapses();

	std::map<synapse_id, Vertex*> mSynVrtMap = m_spHNC->synapse_3dVertex_map();
	for(size_t i=0; i<vActiveSynapsesId.size(); ++i) {
		synapse_id sid = vActiveSynapsesId[i];
		if(mSynVrtMap.find(sid) != mSynVrtMap.end() ) {
			vActiveList.push_back(mSynVrtMap[sid]);
			vCurrent.push_back( m_spHNC->synapse_handler()->current(sid) );
		}
	}


	size_t sz = vActiveList.size();
	for (size_t i = 0; i < sz; ++i)
	{
		Vertex* v = vActiveList[i];

		// get the DoFIndex for this vertex
		std::vector<DoFIndex> vIndex;
		dd->inner_dof_indices(v, m_fctInd, vIndex, false); // we don't need hanging indices as vertices are all from base level

		UG_COND_THROW(!vIndex.size(), "Function given by 'set_flowing_substance_name' is not defined for "
			<< ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), v) << ".")

		UG_ASSERT(vIndex.size() == 1, "Apparently, you are using shape functions different from P1, this is not supported.");

		const DoFIndex& dofInd = vIndex[0];

		// evaluate current (use the same current for all time points)
		// this way, we enforce explicit Euler scheme for current discretization

		// if this happens, get dt from elsewhere (to be set and updated by user)
		UG_COND_THROW(!vScaleStiff, "No stiffness scales given.");

		// sum of stiffness factors should be dt (shouldn't it!?)
		number dt = 0.0;
		size_t ntp = vScaleStiff->size();
		for (size_t tp = 0; tp < ntp; ++tp)
			dt += (*vScaleStiff)[tp];

		// TODO: The currents will have to be scaled properly in case the units
		// of the 1d and 3d simulations do not match!
		// add synaptic current * dt to defect
		// careful, inward current means negative sign in the defect!
		DoFRef(d, dofInd) -= dt * vCurrent[i] / (m_valency*m_F);
	}
}


/*template <typename TDomain, typename TAlgebra>
void HybridSynapseCurrentAssembler<TDomain, TAlgebra>::set_flowing_substance_name(const std::string& fct)
{
	// get function index of whatever it is that the current carries (in our case: calcium)
	FunctionGroup fctGrp(this->m_spApproxSpace->function_pattern());
	try {fctGrp.add(fct);}
	UG_CATCH_THROW("Function " << fct << " could not be identified in given approximation space.");
	m_fctInd = fctGrp.unique_id(0);
}*/






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
