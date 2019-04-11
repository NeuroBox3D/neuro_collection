/*
 * hybrid_synapse_current_assembler_impl.h
 *
 *  Created on: 20.12.2016
 *      Author: mbreit, lreinhardt
 */

#include "lib_grid/algorithms/debug_util.h"  // for ElementDebugInfo
#include "lib_grid/algorithms/volume_calculation.h"  // for CalculateVolume
#include "../cable_neuron/util/functions.h"  // for neuron_identification

#include <algorithm>  // for std::find
#include <vector>

namespace ug {
namespace neuro_collection {


template <typename TDomain, typename TAlgebra>
HybridSynapseCurrentAssembler<TDomain, TAlgebra>::
HybridSynapseCurrentAssembler
(
	SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
	SmartPtr<ApproximationSpace<TDomain> > spApprox1d,
	SmartPtr<cable_neuron::synapse_handler::SynapseHandler<TDomain> > spSH,
	const std::vector<std::string>& plasmaMembraneSubsetName,
	const std::string& fct,
	const std::string& fct_ip3
)
: m_fctInd(0), m_fctInd_ip3(0), m_ip3_set(true), m_F(96485.309), m_valency(2), m_current_percentage(0.1),
  m_spHNC(new hnc_type(spApprox3d, spApprox1d)), m_spDom(spApprox3d->domain()), m_sqSynRadius(0.04),
  m_scaling_3d_to_1d_amount_of_substance(1e-15), m_scaling_3d_to_1d_electric_charge(1.0), m_scaling_3d_to_1d_coordinates(1e-6), m_scaling_3d_to_1d_ip3(1e-15),
  m_j_ip3_max(6e-19), m_j_ip3_decayRate(1.188), m_j_ip3_duration(3.0 / m_j_ip3_decayRate)
{
	// get function index of whatever it is that the current carries (in our case: calcium)
	FunctionGroup fctGrp(spApprox3d->function_pattern());
	try {fctGrp.add(fct);}
	UG_CATCH_THROW("Function " << fct << " could not be identified in given approximation space.");
	m_fctInd = fctGrp.unique_id(0);

	// get function index if use of ip3 is enabled
	fctGrp.remove(fct);
	try {fctGrp.add(fct_ip3);}
	UG_CATCH_THROW("Function " << fct_ip3 << " could not be identified in given approximation space.");
	m_fctInd_ip3 = fctGrp.unique_id(0);

	// init HNC
	m_spHNC->set_synapse_handler(spSH);
	m_spHNC->set_coordinate_scale_factor_3d_to_1d(m_scaling_3d_to_1d_coordinates);
	m_spHNC->set_current_subsets(plasmaMembraneSubsetName);

	// also save plasma membrane subset indices
    SubsetGroup ssGrp;
    try {ssGrp = SubsetGroup(spApprox3d->domain()->subset_handler(), plasmaMembraneSubsetName);}
    UG_CATCH_THROW("Subset group creation failed.");

    const size_t nSs = ssGrp.size();
    for (size_t si = 0; si < nSs; ++si)
    	m_vMembraneSI.push_back(ssGrp[si]);
}


template <typename TDomain, typename TAlgebra>
HybridSynapseCurrentAssembler<TDomain, TAlgebra>::
HybridSynapseCurrentAssembler
(
	SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
	SmartPtr<ApproximationSpace<TDomain> > spApprox1d,
	SmartPtr<cable_neuron::synapse_handler::SynapseHandler<TDomain> > spSH,
	const std::vector<std::string>& plasmaMembraneSubsetName,
	const std::string& fct
)
: m_fctInd(0), m_fctInd_ip3(0), m_ip3_set(false), m_F(96485.309), m_valency(2), m_current_percentage(0.1),
  m_spHNC(new hnc_type(spApprox3d, spApprox1d)), m_spDom(spApprox3d->domain()), m_sqSynRadius(0.04),
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
	m_spHNC->set_current_subsets(plasmaMembraneSubsetName);
}


template <typename TDomain, typename TAlgebra>
number HybridSynapseCurrentAssembler<TDomain, TAlgebra>::get_ip3(synapse_id sid, number time)
{
	// check whether vrt is newly active
	typename std::map<synapse_id, number>::iterator it = m_mSynapseActivationTime.find(sid);
	if (it == m_mSynapseActivationTime.end())
	{
		m_mSynapseActivationTime[sid] = time;
		return m_j_ip3_max;
	}

	// check whether v is still active
	if (time >= it->second + m_j_ip3_duration)
	{
		// erase entry
		m_mSynapseActivationTime.erase(sid);
		return 0.0;
	}

	return m_j_ip3_max * std::exp(m_j_ip3_decayRate*(it->second - time));
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

	// Using the m_spHNC, get a list of all synapse positions
	// that are mapped to a (locally) active (1d) synapse;
	// also get the corresponding current values at the same time.
	typename std::vector<MathVector<dim> > vActiveList;
	std::vector<number> vSynCurrent;
	std::vector<synapse_id> vSynID;

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		typename std::vector<MathVector<dim> > vLocActiveList;
		std::vector<number> vLocSynCurrent;
		std::vector<synapse_id> vLocSynID;
		m_spHNC->gather_synaptic_currents(vLocActiveList, vLocSynCurrent, vLocSynID, time);

		// make all (3d) active synapse positions known to all procs
		pcl::ProcessCommunicator com;
		com.allgatherv(vActiveList, vLocActiveList, NULL, NULL);
		com.allgatherv(vSynCurrent, vLocSynCurrent, NULL, NULL);
		com.allgatherv(vSynID, vLocSynID, NULL, NULL);
	}
	else
#endif
		m_spHNC->gather_synaptic_currents(vActiveList, vSynCurrent, vSynID, time);



	// calculate dt
	// if the following happens, get dt from elsewhere (to be set and updated by user)
	UG_COND_THROW(!vScaleStiff, "No stiffness scales given.");

	// sum of stiffness factors should be dt (shouldn't it!?)
	number dt = 0.0;
	size_t ntp = vScaleStiff->size();
	for (size_t tp = 0; tp < ntp; ++tp)
		dt += (*vScaleStiff)[tp];


	// loop all surface sides in membrane subsets
	typename TDomain::position_accessor_type& aaPos = m_spDom->position_accessor();
	SmartPtr<ISubsetHandler> sh = m_spDom->subset_handler();
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename DoFDistribution::traits<side_type>::const_iterator const_side_iter;

	const size_t nSyn = vActiveList.size();
	std::vector<number> totalSynAreaLocal(nSyn, 0.0);
	std::vector<std::vector<side_type*> > elemsForSyn(nSyn);
	const size_t nSs = m_vMembraneSI.size();
	for (size_t ss = 0; ss < nSs; ++ss)
	{
		const int si = m_vMembraneSI[ss];
		const_side_iter it = dd->begin<side_type>(si);
		const_side_iter itEnd = dd->end<side_type>(si);
		for (; it != itEnd; ++it)
		{
			side_type* side = *it;

			// loop all active synapse positions and find out whether
			// the current element is in their range
			for (size_t s = 0; s < nSyn; ++s)
			{
				if (VecDistanceSq(vActiveList[s], CalculateCenter(side, aaPos)) < m_sqSynRadius)
				{
					totalSynAreaLocal[s] += CalculateVolume(side, aaPos);
					elemsForSyn[s].push_back(side);
				}
				// at least directly adjacent elements need to be used
				else
				{
					const size_t nVrt = side->num_vertices();
					for (size_t v = 0; v < nVrt; ++v)
					{
						if (VecDistanceSq(aaPos[side->vertex(v)], vActiveList[s]) < 1e-10*m_sqSynRadius)
						{
							totalSynAreaLocal[s] += CalculateVolume(side, aaPos);
							elemsForSyn[s].push_back(side);
							break;
						}
					}
				}
			}
		}
	}

	std::vector<number> totalSynArea = totalSynAreaLocal;
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		pc.allreduce(totalSynAreaLocal, totalSynArea, PCL_RO_SUM);
	}
#endif

	// now treat all synapses
	std::vector<DoFIndex> vDoFIndex;
	for (size_t s = 0; s < nSyn; ++s)
	{
		const size_t nElems = elemsForSyn[s].size();
		if (!nElems)
			continue;

		const std::vector<side_type*>& vElems = elemsForSyn[s];

		// if the potential rises high, synaptic currents are reversed;
		// we need to exclude calcium from this effect
		// TODO: this is a bit awkward, better use proper Ca2+ entry modeling
		if (vSynCurrent[s] < 0.0)
		{
			const number substanceCurrent = vSynCurrent[s] * m_current_percentage / (m_valency*m_F) / m_scaling_3d_to_1d_amount_of_substance;
			const number fluxDensity = substanceCurrent / totalSynArea[s];

			// loop all elems participating in that synapse
			for (size_t e = 0; e < nElems; ++e)
			{
				side_type* elem = vElems[e];
				const size_t nSynElemVrts = elem->num_vertices();

				// each node gets an equal part of this synElem's current
				const number currentPerNode = fluxDensity * CalculateVolume(elem, aaPos) / nSynElemVrts;
				for (size_t n = 0; n < nSynElemVrts; ++n)
				{
					Vertex* node = elem->vertex(n);

					// get the DoFIndex for this vertex
					dd->inner_dof_indices(node, m_fctInd, vDoFIndex, true);
					UG_COND_THROW(!vDoFIndex.size(), "Function for flowing substance is not defined for "
						<< ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), node) << ".")

					UG_ASSERT(vDoFIndex.size() == 1, "Apparently, you are using shape functions different from P1,"
						" this is not supported.");
					const DoFIndex& dofInd = vDoFIndex[0];

					// currents are outward in the synapse handler, so we _add_ to defect
					DoFRef(d, dofInd) += dt * currentPerNode;
				}
			}
		}

		// same for IP3 currents
		if (!m_ip3_set)
			return;

		const synapse_id sid = vSynID[s];

		const number substanceCurrent = get_ip3(sid, time) / m_scaling_3d_to_1d_ip3;
		const number fluxDensity = substanceCurrent / totalSynArea[s];

		// loop all elems participating in that synapse
		for (size_t e = 0; e < nElems; ++e)
		{
			side_type* elem = vElems[e];
			const size_t nSynElemVrts = elem->num_vertices();

			// each node gets an equal part of this synElem's current
			const number currentPerNode = fluxDensity * CalculateVolume(elem, aaPos) / nSynElemVrts;
			for (size_t n = 0; n < nSynElemVrts; ++n)
			{
				Vertex* node = elem->vertex(n);

				// get the DoFIndex for this vertex
				dd->inner_dof_indices(node, m_fctInd_ip3, vDoFIndex, true);
				UG_COND_THROW(!vDoFIndex.size(), "Function for IP3 is not defined for "
					<< ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), node) << ".")

				UG_ASSERT(vDoFIndex.size() == 1, "Apparently, you are using shape functions different from P1,"
					" this is not supported.");
				const DoFIndex& dofInd = vDoFIndex[0];

				// currents are inward here, so we _subtract_ from defect
				DoFRef(d, dofInd) -= dt * currentPerNode;
			}
		}
	}
}


template <typename TDomain, typename TAlgebra>
void HybridSynapseCurrentAssembler<TDomain, TAlgebra>::
adjust_error
(
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const std::vector<number>* vScaleMass,
	const std::vector<number>* vScaleStiff
)
{
	//	get the error estimator data object and check that it is of the right type
	UG_COND_THROW(!this->m_spErrEstData.get(),
		"No ErrEstData object has been given to HybridSynapseCurrentAssembler.");

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());
	UG_COND_THROW(!err_est_data, "Dynamic cast to MultipleSideAndElemErrEstData failed."
		<< std::endl << "Make sure you handed the correct type of ErrEstData to this discretization.");


	// Using the m_spHNC, get a list of all synapse positions
	// that are mapped to a (locally) active (1d) synapse;
	// also get the corresponding current values at the same time.
	typename std::vector<MathVector<dim> > vActiveList;
	std::vector<number> vSynCurrent;
	std::vector<synapse_id> vSynID;

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		typename std::vector<MathVector<dim> > vLocActiveList;
		std::vector<number> vLocSynCurrent;
		std::vector<synapse_id> vLocSynID;
		m_spHNC->gather_synaptic_currents(vLocActiveList, vLocSynCurrent, vLocSynID, time);

		// make all (3d) active synapse positions known to all procs
		pcl::ProcessCommunicator com;
		com.allgatherv(vActiveList, vLocActiveList, NULL, NULL);
		com.allgatherv(vSynCurrent, vLocSynCurrent, NULL, NULL);
		com.allgatherv(vSynID, vLocSynID, NULL, NULL);
	}
	else
#endif
		m_spHNC->gather_synaptic_currents(vActiveList, vSynCurrent, vSynID, time);


	// calculate dt
	// if the following happens, get dt from elsewhere (to be set and updated by user)
	UG_COND_THROW(!vScaleStiff, "No stiffness scales given.");

	// sum of stiffness factors should be dt (shouldn't it!?)
	number dt = 0.0;
	size_t ntp = vScaleStiff->size();
	for (size_t tp = 0; tp < ntp; ++tp)
		dt += (*vScaleStiff)[tp];


	// loop all surface sides in membrane subsets
	typename TDomain::position_accessor_type& aaPos = m_spDom->position_accessor();
	SmartPtr<ISubsetHandler> sh = m_spDom->subset_handler();
	typedef typename domain_traits<TDomain::dim>::side_type side_type;
	typedef typename DoFDistribution::traits<side_type>::const_iterator const_side_iter;

	const size_t nSyn = vActiveList.size();
	std::vector<number> totalSynAreaLocal(nSyn, 0.0);
	std::vector<std::vector<side_type*> > elemsForSyn(nSyn);
	const size_t nSs = m_vMembraneSI.size();
	for (size_t ss = 0; ss < nSs; ++ss)
	{
		const int si = m_vMembraneSI[ss];
		const_side_iter it = dd->begin<side_type>(si);
		const_side_iter itEnd = dd->end<side_type>(si);
		for (; it != itEnd; ++it)
		{
			side_type* side = *it;

			// loop all active synapse positions and find out whether
			// the current element is in their range
			for (size_t s = 0; s < nSyn; ++s)
			{
				if (VecDistanceSq(vActiveList[s], CalculateCenter(side, aaPos)) < m_sqSynRadius)
				{
					totalSynAreaLocal[s] += CalculateVolume(side, aaPos);
					elemsForSyn[s].push_back(side);
				}
				// at least directly adjacent elements need to be used
				else
				{
					const size_t nVrt = side->num_vertices();
					for (size_t v = 0; v < nVrt; ++v)
					{
						if (VecDistanceSq(aaPos[side->vertex(v)], vActiveList[s]) < 1e-10*m_sqSynRadius)
						{
							totalSynAreaLocal[s] += CalculateVolume(side, aaPos);
							elemsForSyn[s].push_back(side);
							break;
						}
					}
				}
			}
		}
	}

	std::vector<number> totalSynArea = totalSynAreaLocal;
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator pc;
		pc.allreduce(totalSynAreaLocal, totalSynArea, PCL_RO_SUM);
	}
#endif


	// now treat all synapses
	for (size_t s = 0; s < nSyn; ++s)
	{
		const size_t nElems = elemsForSyn[s].size();
		if (!nElems)
			continue;

		const std::vector<side_type*>& vElems = elemsForSyn[s];

		// if the potential rises high, synaptic currents are reversed;
		// we need to exclude calcium from this effect
		// TODO: this is a bit awkward, better use proper Ca2+ entry modeling
		if (vSynCurrent[s] < 0.0)
		{
			const number substanceCurrent = vSynCurrent[s] * m_current_percentage / (m_valency*m_F) / m_scaling_3d_to_1d_amount_of_substance;
			const number fluxDensity = substanceCurrent / totalSynArea[s];

			// loop all elems participating in that synapse
			for (size_t e = 0; e < nElems; ++e)
			{
				side_type* elem = vElems[e];

				// get reference object id
				ReferenceObjectID roid = elem->reference_object_id();

				// get corner coords (for later use in calculating global IPs)
				std::vector<typename TDomain::position_type> vCoCo;
				CollectCornerCoordinates(vCoCo, elem, *m_spDom, false);

				// substract constant flux density value from every IP on the side
				size_t numSideIPs;
				const MathVector<side_type::dim>* sideLocIPs;
				const MathVector<TDomain::dim>* sideGlobIPs;
				try
				{
					numSideIPs = err_est_data->get(m_fctInd)->num_side_ips(roid);
					sideLocIPs = err_est_data->get(m_fctInd)->template side_local_ips<side_type::dim>(roid);
					sideGlobIPs = err_est_data->get(m_fctInd)->side_global_ips(elem, &vCoCo[0]);
				}
				UG_CATCH_THROW("Global integration points for error estimator cannot be determined.");

				for (size_t ip = 0; ip < numSideIPs; ++ip)
					(*err_est_data->get(m_fctInd))(elem, ip) -= dt * fluxDensity;
			}
		}

		// same for IP3 currents
		if (!m_ip3_set)
			return;

		const synapse_id sid = vSynID[s];

		const number substanceCurrent = get_ip3(sid, time) / m_scaling_3d_to_1d_ip3;
		const number fluxDensity = substanceCurrent / totalSynArea[s];

		// loop all elems participating in that synapse
		for (size_t e = 0; e < nElems; ++e)
		{
			side_type* elem = vElems[e];

			// get reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			// get corner coords (for later use in calculating global IPs)
			std::vector<typename TDomain::position_type> vCoCo;
			CollectCornerCoordinates(vCoCo, elem, *m_spDom, false);

			size_t numSideIPs;
			const MathVector<side_type::dim>* sideLocIPs;
			const MathVector<TDomain::dim>* sideGlobIPs;
			try
			{
				numSideIPs = err_est_data->get(m_fctInd_ip3)->num_side_ips(roid);
				sideLocIPs = err_est_data->get(m_fctInd_ip3)->template side_local_ips<side_type::dim>(roid);
				sideGlobIPs = err_est_data->get(m_fctInd_ip3)->side_global_ips(elem, &vCoCo[0]);
			}
			UG_CATCH_THROW("Global integration points for error estimator cannot be determined.");

			for (size_t ip = 0; ip < numSideIPs; ++ip)
				(*err_est_data->get(m_fctInd_ip3))(elem, ip) += dt * fluxDensity;
		}
	}
}


} // namespace neuro_collection
} // namespace ug
