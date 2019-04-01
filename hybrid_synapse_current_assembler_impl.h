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
number HybridSynapseCurrentAssembler<TDomain, TAlgebra>::get_ip3(Vertex* const vrt, number time)
{
	// check whether vrt is newly active
	typename std::map<Vertex*, IP3Timing>::iterator it = m_mSynapseActivationTime.find(vrt);
	if (it == m_mSynapseActivationTime.end())
	{
		IP3Timing timing;
		timing.t_start = time;
		timing.t_end = time + m_j_ip3_duration;
		m_mSynapseActivationTime[vrt] = timing;

		return m_j_ip3_max;

	}

	// check whether v is still active
	if (time >= it->second.t_end)
	{
		// erase entry
		m_mSynapseActivationTime.erase(vrt);
		return 0.0;
	}

	return m_j_ip3_max * std::exp(m_j_ip3_decayRate*(it->second.t_start - time));
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


	SmartPtr<MultiGrid> mg = m_spDom->grid();
	typename TDomain::position_accessor_type& aaPos = m_spDom->position_accessor();
	SmartPtr<ISubsetHandler> sh = m_spDom->subset_handler();
	ConstSmartPtr<SurfaceView> sv = dd->surface_view();

	typedef typename domain_traits<TDomain::dim>::side_type side_type;

	typedef typename Grid::traits<side_type>::secure_container side_list_type;
	side_list_type sl;

	std::queue<side_type*> q;

	std::vector<side_type*> vSynElems;
	std::vector<number> vSynElemArea;
	std::vector<DoFIndex> vDoFIndex;

	const size_t nSyn = vActiveList.size();
	for (size_t s = 0; s < nSyn; ++s)
	{
		Vertex* synVrt = vActiveList[s];

		// For each synaptic vertex, find all (local) plasma membrane elements within a specific range.
		// The synaptic current will be distributed evenly among them.
		// (A more correct approach would be to find all GLOBAL elements, but this is more complicated
		// and probably not really necessary.)
		vSynElems.clear();
		vSynElemArea.clear();
		number totalSynArea = 0.0;

		const typename TDomain::position_type& synPos = aaPos[synVrt];

		// initialize queue with all neighboring elements (they at least need to be used)
		mg->begin_marking();
		mg->associated_elements(sl, synVrt);
		const size_t slSz = sl.size();
		for (size_t e = 0; e < slSz; ++e)
		{
			side_type* side = sl[e];
			if (std::find(m_vMembraneSI.begin(), m_vMembraneSI.end(), sh->get_subset_index(side))
				!= m_vMembraneSI.end())
			{
				q.push(side);
				mg->mark(side);
			}
		}

		// find all near neighbors
		while (!q.empty())
		{
			side_type* side = q.front();
			q.pop();

			// take care only to adjust the surface
			const size_t numChildren = mg->num_children<side_type>(side);
			if (numChildren)
			{
				for (size_t c = 0; c < numChildren; ++c)
				{
					side_type* childSide = mg->get_child<side_type>(side, c);

					if (!mg->is_marked(childSide)
						&& VecDistanceSq(synPos, CalculateCenter(childSide, aaPos)) < m_sqSynRadius)
					{
						q.push(childSide);
					}

					mg->mark(childSide);
				}
			}
			else if (sv->is_contained(side, GridLevel()))
			{
				vSynElems.push_back(side);
				vSynElemArea.push_back(CalculateVolume(side, aaPos));
				totalSynArea += vSynElemArea.back();
			}
			const size_t nSideVrt = side->num_vertices();
			for (size_t v = 0; v < nSideVrt; ++v)
			{
				Vertex* sideVrt = side->vertex(v);

				mg->associated_elements(sl, sideVrt);
				const size_t slSz = sl.size();
				for (size_t e = 0; e < slSz; ++e)
				{
					side_type* connSide = sl[e];
					if (!mg->is_marked(connSide)
						&& std::find(m_vMembraneSI.begin(), m_vMembraneSI.end(), sh->get_subset_index(connSide))
							!= m_vMembraneSI.end()
						&& VecDistanceSq(synPos, CalculateCenter(connSide, aaPos)) < m_sqSynRadius)
					{
						q.push(connSide);
					}

					mg->mark(connSide);
				}
			}
		}
		mg->end_marking();
		const size_t nSynElems = vSynElems.size();


		// if the potential rises high, synaptic currents are reversed;
		// we need to exclude calcium from this effect
		// TODO: this is a bit awkward, better use proper Ca2+ entry modeling
		if (vCurrent[s] < 0.0)
		{
			const number substanceCurrent = vCurrent[s] * m_current_percentage / (m_valency*m_F) / m_scaling_3d_to_1d_amount_of_substance;
			const number fluxDensity = substanceCurrent / totalSynArea;

			// loop all elems participating in that synapse
			for (size_t e = 0; e < nSynElems; ++e)
			{
				side_type* elem = vSynElems[e];
				const size_t nSynElemVrts = elem->num_vertices();

				// each node gets an equal part of this synElem's current
				const number currentPerNode = fluxDensity * vSynElemArea[e] / nSynElemVrts;
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
		if (!m_ip3_set) return;

		const number substanceCurrent = get_ip3(synVrt, time) / m_scaling_3d_to_1d_ip3;
		const number fluxDensity = substanceCurrent / totalSynArea;

		// loop all elems participating in that synapse
		for (size_t e = 0; e < nSynElems; ++e)
		{
			side_type* elem = vSynElems[e];
			const size_t nSynElemVrts = elem->num_vertices();

			// each node gets an equal part of this synElem's current
			const number currentPerNode = fluxDensity * vSynElemArea[e] / nSynElemVrts;
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


	SmartPtr<MultiGrid> mg = m_spDom->grid();
	typename TDomain::position_accessor_type& aaPos = m_spDom->position_accessor();
	SmartPtr<ISubsetHandler> sh = m_spDom->subset_handler();
	ConstSmartPtr<SurfaceView> sv = dd->surface_view();

	typedef typename domain_traits<TDomain::dim>::side_type side_type;

	typedef typename Grid::traits<side_type>::secure_container side_list_type;
	side_list_type sl;

	std::queue<side_type*> q;

	std::vector<side_type*> vSynElems;
	std::vector<number> vSynElemArea;
	std::vector<DoFIndex> vDoFIndex;

	const size_t nSyn = vActiveList.size();
	for (size_t s = 0; s < nSyn; ++s)
	{
		Vertex* synVrt = vActiveList[s];

		// For each synaptic vertex, find all (local) plasma membrane elements within a specific range.
		// The synaptic current will be distributed evenly among them.
		// (A more correct approach would be to find all GLOBAL elements, but this is more complicated
		// and probably not really necessary.)
		vSynElems.clear();
		vSynElemArea.clear();
		number totalSynArea = 0.0;

		const typename TDomain::position_type& synPos = aaPos[synVrt];

		// initialize queue with all neighboring elements (they at least need to be used)
		mg->begin_marking();
		mg->associated_elements(sl, synVrt);
		const size_t slSz = sl.size();
		for (size_t e = 0; e < slSz; ++e)
		{
			side_type* side = sl[e];
			if (std::find(m_vMembraneSI.begin(), m_vMembraneSI.end(), sh->get_subset_index(side))
				!= m_vMembraneSI.end())
			{
				q.push(side);
				mg->mark(side);
			}
		}

		// find all near neighbors
		while (!q.empty())
		{
			side_type* side = q.front();
			q.pop();

			// take care only to adjust the surface
			const size_t numChildren = mg->num_children<side_type>(side);
			if (numChildren)
			{
				for (size_t c = 0; c < numChildren; ++c)
				{
					side_type* childSide = mg->get_child<side_type>(side, c);

					if (!mg->is_marked(childSide)
						&& VecDistanceSq(synPos, CalculateCenter(childSide, aaPos)) < m_sqSynRadius)
					{
						q.push(childSide);
					}

					mg->mark(childSide);
				}
			}
			else if (sv->is_contained(side, GridLevel()))
			{
				vSynElems.push_back(side);
				vSynElemArea.push_back(CalculateVolume(side, aaPos));
				totalSynArea += vSynElemArea.back();
			}

			const size_t nSideVrt = side->num_vertices();
			for (size_t v = 0; v < nSideVrt; ++v)
			{
				Vertex* sideVrt = side->vertex(v);

				mg->associated_elements(sl, sideVrt);
				const size_t slSz = sl.size();
				for (size_t e = 0; e < slSz; ++e)
				{
					side_type* connSide = sl[e];
					if (!mg->is_marked(connSide)
						&& std::find(m_vMembraneSI.begin(), m_vMembraneSI.end(), sh->get_subset_index(connSide))
							!= m_vMembraneSI.end()
						&& VecDistanceSq(synPos, CalculateCenter(connSide, aaPos)) < m_sqSynRadius)
					{
						q.push(connSide);
					}

					mg->mark(connSide);
				}
			}
		}
		mg->end_marking();


		// if the potential rises high, synaptic currents are reversed;
		// we need to exclude calcium from this effect
		// TODO: this is a bit awkward, better use proper Ca2+ entry modeling
		if (vCurrent[s] < 0.0)
		{
			const number substanceCurrent = vCurrent[s] * m_current_percentage / (m_valency*m_F) / m_scaling_3d_to_1d_amount_of_substance;
			const number fluxDensity = substanceCurrent / totalSynArea;

			// loop all elems participating in that synapse
			const size_t nSynElems = vSynElems.size();
			for (size_t e = 0; e < nSynElems; ++e)
			{
				side_type* elem = vSynElems[e];

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
		if (!m_ip3_set) return;

		const number substanceCurrent = get_ip3(synVrt, time) / m_scaling_3d_to_1d_ip3;
		const number fluxDensity = substanceCurrent / totalSynArea;

		// loop all elems participating in that synapse
		const size_t nSynElems = vSynElems.size();
		for (size_t e = 0; e < nSynElems; ++e)
		{
			side_type* elem = vSynElems[e];

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
