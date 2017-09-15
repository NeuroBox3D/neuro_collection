/*
 * hybrid_synapse_current_assembler_impl.h
 *
 *  Created on: 20.12.2016
 *      Author: mbreit, lreinhardt
 */

#include "lib_grid/algorithms/debug_util.h" // ElementDebugInfo
#include "../cable_neuron/util/functions.h"	// neuron_identification

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
: m_fctInd(0), m_fctInd_ip3(0), m_ip3_set(true), m_F(96485.309), m_valency(2), m_current_percentage(0.1), m_spHNC(new hnc_type(spApprox3d, spApprox1d)),
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


	// loop active synapses
	size_t sz = vActiveList.size();
	for (size_t i = 0; i < sz; ++i)
	{
		Vertex* vrt = vActiveList[i];

		// get the DoFIndex for this vertex
		std::vector<DoFIndex> vIndex;
		dd->inner_dof_indices(vrt, m_fctInd, vIndex, false);

		UG_COND_THROW(!vIndex.size(), "Function given by 'set_flowing_substance_name' is not defined for "
			<< ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), vrt) << ".")

		UG_ASSERT(vIndex.size() == 1, "Apparently, you are using shape functions different from P1, this is not supported.");

		const DoFIndex& dofInd = vIndex[0];

		// Evaluate current (use the same current for all time points).
		// This way, we enforce explicit Euler scheme for current discretization.
		// The currents have to be scaled properly in case the units
		// of the 1d and 3d simulations do not match!
		// Add synaptic current * dt to defect.
		// Currents are outward in the synapse handler, so we _add_ to defect.

		// if the potential rises high, synaptic currents are reversed;
		// we need to exclude calcium from this effect
		// TODO: this is a bit awkward, better use proper Ca2+ entry modeling
		if (vCurrent[i] > 0.0)
			vCurrent[i] = 0.0;
		DoFRef(d, dofInd) += dt * vCurrent[i] * m_current_percentage / (m_valency*m_F) / m_scaling_3d_to_1d_amount_of_substance;

		if (!m_ip3_set) return;
		//get DoFIndex for this vertex (ip3)
		std::vector<DoFIndex> vIndex_ip3;
		dd->inner_dof_indices(vrt, m_fctInd_ip3, vIndex_ip3, false);

		UG_COND_THROW(!vIndex_ip3.size(), "Function for IP3 is not defined for "
			<< ElementDebugInfo(*this->m_spApproxSpace->domain()->grid(), vrt) << ".")

		UG_ASSERT(vIndex_ip3.size() == 1, "Apparently, you are using shape functions different from P1, this is not supported.");

		const DoFIndex& dofInd_ip3 = vIndex_ip3[0];

		// inward currents need to be subtracted!
		DoFRef(d, dofInd_ip3) -= get_ip3(vrt, time) * dt / m_scaling_3d_to_1d_ip3;
	}
}



} // namespace neuro_collection
} // namespace ug
