/*
 * hybrid_synapse_current_assembler.h
 *
 *  Created on: 20.12.2016
 *      Author: mbreit, lreinhardt
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__HYBRID_SYNAPSE_CURRENT_ASSEMBLER
#define UG__PLUGINS__NEURO_COLLECTION__HYBRID_SYNAPSE_CURRENT_ASSEMBLER

#include "../cable_neuron/synapse_handling/synapse_handler.h"
#include "../cable_neuron/synapse_handling/synapses/base_synapse.h"
#include "../cable_neuron/synapse_handling/synapses/post_synapse.h"
#include "../cable_neuron/synapse_handling/synapses/pre_synapse.h"
#include "hybrid_neuron_communicator.h"

namespace ug {
namespace neuro_collection {


template <typename TDomain, typename TAlgebra>
class HybridSynapseCurrentAssembler : public IDomainConstraint<TDomain, TAlgebra>
{

public:
    typedef HybridNeuronCommunicator<TDomain> hnc_type;
    typedef TDomain domain_type;
    typedef TAlgebra algebra_type;
    typedef typename algebra_type::matrix_type matrix_type;
    typedef typename algebra_type::vector_type vector_type;


    /// construcor with IP3
	HybridSynapseCurrentAssembler(
		SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
		SmartPtr<ApproximationSpace<TDomain> > spApprox1d,
		SmartPtr<cable_neuron::synapse_handler::SynapseHandler<TDomain> > spSH,
		const std::vector<std::string>& PlasmaMembraneSubsetName,
		const std::string& fct,
		const std::string& fct_ip3
	);

	/// constructor without IP3
	HybridSynapseCurrentAssembler(
		SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
		SmartPtr<ApproximationSpace<TDomain> > spApprox1d,
		SmartPtr<cable_neuron::synapse_handler::SynapseHandler<TDomain> > spSH,
		const std::vector<std::string>& PlasmaMembraneSubsetName,
		const std::string& fct
	);

	virtual ~HybridSynapseCurrentAssembler(){}

	void adjust_jacobian(matrix_type& J, const vector_type& u,
			                             ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
			                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
										 const number s_a0 = 1.0) {}

	void adjust_defect(vector_type& d, const vector_type& u,
									   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
									   ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
									   const std::vector<number>* vScaleMass = NULL,
									   const std::vector<number>* vScaleStiff = NULL);

	void adjust_linear(matrix_type& mat, vector_type& rhs,
			                           ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0){}

	void adjust_rhs(vector_type& rhs, const vector_type& u,
			                        ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0){}

	void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd, int type,
										 number time = 0.0){}

	int type() const {return CT_CONSTRAINTS;};

	void set_valency(int val) {m_valency = val;}

	void set_current_percentage(number val) {m_current_percentage = val;}

	void set_ip3_production_params(number j_max, number decayRate)
	{
		m_j_ip3_max = j_max;
		m_j_ip3_decayRate = decayRate;

		// We "turn off" the IP3 production when 95% of the total amount have been produced.
		// This amounts to a production time of log(1/(1-0.95)) / decayRate and as log(20)
		// is 3.0 pretty exactly, we simply choose 3.0. :)
		m_j_ip3_duration = 3.0 / m_j_ip3_decayRate;
	}
	//void set_j_ip3_max(number val) {m_j_ip3_max = val; m_J_ip3_max = m_j_ip3_max/m_j_ip3_decayRate;}
	//void set_j_ip3_decayRate(number val) {m_j_ip3_decayRate = val; m_J_ip3_max = m_j_ip3_max/m_j_ip3_decayRate;}

	/**
	 * change scaling factors, that have to be applied to 3d values, so that 1d and 3d values
	 * have equal units
	 */
	void set_scaling_factors(number scaling_3d_to_1d_amount_of_substance = 1.0,
			 	 	 	 	 number scaling_3d_to_1d_coordinates = 1.0,
							 number scaling_3d_to_1d_electric_charge = 1.0,
							 number scaling_3d_to_1d_ip3 = 1.0)
	{
		m_scaling_3d_to_1d_amount_of_substance = scaling_3d_to_1d_amount_of_substance;
		m_scaling_3d_to_1d_electric_charge = scaling_3d_to_1d_electric_charge;
		m_scaling_3d_to_1d_coordinates = scaling_3d_to_1d_coordinates;
		m_scaling_3d_to_1d_ip3 = scaling_3d_to_1d_ip3;

		m_spHNC->set_coordinate_scale_factor_3d_to_1d(m_scaling_3d_to_1d_coordinates);
	}

	void set_3d_neuron_ids(const std::vector<size_t>& ids)
	{
		// uint is not registered, we therefore use size_t as param type
		// but we need to convert this here
		std::vector<uint> vID(ids.size());
		for (size_t i = 0; i < ids.size(); ++i)
			vID[i] = (uint) ids[i];

		m_spHNC->set_neuron_ids(vID);
	}

	//void set_ip3_duration(const number& dur) {m_j_ip3_duration = dur;}

protected:
	number get_ip3(Vertex* const vrt, number time);

private:
	struct IP3Timing {
		number t_start;
		number t_end;
	};

	/// function index of the carried ion species
	size_t m_fctInd;
	size_t m_fctInd_ip3;
	bool m_ip3_set;

	/// Faraday constant
	const number m_F; //in C/mol

	/// valency of the carried ion species
	int m_valency;

	/// which fraction of the current is carried by the ion species in question
	number m_current_percentage;

	SmartPtr<hnc_type> m_spHNC;

	/// scaling factors
	number m_scaling_3d_to_1d_amount_of_substance;
	number m_scaling_3d_to_1d_electric_charge;
	number m_scaling_3d_to_1d_coordinates;
	number m_scaling_3d_to_1d_ip3;					//1d units: (mol um)/(dm^3 s); scaling: u = 1e-6

	// IP3-related
	number m_j_ip3_max; // maximal ip3 "current" (in  mol/s) (adapted from Fink et al.)
	number m_j_ip3_decayRate; //in 1/s
	number m_j_ip3_duration; //duration of ip3 influx until specified ip3 fraction of total is reached (in s)

	std::map<Vertex*, IP3Timing> m_mSynapseActivationTime; //maps a 3d vertex mapped synapse to their activation time for IP3 generation
};

} // namespace neuro_collection
} // namespace ug

#include "hybrid_synapse_current_assembler_impl.h"

#endif // UG__PLUGINS__NEURO_COLLECTION__HYBRID_SYNAPSE_CURRENT_ASSEMBLER
