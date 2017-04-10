/*
 * hybrid_neuron_communicator.h
 *
 *  Created on: 20.12.2016
 *      Author: mbreit, lreinhardt
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR
#define UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR

#include "../cable_neuron/synapse_handling/synapse_handler.h"
#include "../cable_neuron/synapse_handling/synapses/base_synapse.h"
#include "../cable_neuron/synapse_handling/synapses/post_synapse.h"
#include "../cable_neuron/synapse_handling/synapses/pre_synapse.h"


namespace ug {
namespace neuro_collection {


// TODO: register this object as refinement listener
//       and as redistribution listener
//       it has to react to both events with a call to reinit_potential_mappings()
//       as well as Mapping3d()

template <typename TDomain>
class HybridNeuronCommunicator
{
    protected:
        static const int dim = TDomain::dim;
        typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
        typedef typename elem_t::side side_t;
    	typedef typename TDomain::position_type posType;
    	typedef typename TDomain::position_accessor_type aaPos_type;
        typedef CPUAlgebra algebra_t;

    	typedef typename cable_neuron::synapse_handler::SynapseHandler<TDomain> synh_type;
    	typedef Attachment<uint> ANeuronID;

    public:
        /// constructor
        HybridNeuronCommunicator
        (
            SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
            SmartPtr<ApproximationSpace<TDomain> > spApprox1d
        );

        /// destructor
        ~HybridNeuronCommunicator();

        /// set synapse handler object
        void set_synapse_handler(SmartPtr<cable_neuron::synapse_handler::SynapseHandler<TDomain> > spSH);

        /// set subsets on which to communicate the potential values from 1d to 3d
        void set_potential_subsets(const std::vector<std::string>& vSubset);

        /// set subsets on which to communicate the current values from 1d to 3d
        void set_current_subsets(const std::vector<std::string>& vSubset);

        /// set current solution vector and function index of potential
        void set_solution_and_potential_index(ConstSmartPtr<GridFunction<TDomain, algebra_t> > u, size_t fctInd);

        /**
         * @brief Set scale factor from 3d to 1d.
         * This is required as 3d and 1d morphology might be differently scaled.
         * For example, 1d geometries used in the cable_neuron code should be scaled to length unit [m],
         * while 3d simulation length scale is often [um]. In that case, we need a scaling factor of 1e-6.
         */
        void set_coordinate_scale_factor_3d_to_1d(number scale);

        void set_neuron_ids(const std::vector<uint>& vNid);

        ConstSmartPtr<synh_type> synapse_handler() const {return m_spSynHandler;}

        /// communicate potential values
        void coordinate_potential_values();

        /// coordinate synaptic current values (only active post-synapses)
        void gather_synaptic_currents(std::vector<Vertex*>& vActSynOut, std::vector<number>& vSynCurrOut, number time);

        /// get potential value for high-dim side element
        number potential(side_t* elem) const;

    private:
        struct VecIndexCompare
        {
            VecIndexCompare(const std::vector<MathVector<dim> >& _vec, size_t _cmp)
            : vec(_vec), cmp(_cmp) {};

            bool operator()(const size_t& a, const size_t& b)
            {return vec[a][cmp] < vec[b][cmp];}

            private:
                const std::vector<MathVector<dim> >& vec;
                size_t cmp;
        };

        int nearest_neighbor_search
		(
			const std::vector<posType>& queryPts,
			const std::vector<posType>& dataPts,
			std::vector<size_t>& vNNout,
			std::vector<typename posType::value_type>& vDistOut
		) const;

    public:
    	/**
    	 * Calculates real coordinates of given Synapse (by ID) and returns in vCoords
    	 * If SynapseId not found vCoords is empty.
    	 */
    	void get_postsyn_coordinates(synapse_id id, MathVector<dim>& vCoords);
    	uint get_postsyn_neuron_id(synapse_id id);

    	const std::map<synapse_id, Vertex*>& synapse_3dVertex_map() const {return m_mSynapse3dVertex;}


    protected:
    	/**
    	 * reinitialize internal mappings
         * This method needs to be called whenever there are any changes in one of the geometries
         * or when one of the geometries has been redistributed.
        **/
    	void reinit();

        ///reinitialize mappings for 3d elem -> 1d vertex potential value mapping
        void reinit_potential_mapping();


        /// reinitialize mappings for 1d syn -> 3d vertex mapping
        void reinit_synapse_mapping();

    private:
        SmartPtr<synh_type> m_spSynHandler;

        /// memory for side element potential values
        std::map<side_t*, number> m_mElemPot;

        //Synapse to 3d-Vertex mapping
        std::map<synapse_id, Vertex*> m_mSynapse3dVertex;

#ifdef UG_PARALLEL
        /// list of 1d sender vertices on this proc and who they send to
        std::map<int, std::vector<Vertex*> > m_mSendInfo;

        /// list of 3d receiver elems and who they receive from
        std::map<int, std::vector<side_t*> > m_mReceiveInfo;

        int* rcvSize;
        int* rcvFrom;
        void* rcvBuf;

        int* sendSize;
        int* sendTo;
        void* sendBuf;
#endif

        std::map<side_t*, Vertex*> m_mPotElemToVertex;  // direct map (only for the serial case)

        SmartPtr<ApproximationSpace<TDomain> > m_spApprox1d;
        SmartPtr<ApproximationSpace<TDomain> > m_spApprox3d;

        ConstSmartPtr<GridFunction<TDomain, algebra_t> > m_spU;
        size_t m_potFctInd;

        SmartPtr<MultiGrid> m_spGrid1d;
        SmartPtr<MultiGrid> m_spGrid3d;
        SmartPtr<MultiGridSubsetHandler> m_spMGSSH3d;


        number m_scale_factor_from_3d_to_1d;

        aaPos_type m_aaPos1d;
        aaPos_type m_aaPos3d;

    	ANeuronID m_aNID;
    	Grid::VertexAttachmentAccessor<ANeuronID> m_aaNID;

    	std::vector<uint> m_vNid;

        std::vector<int> m_vPotSubset3d;
        std::vector<int> m_vCurrentSubset3d;

        bool m_bInited;
};


template <typename TDomain, typename TAlgebra>
class HybridSynapseCurrentAssembler : public IDomainConstraint<TDomain, TAlgebra>
{

public:
    typedef HybridNeuronCommunicator<TDomain> hnc_type;
    typedef TDomain domain_type;
    typedef TAlgebra algebra_type;
    typedef typename algebra_type::matrix_type matrix_type;
    typedef typename algebra_type::vector_type vector_type;


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

	/**
	 * change scaling factors, that have to be applied to 3d values, so that 1d and 3d values
	 * have equal units
	 */
	void set_scaling_factors(number scaling_3d_to_1d_amount_of_substance = 1.0,
			 	 	 	 	 number scaling_3d_to_1d_coordinates = 1.0,
							 number scaling_3d_to_1d_electric_charge = 1.0)
	{
		m_scaling_3d_to_1d_amount_of_substance = scaling_3d_to_1d_amount_of_substance;
		m_scaling_3d_to_1d_electric_charge = scaling_3d_to_1d_electric_charge;
		m_scaling_3d_to_1d_coordinates = scaling_3d_to_1d_coordinates;

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

private:
	/// function index of the carried ion species
	size_t m_fctInd;

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
};

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR
