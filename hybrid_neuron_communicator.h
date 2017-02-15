/*
 * hybrid_neuron_communicator.h
 *
 *  Created on: 20.12.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR
#define UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR

#include "../cable_neuron/cable_disc/cable_equation.h"


namespace ug {
namespace neuro_collection {


template <typename TDomain>
class HybridNeuronCommunicator
{
    protected:
        static const int dim = TDomain::dim;
        typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
        typedef typename elem_t::side side_t;
    	typedef typename TDomain::position_type posType;
    	typedef typename TDomain::position_accessor_type aaPos_type;

    	typedef typename cable_neuron::synapse_handler::SplitSynapseHandler<TDomain> synh_type;
    	typedef Attachment<int> ANeuronID;

    public:
        /// constructor
        HybridNeuronCommunicator
        (
            SmartPtr<ApproximationSpace<TDomain> > spApprox3d,
            SmartPtr<ApproximationSpace<TDomain> > spApprox1d
        );

        /// destructor
        ~HybridNeuronCommunicator();

        /// set cable equation object
        void set_ce_object(ConstSmartPtr<cable_neuron::CableEquation<TDomain> > spCEDisc);

        /// set subsets on which to communicate the potential values from 1d to 3d
        void set_potential_subsets(const std::vector<std::string>& vSubset);

        /// reinitialize communication
        /**
         * This method needs to be called whenever there are any changes in one of the geometries
         * or when one of the geometries has been redistributed.
         */
        void reinit_potential_mappings();

        /// communicate potential values
        void coordinate_potential_values();

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
    	void get_coordinates(SYNAPSE_ID id, MathVector<dim>& vCoords);
    	int get_neuron_id(SYNAPSE_ID id);
    	void prep_timestep(const number& t, const int& id, std::vector<number>& vCurr, std::vector<SYNAPSE_ID>& vSid);

    protected:
        void neuron_identification();
		int deep_first_search(Vertex* v, int id);
		int Mapping3d(int neuron_id, std::vector<Vertex*>& vMinimizing3dVertices, std::vector<MathVector<dim> >& vMinimizingSynapseCoords);

		/**
		 * Takes a std::vector of Synapses by coordinates and a std::vector of Vertexcoordinates and writes the nearest neighbor of each
		 * synapse out in vMap, so that the order of v1dSynapses equals the order of vMap
		 */
		int nearest_neighbor(	const std::vector<MathVector<dim> >& v1dCoords,
								const std::vector<Vertex*>& v3dVertices,
								std::vector<Vertex*>& vMap,
								std::vector<number>& vDistances);


    private:
        /// access to 1d cable discretization
        ConstSmartPtr<cable_neuron::CableEquation<TDomain> > m_spCE;

        SmartPtr<synh_type> m_spSynHandler;

        /// memory for side element potential values
        std::map<side_t*, number> m_mElemPot;

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
        SmartPtr<MultiGrid> m_spGrid1d;

        aaPos_type m_aaPos1d;
        aaPos_type m_aaPos3d;

    	ANeuronID m_aNID;
    	Grid::VertexAttachmentAccessor<ANeuronID> m_aaNID;

        std::vector<int> m_vPotSubset3d;
};


template <typename TDomain, typename TAlgebra>
class HybridSynapseCurrentAssembler : public IDomainConstraint<TDomain, TAlgebra>
{

private:
	ConstSmartPtr<hnc_type> m_spHNC;
public:

	typedef typename HybridNeuronCommunicator<TDomain> hnc_type;
	typedef TDomain domain_type;
	typedef TAlgebra algebra_type;
	typedef typename algebra_type::matrix_type matrix_type;
	typedef typename algebra_type::vector_type vector_type;




	HybridSynapseCurrentAssembler();
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

	int type() const;

	void set_flowing_substance_name(const std::string& fct); // todo: better name!

private:
	/// function index of the carried ion species
	size_t m_fctInd;

	/// Faraday constant
	const number m_F; // todo: set in constructor

	/// valency of the carried ion species
	number m_valency; // todo: implement a setter and default value for constructor!

	/// which fraction of the current is carried by the ion species in question
	number m_current_percentage; // todo: implement a setter and default value for constructor

};

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR
