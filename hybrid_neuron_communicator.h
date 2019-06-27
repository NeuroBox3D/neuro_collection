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

#ifndef UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR_H
#define UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR_H

#include "../cable_neuron/synapse_handling/synapse_handler.h"
#include "../cable_neuron/synapse_handling/synapses/base_synapse.h"
#include "../cable_neuron/synapse_handling/synapses/post_synapse.h"
#include "../cable_neuron/synapse_handling/synapses/pre_synapse.h"


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
        void gather_synaptic_currents
		(
			std::vector<MathVector<dim> >& vActSynPosOut,
			std::vector<number>& vSynCurrOut,
			std::vector<synapse_id>& vSynIDOut,
			number time
		);

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


    protected:
        ///reinitialize mappings for 3d elem -> 1d vertex potential value mapping
        void reinit_potential_mapping();


        /// reinitialize mappings for 1d syn -> 3d vertex mapping
        void reinit_synapse_mapping();


		MessageHub::SPCallbackId 	m_spGridAdaptionCallbackID;
		void grid_adaption_callback(const GridMessage_Adaption& msg);

		MessageHub::SPCallbackId m_spGridDistributionCallbackID;
		void grid_distribution_callback(const GridMessage_Distribution& gmd);

    private:
        SmartPtr<synh_type> m_spSynHandler;

        /// memory for side element potential values
        std::map<side_t*, number> m_mElemPot;

        /// synapse to 3d coordinate vertex mapping
        std::map<synapse_id, MathVector<dim> > m_mSynapse3dCoords;

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

        bool m_bPotentialMappingNeedsUpdate;
        bool m_bSynapseMappingNeedsUpdate;
};


} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR_H
