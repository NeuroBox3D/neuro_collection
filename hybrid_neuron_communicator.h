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


    private:
        /// access to 1d cable discretization
        ConstSmartPtr<cable_neuron::CableEquation<TDomain> > m_spCE;

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
        std::vector<int> m_vPotSubset3d;
};

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__HYBRID_NEURON_COMMUNICATOR
