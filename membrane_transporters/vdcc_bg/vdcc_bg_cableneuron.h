/*
 * vdcc_bg_cableneuron.h
 *
 *  Created on: 20.12.2016
 *      Author: markus
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_CN__
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_CN__

#include "vdcc_bg.h"
#include "../../hybrid_neuron_communicator.h"

namespace ug {
namespace neuro_collection {


/// @addtogroup plugin_neuro_collection
/// @{

/// Borg Graham type VGCCs with membrane potential supply by UG's cable_neuron project.
/** This class is a specialization of the Borg-Graham interface.
 *  It supplies the channel with the necessary membrane potential values by an object of
 *  cable_neuron::CableEquation object discretizing a 1d representation of a 3d geometry.
**/

template <typename TDomain>
class VDCC_BG_CN
: public VDCC_BG<TDomain>
{
    protected:
        typedef typename GeomObjBaseTypeByDim<TDomain::dim>::base_obj_type elem_t;
        typedef typename elem_t::side side_t;

    public:
        /// constructor
        VDCC_BG_CN
        (
            const std::vector<std::string>& fcts,
            const std::vector<std::string>& subsets,
            SmartPtr<ApproximationSpace<TDomain> > approx
        );

        /// destructor
        virtual ~VDCC_BG_CN();

        /// setting a communicator object for hybrid neuron treatment
        void set_hybrid_neuron_communicator(SmartPtr<HybridNeuronCommunicator<TDomain> > spHNC);

        /// @copydoc IMembraneTransporter::prep_timestep()
        virtual void prep_timestep(const number time, VectorProxyBase* upb);

        /// @copydoc VDCC_BG::update_potential()
        virtual void update_potential(side_t* elem);

        /// @copydoc IMembraneTransporter::name()
        virtual const std::string name() const;

        /// @copydoc IMembraneTransporter::print_units()
        virtual void print_units() const;

    private:
        SmartPtr<HybridNeuronCommunicator<TDomain> > m_spHNC;
};

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_CN__
