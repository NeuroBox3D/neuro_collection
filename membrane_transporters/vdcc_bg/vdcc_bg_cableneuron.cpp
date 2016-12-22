/*
 * vdcc_bg_cableneuron.cpp
 *
 *  Created on: 20.12.2016
 *      Author: mbreit
 */

#include "vdcc_bg_cableneuron.h"

namespace ug {
namespace neuro_collection {


template <typename TDomain>
VDCC_BG_CN<TDomain>::VDCC_BG_CN
(
    const std::vector<std::string>& fcts,
    const std::vector<std::string>& subsets,
    SmartPtr<ApproximationSpace<TDomain> > approx
)
: VDCC_BG<TDomain>(fcts, subsets, approx) {}


template <typename TDomain>
VDCC_BG_CN<TDomain>::~VDCC_BG_CN()
{}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::set_hybrid_neuron_communicator(SmartPtr<HybridNeuronCommunicator<TDomain> > spHNC)
{
    m_spHNC = spHNC;
    m_spHNC->set_potential_subsets(this->m_vSubset);
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::prep_timestep(const number time, VectorProxyBase* upb)
{
    // first: communicate current membrane potential values
    m_spHNC->communicate_potential_values();

    // call regular prep timestep from base class
    VDCC_BG<TDomain>::prep_timestep(time, upb);
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::update_potential(side_t* elem)
{
    // fill attachments with renewed values
    this->m_aaVm[elem] = m_spHNC->potential(elem);
}


template <typename TDomain>
const std::string VDCC_BG_CN<TDomain>::name() const
{
    return std::string("VDCC_BG_CN");
}


template <typename TDomain>
void VDCC_BG_CN<TDomain>::print_units() const
{
    std::string nm = name();
    size_t n = nm.size();
    UG_LOG(std::endl);
    UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
    UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
    UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
    UG_LOG("|    Input                                                                     |"<< std::endl);
    UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
    UG_LOG("|      [Ca_ext]  mM (= mol/m^3)                                                |"<< std::endl);
    UG_LOG("|      V_m       V                                                             |"<< std::endl);
    UG_LOG("|                                                                              |"<< std::endl);
    UG_LOG("|    Output                                                                    |"<< std::endl);
    UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
    UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
    UG_LOG(std::endl);
}


// explicit template specializations
#ifdef UG_DIM_1
    template class VDCC_BG_CN<Domain1d>;
#endif
#ifdef UG_DIM_2
    template class VDCC_BG_CN<Domain2d>;
#endif
#ifdef UG_DIM_3
    template class VDCC_BG_CN<Domain3d>;
#endif


} // namespace neuro_collection
} // namespace ug
