/*
 * vdcc_bg_neuron.cpp
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#include "vdcc_bg_neuron.h"

namespace ug{
namespace neuro_collection{


template<typename TDomain>
VDCC_BG_VM2UG_NEURON<TDomain>::VDCC_BG_VM2UG_NEURON
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	SmartPtr<membrane_potential_mapping::Transformator> transformator,
	const std::string baseName,
	const char* timeFmt,
	const std::string ext,
	const bool posCanChange
)
: VDCC_BG<TDomain>(fcts, subsets, approx),
  m_NrnInterpreter(transformator), m_tFmt(timeFmt), m_vmTime(0.0)
{
	// nothing more to do
};

template<typename TDomain>
VDCC_BG_VM2UG_NEURON<TDomain>::VDCC_BG_VM2UG_NEURON
(
	const char* fcts,
	const char* subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	SmartPtr<membrane_potential_mapping::Transformator> transformator,
	const std::string baseName,
	const char* timeFmt,
	const std::string ext,
	const bool posCanChange
)
: VDCC_BG<TDomain>(fcts, subsets, approx),
  m_NrnInterpreter(transformator), m_tFmt(timeFmt), m_vmTime(0.0)
{
	// nothing more to do
};



template <typename TDomain>
VDCC_BG_VM2UG_NEURON<TDomain>::~VDCC_BG_VM2UG_NEURON()
{
	// nothing to do
}

template<typename TDomain>
void VDCC_BG_VM2UG_NEURON<TDomain>::init(number time)
{
	try
	{
		/// init to -65 mV (or to the finalize value supplied to the transformator instance)
		this->m_time = time;
		this->m_mapper.get()->get_transformator()->extract_vms(1, 0);

		/// with the extracted vms we build the tree then
		std::cout << "vms could be extracted" << std::endl;
		std::cout << "Our NEURON setup: " << std::endl;
		this->m_mapper.get()->get_transformator()->print_setup(true);
		// std::cout << "Our potential: " << this->m_vmProvider.get()->get_potential(0, 0, 0);
		this->m_mapper.get()->build_tree();
	}
	UG_CATCH_THROW("NEURON interpreter could not advance, vmProvider could not be created.");

	typedef typename DoFDistribution::traits<side_t>::const_iterator itType;
	SubsetGroup ssGrp;
	try {ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	const typename TDomain::position_accessor_type& aaPos = this->m_dom->position_accessor();
	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<side_t>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<side_t>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			/// retrieve membrane potential via vm2ug
			number vm;
			try
			{
				const typename TDomain::position_type& coords = CalculateCenter(*iter, aaPos);
				vm = this->m_mapper.get()->get_vm(coords[0], coords[1], coords[2]);
				std::cout << "VM: " << vm << std::endl;
				//vm = this->m_vmProvider.get()->get_data_from_nearest_neighbor(coords);
			}
			UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

			this->m_aaMGate[*iter] = this->calc_gating_start(this->m_gpMGate, vm);
			if (has_hGate()) this->m_aaHGate[*iter] = this->calc_gating_start(this->m_gpHGate, vm);
			this->m_aaVm[*iter] = 0.001 * vm; // mV -> V
		}
	}

	this->m_time = time;
	this->m_initiated = true;
}


template<typename TDomain>
void VDCC_BG_VM2UG_NEURON<TDomain>::update_time(number newTime) {
	/// only work if really necessary
	if (newTime == this->m_time) return;

	/// get dt from NEURON
	number dt = m_NrnInterpreter.get()->get_dt();

	/// set dt for NEURON explicit again
	std::stringstream ss;
	ss << "dt = " << dt;
	m_NrnInterpreter.get()->execute_hoc_stmt(ss.str());
	ss.str(""); ss.clear();

	/// advance the NEURON simulation by the UG dt, then update mapper
	/// with new membrane potential values for given time from NEURON
	size_t num_fadvances = (size_t) floor ((newTime - this->m_time) / dt);
	this->m_mapper.get()->get_transformator()->extract_vms(1, num_fadvances);

	/// update UG time
	this->m_oldTime = this->m_time;
	this->m_time = newTime;
}


template<typename TDomain>
void VDCC_BG_VM2UG_NEURON<TDomain>::update_potential(side_t* elem)
{
	// retrieve membrane potential via NEURON and the VMProvider/Mapper
	number vm;
	try
	{
		const typename TDomain::position_type& coords = CalculateCenter(elem, this->m_aaPos);
		vm = this->m_mapper.get()->get_vm(coords[0], coords[1], coords[2]);
		// vm = this->m_vmProvider.get()->get_data_from_nearest_neighbor(coords);
	}
	UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

	// set membrane potential value
	this->m_aaVm[elem] = 0.001 * vm; /// v -> mV
}



template<typename TDomain>
void VDCC_BG_VM2UG_NEURON<TDomain>::print_units() const
{
	std::string nm = this->name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_ext]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      V_m       mV                                                            |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);

	/// print NEURON setup
	this->m_mapper.get()->get_transformator()->print_setup(true);
	UG_LOG(std::endl);
}



// explicit template specializations
#ifdef UG_DIM_1
		template class VDCC_BG_VM2UG_NEURON<Domain1d>;
#endif
#ifdef UG_DIM_2
		template class VDCC_BG_VM2UG_NEURON<Domain2d>;
#endif
#ifdef UG_DIM_3
		template class VDCC_BG_VM2UG_NEURON<Domain3d>;
#endif

} // neuro_collection
} // namespace ug
