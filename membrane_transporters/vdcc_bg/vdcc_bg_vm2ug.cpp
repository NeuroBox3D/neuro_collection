/*
 * vdcc_bg_vm2ug.cpp
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#include "vdcc_bg_vm2ug.h"

#include <locale>	// for control over the decimal separator (point instead of comma, please!)


namespace ug{
namespace neuro_collection{


///////////////////////////////////////////////////////////
///////////   BorgGrahamWithVM2UG   ///////////////////////
///////////////////////////////////////////////////////////

template <typename TDomain>
VDCC_BG_VM2UG<TDomain>::VDCC_BG_VM2UG
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const std::string baseName,
	const char* timeFmt,
	const std::string ext,
	const bool posCanChange
)
: VDCC_BG<TDomain>(fcts, subsets, approx),
 /// m_vmProvider(baseName, ext, !posCanChange), m_tFmt(timeFmt),
  m_vmProvider(), m_tFmt(timeFmt),
  m_fileInterval(0.0), m_fileOffset(0.0), m_baseName(baseName), m_ext(ext)
{
	// nothing further to do
}

template <typename TDomain>
VDCC_BG_VM2UG<TDomain>::VDCC_BG_VM2UG
(
	const char* fcts,
	const char* subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	const std::string baseName,
	const char* timeFmt,
	const std::string ext,
	const bool posCanChange
)
: VDCC_BG<TDomain>(fcts, subsets, approx),
 /// m_vmProvider(baseName, ext, !posCanChange), m_tFmt(timeFmt),
  m_vmProvider(), m_tFmt(timeFmt),
  m_fileInterval(0.0), m_fileOffset(0.0), m_baseName(baseName), m_ext(ext)
{
	// nothing further to do
}


template <typename TDomain>
VDCC_BG_VM2UG<TDomain>::~VDCC_BG_VM2UG()
{
	// nothing to do
}


template<typename TDomain>
void VDCC_BG_VM2UG<TDomain>::init(number time)
{
	// fill attachments with initial values

	// truncate time to last time that data exists for (if known)
	number vm_time;
	if (m_fileInterval >= 1e-9)
		vm_time = floor((time-m_fileOffset)/m_fileInterval)*m_fileInterval;
	else
		vm_time = time;

	try
	{
		char buffer[100];
		char* oldLocale = setlocale (LC_ALL, NULL);
		setlocale(LC_NUMERIC, "C");		// ensure decimal point is a . (not a ,)
		sprintf(buffer, m_tFmt.c_str(), vm_time);
		setlocale(LC_NUMERIC, oldLocale);
		m_timeAsString = buffer;
	}
	UG_CATCH_THROW("Time format string provided does not meet requirements.\n"
				<< "It must contain exactly one placeholder (which has to convert a floating point number)"
				<< "and must not produce an output string longer than 100 chars.");

	try {m_vmProvider.build_tree(m_baseName + m_timeAsString + m_ext, " ");}
	UG_CATCH_THROW("Underlying Vm2uG object could not build its tree on given file.\n"
		<< "If this is due to an inappropriate point in time, you might consider\n"
		"using set_file_times(fileInterval, fileOffset).");

	typedef typename DoFDistribution::traits<side_t>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	const typename TDomain::position_accessor_type& aaPos = this->m_dom->position_accessor();
	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<side_t>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<side_t>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// retrieve membrane potential via vm2ug
			number vm;
			try
			{
				const typename TDomain::position_type& coords = CalculateCenter(*iter, aaPos);
				vm = m_vmProvider.get_vm(coords);
			}
			UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");


			this->m_aaMGate[*iter] = this->calc_gating_start(this->m_gpMGate, vm);
			if (has_hGate()) this->m_aaHGate[*iter] = this->calc_gating_start(this->m_gpHGate, vm);
			this->m_aaVm[*iter] = 0.001 * vm;
		}
	}

	this->m_time = time;
	this->m_initiated = true;
}


template<typename TDomain>
void VDCC_BG_VM2UG<TDomain>::update_time(number newTime)
{
	// only work if really necessary
	if (newTime == this->m_time)
		return;

	// truncate time to last time that data exists for (if known)
	number vm_time;
	if (m_fileInterval >= 1e-9)
		vm_time = floor((newTime-m_fileOffset)/m_fileInterval)*m_fileInterval;
	else
		vm_time = newTime;

	try
	{
		char buffer[100];
		char* oldLocale = setlocale (LC_ALL, NULL);
		setlocale(LC_NUMERIC, "C");		// ensure decimal point is a . (not a ,)
		sprintf(buffer, m_tFmt.c_str(), vm_time);
		setlocale(LC_NUMERIC, oldLocale);
		m_timeAsString = buffer;
	}
	UG_CATCH_THROW("Time format string provided does not meet requirements.\n"
				<< "It must contain exactly one placeholder (which has to convert a floating point number)"
				<< "and must not produce an output string longer than 100 chars.");

	m_vmProvider.build_tree(m_baseName + m_timeAsString + m_ext, " ");
	this->m_oldTime = this->m_time;
	this->m_time = newTime;
}


template<typename TDomain>
void VDCC_BG_VM2UG<TDomain>::update_potential(side_t* elem)
{
	// retrieve membrane potential via vm2ug
	number vm;
	try
	{
		const typename TDomain::position_type& coords = CalculateCenter(elem, this->m_aaPos);
		vm = m_vmProvider.get_vm(coords);
	}
	UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

	// set membrane potential value
	this->m_aaVm[elem] = 0.001 * vm;
}

template<typename TDomain>
void VDCC_BG_VM2UG<TDomain>::print_units() const
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
}



// explicit template specializations
#ifdef UG_DIM_1
	template class VDCC_BG_VM2UG<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class VDCC_BG_VM2UG<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class VDCC_BG_VM2UG<Domain3d>;
#endif

} // neuro_collection
} // namespace ug
