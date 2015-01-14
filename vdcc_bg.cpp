/*
 * vdcc_bg.cpp
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#include "vdcc_bg.h"

namespace ug{
namespace neuro_collection{


///////////////////////////////////////////////////////////
///////////   BorgGraham   ////////////////////////////////
///////////////////////////////////////////////////////////
template<typename TDomain>
VDCC_BG<TDomain>::VDCC_BG
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: IMembraneTransporter(fcts),
  R(8.314), T(310.0), F(96485.0),
  m_dom(approx->domain()), m_mg(m_dom->grid()), m_dd(approx->dof_distribution(GridLevel::TOP)),
  m_sh(m_dom->subset_handler()), m_aaPos(m_dom->position_accessor()), m_vSubset(subsets),
  m_gpMGate(3.4, -21.0, 1.5), m_gpHGate(-2.0, -40.0, 75.0), m_time(0.0), m_oldTime(0.0),
  m_perm(2.4e-19), m_mp(2), m_hp(1), m_channelType(BG_Ntype), m_initiated(false)
{
	// process subsets

	//	remove white space
	for(size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);

	//	if no subset passed, clear subsets
	if (m_vSubset.size() == 1 && m_vSubset[0].empty())
		m_vSubset.clear();

	//	if subsets passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < m_vSubset.size(); ++i)
	{
		if (m_vSubset.empty())
		{
			UG_THROW("Error while setting subsets in an ElemDisc: passed "
			 		 "subset string lacks a subset specification at position "
					 << i << "(of " << m_vSubset.size()-1 << ")");
		}
	}

	// attach attachments
	if (m_mg->template has_attachment<side_t>(this->m_MGate))
		UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_MGate);

	if (has_hGate())
	{
		if (m_mg->template has_attachment<side_t>(this->m_HGate))
			UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
					 "could not be made, since it already exists.");
		m_mg->template attach_to<side_t>(this->m_HGate);
	}

	if (m_mg->template has_attachment<side_t>(this->m_Vm))
		UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
				 "could not be made, since it already exists.");
	m_mg->template attach_to<side_t>(this->m_Vm);


	// create attachment accessors
	m_aaMGate = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_MGate);
	if (has_hGate()) m_aaHGate = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_HGate);
	m_aaVm = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_Vm);
}


template<typename TDomain>
VDCC_BG<TDomain>::~VDCC_BG()
{
	m_mg->template detach_from<side_t>(this->m_MGate);
	m_mg->template detach_from<side_t>(this->m_HGate);
	m_mg->template detach_from<side_t>(this->m_Vm);
}


template<typename TDomain>
number VDCC_BG<TDomain>::calc_gating_start(GatingParams& gp, number Vm)
{
	return 1.0 / (1.0 + exp(-gp.z * (Vm - gp.V_12)/1000.0 * F/(R*T)));
}


template<typename TDomain>
void VDCC_BG<TDomain>::calc_gating_step(GatingParams& gp, number Vm, number dt, number& currVal)
{
	// forward step: implicit
	if (dt>=0) currVal = (currVal + dt/gp.tau_0 * calc_gating_start(gp,Vm)) / (1.0 + dt/gp.tau_0);
	// backward step: explicit
	else currVal += dt/gp.tau_0 * (calc_gating_start(gp,Vm) - currVal);
}


template<typename TDomain>
void VDCC_BG<TDomain>::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	side_t* elem = dynamic_cast<side_t*>(e);
	if (!elem) UG_THROW("OneSidedBorgGrahamFV1 fluxDensityFunction called with the wrong type of element.");

	number gating = pow(m_aaMGate[elem], m_mp);
	if (has_hGate()) gating *= pow(m_aaHGate[elem], m_hp);

	// flux derived from Goldman-Hodgkin-Katz equation,
	number maxFlux;
	number vm = m_aaVm[elem];
	number caCyt = u[_CCYT_];		// cytosolic Ca2+ concentration
	number caExt = u[_CEXT_];		// extracellular Ca2+ concentration

	// near V_m == 0: approximate by first order Taylor to avoid relative errors and div-by-0
	if (fabs(vm) < 1e-8) maxFlux = m_perm * ((caExt - caCyt) - F/(R*T) * (caExt + caCyt)*vm);
	else maxFlux = -m_perm * 2*F/(R*T) * vm * (caExt - caCyt*exp(2*F/(R*T)*vm)) / (1.0 - exp(2*F/(R*T)*vm));

	flux[0] = gating * maxFlux;
}


template<typename TDomain>
void VDCC_BG<TDomain>::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	side_t* elem = dynamic_cast<side_t*>(e);
	if (!elem) UG_THROW("OneSidedBorgGrahamFV1 fluxDensityFunction called with the wrong type of element.");

	number gating = pow(m_aaMGate[elem], m_mp);
	if (has_hGate()) gating *= pow(m_aaHGate[elem], m_hp);

	number dMaxFlux_dCyt, dMaxFlux_dExt;
	number vm = m_aaVm[elem];

	// near V_m == 0: approximate by first order Taylor to avoid relative errors and div-by-0
	if (fabs(vm) < 1e-8)
	{
		dMaxFlux_dCyt =  m_perm * (-1.0 - F/(R*T) * vm);
		dMaxFlux_dExt =  m_perm * (1.0 - F/(R*T) * vm);
	}
	else
	{
		dMaxFlux_dCyt = m_perm * 2*F/(R*T) * vm / (exp(-2*F/(R*T)*vm) - 1.0);
		dMaxFlux_dExt = m_perm * 2*F/(R*T) * vm / (exp(2*F/(R*T)*vm) - 1.0);
	}

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = gating * dMaxFlux_dCyt;
		i++;
	}
	if (!has_constant_value(_CEXT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CEXT_);
		flux_derivs[0][i].second = gating * dMaxFlux_dExt;
	}
}


template<typename TDomain>
const size_t VDCC_BG<TDomain>::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_CCYT_)) n--;
	if (has_constant_value(_CEXT_)) n--;

	return n;
}


template<typename TDomain>
size_t VDCC_BG<TDomain>::n_fluxes() const
{
	return 1;
}


template<typename TDomain>
const std::pair<size_t,size_t> VDCC_BG<TDomain>::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CEXT_)) from = local_fct_index(_CEXT_); else from = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


template<typename TDomain>
const std::string VDCC_BG<TDomain>::name() const
{
	return std::string("VDCC_BG");
}


template<typename TDomain>
void VDCC_BG<TDomain>::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CEXT_))
	{
		UG_THROW("Supplying neither cytosolic nor extracellular calcium concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}
}


template<typename TDomain>
void VDCC_BG<TDomain>::print_units() const
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


template<typename TDomain>
void VDCC_BG<TDomain>::set_permeability(const number perm)
{
	m_perm = perm;
}


template<typename TDomain>
void VDCC_BG<TDomain>::init(number time)
{
	this->m_time = time;

	typedef typename DoFDistribution::traits<side_t>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = m_dd->template begin<side_t>(ssGrp[si]);
		itType iterEnd = m_dd->template end<side_t>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// get potential data
			update_potential(*iter);
			number vm = m_aaVm[*iter];

			// calculate corresponding start condition for gates
			m_aaMGate[*iter] = calc_gating_start(m_gpMGate, vm);
			if (has_hGate()) m_aaHGate[*iter] = calc_gating_start(m_gpHGate, vm);
		}
	}

	this->m_initiated = true;
}


template<typename TDomain>
void VDCC_BG<TDomain>::update_gating(side_t* elem)
{
	if (!this->m_initiated)
		UG_THROW("Borg-Graham not initialized.\n"
				  << "Do not forget to do so before any updates by calling init(initTime).");

	// set new gating particle values
	number dt = 1000.0*(m_time - m_oldTime);	// calculating in ms
	calc_gating_step(m_gpMGate, 1000.0*m_aaVm[elem], dt, m_aaMGate[elem]);
	if (has_hGate()) calc_gating_step(m_gpHGate, 1000.0*m_aaVm[elem], dt, m_aaHGate[elem]);
}

template<typename TDomain>
void VDCC_BG<TDomain>::prep_timestep_elem
(
	const number time,
	const LocalVector& u,
	GridObject* elem
)
{
	update_time(time);

	side_t* side = dynamic_cast<side_t*>(elem);
	if (!side) UG_THROW("OneSidedBorgGrahamFV1::prep_timestep_elem() called with improper element type.");

	update_potential(side);
	update_gating(side);
}


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
  m_vmProvider(baseName, ext, !posCanChange), m_tFmt(timeFmt),
  m_fileInterval(0.0), m_fileOffset(0.0)
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

	std::string timeAsString;
	try
	{
		char buffer[100];
		char* oldLocale = setlocale (LC_ALL, NULL);
		setlocale(LC_NUMERIC, "C");		// ensure decimal point is a . (not a ,)
		sprintf(buffer, m_tFmt.c_str(), vm_time);
		setlocale(LC_NUMERIC, oldLocale);
		timeAsString = buffer;
	}
	UG_CATCH_THROW("Time format string provided does not meet requirements.\n"
				<< "It must contain exactly one placeholder (which has to convert a floating point number)"
				<< "and must not produce an output string longer than 100 chars.");

	try {m_vmProvider.buildTree(timeAsString);}
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
				if (coords.size() == 2)
					vm = m_vmProvider.get_potential(coords[0], coords[1], 0.0, timeAsString);
				else if (coords.size() == 3)
					vm = m_vmProvider.get_potential(coords[0], coords[1], coords[2], timeAsString);
				else
					UG_THROW("Coordinates of vertex are not 2d or 3d"
						  << "(which are the only valid dimensionalities for this discretization).");
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

	if (!m_vmProvider.treeBuild())
	UG_THROW("Underlying Vm2uG object's tree is not yet built.\n"
		  << "Do not forget to initialize the Borg-Graham object first by calling init(initTime).");

	// set new timestep file in vm2ug object (not strictly necessary)
	const std::string ts = m_timeAsString;
	m_vmProvider.setTimestep(ts);

	this->m_oldTime = this->m_time;
	this->m_time = newTime;
}


template<typename TDomain>
void VDCC_BG_VM2UG<TDomain>::update_potential(side_t* elem)
{
	if (!m_vmProvider.treeBuild())
		UG_THROW("Underlying Vm2uG object's tree is not yet built.\n"
			  << "Do not forget to initialize the Borg-Graham object first by calling init(initTime).");

	// retrieve membrane potential via vm2ug
	number vm;
	try
	{
		const typename TDomain::position_type& coords = CalculateCenter(elem, this->m_aaPos);
		if (coords.size() == 2)
			vm = m_vmProvider.get_potential(coords[0], coords[1], 0.0, m_timeAsString);
		else if (coords.size() == 3)
			vm = m_vmProvider.get_potential(coords[0], coords[1], coords[2], m_timeAsString);
		else
			UG_THROW("Coordinates of vertex are not 2d or 3d"
				  << "(which are the only valid dimensionalities for this discretization).");
	}
	UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

	// set membrane potential value
	this->m_aaVm[elem] = 0.001 * vm;
}



///////////////////////////////////////////////////////////
///////////   BorgGrahamWithNEURON   ///////////////////////
///////////////////////////////////////////////////////////
#ifdef MPMNEURON

template<typename TDomain>
VDCC_BG_VM2UG_NEURON<TDomain>::VDCC_BG_VM2UG_NEURON
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx,
	SmartPtr<Transformator> transformator,
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
		// make zero steps (e. g. init) and extract the membrane potentials
		this->m_NrnInterpreter.get()->extract_vms(1, 3);
		// with the extracted vms we build the tree then
		std::cout << "vms could be extracted" << std::endl;
		std::cout << "Our NEURON setup: " << std::endl;
		this->m_NrnInterpreter.get()->print_setup(true);
		std::cout << "Our potential: " << this->m_vmProvider.get()->get_potential(0, 0, 0);
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
			// retrieve membrane potential via vm2ug
			number vm;
			try
			{
				const typename TDomain::position_type& coords = CalculateCenter(*iter, aaPos);
				if (coords.size() == 2)
					vm = this->m_vmProvider.get()->get_potential(coords[0], coords[1], 0.0);
				else if (coords.size() == 3)
					vm = this->m_vmProvider.get()->get_potential(coords[0], coords[1], coords[2]);
				else
					UG_THROW("Coordinates of vertex are not 2d or 3d"
						  << "(which are the only valid dimensionalities for this discretization).");
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
void VDCC_BG_VM2UG_NEURON<TDomain>::update_time(number newTime)
{
	/// set dt for neuron interpreter
	number dt = this->m_time - this->m_oldTime;
	this->m_oldTime = this->m_time;
	this->m_time = newTime;
	std::stringstream ss;
	ss << "dt = " << dt;
	m_NrnInterpreter.get()->execute_hoc_stmt(ss.str());
	ss.str(""); ss.clear();

	// only work if really necessary
	if (newTime == this->m_time) {
		return;
	}

	// advance by dt (set above)
	m_NrnInterpreter.get()->fadvance();
	this->m_oldTime = this->m_time;
	this->m_time = newTime;
}


template<typename TDomain>
void VDCC_BG_VM2UG_NEURON<TDomain>::update_potential(side_t* elem)
{
	if (!this->m_vmProvider.get()->treeBuild())
	UG_THROW("Underlying Vm2uG object's tree is not yet built.\n"
		  << "Do not forget to initialize the Borg-Graham object first by calling init(initTime).");

	// retrieve membrane potential via vm2ug
	number vm;
	try
	{
		const typename TDomain::position_type& coords = CalculateCenter(elem, this->m_aaPos);
		if (coords.size() == 2)
			vm = this->m_vmProvider.get()->get_potential(coords[0], coords[1], 0.0);
		else if (coords.size() == 3)
			vm = this->m_vmProvider.get()->get_potential(coords[0], coords[1], coords[2]);
		else
		{
			UG_THROW("Coordinates of vertex are not 2d or 3d"
				  << "(which are the only valid dimensionalities for this discretization).");
		}
	}
	UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

	// set membrane potential value
	this->m_aaVm[elem] = 0.001 * vm;
}

#endif




///////////////////////////////////////////////////////////
///////////   BorgGrahamWithUserData   ////////////////////
///////////////////////////////////////////////////////////
template <typename TDomain>
VDCC_BG_UserData<TDomain>::VDCC_BG_UserData
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	SmartPtr<ApproximationSpace<TDomain> > approx
)
: VDCC_BG<TDomain>(fcts, subsets, approx), m_bIsConstData(false)
{
	// nothing to do
}

template <typename TDomain>
VDCC_BG_UserData<TDomain>::~VDCC_BG_UserData()
{
	// nothing to do
}

template<typename TDomain>
void VDCC_BG_UserData<TDomain>::set_potential_function(const char* name)
{
	// name must be a valid lua function name conforming to LuaUserNumber specs
	if (LuaUserData<number, dim>::check_callback_returns(name))
	{
		set_potential_function(LuaUserDataFactory<number, dim>::create(name));
		return;
	}

	// no match found
	if (!CheckLuaCallbackName(name))
		UG_THROW("Lua-Callback with name '" << name << "' does not exist.");

	// name exists, but wrong signature
	UG_THROW("Cannot find matching callback signature. Use:\n"
			"Number - Callback\n" << (LuaUserData<number, dim>::signature()) << "\n");
}

template<typename TDomain>
void VDCC_BG_UserData<TDomain>::set_potential_function(const number value)
{
	// if value is null, smart pointer will point to null
	if (value == 0.0) set_potential_function(SmartPtr<CplUserData<number, dim> >());
	else set_potential_function(make_sp(new ConstUserNumber<dim>(value)));

	m_bIsConstData = true;
}

template<typename TDomain>
void VDCC_BG_UserData<TDomain>::set_potential_function(SmartPtr<CplUserData<number, dim> > spPotFct)
{
	m_spPotential = spPotFct;
	m_bIsConstData = false;
}


template<typename TDomain>
void VDCC_BG_UserData<TDomain>::update_potential(side_t* elem)
{
	// only work if really necessary
	if (m_bIsConstData) return;

	// fill attachments with renewed values
	const typename TDomain::position_type& coords = CalculateCenter(elem, this->m_aaPos);
	number vm = -65.0;
	if (this->m_spPotential.valid())
		(*this->m_spPotential)(vm, coords, this->m_time, this->m_sh->get_subset_index(elem));

	// set membrane potential value
	this->m_aaVm[elem] = 0.001 * vm;
}


// explicit template specializations
#ifdef UG_DIM_1
	template class VDCC_BG<Domain1d>;
	template class VDCC_BG_VM2UG<Domain1d>;
	#ifdef MPMNEURON
		template class VDCC_BG_VM2UG_NEURON<Domain1d>;
	#endif
	template class VDCC_BG_UserData<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class VDCC_BG<Domain2d>;
	template class VDCC_BG_VM2UG<Domain2d>;
	#ifdef MPMNEURON
		template class VDCC_BG_VM2UG_NEURON<Domain2d>;
	#endif
	template class VDCC_BG_UserData<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class VDCC_BG<Domain3d>;
	template class VDCC_BG_VM2UG<Domain3d>;
	#ifdef MPMNEURON
		template class VDCC_BG_VM2UG_NEURON<Domain1d>;
	#endif
	template class VDCC_BG_UserData<Domain3d>;
#endif

} // neuro_collection
} // namespace ug
