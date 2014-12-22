/*
 * BorgGraham_impl.h
 *
 *  Created on: 05.02.2013
 *      Author: mbreit
 */

#include "one_sided_borg_graham_fv1.h"

namespace ug
{
namespace neuro_collection
{

///////////////////////////////////////////////////////////
///////////   BorgGraham   ////////////////////////////////
///////////////////////////////////////////////////////////
template<typename TDomain>
OneSidedBorgGrahamFV1<TDomain>::OneSidedBorgGrahamFV1
(
	const char* functions,
	const char* subsets,
	ApproximationSpace<TDomain>& approx
)
	: OneSidedMembraneTransportFV1<TDomain>(functions, subsets), RHO_BG(5.0),
	  m_dom(approx.domain()), m_mg(approx.domain()->grid()), m_dd(approx.dof_distribution(GridLevel::TOP)),
	  m_sh(m_dom->subset_handler()), m_aaPos(m_dom->position_accessor()),
	  m_gpMGate(3.4, -21.0, 1.5), m_gpHGate(-2.0, -40.0, 75.0), m_time(0.0), m_oldTime(0.0),
	  m_lambda(2.0e-11), m_mp(2), m_hp(1), m_channelType(BG_Ntype), m_initiated(false)
{
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
void OneSidedBorgGrahamFV1<TDomain>::prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	// remember
	m_bNonRegularGrid = bNonRegularGrid;

	// set assembling functions from base class first
	this->DependentNeumannBoundaryFV1<TDomain>::prepare_setting(vLfeID, bNonRegularGrid);

	// update assemble functions
	register_all_fv1_funcs();
}


template<typename TDomain>
template<int TType>
void OneSidedBorgGrahamFV1<TDomain>::set_channel_type()
{
	if (m_initiated)
	{
		UG_THROW("Borg-Graham channel type can not be changed after initialization.\n"
				 "Call set_channel_type() BEFORE init().");
	}

	switch (TType)
	{
		// - all gating params according to table 5, page 98 of the Borg-Graham article
		//   "Interpretations of data and mechanisms for hippocampal pyramidal cell models"
		// - conductance values from Fisher et al. - "Properties and distribution of single
		//   voltage-gated calcium channels in adult hippocampal neurons"
		case BG_Ntype:
			m_gpMGate = GatingParams(3.4, -21.0, 1.5);
			m_gpHGate = GatingParams(-2.0, -40.0, 75.0);
			m_lambda = 0.2e-12;	// reference concentration: 1mM
			m_mp = 2;
			m_hp = 1;
			break;
		case BG_Ltype:
			m_gpMGate = GatingParams(4.6, -1.0, 1.5);
			m_gpHGate = GatingParams(0.0, 0.0, 0.0);
			m_lambda = 0.3e-12;	// reference concentration: 1mM
			m_mp = 2;
			m_hp = 0;
			break;
		case BG_Ttype:
			m_gpMGate = GatingParams(3.0, -36.0, 1.5);
			m_gpHGate = GatingParams(-5.2, -68.0, 10.0);
			m_lambda = 0.1e-12;	// reference concentration: 1mM
			m_mp = 2;
			m_hp = 1;
			break;
		default:
			UG_THROW("Type of Borg-Graham channel does not match any of the pre-implemented.");
			break;
	}
}


template<typename TDomain>
number OneSidedBorgGrahamFV1<TDomain>::calc_gating_start(GatingParams& gp, number Vm)
{
	return 1.0 / (1.0 + exp(-gp.z * (Vm - gp.V_12)/1000.0 * F/(R*T)));
}


template<typename TDomain>
void OneSidedBorgGrahamFV1<TDomain>::calc_gating_step(GatingParams& gp, number Vm, number dt, number& currVal)
{
	// forward step: implicit
	if (dt>=0) currVal = (currVal + dt/gp.tau_0 * calc_gating_start(gp,Vm)) / (1.0 + dt/gp.tau_0);
	// backward step: explicit
	else currVal += dt/gp.tau_0 * (calc_gating_start(gp,Vm) - currVal);
}


template<typename TDomain>
template<typename TElem>
void OneSidedBorgGrahamFV1<TDomain>::prep_timestep_elem
(
	const number time,
	const LocalVector& u,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	update_time(time);

	side_t* side = dynamic_cast<side_t*>(elem);
	if (!side) UG_THROW("OneSidedBorgGrahamFV1::prep_timestep_elem() called with improper element type.");

	update_potential(side);
	update_gating(side);
}


template<typename TDomain>
bool OneSidedBorgGrahamFV1<TDomain>::fluxDensityFct(const std::vector<LocalVector::value_type>& u, GridObject* e, const MathVector<dim>& coords, int si, NFluxCond& fc)
{
	side_t* elem = dynamic_cast<side_t*>(e);
	if (!elem) UG_THROW("OneSidedBorgGrahamFV1 fluxDensityFunction called with the wrong type of element.");

	number current = ionic_current(elem);

	number density;
	if (this->m_spDensityFct.valid())
		(*this->m_spDensityFct)(density, coords, this->time(), si);
	else
		density = RHO_BG;

	number flux = density * current;

	// dimensional correction: concentrations are mol/dm^3, but length unit is um
	flux *= 1e15;

	fc.flux.resize(1, 0.0);	fc.flux[0] = flux;
	fc.to.resize(1);	fc.to[0] = 0;

	return true;
}


template<typename TDomain>
bool OneSidedBorgGrahamFV1<TDomain>::fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, GridObject* e, const MathVector<dim>& coords, int si, NFluxDerivCond& fdc)
{
	// nothing, since not dependent on any unknowns in the current implementation

	// add to Jacobian
	fdc.fluxDeriv.resize(1);
	fdc.fluxDeriv[0].resize(1);

	fdc.fluxDeriv[0][0] = 0.0;

	fdc.to.resize(1);	fdc.to[0] = 0;

	return true;
}


template<typename TDomain>
void OneSidedBorgGrahamFV1<TDomain>::init(number time)
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
void OneSidedBorgGrahamFV1<TDomain>::update_gating(side_t* elem)
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
number OneSidedBorgGrahamFV1<TDomain>::ionic_current(side_t* e)
{
	number gating = pow(m_aaMGate[e], m_mp);
	if (has_hGate()) gating *= pow(m_aaHGate[e], m_hp);

	// simplified form of the correct flux derived from Goldman-Hodgkin-Katz equation,
	// very accurate for Vm < 0mV and still reasonably accurate for Vm < 50mV
	number maxFlux;
	// just to be on the safe side near V_m == 0:
	if (fabs(m_aaVm[e]) < 1e-8) maxFlux = m_lambda * R*T/(2*F);
	else maxFlux = - m_lambda * m_aaVm[e] / (1.0 - exp(2*F/(R*T) * m_aaVm[e]));

	return gating * maxFlux / (2*F);
}



///////////////////////////////////////////////////////////
///////////   BorgGrahamWithVM2UG   ///////////////////////
///////////////////////////////////////////////////////////

template<typename TDomain>
void OneSidedBorgGrahamFV1WithVM2UG<TDomain>::init(number time)
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
void OneSidedBorgGrahamFV1WithVM2UG<TDomain>::update_time(number newTime)
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
void OneSidedBorgGrahamFV1WithVM2UG<TDomain>::update_potential(side_t* elem)
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
void OneSidedBorgGrahamFV1WithVM2UGNEURON<TDomain>::init(number time) {
	try {
		/// make zero steps (e. g. init) and extract the membrane potentials
		this->m_NrnInterpreter.get()->extract_vms(1, 3);
		/// with the extracted vms we build the tree then
		std::cout << "vms could be extracted" << std::endl;
		std::cout << "Our NEURON setup: " << std::endl;
		this->m_NrnInterpreter.get()->print_setup(true);
		std::cout << "Our potential: " << this->m_vmProvider.get()->get_potential(0, 0, 0);
	} UG_CATCH_THROW("NEURON interpreter could not advance, vmProvider could not be created.");

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
void OneSidedBorgGrahamFV1WithVM2UGNEURON<TDomain>::update_time(number newTime)
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
void OneSidedBorgGrahamFV1WithVM2UGNEURON<TDomain>::update_potential(side_t* elem)
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
			UG_THROW("Coordinates of vertex are not 2d or 3d"
				  << "(which are the only valid dimensionalities for this discretization).");
	}
	UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

	// set membrane potential value
	this->m_aaVm[elem] = 0.001 * vm; // todo: next step by dt * vm not by 0.001?
									// -> no, this is not the step size, but correction
									// for units being [V] in the attachments, but [mV] in Vm2uG.
}


#endif




///////////////////////////////////////////////////////////
///////////   BorgGrahamWithUserData   ////////////////////
///////////////////////////////////////////////////////////

template<typename TDomain>
void OneSidedBorgGrahamFV1WithUserData<TDomain>::set_potential_function(const char* name)
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
void OneSidedBorgGrahamFV1WithUserData<TDomain>::set_potential_function(const number value)
{
	// if value is null, smart pointer will point to null
	if (value == 0.0) set_potential_function(SmartPtr<CplUserData<number, dim> >());
	else set_potential_function(make_sp(new ConstUserNumber<dim>(value)));

	m_bIsConstData = true;
}

template<typename TDomain>
void OneSidedBorgGrahamFV1WithUserData<TDomain>::set_potential_function(SmartPtr<CplUserData<number, dim> > spPotFct)
{
	m_spPotential = spPotFct;
	m_bIsConstData = false;
}


template<typename TDomain>
void OneSidedBorgGrahamFV1WithUserData<TDomain>::update_potential(side_t* elem)
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


// register for 1D
template<typename TDomain>
void OneSidedBorgGrahamFV1<TDomain>::register_all_fv1_funcs()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::ManifoldElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList>(RegisterFV1(this));
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void OneSidedBorgGrahamFV1<TDomain>::register_fv1_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->set_prep_timestep_elem_fct(id, &T::template prep_timestep_elem<TElem>);
}
} // neuro_collection
} // namespace ug
