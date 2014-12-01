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
bool OneSidedBorgGrahamFV1<TDomain>::fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, NFluxCond& fc)
{
	Vertex* vert = this->m_currVertex;
	number current = ionic_current(vert);

	number density;
	if (this->m_spDensityFct.valid())
		(*this->m_spDensityFct)(density, coords, this->time(), si);
	else
		density = RHO_BG;

	number flux = density * current;

	// dimensional correction: concentrations are mol/dm^3, but length unit is um
	flux *= 1e15;

	fc.flux.resize(1, 0.0);	fc.flux[0] = flux;
	fc.to.resize(1);	fc.to[0] = 0;		// cytosol

	return true;
}


template<typename TDomain>
bool OneSidedBorgGrahamFV1<TDomain>::fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& coords, int si, NFluxDerivCond& fdc)
{
	// nothing, since not dependent on any unknowns in the current implementation

	// add to Jacobian
	fdc.fluxDeriv.resize(1);
	fdc.fluxDeriv[0].resize(1);

	fdc.fluxDeriv[0][0] = 0.0;

	fdc.to.resize(1);	fdc.to[0] = 0;		// cytosol

	return true;
}


///////////////////////////////////////////////////////////
///////////   BorgGrahamWithVM2UG   ///////////////////////
///////////////////////////////////////////////////////////

template<typename TDomain>
void OneSidedBorgGrahamFV1WithVM2UG<TDomain>::init(number time)
{
	// attach attachments
	if (this->m_mg->has_vertex_attachment(this->m_MGate))
		UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_MGate);

	if (has_hGate())
	{
		if (this->m_mg->has_vertex_attachment(this->m_HGate))
			UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
					 "could not be made, since it already exists.");
		this->m_mg->attach_to_vertices(this->m_HGate);
	}

	if (this->m_mg->has_vertex_attachment(this->m_Vm))
		UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_Vm);


	// create attachment accessors
	this->m_aaMGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_MGate);
	if (has_hGate()) this->m_aaHGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_HGate);
	this->m_aaVm = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_Vm);


	// fill attachments with initial values
	std::string timeAsString;
	try
	{
		char buffer[100];
		char* oldLocale = setlocale (LC_ALL, NULL);
		setlocale(LC_NUMERIC, "C");		// ensure decimal point is a . (not a ,)
		sprintf(buffer, m_tFmt.c_str(), time);
		setlocale(LC_NUMERIC, oldLocale);
		timeAsString = buffer;
	}
	UG_CATCH_THROW("Time format string provided does not meet requirements.\n"
				<< "It must contain exactly one placeholder (which has to convert a floating point number)"
				<< "and must not produce an output string longer than 100 chars.");

	try {m_vmProvider.buildTree(timeAsString);}
	UG_CATCH_THROW("Underlying Vm2uG object could not build its tree on given file.");

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	const typename TDomain::position_accessor_type& aaPos = this->m_dom->position_accessor();
	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// retrieve membrane potential via vm2ug
			number vm;
			try
			{
				const typename TDomain::position_type& coords = aaPos[*iter];
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
void OneSidedBorgGrahamFV1WithVM2UG<TDomain>::update_potential(number newTime)
{
	// only work if really necessary
	if (newTime == this->m_vmTime) return;

	// fill attachments with renewed values
	std::string timeAsString;
	try
	{
		char buffer[100];
		char* oldLocale = setlocale (LC_ALL, NULL);
		setlocale(LC_NUMERIC, "C");		// ensure decimal point is a . (not a ,)
		sprintf(buffer, m_tFmt.c_str(), newTime);
		setlocale(LC_NUMERIC, oldLocale);
		timeAsString = buffer;
	}
	UG_CATCH_THROW("Time format string provided does not meet requirements.\n"
				<< "It must contain exactly one placeholder (which has to convert a floating point number)"
				<< "and must not produce an output string longer than 100 chars.");

	if (!m_vmProvider.treeBuild())
	UG_THROW("Underlying Vm2uG object's tree is not yet built.\n"
		  << "Do not forget to initialize the Borg-Graham object first by calling init(initTime).");

	// set new timestep file in vm2ug object (not strictly necessary)
	const std::string ts = timeAsString;
	m_vmProvider.setTimestep(ts);

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	const typename TDomain::position_accessor_type& aaPos = this->m_dom->position_accessor();
	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// retrieve membrane potential via vm2ug
			number vm;
			try
			{
				const typename TDomain::position_type& coords = aaPos[*iter];
				if (coords.size() == 2)
					vm = m_vmProvider.get_potential(coords[0], coords[1], 0.0, timeAsString);
				else if (coords.size() == 3)
					vm = m_vmProvider.get_potential(coords[0], coords[1], coords[2], timeAsString);
				else
					UG_THROW("Coordinates of vertex are not 2d or 3d"
						  << "(which are the only valid dimensionalities for this discretization).");
			}
			UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

			// set membrane potential value
			this->m_aaVm[*iter] = 0.001 * vm;
		}
	}

	this->m_vmTime = newTime;
}


template<typename TDomain>
void OneSidedBorgGrahamFV1WithVM2UG<TDomain>::update_gating(number newTime)
{
	if (!this->m_initiated)
		UG_THROW("Borg-Graham not initialized.\n"
				  << "Do not forget to do so before any updates by calling init(initTime).");

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// set new gating particle values
			number dt = 1000.0*(newTime - this->m_time);	// calculating in ms
			this->calc_gating_step(this->m_gpMGate, 1000.0*this->m_aaVm[*iter], dt, this->m_aaMGate[*iter]);
			if (has_hGate()) this->calc_gating_step(this->m_gpHGate, 1000.0*this->m_aaVm[*iter], dt, this->m_aaHGate[*iter]);
		}
	}

	this->m_time = newTime;
}


template<typename TDomain>
number OneSidedBorgGrahamFV1WithVM2UG<TDomain>::ionic_current(Vertex* v)
{
	number gating = pow(this->m_aaMGate[v], this->m_mp);
	if (has_hGate()) gating *= pow(this->m_aaHGate[v], this->m_hp);

	// simplified form of the correct flux derived from Goldman-Hodgkin-Katz equation,
	// very accurate for Vm < 0mV and still reasonably accurate for Vm < 50mV
	number maxFlux;
	// just to be on the safe side near V_m == 0:
	if (fabs(this->m_aaVm[v]) < 1e-8) maxFlux = this->m_lambda * R*T/(2*F);
	else maxFlux = - this->m_lambda * this->m_aaVm[v] / (1.0 - exp(2*F/(R*T) * this->m_aaVm[v]));

	return gating * maxFlux / (2*F);
}


///////////////////////////////////////////////////////////
///////////   BorgGrahamWithNEURON   ///////////////////////
///////////////////////////////////////////////////////////
#ifdef MPMNEURON
template<typename TDomain>
void OneSidedBorgGrahamFV1WithVM2UGNEURON<TDomain>::init(number time)
{
	// attach attachments
	if (this->m_mg->has_vertex_attachment(this->m_MGate))
		UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_MGate);

	if (has_hGate())
	{
		if (this->m_mg->has_vertex_attachment(this->m_HGate))
			UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
					 "could not be made, since it already exists.");
		this->m_mg->attach_to_vertices(this->m_HGate);
	}

	if (this->m_mg->has_vertex_attachment(this->m_Vm))
		UG_THROW("Attachment necessary for Borg-Graham channel dynamics "
				 "could not be made, since it already exists.");
	this->m_mg->attach_to_vertices(this->m_Vm);


	// create attachment accessors
	this->m_aaMGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_MGate);
	if (has_hGate()) this->m_aaHGate = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_HGate);
	this->m_aaVm = Grid::AttachmentAccessor<Vertex, ADouble>(*this->m_mg, this->m_Vm);

	try {
		// this can be improved todo
		m_NrnInterpreter.get()->setup_hoc(time, 10000, 0.001, -75);
	}
	UG_CATCH_THROW("Either underlying NEURON interpreter not available or not constructable or vm2ug tree could not be built.");

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	const typename TDomain::position_accessor_type& aaPos = this->m_dom->position_accessor();
	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// retrieve membrane potential via vm2ug
			number vm;
			try
			{
				const typename TDomain::position_type& coords = aaPos[*iter];
				vm = m_vmProvider.get_potential(coords[0], coords[1], coords[2]);
			}
			UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

			this->m_aaMGate[*iter] = this->calc_gating_start(this->m_gpMGate, vm);
			if (has_hGate()) this->m_aaHGate[*iter] = this->calc_gating_start(this->m_gpHGate, vm);
			this->m_aaVm[*iter] = 0.001 * vm; // what is this? 0.001 -> should this not be dt todo
		}
	}

	this->m_time = time;
	this->m_initiated = true;
}


template<typename TDomain>
void OneSidedBorgGrahamFV1WithVM2UGNEURON<TDomain>::update_potential(number newTime)
{
	// only work if really necessary
	if (newTime == this->m_vmTime) return;

	if (!m_vmProvider.treeBuild())
	UG_THROW("Underlying Vm2uG object's tree is not yet built.\n"
		  << "Do not forget to initialize the Borg-Graham object first by calling init(initTime).");

	// set new timestep file in vm2ug object (not strictly necessary)
	const std::string ts = this->timeAsString;
	// todo advance by predefined dt in init method of OneSidedBorgGrahamFV1WithVM2UG above... (see above)
	m_NrnInterpreter.get()->fadvance();

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	const typename TDomain::position_accessor_type& aaPos = this->m_dom->position_accessor();
	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// retrieve membrane potential via vm2ug
			number vm;
			try
			{
				const typename TDomain::position_type& coords = aaPos[*iter];
				vm = m_vmProvider.get_potential(coords[0], coords[1], coords[2]);
			}
			UG_CATCH_THROW("Vm2uG object failed to retrieve a membrane potential for the vertex.");

			// set membrane potential value
			this->m_aaVm[*iter] = 0.001 * vm; // todo: next step by dt * vm not by 0.001?
		}
	}

	this->m_vmTime = newTime;
}


template<typename TDomain>
void OneSidedBorgGrahamFV1WithVM2UGNEURON<TDomain>::update_gating(number newTime)
{
	if (!this->m_initiated)
		UG_THROW("Borg-Graham not initialized.\n"
				  << "Do not forget to do so before any updates by calling init(initTime).");

	typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;
	SubsetGroup ssGrp;
	try { ssGrp = SubsetGroup(this->m_dom->subset_handler(), this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	for (std::size_t si = 0; si < ssGrp.size(); si++)
	{
		itType iterBegin = this->m_dd->template begin<Vertex>(ssGrp[si]);
		itType iterEnd = this->m_dd->template end<Vertex>(ssGrp[si]);

		for (itType iter = iterBegin; iter != iterEnd; ++iter)
		{
			// set new gating particle values
			number dt = 1000.0*(newTime - this->m_time);	// calculating in ms
			this->calc_gating_step(this->m_gpMGate, 1000.0*this->m_aaVm[*iter], dt, this->m_aaMGate[*iter]);
			if (has_hGate()) this->calc_gating_step(this->m_gpHGate, 1000.0*this->m_aaVm[*iter], dt, this->m_aaHGate[*iter]);
		}
	}

	this->m_time = newTime;
}


template<typename TDomain>
number OneSidedBorgGrahamFV1WithVM2UGNEURON<TDomain>::ionic_current(Vertex* v)
{
	number gating = pow(this->m_aaMGate[v], this->m_mp);
	if (has_hGate()) gating *= pow(this->m_aaHGate[v], this->m_hp);

	// simplified form of the correct flux derived from Goldman-Hodgkin-Katz equation,
	// very accurate for Vm < 0mV and still reasonably accurate for Vm < 50mV
	number maxFlux;
	// just to be on the safe side near V_m == 0:
	if (fabs(this->m_aaVm[v]) < 1e-8) maxFlux = this->m_lambda * R*T/(2*F);
	else maxFlux = - this->m_lambda * this->m_aaVm[v] / (1.0 - exp(2*F/(R*T) * this->m_aaVm[v]));

	return gating * maxFlux / (2*F);
}
#endif

} // neuro_collection
} // namespace ug
