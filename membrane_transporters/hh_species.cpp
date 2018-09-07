/*
 *  Discretization of Hodgkin-Huxley voltage-dependent channels in the plasma membrane
 *
 *  Created on: 2018-09-05
 *      Author: mbreit
 */

#include "hh_species.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"  // for GeomProvider
#include "lib_disc/function_spaces/grid_function.h"  // for GridFunction

namespace ug {
namespace neuro_collection {


template <typename TDomain>
void HHSpecies<TDomain>::set_conductances(number gk, number gna)
{
	m_gK = gk;
	m_gNa = gna;
}


template <typename TDomain>
void HHSpecies<TDomain>::set_reversal_potentials(number ek, number ena)
{
	m_eK = ek;
	m_eNa = ena;
}


template <typename TDomain>
void HHSpecies<TDomain>::set_reference_time(number refTime)
{
	m_refTime = refTime;
}


template <typename TDomain>
void HHSpecies<TDomain>::use_exact_gating_mode()
{
	m_bVoltageExplicitDiscMode = true;
}


template <typename TDomain>
HHSpecies<TDomain>::HHSpecies
(
	const std::vector<std::string>& fcts,
	const std::vector<std::string>& subsets,
	ConstSmartPtr<ISubsetHandler> spSH
)
: IMembraneTransporter(fcts),
  m_gK(2e-11), m_gNa(2e-11),
  m_eK(-0.077), m_eNa(0.05),
  m_spSH(spSH),
  m_vSubset(subsets),
  m_refTime(1.0),
  m_time(0.0),
  m_initTime(0.0),
  m_oldTime(0.0),
  m_bInitiated(false),
  m_bVoltageExplicitDiscMode(false)
{
	// nothing to do
}


template <typename TDomain>
HHSpecies<TDomain>::HHSpecies(const char* fcts, const char* subsets, ConstSmartPtr<ISubsetHandler> spSH)
: IMembraneTransporter(fcts),
  m_gK(2e-11), m_gNa(2e-11),
  m_eK(-0.077), m_eNa(0.05),
  m_spSH(spSH),
  m_vSubset(TokenizeString(subsets, ',')),
  m_refTime(1.0),
  m_time(0.0),
  m_initTime(0.0),
  m_oldTime(0.0),
  m_bInitiated(false),
  m_bVoltageExplicitDiscMode(false)
{
	// nothing to do
}


template <typename TDomain>
HHSpecies<TDomain>::~HHSpecies()
{
	// nothing to do
}


template <typename TDomain>
void HHSpecies<TDomain>::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	const number vm = u[_PHII_] - u[_PHIO_]; // membrane potential
	typename GatingMap::const_iterator it = m_mGating.find(e);
	UG_COND_THROW(it == m_mGating.end(), "No gating value available on element.");
	const number n = it->second.n;  // gating state m
	const number m = it->second.m;  // gating state m
	const number h = it->second.h;  // gating state h

	const number currentK = m_gK * n*n*n*n * (vm - m_eK);
	const number currentNa = m_gNa * m*m*m*h * (vm - m_eNa);

	flux[0] = currentK;
	flux[1] = currentNa;
//UG_LOGN(std::setprecision(16) << "Kleak: " << currentK << "   Naleak: " << currentNa);
}


template <typename TDomain>
void HHSpecies<TDomain>::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	typename GatingMap::const_iterator it = m_mGating.find(e);
	UG_COND_THROW(it == m_mGating.end(), "No gating value available on element.");
	const number n = it->second.n;  // gating state n
	const number m = it->second.m;  // gating state m
	const number h = it->second.h;  // gating state h

	size_t i = 0;
	if (!has_constant_value(_PHII_))
	{
		flux_derivs[0][i].first = local_fct_index(_PHII_);
		flux_derivs[0][i].second = m_gK * n*n*n*n;
		flux_derivs[1][i].first = local_fct_index(_PHII_);
		flux_derivs[1][i].second = m_gNa * m*m*m*h;
		++i;
	}
	if (!has_constant_value(_PHIO_))
	{
		flux_derivs[0][i].first = local_fct_index(_PHIO_);
		flux_derivs[0][i].second = - m_gK * n*n*n*n;
		flux_derivs[1][i].first = local_fct_index(_PHIO_);
		flux_derivs[1][i].second = - m_gNa * m*m*m*h;
		++i;
	}
}


template <typename TDomain>
size_t HHSpecies<TDomain>::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_PHII_))
		--n;
	if (has_constant_value(_PHIO_))
		--n;

	return n;
}


template <typename TDomain>
size_t HHSpecies<TDomain>::n_fluxes() const
{
	return 2;
}


template <typename TDomain>
const std::pair<size_t,size_t> HHSpecies<TDomain>::flux_from_to(size_t flux_i) const
{
	// current goes from the inside charge density to the outside charge density
	size_t from, to;
	if (flux_i == 0)
	{
		to = allows_flux(_KO_) ? local_fct_index(_KO_) : InnerBoundaryConstants::_IGNORE_;
		from = allows_flux(_KI_) ? local_fct_index(_KI_) : InnerBoundaryConstants::_IGNORE_;
	}
	else if (flux_i == 1)
	{
		to = allows_flux(_NAO_) ? local_fct_index(_NAO_) : InnerBoundaryConstants::_IGNORE_;
		from = allows_flux(_NAI_) ? local_fct_index(_NAI_) : InnerBoundaryConstants::_IGNORE_;
	}
	else UG_THROW("Flux only has 2 components, but component " << flux_i << " was queried for.");

	return std::pair<size_t, size_t>(from, to);
}


template <typename TDomain>
const std::string HHSpecies<TDomain>::name() const
{
	return std::string("HH channel");
}


template <typename TDomain>
void HHSpecies<TDomain>::check_supplied_functions() const
{
	// Check that not both, inner and outer K+/Na+ concentrations are not supplied;
	// in that case, calculation of a current would be of no consequence.
	if (!allows_flux(_KI_) && !allows_flux(_KO_))
	{
		UG_THROW("Supplying neither inner nor outer Na concentration is not allowed.\n"
				"This would mean that the current calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}
	if (!allows_flux(_NAI_) && !allows_flux(_NAO_))
	{
		UG_THROW("Supplying neither inner nor outer Na concentration is not allowed.\n"
				"This would mean that the current calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}
}


template <typename TDomain>
void HHSpecies<TDomain>::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      Phi_i    inner potential  V                                             |"<< std::endl);
	UG_LOG("|      Phi_o    outer potential  V                                             |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|      E_K      K reversal potential    V                                      |"<< std::endl);
	UG_LOG("|      E_Na     Na reversal potential   V                                      |"<< std::endl);
	UG_LOG("|      g_K      K channel conductance   C/(Vs)                                 |"<< std::endl);
	UG_LOG("|      g_Na     Na channel conductance  C/(Vs)                                 |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      current  C/s                                                            |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}



static number alpha_n(number u)
{
	const number x = -(u + 0.055);
	if (fabs(x) > 1e-10)
		return 1e4*x / (exp(100.0*x) - 1.0);

	return 1e4 * (0.01 - 0.5*x);
}

static number beta_n(number u)
{
	return 125.0 * exp(-(u + 0.065) / 0.08);
}

static number n_infty(number u)
{
	return alpha_n(u) / (alpha_n(u) + beta_n(u));
}

static number tau_n(number u)
{
	return 1.0 / (alpha_n(u) + beta_n(u));
}


static number alpha_m(number u)
{
	const number x = -(u + 0.04);
	if (fabs(x) > 1e-10)
		return  1e5*x / (exp(100.0*x) - 1.0);

	return 1e5 * (0.01 - 0.5*x);
}

static number beta_m(number u)
{
	return 4e3 * exp(-(u + 0.065) / 0.018);
}

static number m_infty(number u)
{
	return alpha_m(u) / (alpha_m(u) + beta_m(u));
}

static number tau_m(number u)
{
	return 1.0 / (alpha_m(u) + beta_m(u));
}


static number alpha_h(number u)
{
	return 70.0 * exp(-(u + 0.065) / 0.02);
}

static number beta_h(number u)
{
	return 1e3 / (exp(-(u + 0.035) / 0.01) + 1.0);
}

static number h_infty(number u)
{
	return alpha_h(u) / (alpha_h(u) + beta_h(u));
}

static number tau_h(number u)
{
	return 1.0 / (alpha_h(u) + beta_h(u));
}



template<typename TDomain>
template <typename TAlgebra, int locDim>
void HHSpecies<TDomain>::update_potential
(
	ConstSmartPtr<DoFDistribution> dd,
	int si,
	const typename TAlgebra::vector_type& u
)
{
	size_t globFctIndPhiIn;
	size_t globFctIndPhiOut;
	if (!has_constant_value(_PHII_))
	{
		size_t locFctIndI = local_fct_index(_PHII_);
		try {globFctIndPhiIn = dd->fct_id_by_name(this->m_vFct[locFctIndI].c_str());}
		UG_CATCH_THROW("Function '" << this->m_vFct[locFctIndI] << "' not found in function pattern.");
	}
	if (!has_constant_value(_PHIO_))
	{
		size_t locFctIndO = local_fct_index(_PHII_);
		try {globFctIndPhiOut = dd->fct_id_by_name(this->m_vFct[locFctIndO].c_str());}
		UG_CATCH_THROW("Function '" << this->m_vFct[locFctIndO] << "' not found in function pattern.");
	}


	typedef typename GeomObjBaseTypeByDim<locDim>::base_obj_type elem_t;
	typedef typename DoFDistribution::traits<elem_t>::const_iterator it_type;

	std::vector<DoFIndex> vDI;

	it_type it = dd->begin<elem_t>(si);
	it_type it_end = dd->end<elem_t>(si);
	for (; it != it_end; ++it)
	{
		elem_t* elem = *it;

		// get inner potential
		number phiI = 0.0;
		if (!has_constant_value(_PHII_, phiI))
		{
			// get function index for inner potential
			dd->dof_indices(elem, globFctIndPhiIn, vDI, false, true);
			const size_t nDI = vDI.size();
			UG_COND_THROW(nDI < 1, "No DoF indices for function " << globFctIndPhiIn << " on element.");

			// calculate value in middle of elem
			for (size_t di = 0; di < nDI; ++di)
				phiI += DoFRef(u, vDI[di]);
			phiI /= nDI;
		}

		// get outer potential
		number phiO = 0.0;
		if (!has_constant_value(_PHIO_, phiO))
		{
			// get function index for outer potential
			dd->dof_indices(elem, globFctIndPhiOut, vDI, false, true);
			const size_t nDI = vDI.size();
			UG_COND_THROW(nDI < 1, "No DoF indices for function " << globFctIndPhiOut << " on element.");

			// calculate value in middle of elem
			for (size_t di = 0; di < nDI; ++di)
				phiO += DoFRef(u, vDI[di]);
			phiO /= nDI;
		}

		// save (scaled) membrane potential
		m_mGating[elem].vm = (phiI*scale_input(_PHII_) - phiO*scale_input(_PHIO_));
	}
}


template<typename TDomain>
void HHSpecies<TDomain>::init_gating(GridObject* elem)
{
	// calculate corresponding start condition for gates and save them
	typename GatingMap::iterator it = m_mGating.find(elem);
	UG_COND_THROW(it == m_mGating.end(), "No gating value available on element.");
	GatingInfo& gatings = it->second;
	const number& vm = gatings.vm;
	number& n = gatings.n;
	number& m = gatings.m;
	number& h = gatings.h;

	n = n_infty(vm);
	m = m_infty(vm);
	h = h_infty(vm);
}


template<typename TDomain>
void HHSpecies<TDomain>::update_gating(GridObject* elem)
{
	// time step
	number dt = m_refTime * (m_time - m_oldTime);
	typename GatingMap::iterator it = m_mGating.find(elem);
	UG_COND_THROW(it == m_mGating.end(), "No gating value available on element.");
	GatingInfo& gatings = it->second;
	const number& vm = gatings.vm;
	number& n = gatings.n;
	number& m = gatings.m;
	number& h = gatings.h;

	const number ninf = n_infty(vm);
	const number minf = m_infty(vm);
	const number hinf = h_infty(vm);
	const number taun = tau_n(vm);
	const number taum = tau_m(vm);
	const number tauh = tau_h(vm);

	if (m_bVoltageExplicitDiscMode)
	{
		// solve gating ODEs directly
		n = ninf + (n - ninf) * exp(-dt/taun);
		m = minf + (m - minf) * exp(-dt/taum);
		h = hinf + (h - hinf) * exp(-dt/tauh);
	}
	else
	{
		// solve gating ODEs using implicit Euler
		n = (taun * n + dt * ninf) / (dt + taun);
		m = (taum * m + dt * minf) / (dt + taum);
		h = (tauh * h + dt * hinf) / (dt + tauh);
	}
}

template<typename TDomain>
void HHSpecies<TDomain>::update_time(const number newTime)
{
	if (newTime != m_time)
	{
		m_oldTime = m_time;
		m_time = newTime;
	}
};


template <typename TDomain>
template <typename TAlgebra, int locDim>
void HHSpecies<TDomain>::prep_timestep_with_algebra_type_and_dim
(
	number future_time,
	const number time,
	const typename TAlgebra::vector_type& u,
	int si
)
{
	typedef typename GeomObjBaseTypeByDim<locDim>::base_obj_type elem_t;
	typedef typename DoFDistribution::traits<elem_t>::const_iterator it_type;

	// get DoFDistro from solution
	const GridFunction<TDomain, TAlgebra>* gfU = dynamic_cast<const GridFunction<TDomain, TAlgebra>*>(&u);
	UG_COND_THROW(!gfU, "Passed solution vector is not a grid function.");
	ConstSmartPtr<DoFDistribution> dd = gfU->dd();

	// update or init potential values
	update_potential<TAlgebra, locDim>(dd, si, u);

	// initiate gatings if this has not already been done (or init again; stationary case)
	if (!m_bInitiated || future_time == m_initTime)
	{
		m_time = time;
		m_initTime = time;

		it_type it = dd->begin<elem_t>(si);
		it_type it_end = dd->end<elem_t>(si);
		for (; it != it_end; ++it)
			init_gating(*it);
	}

	update_time(future_time);

	// loop sides and update gatings
	it_type it = dd->begin<elem_t>(si);
	it_type it_end = dd->end<elem_t>(si);
	for (; it != it_end; ++it)
		update_gating(*it);
}


template <typename TDomain>
template <typename TAlgebra>
void HHSpecies<TDomain>::prep_timestep_with_algebra_type
(
	number future_time,
	const number time,
	const typename TAlgebra::vector_type& u
)
{
	SubsetGroup ssGrp;
	try {ssGrp = SubsetGroup(m_spSH, this->m_vSubset);}
	UG_CATCH_THROW("Subset group creation failed.");

	const size_t nSs = ssGrp.size();
	for (std::size_t si = 0; si < nSs; ++si)
	{
		int ssDim = DimensionOfSubset(*m_spSH, ssGrp[si]);
		if (ssDim == 3)
			prep_timestep_with_algebra_type_and_dim<TAlgebra, 3>(future_time, time, u, ssGrp[si]);
		else if (ssDim == 2)
			prep_timestep_with_algebra_type_and_dim<TAlgebra, 2>(future_time, time, u, ssGrp[si]);
		else if (ssDim == 1)
			prep_timestep_with_algebra_type_and_dim<TAlgebra, 1>(future_time, time, u, ssGrp[si]);
		else if (ssDim == 0)
			prep_timestep_with_algebra_type_and_dim<TAlgebra, 0>(future_time, time, u, ssGrp[si]);
		else UG_THROW("Subset dimension " << ssDim << " is not supported.");
	}

	m_bInitiated = true;
}


template <typename TAlgebra>
static bool vector_from_algebra_type(VectorProxyBase* upb)
{
	typedef typename TAlgebra::vector_type v_type;
	typedef VectorProxy<v_type> vp_type;
	vp_type* up = dynamic_cast<vp_type*>(upb);
	return up != NULL;
}


template <typename TDomain>
void HHSpecies<TDomain>::prep_timestep(number future_time, const number time, VectorProxyBase* upb)
{
	// as long as we cannot template this method, we have to find out the algebra type ourselves
#ifdef UG_CPU_1
	if (vector_from_algebra_type<CPUAlgebra>(upb))
	{
		prep_timestep_with_algebra_type<CPUAlgebra>(future_time, time,
			(dynamic_cast<VectorProxy<typename CPUAlgebra::vector_type>*>(upb))->m_v);
		return;
	}
#endif
#ifdef UG_CPU_5
	if (vector_from_algebra_type<CPUBlockAlgebra<5> >(upb))
	{
		prep_timestep_with_algebra_type<CPUBlockAlgebra<5> >(future_time, time,
			(dynamic_cast<VectorProxy<typename CPUBlockAlgebra<5>::vector_type>*>(upb))->m_v);
		return;
	}
#endif
	UG_THROW("Given algebra type is not treated by this class.");
}



// explicit template specializations
#ifdef UG_DIM_1
	template class HHSpecies<Domain1d>;
#endif
#ifdef UG_DIM_2
	template class HHSpecies<Domain2d>;
#endif
#ifdef UG_DIM_3
	template class HHSpecies<Domain3d>;
#endif



} // namespace neuro_collection
} // namespace ug

