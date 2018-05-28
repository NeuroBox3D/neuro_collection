/*
 *	Discretization for a simplicstic NMDA receptor channel.
 *
 *  Created on: 2018-04-19
 *      Author: mbreit
 */

#include "nmdar.h"


namespace ug {
namespace neuro_collection {


NMDAR::NMDAR(const std::vector<std::string>& fcts)
: IMembraneTransporter(fcts),
  m_z(2.0), m_F(96485.3329), m_R(8.31446),
  m_T(310.0), m_t0(0.0), m_tau(0.15), m_perm(7.61e-21), m_vm(-0.07),
  m_time(0.0)
{
	// nothing to do
}

NMDAR::NMDAR(const char* fcts)
: IMembraneTransporter(fcts),
  m_z(2.0), m_F(96485.3329), m_R(8.31446),
  m_T(310.0), m_t0(0.0),
  m_tau(0.15), m_perm(7.61e-21), m_vm(-0.07),
  m_time(0.0)
{
	// nothing to do
}

NMDAR::~NMDAR()
{
	// nothing to do
}


void NMDAR::set_activation_time(number t0)
{
	m_t0 = t0;
}


void NMDAR::set_decay_time(number tau)
{
	m_tau = tau;
}


void NMDAR::set_permeability(number perm)
{
	m_perm = perm;
}


void NMDAR::set_membrane_potential(number vm)
{
	m_vm = vm;
}


void NMDAR::set_temperature(number temp)
{
	m_T = temp;
}



void NMDAR::prep_timestep(number future_time, const number time, VectorProxyBase* upb)
{
	m_time = future_time;
}


void NMDAR::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	const number c_ext = u[_EXT_];	// source concentration
	const number c_int = u[_CYT_];	// target concentration

	const number p_open = (m_time >= m_t0) ? exp(-(m_time - m_t0) / m_tau) : 0.0;
	const number current = (fabs(m_vm) < 1e-8) ?
		m_perm * ((c_ext - c_int) - m_z*m_F/(2*m_R*m_T) * (c_ext + c_int) * m_vm)
		: -m_perm * m_z*m_F/(m_R*m_T) * m_vm * (c_ext - c_int*exp(m_z*m_F/(m_R*m_T)*m_vm)) / (1.0 - exp(m_z*m_F/(m_R*m_T)*m_vm));
	flux[0] = p_open * current;
}


void NMDAR::calc_flux_deriv
(
	const std::vector<number>& u,
	GridObject* e,
	std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs
) const
{
	number deriv_ext;
	number deriv_cyt;
	const number p_open = (m_time >= m_t0) ? exp(-(m_time - m_t0) / m_tau) : 0.0;
	if (fabs(m_vm) < 1e-8)
	{
		deriv_ext = -p_open * m_perm * (-1.0 + m_z*m_F/(2*m_R*m_T) * m_vm);
		deriv_cyt = -p_open * m_perm * (1.0 + m_z*m_F/(2*m_R*m_T) * m_vm);
	}
	else
	{
		const number zFRTVm = m_z*m_F/(m_R*m_T) * m_vm;
		deriv_ext = -p_open * m_perm * zFRTVm / (1.0 - exp(zFRTVm));
		deriv_cyt = -p_open * m_perm * zFRTVm / (1.0 - exp(-zFRTVm));
	}

	size_t i = 0;
	if (!has_constant_value(_EXT_))
	{
		flux_derivs[0][i].first = local_fct_index(_EXT_);
		flux_derivs[0][i].second = deriv_ext;
		i++;
	}
	if (!has_constant_value(_CYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CYT_);
		flux_derivs[0][i].second = deriv_cyt;
		i++;
	}
}


// return number of unknowns this transport mechanism depends on
size_t NMDAR::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_EXT_))
		n--;
	if (has_constant_value(_CYT_))
		n--;

	return n;
}

// return number of fluxes calculated by this machanism
size_t NMDAR::n_fluxes() const
{
	return 1;
};

// from where to where do the fluxes occur
const std::pair<size_t,size_t> NMDAR::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_EXT_)) from = local_fct_index(_EXT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CYT_)) to = local_fct_index(_CYT_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}

const std::string NMDAR::name() const
{
	return std::string("NMDAR");
}

void NMDAR::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a current would be of no consequence.
	if (!allows_flux(_EXT_) && !allows_flux(_CYT_))
	{
		UG_THROW("Supplying neither source nor target concentrations is not allowed.\n"
				"This would mean that the current calculation would be of no consequence\n"
				"and this NMDAR would not do anything.");
	}
}


void NMDAR::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      Concentrations     mM (= mol/m^3)                                       |"<< std::endl);
	UG_LOG("|      Temperature        K                                                    |"<< std::endl);
	UG_LOG("|      Time               s                                                    |"<< std::endl);
	UG_LOG("|      Permeability       m^3/s                                                |"<< std::endl);
	UG_LOG("|      Membrane voltage   V                                                    |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output (in MembraneTransport)                                             |"<< std::endl);
	UG_LOG("|      Flux DENSITY       mol/(m^2 s)                                          |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}



} // namespace neuro_collection
} // namespace ug

