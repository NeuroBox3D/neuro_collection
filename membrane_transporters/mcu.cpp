/*
 *	 Discretization for the mitochondrial uniporter MCU in the mitochondrial membrane
 *
 *  Created on: 30.06.2015
 *      Author: mstepnie
 */

#include "mcu.h"

namespace ug{
namespace neuro_collection{


MCU::MCU(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts),
F(0.096484), RT(2.5775),
K_C(3.965e-6), K_M(0.655e-3), K_Pi(0.2e-3), g(3.9), nH(2.65)
{
	m_k = 1.27e-3;

	m_pi_cyt  = 0.0;
    m_mg_cyt  = 0.0;
    m_mg_mit  = 0.0;
	m_psi 	  = 0.0;

	K_CC 	= K_C * (1 + m_pi_cyt/(K_Pi+m_pi_cyt));
	K_MM 	= K_M / (1 + m_pi_cyt/(K_Pi+m_pi_cyt));

	m_mit_volume  = 0.0;
	m_mit_surface = 0.0;
}

MCU::MCU(const char* fcts) : IMembraneTransporter(fcts),
F(0.096484), RT(2.5775),
K_C(3.965e-6), K_M(0.655e-3), K_Pi(0.2e-3), g(3.9), nH(2.65)
{
	m_k = 1.27e-3;

	m_pi_cyt  = 0.0;
    m_mg_cyt  = 0.0;
    m_mg_mit  = 0.0;
	m_psi 	  = 0.0;

	K_CC 	= K_C * (1 + m_pi_cyt/(K_Pi+m_pi_cyt));
	K_MM 	= K_M / (1 + m_pi_cyt/(K_Pi+m_pi_cyt));

	m_mit_volume  = 0.0;
	m_mit_surface = 0.0;
}


MCU::~MCU()
{
	// nothing to do
}


void MCU::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	if(m_mit_volume == 0.0 || m_mit_surface == 0.0)
		UG_THROW("ERROR in MCU membrane transport: mitochondrial volume or surface not specified.");

// 	get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caMit = u[_CMIT_];	// mitochondrial Ca2+ concentration

	number mgCyt = m_mg_cyt;
	number mgMit = m_mg_mit;

	double phi = 2.0 * F * m_psi / RT;

	double beta_e = 0.5 * (1 + nH/phi * log((phi/nH) / (sinh(phi/nH))));
	double beta_x = 0.5 * (1 - nH/phi * log((phi/nH) / (sinh(phi/nH))));

	double k_i = m_k * exp(2*beta_e*phi);
	double k_o = m_k * exp(-2*beta_x*phi);

	double D = 	1 + ((caCyt*caCyt)/(K_CC*K_CC)) + ((caMit*caMit)/(K_CC*K_CC)) +
					((mgCyt*mgCyt)/(K_MM*K_MM)) + ((mgMit*mgMit)/(K_MM*K_MM)) +
					((caCyt*caCyt*mgCyt*mgCyt)/(pow(g, 4)*K_CC*K_CC*K_MM*K_MM)) +
					((caMit*caMit*mgMit*mgMit)/(pow(g, 4)*K_CC*K_CC*K_MM*K_MM));

	flux[0] = 1/D * (k_i*((caCyt*caCyt)/(K_CC*K_CC)) - k_o*((caMit*caMit)/(K_CC*K_CC)));	// in nmol/mg/s

	UG_LOG("MCU::calc_flux: " << flux[0] << std::endl);

//	Transform original flux nmol/mg/s to mitochondrial flux nmol/s
//	1um^3 mitochondrial volume = 1e-9mg mitochondrial protein
	flux[0] *= 1e-9 * m_mit_volume;

//  Transform mitochondrial flux nmol/s to flux nmol/um^2/s
	flux[0] *= 1.0 / m_mit_surface;

//  Transform mitochondrial flux nmol/um^2/s to mol/um^2/s
	flux[0] *= 1e-9;
}


void MCU::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	if(m_mit_volume == 0.0 || m_mit_surface == 0.0)
		UG_THROW("ERROR in MCU membrane transport: mitochondrial volume or surface not specified.");

// 	get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caMit = u[_CMIT_];	// mitochondrial Ca2+ concentration

	number mgCyt = m_mg_cyt;
	number mgMit = m_mg_mit;

	double phi = 2.0 * F * m_psi / RT;

	double beta_e = 0.5 * (1 + nH/phi * log((phi/nH) / (sinh(phi/nH))));
	double beta_x = 0.5 * (1 - nH/phi * log((phi/nH) / (sinh(phi/nH))));

	double k_o = m_k * exp(-2*beta_x*phi);
	double k_i = m_k * exp(2*beta_e*phi);

	double D = 	1 + ((caCyt*caCyt)/(K_CC*K_CC)) + ((caMit*caMit)/(K_CC*K_CC)) +
					((mgCyt*mgCyt)/(K_MM*K_MM)) + ((mgMit*mgMit)/(K_MM*K_MM)) +
					((caCyt*caCyt*mgCyt*mgCyt)/(pow(g, 4)*K_CC*K_CC*K_MM*K_MM)) +
					((caMit*caMit*mgMit*mgMit)/(pow(g, 4)*K_CC*K_CC*K_MM*K_MM));

	double dD_caCyt = 2*caCyt/(K_C*K_C) + 2*caCyt*mgCyt*mgCyt/(pow(g, 4)*K_CC*K_CC*K_MM*K_MM);
	double dFlux_caCyt = 2*k_i*caCyt/(K_CC*K_CC)/D - dD_caCyt*(k_i*caCyt*caCyt/(K_CC*K_CC) - k_o*caMit*caMit/(K_CC*K_CC))/(D*D);
	dFlux_caCyt *= 1e-9 * m_mit_volume / m_mit_surface * 1e-9;

	double dD_caMit = 2*caMit/(K_C*K_C) + 2*caMit*mgMit*mgMit/(pow(g, 4)*K_CC*K_CC*K_MM*K_MM);
	double dFlux_caMit = -2*k_o*caMit/(K_CC*K_CC)/D - dD_caMit*(k_i*caCyt*caCyt/(K_CC*K_CC) - k_o*caMit*caMit/(K_CC*K_CC))/(D*D);
	dFlux_caMit *= 1e-9 * m_mit_volume / m_mit_surface * 1e-9;

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = dFlux_caCyt;
		i++;
	}
	if (!has_constant_value(_CMIT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CMIT_);
		flux_derivs[0][i].second = dFlux_caMit;
		i++;
	}
}


const size_t MCU::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_CCYT_)) n--;
	if (has_constant_value(_CMIT_)) n--;

	return n;
}


size_t MCU::n_fluxes() const
{
	return 1;
};


const std::pair<size_t,size_t> MCU::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) from = local_fct_index(_CCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CMIT_)) to = local_fct_index(_CMIT_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


const std::string MCU::name() const
{
	return std::string("MCU");
};


void MCU::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CMIT_))
	{
		UG_THROW("Supplying neither cytosolic nor mitochondrial calcium concentrations is not allowed.\n"
			 	 "This would mean that the flux calculation would be of no consequence\n"
				 "and this pump mechanism would not do anything.");
	}
}


void MCU::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      [Ca_cyt]  M (= mol/dm^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_mit]  M (= mol/dm^3)                                                |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


void MCU::set_mit_volume(number mit_volume)
{
	m_mit_volume  = mit_volume;
}


void MCU::set_mit_surface(number mit_surface)
{
	m_mit_surface = mit_surface;
}


void MCU::set_pi_cyt(number pi_cyt)
{
	m_pi_cyt = pi_cyt;

	K_CC 	= K_C * (1 + m_pi_cyt/(K_Pi+m_pi_cyt));
	K_MM 	= K_M / (1 + m_pi_cyt/(K_Pi+m_pi_cyt));
}


void MCU::set_psi(number psi)
{
	m_psi = psi;
}


void MCU::set_mg_cyt(number mg_cyt)
{
	m_mg_cyt = mg_cyt;
}


void MCU::set_mg_mit(number mg_mit)
{
	m_mg_mit = mg_mit;
}

void MCU::set_rate_constant(number k)
{
	m_k = k;
}


} // namespace neuro_collection
} // namespace ug


