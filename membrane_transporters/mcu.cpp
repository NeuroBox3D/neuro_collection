/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Martin Stepniewski
 * Creation date: 2015-06-30
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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


size_t MCU::n_dependencies() const
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

number MCU::get_flux(number ca_cyt, number ca_mit, number pi, number mg_cyt, number mg_mit, number psi)
{
// 	get values of the unknowns in associated node
	number caCyt = ca_cyt;	// cytosolic Ca2+ concentration
	number caMit = ca_mit;	// mitochondrial Ca2+ concentration

	number mgCyt = mg_cyt;
	number mgMit = mg_mit;

	double KC 	= K_C * (1 + pi/(K_Pi+pi));
	double KM 	= K_M / (1 + pi/(K_Pi+pi));

	double phi = 2.0 * F * psi / RT;

	double beta_e = 0.5 * (1 + nH/phi * log((phi/nH) / (sinh(phi/nH))));
	double beta_x = 0.5 * (1 - nH/phi * log((phi/nH) / (sinh(phi/nH))));

	double k_i = m_k * exp(2*beta_e*phi);
	double k_o = m_k * exp(-2*beta_x*phi);

	double D = 	1 + ((caCyt*caCyt)/(KC*KC)) + ((caMit*caMit)/(KC*KC)) +
					((mgCyt*mgCyt)/(KM*KM)) + ((mgMit*mgMit)/(KM*KM)) +
					((caCyt*caCyt*mgCyt*mgCyt)/(pow(g, 4)*KC*KC*KM*KM)) +
					((caMit*caMit*mgMit*mgMit)/(pow(g, 4)*KC*KC*KM*KM));

	double flux = 1/D * (k_i*((caCyt*caCyt)/(KC*KC)) - k_o*((caMit*caMit)/(KC*KC)));	// in nmol/mg/s
	return flux;
}


} // namespace neuro_collection
} // namespace ug


