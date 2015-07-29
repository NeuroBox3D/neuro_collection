/*
 *	 Discretization for the mitochondrial NCX in the mitochondrial membrane
 *
 *  Created on: 30.06.2015
 *      Author: mstepnie
 */

/*
 * Complex NCX model from Paper: "A Biophysically Based Mathematical Model for the Kinetics of Mitochondrial
 * Na-Ca Antiporter" (using Model 1 and kinectic values of reference 1), Pradhan et al. 2010
 *
 */

#include "mncx.h"

namespace ug{
namespace neuro_collection{


MNCX::MNCX(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts),
F(0.096484), RT(2.5775),
K_C(2.28e-9), K_N(9.14e-3), k(4.9)
//K_C(2.1e-6), K_N(8.2e-3), k(0.0846) // Pradhan 2008
{
	m_psi 	= 0.0;

	m_mit_volume  = 0.0;
	m_mit_surface = 0.0;
}

MNCX::MNCX(const char* fcts) : IMembraneTransporter(fcts),
F(0.096484), RT(2.5775),
K_C(2.28e-9), K_N(9.14e-3), k(4.9)
//K_C(2.1e-6), K_N(8.2e-3), k(0.0846) // // Pradhan 2008
{
	m_psi 	= 0.0;

	m_mit_volume  = 0.0;
	m_mit_surface = 0.0;
}


MNCX::~MNCX()
{
	// nothing to do
}


void MNCX::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	if(m_mit_volume == 0.0 || m_mit_surface == 0.0)
		UG_THROW("ERROR in MNCX membrane transport: mitochondrial volume or surface not specified.");

// 	get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caMit = u[_CMIT_];	// mitochondrial Ca2+ concentration

	number naCyt = u[_NCYT_];	// cytosolic Na+ concentration
	number naMit = u[_NMIT_];	// mitochondrial Na2+ concentration

	double phi = F * m_psi / RT;

	double k_i  = k * exp(0.5*phi);
	double k_o  = k * exp(-0.5*phi);

	double in_term  = pow(naCyt, 3)/pow(K_N, 3) + caMit/K_C + caMit*pow(naCyt, 3)/(K_C*(pow(K_N, 3)));
	double out_term = pow(naMit, 3)/pow(K_N, 3) + caCyt/K_C + caCyt*pow(naMit, 3)/(K_C*(pow(K_N, 3)));

	double D = 1 + in_term + out_term;

	double flux_in =  caMit*pow(naCyt, 3)/(K_C*pow(K_N, 3)) * k_i;
	double flux_out = caCyt*pow(naMit, 3)/(K_C*pow(K_N, 3)) * k_o;

	double fluxes = 1/D*(flux_in - flux_out);

	UG_LOG("fluxes: " << fluxes << std::endl);

//	Transform original flux umol/mg/min to mitochondrial flux umol/s
//	1um^3 mitochondrial volume = 1e-9mg mitochondrial protein
	fluxes *= 1e-9 * m_mit_volume;

//  Transform mitochondrial flux umol/min to flux umol/um^2/s
	fluxes *= 1.0 / m_mit_surface;

//  Transform mitochondrial flux umol/um^2/min to mol/um^2/min
	fluxes *= 1e-6;

//  Transform mitochondrial flux mol/um^2/min to mol/um^2/s
	fluxes *= (1.0/60.0);

	flux[0] = fluxes;
	flux[1] = fluxes*3;
}


void MNCX::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	if(m_mit_volume == 0.0 || m_mit_surface == 0.0)
		UG_THROW("ERROR in MNCX membrane transport: mitochondrial volume or surface not specified.");

// 	get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caMit = u[_CMIT_];	// mitochondrial Ca2+ concentration

	number naCyt = u[_NCYT_];	// cytosolic Na+ concentration
	number naMit = u[_NMIT_];	// mitochondrial Na2+ concentration

	double phi = F * m_psi / RT;

	double k_i  = k * exp(0.5*phi);
	double k_o  = k * exp(-0.5*phi);

	double in_term  = pow(naCyt, 3)/pow(K_N, 3) + caMit/K_C + caMit*pow(naCyt, 3)/(K_C*(pow(K_N, 3)));
	double out_term = pow(naMit, 3)/pow(K_N, 3) + caCyt/K_C + caCyt*pow(naMit, 3)/(K_C*(pow(K_N, 3)));

	double D = 1 + in_term + out_term;

	double flux_in =  caMit*pow(naCyt, 3)/(K_C*pow(K_N, 3)) * k_i;
	double flux_out = caCyt*pow(naMit, 3)/(K_C*pow(K_N, 3)) * k_o;

//  derivatives of caCyt
	double dD_caCyt = 1.0/K_C + pow(naMit, 3)/(K_C*pow(K_N, 3));
	double dFlux_caCyt = -k_o*pow(naMit, 3)/(K_C*pow(K_N, 3))/D - ((flux_in-flux_out)*dD_caCyt)/(D*D);
	dFlux_caCyt *= 1e-9 * m_mit_volume / m_mit_surface * 1e-6 / 60.0;

//  derivatives of caMit
	double dD_caMit = 1.0/K_C + pow(naCyt, 3)/(K_C*pow(K_N, 3));
	double dFlux_caMit = k_i*pow(naCyt, 3)/(K_C*pow(K_N, 3))/D - ((flux_in-flux_out)*dD_caMit)/(D*D);
	dFlux_caMit *= 1e-9 * m_mit_volume / m_mit_surface * 1e-6 / 60.0;

//  derivatives of naCyt
	double dD_naCyt = 1.0/pow(K_N, 3) + 3*pow(naCyt, 2)*caMit/(K_C*pow(K_N, 3));
	double dFlux_naCyt = k_i*caMit*3*pow(naCyt, 2)/(K_C*pow(K_N, 3))/D - ((flux_in-flux_out)*dD_naCyt)/(D*D);
	dFlux_naCyt *= 3 * 1e-9 * m_mit_volume / m_mit_surface * 1e-6 / 60.0;

//  derivatives of naMit
	double dD_naMit = 1.0/pow(K_N, 3) + 3*pow(naMit, 2)*caCyt/(K_C*pow(K_N, 3));
	double dFlux_naMit = -k_o*caCyt*3*pow(naMit, 2)/(K_C*pow(K_N, 3))/D - ((flux_in-flux_out)*dD_naMit)/(D*D);
	dFlux_naMit *= 3 * 1e-9 * m_mit_volume / m_mit_surface * 1e-6 / 60.0;

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
	if (!has_constant_value(_NCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_NCYT_);
		flux_derivs[0][i].second = dFlux_naCyt;
		i++;
	}
	if (!has_constant_value(_NMIT_))
	{
		flux_derivs[0][i].first = local_fct_index(_NMIT_);
		flux_derivs[0][i].second = dFlux_naMit;
		i++;
	}
}


const size_t MNCX::n_dependencies() const
{
	size_t n = 4;
	if (has_constant_value(_CCYT_)) n--;
	if (has_constant_value(_CMIT_)) n--;
	if (has_constant_value(_NCYT_)) n--;
	if (has_constant_value(_NMIT_)) n--;

	return n;
}


size_t MNCX::n_fluxes() const
{
	return 2;
};


const std::pair<size_t,size_t> MNCX::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CMIT_)) from = local_fct_index(_CMIT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;

	if (allows_flux(_NCYT_)) from = local_fct_index(_NCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_NMIT_)) to = local_fct_index(_NMIT_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


const std::string MNCX::name() const
{
	return std::string("MNCX");
};


void MNCX::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CMIT_))
	{
		UG_THROW("Supplying neither cytosolic nor mitochondrial calcium concentrations is not allowed.\n"
			 	 "This would mean that the flux calculation would be of no consequence\n"
				 "and this pump mechanism would not do anything.");
	}

	if (!allows_flux(_NCYT_) && !allows_flux(_NMIT_))
	{
		UG_THROW("Supplying neither cytosolic nor mitochondrial sodium concentrations is not allowed.\n"
			 	 "This would mean that the flux calculation would be of no consequence\n"
				 "and this pump mechanism would not do anything.");
	}
}


void MNCX::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      [Ca_cyt]  M (= mol/dm^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_out]  M (= mol/dm^3)                                                |"<< std::endl);
	UG_LOG("|      [Na_cyt]  M (= mol/dm^3)                                                |"<< std::endl);
	UG_LOG("|      [Na_out]  M (= mol/dm^3)                                                |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("|      Na flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


void MNCX::set_mit_volume(number mit_volume)
{
	m_mit_volume  = mit_volume;
}


void MNCX::set_mit_surface(number mit_surface)
{
	m_mit_surface = mit_surface;
}


void MNCX::set_psi(number psi)
{
	m_psi = psi;
}

number MNCX::get_flux(number ca_cyt, number ca_mit, number na_cyt, number na_mit, number psi)
{
// 	get values of the unknowns in associated node
	number caCyt = ca_cyt;	// cytosolic Ca2+ concentration
	number caMit = ca_mit;	// mitochondrial Ca2+ concentration

	number naCyt = na_cyt;	// cytosolic Na+ concentration
	number naMit = na_mit;	// mitochondrial Na2+ concentration

	double phi = F * psi / RT;

	double k_i  = k * exp(0.5*phi);
	double k_o  = k * exp(-0.5*phi);

	double in_term  = pow(naCyt, 3)/pow(K_N, 3) + caMit/K_C + caMit*pow(naCyt, 3)/(K_C*(pow(K_N, 3)));
	double out_term = pow(naMit, 3)/pow(K_N, 3) + caCyt/K_C + caCyt*pow(naMit, 3)/(K_C*(pow(K_N, 3)));

	double D = 1 + in_term + out_term;

	double flux_in =  caMit*pow(naCyt, 3)/(K_C*pow(K_N, 3)) * k_i;
	double flux_out = caCyt*pow(naMit, 3)/(K_C*pow(K_N, 3)) * k_o;

	double fluxes = 1/D*(flux_in - flux_out);

	return fluxes;
}


} // namespace neuro_collection
} // namespace ug


