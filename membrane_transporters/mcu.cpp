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
F(0.096484), RT(2.5775), ca_valence(2.0)
{
//	IMPLEMENT INITIALIZATION

	K_C 	= 3.965e-6 * 1e3;     // in mMol
	K_M 	= 0.655e-3 * 1e3; 	  // in mMol
	gamma 	= 3.9;
	//k = 16.51e-3 * 1e-6;		  // Reference (14) Vinogadrov and Scarpa in mmol/mg/s
	k 		= 1.27e-3 * 1e-6;	  // Reference  (9) Crompton et al. in mmol/mg/s
	nH 		= 2.65;
	K_Pi 	= 0.2e-3 * 1e3;       // in mMol
	K_CC 	= K_C; // * (1 + pi_cyt/(K_Pi+pi_cyt));
	K_MM 	= K_M; // / (1 + pi_cyt/(K_Pi+pi_cyt));
	m_psi 	= 0.0; //psi;

	mg_int = 0.0;
	mg_ext = 0.0;
}

MCU::MCU(const char* fcts) : IMembraneTransporter(fcts),
F(0.096484), RT(2.5775), ca_valence(2.0)
{
//	IMPLEMENT INITIALIZATION

	K_C 	= 3.965e-6 * 1e3;     // in mMol
	K_M 	= 0.655e-3 * 1e3; 	  // in mMol
	gamma 	= 3.9;
	//k = 16.51e-3 * 1e-6;		  // Reference (14) Vinogadrov and Scarpa in mmol/mg/s
	k 		= 1.27e-3 * 1e-6;	  // Reference  (9) Crompton et al. in mmol/mg/s
	nH 		= 2.65;
	K_Pi 	= 0.2e-3 * 1e3;       // in mMol
	K_CC 	= K_C; // * (1 + pi_cyt/(K_Pi+pi_cyt));
	K_MM 	= K_M; // / (1 + pi_cyt/(K_Pi+pi_cyt));
	m_psi 	= 0.0; //psi;

	mg_int = 0.0;
	mg_ext = 0.0;
}


MCU::~MCU()
{
	// nothing to do
}

// TODO right unit of flux
void MCU::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	// 	get values of the unknowns in associated node
	number caInt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caExt = u[_CEXT_];	// extracellular Ca2+ concentration

	// flux needs to be in M/s
	double phi = ca_valence * F * m_psi / RT;

	double beta_e = 0.5 * (1 + nH/phi * log((phi/nH) / (sinh(phi/nH))));
	double beta_x = 0.5 * (1 - nH/phi * log((phi/nH) / (sinh(phi/nH))));

	double k_o = k * exp(-2*beta_x*phi);
	double k_i = k * exp(2*beta_e*phi);

	double D = 	1 + ((caExt*caExt)/(K_CC*K_CC)) + ((caInt*caInt)/(K_CC*K_CC)) +
					((mg_ext*mg_ext)/(K_MM*K_MM)) + ((mg_int*mg_int)/(K_MM*K_MM)) +
					((caExt*caExt*mg_ext*mg_ext)/(pow(gamma, 4)*K_CC*K_CC*K_MM*K_MM)) +
					((caInt*caInt*mg_int*mg_int)/(pow(gamma, 4)*K_CC*K_CC*K_MM*K_MM));


	// flux in
	double fluxes = 1/D * (k_i*((caExt*caExt)/(K_CC*K_CC)) - k_o*((caInt*caInt)/(K_CC*K_CC)));

	flux[0]=(fluxes);
}


void MCU::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
// 	get values of the unknowns in associated node
	number caInt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caExt = u[_CEXT_];	// extracellular Ca2+ concentration

	//values neede for both derivations
	double phi = ca_valence * F * m_psi / RT;

	double beta_e = 0.5 * (1 + nH/phi * log((phi/nH) / (sinh(phi/nH))));
	double beta_x = 0.5 * (1 - nH/phi * log((phi/nH) / (sinh(phi/nH))));

	double k_o = k * exp(-2*beta_x*phi);
	double k_i = k * exp(2*beta_e*phi);

	double D = 	1 + ((caExt*caExt)/(K_CC*K_CC)) + ((caInt*caInt)/(K_CC*K_CC)) +
					((mg_ext*mg_ext)/(K_MM*K_MM)) + ((mg_int*mg_int)/(K_MM*K_MM)) +
					((caExt*caExt*mg_ext*mg_ext)/(pow(gamma, 4)*K_CC*K_CC*K_MM*K_MM)) +
					((caInt*caInt*mg_int*mg_int)/(pow(gamma, 4)*K_CC*K_CC*K_MM*K_MM));

	//deriv of caInt
	double dD_caInt = 2*caInt/(K_C*K_C) + 2*caInt*mg_int*mg_int/(pow(gamma, 4)*K_CC*K_CC*K_MM*K_MM);

	double dFlux_caInt = -2*k_o*caInt/(K_CC*K_CC)/D - dD_caInt*(k_i*caExt*caExt/(K_CC*K_CC) - k_o*caInt*caInt/(K_CC*K_CC))/(D*D);

	flux_derivs[0][0].first = local_fct_index(_CCYT_);
	flux_derivs[0][0].second = dFlux_caInt; //needs right flux values


	//deriv of caExt
	double dD_caExt = 2*caExt/(K_C*K_C) + 2*caExt*mg_ext*mg_ext/(pow(gamma, 4)*K_CC*K_CC*K_MM*K_MM);

	double dFlux_caExt = 2*k_i*caExt/(K_CC*K_CC)/D - dD_caExt*(k_i*caExt*caExt/(K_CC*K_CC) - k_o*caInt*caInt/(K_CC*K_CC))/(D*D);

	flux_derivs[0][1].first = local_fct_index(_CEXT_);
	flux_derivs[0][1].second = dFlux_caExt; //needs right flux values

}


const size_t MCU::n_dependencies() const
{
//	Todo
	return 2;
}


size_t MCU::n_fluxes() const
{
//	Todo
	return 1;
};


const std::pair<size_t,size_t> MCU::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) from = local_fct_index(_CCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CEXT_)) to = local_fct_index(_CEXT_); else to = InnerBoundaryConstants::_IGNORE_;

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
	if (!allows_flux(_CCYT_) && !allows_flux(_CEXT_))
	{
		UG_THROW("Supplying neither cytosolic nor extracellular calcium concentrations is not allowed.\n"
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
	UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_out]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


} // namespace neuro_collection
} // namespace ug


