/*
 *	 Discretization for the mitochondrial uniporter MCU in the mitochondrial membrane
 *
 *  Created on: 30.06.2015
 *      Author: mstepnie
 */

#include "mcu.h"

namespace ug{
namespace neuro_collection{


MCU::MCU(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts)
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
}

MCU::MCU(const char* fcts) : IMembraneTransporter(fcts)
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
}


MCU::~MCU()
{
	// nothing to do
}


void MCU::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
// 	get values of the unknowns in associated node
//	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
//	number caExt = u[_CEXT_];	// extracellular Ca2+ concentration

//	IMPLEMENT ME
}


void MCU::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
// 	get values of the unknowns in associated node
//	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
//	number caExt = u[_CEXT_];	// extracellular Ca2+ concentration

//	IMPLEMENT ME
}


const size_t MCU::n_dependencies() const
{
//	Todo
	return 1;
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


