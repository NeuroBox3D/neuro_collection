/*
 *	Discretization for the SERCA calcium pump in the ER membrane
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#include "serca.h"

namespace ug{
namespace neuro_collection{
    
    
SERCA::SERCA(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts),
VS(6.5e-24), KS(1.8e-4)
{
	// nothing to do
}


SERCA::SERCA(const char* fcts) : IMembraneTransporter(fcts),
VS(6.5e-24), KS(1.8e-4)
{
	// nothing to do
}


SERCA::~SERCA()
{
	// nothing to do
}


void SERCA::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	// get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration

	flux[0] = VS*caCyt / ((KS + caCyt) * caER);
}


void SERCA::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = (VS*KS) / (pow(KS+caCyt,2)*caER);
		i++;
	}
	if (!has_constant_value(_CER_))
	{
		flux_derivs[0][i].first = local_fct_index(_CER_);
		flux_derivs[0][i].second = - VS*caCyt / ((KS+caCyt)*caER*caER);
		i++;
	}
}


size_t SERCA::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_CCYT_))
		n--;
	if (has_constant_value(_CER_))
		n--;

	return n;
}


size_t SERCA::n_fluxes() const
{
	return 1;
}


const std::pair<size_t,size_t> SERCA::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) from = local_fct_index(_CCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CER_)) to = local_fct_index(_CER_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


const std::string SERCA::name() const
{
	return std::string("SERCA");
}


void SERCA::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CER_))
	{
		UG_THROW("Supplying neither cytosolic nor endoplasmic calcium concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this pump would not do anything.");
	}
}


void SERCA::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      [Ca_cyt]  mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|      [Ca_er]   mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}
    
    
} // namespace neuro_collection
} // namespace ug
