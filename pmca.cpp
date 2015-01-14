/*
 *	 Discretization for the PMCA pump in the PM
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#include "pmca.h"

namespace ug{
namespace neuro_collection{


PMCA::PMCA(std::vector<std::string> fcts) : IMembraneTransporter(fcts),
KD_P(6.0e-5), IMAX_P(1.7e-23)
{
	// nothing to do
}


PMCA::~PMCA()
{
	// nothing to do
}


void PMCA::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	// get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration

	number gatingFactor = caCyt*caCyt / (KD_P*KD_P + caCyt*caCyt);

	flux[0] = gatingFactor * IMAX_P;
}


void PMCA::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration

	number dGating_dCyt = 2*KD_P*KD_P*caCyt / std::pow(KD_P*KD_P + caCyt*caCyt, 2);

	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][0].first = local_fct_index(_CCYT_);
		flux_derivs[0][0].second = dGating_dCyt * IMAX_P;
	}
}


const size_t PMCA::n_dependencies() const
{
	return 1;
}


size_t PMCA::n_fluxes() const
{
	return 1;
};


const std::pair<size_t,size_t> PMCA::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) from = local_fct_index(_CCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CEXT_)) to = local_fct_index(_CEXT_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


const std::string PMCA::name() const
{
	return std::string("PMCA");
};


void PMCA::check_supplied_functions() const
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


void PMCA::print_units() const
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


