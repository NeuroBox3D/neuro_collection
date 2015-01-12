/*
 *	Discretization for the leakage flux through a membrane
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#include "leak.h"


namespace ug{
namespace neuro_collection{


Leak::Leak(std::vector<std::string> fcts) : IMembraneTransporter(fcts)
{
	// nothing to do
}

Leak::~Leak()
{
	// nothing to do
}

void Leak::calc_flux(const std::vector<number>& u, std::vector<number>& flux) const
{
	number source_conc = u[_S_];	// source concentration
	number target_conc = u[_T_];	// target concentration

	// membrane current corresponding to diffusion pressure
	// cheating a little here: utilizing density for leakage flux "density" constant,
	// since the leakage flux is not really caused by a density, but depends of course
	// on the densities of all the pumps/channels in the membrane and is therefore
	// position/subset dependent

	// the actual flux density is flux[0] * density_fct
	flux[0] = source_conc-target_conc;
}


void Leak::calc_flux_deriv(const std::vector<number>& u, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	size_t i = 0;
	if (!has_constant_value(_S_))
	{
		flux_derivs[0][i].first = local_fct_index(_S_);
		flux_derivs[0][i].second = 1.0;
		i++;
	}
	if (!has_constant_value(_T_))
	{
		flux_derivs[0][i].first = local_fct_index(_T_);
		flux_derivs[0][i].second = -1.0;
		i++;
	}
}


// return number of unknowns this transport mechanism depends on
const size_t Leak::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_S_))
		n--;
	if (has_constant_value(_T_))
		n--;

	return n;
}

// return number of fluxes calculated by this machanism
size_t Leak::n_fluxes() const
{
	return 1;
};

// from where to where do the fluxes occur
const std::pair<size_t,size_t> Leak::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_S_)) from = local_fct_index(_S_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_T_)) to = local_fct_index(_T_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}

const std::string Leak::name() const
{
	return std::string("Leak");
}

void Leak::check_constant_allowed(const size_t i, const number val) const
{
	// Check that not both, inner and outer calcium concentrations are set constant;
	// in that case, calculation of a flux would be of no consequence.
	if ((has_constant_value(_S_) && i == _T_)
		|| (has_constant_value(_T_) && i == _S_))
	{
		UG_THROW("It is not allowed to set both, the source and the\n"
				"target concentrations to a constant value.\n"
				"This would mean that the flux calculation would be of\n"
				"no consequence and this leak would not do anything.");
	}
}


void Leak::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      Concentrations     mM (= mol/m^3)                                       |"<< std::endl);
	UG_LOG("|      Leakage constant   m/s                                                  |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output (in TwoSidedMembraneTransport)                                     |"<< std::endl);
	UG_LOG("|      flux DENSITY       mol/(m^2 s)                                          |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


} // namespace neuro_collection
} // namespace ug

