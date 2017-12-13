/*
 *	Discretization for an ohmic leakage flux through a membrane
 *
 *  Created on: 2017-11-01
 *      Author: mbreit
 */

#include "leakage_ohmic.h"


namespace ug {
namespace neuro_collection {


OhmicLeakage::OhmicLeakage(const std::vector<std::string>& fcts)
: IMembraneTransporter(fcts), m_g(3.0), m_eL(-0.0546)
{
	// nothing to do
}

OhmicLeakage::OhmicLeakage(const char* fcts)
: IMembraneTransporter(fcts), m_g(3.0), m_eL(-0.0546)
{
	// nothing to do
}

OhmicLeakage::~OhmicLeakage()
{
	// nothing to do
}


void OhmicLeakage::set_conductance(number g)
{
	m_g = g;
}


void OhmicLeakage::set_reversal_potential(number el)
{
	m_eL = el;
}


void OhmicLeakage::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	const number vm = u[_PHII_] - u[_PHIO_];	// membrane potential

	// the actual flux density is flux[0] * density_fct
	flux[0] = m_g * (vm - m_eL);
}


void OhmicLeakage::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	size_t i = 0;
	if (!has_constant_value(_PHII_))
	{
		flux_derivs[0][i].first = local_fct_index(_PHII_);
		flux_derivs[0][i].second = m_g;
		i++;
	}
	if (!has_constant_value(_PHIO_))
	{
		flux_derivs[0][i].first = local_fct_index(_PHIO_);
		flux_derivs[0][i].second = -m_g;
	}
}


// return number of unknowns this transport mechanism depends on
size_t OhmicLeakage::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_PHII_))
		n--;
	if (has_constant_value(_PHIO_))
		n--;

	return n;
}

// return number of fluxes calculated by this machanism
size_t OhmicLeakage::n_fluxes() const
{
	return 1;
};

// from where to where do the fluxes occur
const std::pair<size_t,size_t> OhmicLeakage::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_RHOI_)) from = local_fct_index(_RHOI_); else from = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_RHOO_)) to = local_fct_index(_RHOO_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}

const std::string OhmicLeakage::name() const
{
	return std::string("OhmicLeakage");
}

void OhmicLeakage::check_supplied_functions() const
{
	// Check that not both, inner and outer charge densities are not supplied;
	// in that case, calculation of a current would be of no consequence.
	if (!allows_flux(_RHOI_) && !allows_flux(_RHOO_))
	{
		UG_THROW("Supplying neither source nor target charge densities is not allowed.\n"
				"This would mean that the current calculation would be of no consequence\n"
				"and this leak would not do anything.");
	}
}


void OhmicLeakage::print_units() const
{
	std::string nm = name();
	size_t n = nm.size();
	UG_LOG(std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG("|  Units used in the implementation of " << nm << std::string(n>=40?0:40-n, ' ') << "|" << std::endl);
	UG_LOG("|------------------------------------------------------------------------------|"<< std::endl);
	UG_LOG("|    Input                                                                     |"<< std::endl);
	UG_LOG("|      Potentials    V                                                         |"<< std::endl);
	UG_LOG("|      Conductance   C/(Vs) = A/V = S                                          |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output (in TwoSidedMembraneTransport)                                     |"<< std::endl);
	UG_LOG("|      Current       C/s = A                                                   |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}

} // namespace neuro_collection
} // namespace ug

