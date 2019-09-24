/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-11-01
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
	const number vm = u[_PHII_] - u[_PHIO_];

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
	if (is_supplied(_PHII_)) from = local_fct_index(_PHII_); else from = InnerBoundaryConstants::_IGNORE_;
	if (is_supplied(_PHIO_)) to = local_fct_index(_PHIO_); else to = InnerBoundaryConstants::_IGNORE_;

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
	if (!is_supplied(_PHII_) && !is_supplied(_PHIO_))
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






OhmicLeakageCharges::OhmicLeakageCharges(const std::vector<std::string>& fcts)
: IMembraneTransporter(fcts), m_g(3.0), m_eL(-0.0546)
{
	// nothing to do
}

OhmicLeakageCharges::OhmicLeakageCharges(const char* fcts)
: IMembraneTransporter(fcts), m_g(3.0), m_eL(-0.0546)
{
	// nothing to do
}

OhmicLeakageCharges::~OhmicLeakageCharges()
{
	// nothing to do
}


void OhmicLeakageCharges::set_conductance(number g)
{
	m_g = g;
}


void OhmicLeakageCharges::set_reversal_potential(number el)
{
	m_eL = el;
}


void OhmicLeakageCharges::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	const number vm = u[_PHII_] - u[_PHIO_];

	flux[0] = m_g * (vm - m_eL);
}


void OhmicLeakageCharges::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
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
size_t OhmicLeakageCharges::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_PHII_))
		n--;
	if (has_constant_value(_PHIO_))
		n--;

	return n;
}

// return number of fluxes calculated by this machanism
size_t OhmicLeakageCharges::n_fluxes() const
{
	return 1;
};

// from where to where do the fluxes occur
const std::pair<size_t,size_t> OhmicLeakageCharges::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (is_supplied(_RHOI_)) from = local_fct_index(_RHOI_); else from = InnerBoundaryConstants::_IGNORE_;
	if (is_supplied(_RHOO_)) to = local_fct_index(_RHOO_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}

const std::string OhmicLeakageCharges::name() const
{
	return std::string("OhmicLeakageCharges");
}

void OhmicLeakageCharges::check_supplied_functions() const
{
	// Check that not both, inner and outer charge densities are not supplied;
	// in that case, calculation of a current would be of no consequence.
	if (!is_supplied(_RHOI_) && !is_supplied(_RHOO_))
	{
		UG_THROW("Supplying neither source nor target charge densities is not allowed.\n"
				"This would mean that the current calculation would be of no consequence\n"
				"and this leak would not do anything.");
	}
}


void OhmicLeakageCharges::print_units() const
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

