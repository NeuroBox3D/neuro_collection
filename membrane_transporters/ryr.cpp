/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2011-12-20
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

#include "ryr.h"


namespace ug {
namespace neuro_collection {


RyR::RyR(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
KA(5.21e13), KB(3.89e9), KC(17.5), MU_RYR(5.0e-11),
REF_CA_ER(2.5e-1)
{
	// nothing to do
};

RyR::RyR(const char* fcts) : IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
KA(5.21e13), KB(3.89e9), KC(17.5), MU_RYR(5.0e-11),
REF_CA_ER(2.5e-1)
{
	// nothing to do
};


RyR::~RyR()
{
	// nothing to do
};


void RyR::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration

	// membrane current corresponding to diffusion pressure
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);

	// open probability
	number pOpen =  (1.0 + KB*pow(caCyt,3)) / (1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3));

	flux[0] = pOpen * current;
}


void RyR::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration

	// membrane potential
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);
	number dI_dCyt = - R*T/(4*F*F) * MU_RYR/REF_CA_ER;
	number dI_dER = -dI_dCyt;

	// RyR flux derivs
	number schlonz1 = 1.0 + KB*pow(caCyt,3);
	number schlonz2 = 1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3);
	number pOpen = schlonz1 / schlonz2;

	number dOpen_dCyt = (3.0*KB*caCyt*caCyt + schlonz1/schlonz2*(4.0/(KA*pow(caCyt,5)) - 3.0*KB*caCyt*caCyt)) / schlonz2;

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = dOpen_dCyt * current + pOpen * dI_dCyt;
		i++;
	}
	if (!has_constant_value(_CER_))
	{
		flux_derivs[0][i].first = local_fct_index(_CER_);
		flux_derivs[0][i].second = pOpen * dI_dER;
		i++;
	}
}


size_t RyR::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_CCYT_))
		n--;
	if (has_constant_value(_CER_))
		n--;

	return n;
}


size_t RyR::n_fluxes() const
{
	return 1;
};


const std::pair<size_t,size_t> RyR::flux_from_to(size_t flux_i) const
{
    size_t from, to;
    if (is_supplied(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
    if (is_supplied(_CER_)) from = local_fct_index(_CER_); else from = InnerBoundaryConstants::_IGNORE_;
    
    return std::pair<size_t, size_t>(from, to);
}


const std::string RyR::name() const
{
	return std::string("RyR");
};


void RyR::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!is_supplied(_CCYT_) && !is_supplied(_CER_))
	{
		UG_THROW("Supplying neither cytosolic nor endoplasmic calcium concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}
}


void RyR::print_units() const
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
