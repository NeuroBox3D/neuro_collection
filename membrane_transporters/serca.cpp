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

#include "serca.h"

namespace ug {
namespace neuro_collection {


SERCA::SERCA(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts),
VS(6.5e-24), KS(1.8e-4), m_bLinearizedAssembling(false)
{
	// nothing to do
}


SERCA::SERCA(const char* fcts) : IMembraneTransporter(fcts),
VS(6.5e-24), KS(1.8e-4), m_bLinearizedAssembling(false)
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

	//UG_COND_THROW(flux[0] != flux[0],
	//	"SERCA NaN: caCyt = " << caCyt << ", caER = " << caER);
}


void SERCA::calc_flux
(
	const std::vector<number>& u,
	const std::vector<number>& uOld,
	GridObject* e,
	std::vector<number>& flux
) const
{
	flux[0] = VS * u[_CCYT_] / ((KS + uOld[_CCYT_]) * uOld[_CER_]);
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


void SERCA::calc_flux_deriv
(
	const std::vector<number>& u,
	const std::vector<number>& uOld,
	GridObject* e,
	std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs
) const
{
	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][0].first = local_fct_index(_CCYT_);
		flux_derivs[0][0].second = VS / ((KS + uOld[_CCYT_]) * uOld[_CER_]);
		++i;
	}
	if (!has_constant_value(_CER_))
	{
		flux_derivs[0][i].first = local_fct_index(_CER_);
		flux_derivs[0][i].second = 0.0;
		++i;
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
	if (is_supplied(_CCYT_)) from = local_fct_index(_CCYT_); else from = InnerBoundaryConstants::_IGNORE_;
	if (is_supplied(_CER_)) to = local_fct_index(_CER_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


const std::string SERCA::name() const
{
	return std::string("SERCA");
}


bool SERCA::needsPreviousSolution() const
{
	return m_bLinearizedAssembling;
}


void SERCA::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!is_supplied(_CCYT_) && !is_supplied(_CER_))
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


void SERCA::set_linearized_assembling()
{
	m_bLinearizedAssembling = true;
}




} // namespace neuro_collection
} // namespace ug
