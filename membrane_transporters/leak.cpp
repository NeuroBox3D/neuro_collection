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

#include "leak.h"


namespace ug{
namespace neuro_collection{


Leak::Leak(const std::vector<std::string>& fcts)
: IMembraneTransporter(fcts),
  m_perm(1.0), m_temp(298.15), m_z(1), m_bNoVoltage(false)
{
	if (!is_supplied(_PHIS_) && !is_supplied(_PHIT_))
		m_bNoVoltage = true;
}

Leak::Leak(const char* fcts)
: IMembraneTransporter(fcts),
  m_perm(1.0), m_temp(298.15), m_z(1), m_bNoVoltage(false)
{
	if (!is_supplied(_PHIS_) && !is_supplied(_PHIT_))
		m_bNoVoltage = true;
}

Leak::~Leak()
{
	// nothing to do
}

void Leak::set_permeability(number p)
{
	m_perm = p;
}


void Leak::set_temperature(number t)
{
	m_temp = t;
}

void Leak::set_valency(int v)
{
	m_z = v;
}


void Leak::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	number source_conc = u[_S_];	// source concentration
	number target_conc = u[_T_];	// target concentration
	if (m_bNoVoltage)
	{
		flux[0] = m_perm * (source_conc - target_conc);
		return;
	}

	number source_pot = u[_PHIS_];	// source potential
	number target_pot = u[_PHIT_];	// target potential

	const number zfrt = m_z*96485.0 / (8.31451 * m_temp);
	const number v =  target_pot - source_pot;
	if (fabs(v) < 1e-8)
		flux[0] = m_perm * (source_conc - target_conc) * (1.0 - 0.5*v*zfrt);
	else
		flux[0] = - m_perm * zfrt * v * (source_conc - target_conc * exp(zfrt*v)) / (1.0 - exp(zfrt*v));
}


void Leak::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	if (m_bNoVoltage)
	{
		size_t i = 0;
		if (!has_constant_value(_S_))
		{
			flux_derivs[0][i].first = local_fct_index(_S_);
			flux_derivs[0][i].second = m_perm;
			++i;
		}
		if (!has_constant_value(_T_))
		{
			flux_derivs[0][i].first = local_fct_index(_T_);
			flux_derivs[0][i].second = -m_perm;
		}

		return;
	}

	number source_conc = u[_S_];	// source concentration
	number target_conc = u[_T_];	// target concentration
	number source_pot = u[_PHIS_];	// source potential
	number target_pot = u[_PHIT_];	// target potential

	const number zfrt = m_z*96485.0 / (8.31451 * m_temp);
	const number v =  target_pot - source_pot;

	number dcs = 0.0;
	number dct = 0.0;
	number dps = 0.0;
	number dpt = 0.0;
	if (fabs(v) < 1e-8)
	{
		dcs = m_perm * (1.0 - 0.5*v*zfrt);
		dct = - dcs;
		dps = 0.5 * m_perm * (source_conc - target_conc) * zfrt;
		dpt = -dps;
	}
	else
	{
		dcs = - m_perm * zfrt * v / (1.0 - exp(zfrt*v));
		dct = - m_perm * zfrt * v / (1.0 - exp(-zfrt*v));

		number in = zfrt*v;
		number ex = exp(in);
		dps = m_perm * zfrt * (source_conc*(1.0 - (1.0 - in)*ex) + target_conc*ex*(ex - 1.0 - in))
						/ ((1.0 - ex) * (1.0 - ex));
		dpt = - dps;
	}

	size_t i = 0;
	if (!has_constant_value(_S_))
	{
		flux_derivs[0][i].first = local_fct_index(_S_);
		flux_derivs[0][i].second = dcs;
		++i;
	}
	if (!has_constant_value(_T_))
	{
		flux_derivs[0][i].first = local_fct_index(_T_);
		flux_derivs[0][i].second = dct;
		++i;
	}
	if (!has_constant_value(_PHIS_))
	{
		flux_derivs[0][i].first = local_fct_index(_PHIS_);
		flux_derivs[0][i].second = dps;
		++i;
	}
	if (!has_constant_value(_PHIT_))
	{
		flux_derivs[0][i].first = local_fct_index(_PHIT_);
		flux_derivs[0][i].second = dpt;
	}

	return;
}


// return number of unknowns this transport mechanism depends on
size_t Leak::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_S_))
		--n;
	if (has_constant_value(_T_))
		--n;

	if (!m_bNoVoltage)
	{
		n += 2;
		if (has_constant_value(_PHIS_))
			--n;
		if (has_constant_value(_PHIT_))
			--n;
	}

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
	if (is_supplied(_S_)) from = local_fct_index(_S_); else from = InnerBoundaryConstants::_IGNORE_;
	if (is_supplied(_T_)) to = local_fct_index(_T_); else to = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}

const std::string Leak::name() const
{
	return std::string("Leak");
}

void Leak::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!is_supplied(_S_) && !is_supplied(_T_))
	{
		UG_THROW("Supplying neither source nor target concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this leak would not do anything.");
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
	UG_LOG("|      Potentials         V                                                    |"<< std::endl);
	UG_LOG("|      Permeability       m/s                                                  |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output (in MembraneTransportFV1)                                          |"<< std::endl);
	UG_LOG("|      flux DENSITY       mol/(m^2 s)                                          |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


} // namespace neuro_collection
} // namespace ug

