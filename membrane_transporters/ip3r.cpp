/*
 *	Discretization for the IP3R calcium channel in the ER membrane
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#include "ip3r.h"


namespace ug{
namespace neuro_collection{


IP3R::IP3R(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
D1(1.3e-4), D2(1.05e-3), D3(9.4e-4), D5(8.23e-5), MU_IP3R(1.6e-12),
REF_CA_ER(2.5e-1)
{
	// nothing to do
}


IP3R::IP3R(const char* fcts) : IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
D1(1.3e-4), D2(1.05e-3), D3(9.4e-4), D5(8.23e-5), MU_IP3R(1.6e-12),
REF_CA_ER(2.5e-1)
{
	// nothing to do
}


IP3R::~IP3R()
{
	// nothing to do
}


void IP3R::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration
	number ip3 = u[_IP3_];		// IP3 concentration

	// TODO: Use single-channel current directly.
	// A seemingly quite good estimate on how to do this is given in Bezprozvanny and Ehrlich, 1994
	// This would be modeled by:
	// I_ca = 1/(2F) * I_el_ref * (caER - caCyt) / REF_CA_ER

	// membrane current corresponding to diffusion pressure
	number current = R*T/(4*F*F) * MU_IP3R/REF_CA_ER * (caER - caCyt);

	// open probability
	number pOpen = pow(caCyt*ip3*D2 / ((caCyt*ip3 + ip3*D2 + D1*D2 + caCyt*D3) * (caCyt+D5)),3);

	flux[0] = pOpen * current;
}


void IP3R::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration
	number ip3 = u[_IP3_];		// IP3 concentration

	// membrane potential
	number current = R*T/(4*F*F) * MU_IP3R/REF_CA_ER * (caER - caCyt);
	number dI_dCyt = - R*T/(4*F*F) * MU_IP3R/REF_CA_ER;
	number dI_dER = -dI_dCyt;

	// IP3R flux derivatives
	number schlonz1 = caCyt*ip3 + ip3*D2 + D1*D2 + caCyt*D3;
	number schlonz2 = schlonz1 * (caCyt+D5);
	number schlonz3 = (caCyt*ip3*D2) / schlonz2;
	number pOpen = pow(schlonz3,3);

	number dOpen_dCyt = 3.0*schlonz3*schlonz3*ip3*D2 * (1.0 - caCyt/schlonz2 * ( (ip3+D3)*(caCyt+D5) + schlonz1 )) / schlonz2;
	number dOpen_dIP3 = 3.0*schlonz3*schlonz3*caCyt*D2 * (1.0 - ip3/schlonz2 * (caCyt+D2)*(caCyt+D5)) / schlonz2;

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
	if (!has_constant_value(_IP3_))
	{
		flux_derivs[0][i].first = local_fct_index(_IP3_);
		flux_derivs[0][i].second = dOpen_dIP3 * current;
	}
}


size_t IP3R::n_dependencies() const
{
	size_t n = 3;
	for (size_t i = 0; i < 3; i++)
	{
		if (has_constant_value(i))
			n--;
	}

	return n;
}


size_t IP3R::n_fluxes() const
{
	return 1;
}


const std::pair<size_t,size_t> IP3R::flux_from_to(size_t flux_i) const
{
	size_t from, to;
	if (allows_flux(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
	if (allows_flux(_CER_)) from = local_fct_index(_CER_); else from = InnerBoundaryConstants::_IGNORE_;

	return std::pair<size_t, size_t>(from, to);
}


const std::string IP3R::name() const
{
	return std::string("IP3R");
}


void IP3R::check_supplied_functions() const
{
	// Check that not both, inner and outer calcium concentrations are not supplied;
	// in that case, calculation of a flux would be of no consequence.
	if (!allows_flux(_CCYT_) && !allows_flux(_CER_))
	{
		UG_THROW("Supplying neither cytosolic nor endoplasmic calcium concentrations is not allowed.\n"
				"This would mean that the flux calculation would be of no consequence\n"
				"and this channel would not do anything.");
	}
}


void IP3R::print_units() const
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
	UG_LOG("|      [IP3]     mM (= mol/m^3)                                                |"<< std::endl);
	UG_LOG("|                                                                              |"<< std::endl);
	UG_LOG("|    Output                                                                    |"<< std::endl);
	UG_LOG("|      Ca flux   mol/s                                                         |"<< std::endl);
	UG_LOG("+------------------------------------------------------------------------------+"<< std::endl);
	UG_LOG(std::endl);
}


} // namespace neuro_collection
} // namespace ug

