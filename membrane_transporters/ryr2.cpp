/*
 *	Discretization for the RyR calcium channel in the ER membrane
 *
 *  Created on: Nov 12, 2015
 *      Author: marcuskessler
 */

#include "ryr2.h"


namespace ug{
namespace neuro_collection{


RyR2::RyR2(const std::vector<std::string>& fcts) : IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
KAplus(5.21e13), KBplus(3.89e9), KCplus(17.5), MU_RYR(5.0e-11),
KAminus(5.21e13), KBminus(3.89e9), KCminus(17.5), REF_CA_ER(2.5e-1),
m_dom(approx->domain()), m_mg(m_dom->grid()), m_time(0.0), m_oldTime(0.0)
{
	if (m_mg->template has_attachment<side_t>(this->m_O1))
			UG_THROW("Attachment necessary for RyR channel dynamics "
					 "could not be made, since it already exists.");
		m_mg->template attach_to<side_t>(this->m_O1);

	m_aaO1 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_O1);
};

RyR2::RyR2(const char* fcts) : IMembraneTransporter(fcts),
R(8.314), T(310.0), F(96485.0),
KAplus(5.21e13), KBplus(3.89e9), KCplus(17.5), MU_RYR(5.0e-11),
KAminus(5.21e13), KBminus(3.89e9), KCminus(17.5), REF_CA_ER(2.5e-1),
m_dom(approx->domain()), m_time(0.0), m_oldTime(0.0), m_mg(m_dom->grid())
{
	if (m_mg->template has_attachment<side_t>(this->m_O1))
				UG_THROW("Attachment necessary for RyR channel dynamics "
						 "could not be made, since it already exists.");
		m_mg->template attach_to<side_t>(this->m_O1);

	if (m_mg->template has_attachment<side_t>(this->m_O2))
				UG_THROW("Attachment necessary for RyR channel dynamics "
						 "could not be made, since it already exists.");
		m_mg->template attach_to<side_t>(this->m_O2);

	if (m_mg->template has_attachment<side_t>(this->m_C1))
				UG_THROW("Attachment necessary for RyR channel dynamics "
						 "could not be made, since it already exists.");
		m_mg->template attach_to<side_t>(this->m_C1);

	if (m_mg->template has_attachment<side_t>(this->m_C2))
				UG_THROW("Attachment necessary for RyR channel dynamics "
						 "could not be made, since it already exists.");
		m_mg->template attach_to<side_t>(this->m_C2);

	m_aaO1 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_O1);
	m_aaO2 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_O2);
	m_aaC1 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_C1);
	m_aaC2 = Grid::AttachmentAccessor<side_t, ADouble>(*m_mg, m_C2);
};


RyR2::~RyR2()
{
	m_mg->template detach_from<side_t>(this->m_Vm);
};


void RyR2::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
{
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration
	number h = m_time - m_oldTime;

	// membrane current corresponding to diffusion pressure
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);

	// open probability
	number pO2 = pO2 + 1/h * dpO2;
	number pC1 = pC1 + 1/h * dpC1;
	number pC2 = pC1 + 1/h * dpC1;
	number pO1 = 1 - pO2 - pC1 - pC2;
	number pOpen = pO1 + pO2;
	//number pOpen =  (1.0 + KB*pow(caCyt,3)) / (1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3));

	flux[0] = pOpen * current;
}


void RyR2::calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// get values of the unknowns in associated node
	number caCyt = u[_CCYT_];	// cytosolic Ca2+ concentration
	number caER = u[_CER_];		// ER Ca2+ concentration

	// membrane potential
	number current = R*T/(4*F*F) * MU_RYR/REF_CA_ER * (caER - caCyt);
	number dI_dCyt = - R*T/(4*F*F) * MU_RYR/REF_CA_ER;
	number dI_dER = -dI_dCyt;

	// RyR flux derivs
	//number schlonz1 = 1.0 + KB*pow(caCyt,3);
	//number schlonz2 = 1.0 + KC + 1.0/(KA*pow(caCyt,4)) + KB*pow(caCyt,3);
	//number pOpen = schlonz1 / schlonz2;

	//number dOpen_dCyt = (3.0*KB*caCyt*caCyt + schlonz1/schlonz2*(4.0/(KA*pow(caCyt,5)) - 3.0*KB*caCyt*caCyt)) / schlonz2;

	// RyR2 derivs - change to equations from model0 paper
	number dpO1 = KAplus * caCyt * pC1 - KAminus * pO1 - KBplus * caCyt * pO1 + KBminus * pO2 - KCplus * pO1 + KCminus * pC2;
	number dpO2 = KBplus * caCyt * pO1 - KBminus * pO2;
	number dpC1 = KAminus * pO1 - KAplus * caCyt * pC1;
	number dpC2 = KCplus * pO1 - KCminus * pC2;

	size_t i = 0;
	if (!has_constant_value(_CCYT_))
	{
		flux_derivs[0][i].first = local_fct_index(_CCYT_);
		flux_derivs[0][i].second = (dpO1 + dpO2) * current + (pO1 + pO2) * dI_dCyt;
		i++;
	}
	if (!has_constant_value(_CER_))
	{
		flux_derivs[0][i].first = local_fct_index(_CER_);
		flux_derivs[0][i].second = (pO1 + pO2) * dI_dER;
		i++;
	}
}


const size_t RyR2::n_dependencies() const
{
	size_t n = 2;
	if (has_constant_value(_CCYT_))
		n--;
	if (has_constant_value(_CER_))
		n--;

	return n;
}


size_t RyR2::n_fluxes() const
{
	return 1;
};


const std::pair<size_t,size_t> RyR2::flux_from_to(size_t flux_i) const
{
    size_t from, to;
    if (allows_flux(_CCYT_)) to = local_fct_index(_CCYT_); else to = InnerBoundaryConstants::_IGNORE_;
    if (allows_flux(_CER_)) from = local_fct_index(_CER_); else from = InnerBoundaryConstants::_IGNORE_;

    return std::pair<size_t, size_t>(from, to);
}


const std::string RyR2::name() const
{
	return std::string("RyR2");
};


void RyR2::check_supplied_functions() const
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


void RyR2::print_units() const
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

void RyR2::update_time(const number newTime)
{
	if (newTime != m_time)
	{
		m_oldTime = m_time;
		m_time = newTime;
	}
}


} // namespace neuro_collection
} // namespace ug



