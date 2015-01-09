/*
 *  General channel/pump transport interface class implementation
 *
 *  Created on: 07.01.2015
 *     Authors: mbreit, mstepniewski
 */


#include "membrane_transporter_interface.h"

namespace ug{
namespace neuro_collection{


IMembraneTransporter::IMembraneTransporter(std::vector<std::string> vFct)
: n_fct(vFct.size()), m_bLocked(false)
{
	// check all unknowns given
	for (size_t i = 0; i < vFct.size(); i++)
	{
		// the corresponding unknown has been given
		if (vFct[i] != "")
		{
			m_mfInd[i] = m_vFct.size();
			m_vFct.push_back(vFct[i]);
		}
		// else: ignore here
	}

	// resize input scaling vector
	// (fluxes scaling vector cannot be resized at this point, since n_fluxes() is not yet available)
	m_vScaleInputs.resize(n_fct,1.0);
};

IMembraneTransporter::~IMembraneTransporter()
{
	// do nothing
}

void IMembraneTransporter::flux(const std::vector<number>& u, std::vector<number>& flux) const
{
	// construct input vector for flux calculation with constant values
	std::vector<number> u_with_consts(n_fct);
	create_local_vector_with_constants(u,u_with_consts);

	// scale each entry
	for (size_t i = 0; i < n_fct; i++)
		u_with_consts[i] *= m_vScaleInputs[i];

	// calculate fluxes
	calc_flux(u_with_consts, flux);

	// scale each flux
	for (size_t i = 0; i < flux.size(); i++)
		flux[i] *= m_vScaleFluxes[i];
}


void IMembraneTransporter::flux_deriv(const std::vector<number>& u, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const
{
	// construct input vector for flux derivative calculation with constant values
	std::vector<number> u_with_consts(n_fct);
	create_local_vector_with_constants(u,u_with_consts);

	// scale each entry
	for (size_t i = 0; i < n_fct; i++)
		u_with_consts[i] *= m_vScaleInputs[i];

	// calculate flux derivatives
	calc_flux_deriv(u_with_consts, flux_derivs);

	// scale each flux deriv
	for (size_t i = 0; i < flux_derivs.size(); i++)
		for (size_t j = 0; j < flux_derivs[i].size(); j++)
			flux_derivs[i][j].second *= m_vScaleFluxes[i] * m_vScaleInputs[flux_derivs[i][j].first];
}


void IMembraneTransporter::create_local_vector_with_constants(const std::vector<number>& u, std::vector<number>& u_wc) const
{
	for (size_t i = 0; i < n_fct; i++)
	{
		// check if not constant
		if (!has_constant_value(i,u_wc[i]))
		{
			UG_ASSERT(m_mfInd.find(i) != m_mfInd.end(), "Neither function supplied nor constant. "
					  "This is not supposed to happen! Check implementation!");
			u_wc[i] = u[m_mfInd.find(i)->second];
		}
		// else: the constant is already written to u_with_consts by has_constant_value()
	}
}

const std::vector<std::string>& IMembraneTransporter::symb_fcts() const
{
	return m_vFct;
}


void IMembraneTransporter::check_constant_allowed(const size_t i, const number val) const
{
	// standard implementation will not do anything
}


void IMembraneTransporter::set_constant(const size_t i, const number val)
{
	// check that not already locked
	if (m_bLocked)
	{
		UG_THROW("The membrane transport mechanism of type \"" << name() << "\" is locked "
				 "and cannot be altered after that.\n"
				 "Please set any constants before passing to "
				 "an instance of TwoSidedMembraneTransport.");
	}

	// check that this unknown is allowed to be set constant (mechanism specific)
	try
	{
		check_constant_allowed(i,val);
	}
	UG_CATCH_THROW("Setting the unknown with index " << i << " to a constant value\n"
				   "is not allowed by an instance of the membrane transport mechanism "
				   "of type \"" << name() << "\".");

	m_mConstVal[i] = val;
}

const bool IMembraneTransporter::has_constant_value(const size_t i, number& val) const
{
	std::map<size_t, number>::const_iterator mit = m_mConstVal.find(i);
	if (mit == m_mConstVal.end())
		return false;

	val = mit->second;

	return true;
}

const bool IMembraneTransporter::has_constant_value(const size_t i) const
{
	if (m_mConstVal.find(i) == m_mConstVal.end())
		return false;

	return true;
}

const bool IMembraneTransporter::allows_flux(const size_t i) const
{
	if (m_mfInd.find(i) == m_mfInd.end())
		return false;

	return true;
}

void IMembraneTransporter::print_units() const
{
	UG_LOG("The units of the transport mechanism of type \"" << name() << "\"\n"
			"are not (yet) supplied by the implementation.\n");
}

void IMembraneTransporter::set_scale_inputs(const std::vector<number>& scale)
{
	if (scale.size() != n_fct)
	{
		UG_THROW("Scaling factors vector does not have the required length\n"
				 "(" << scale.size() << " instead of " << n_fct <<")"
				 " for transport mechanism of type \"" << name() << "\".\n"
				 "Make sure you set exactly as many scaling factors as the\n"
				 "number of functions that need to be defined in the constructor.");
	}

	for (size_t i = 0; i < n_fct; i++)
		m_vScaleInputs[i] = scale[i];
}

void IMembraneTransporter::set_scale_fluxes(const std::vector<number>& scale)
{
	m_vScaleFluxes.resize(n_fluxes(), 1.0);
	if (scale.size() != n_fluxes())
	{
		UG_THROW("Scaling factors vector does not have the required length\n"
				 "(" << scale.size() << " instead of " << n_fluxes() <<")"
				 " for transport mechanism of type \"" << name() << "\".\n"
				 "Make sure you set exactly as many scaling factors as the\n"
				 "number of fluxes that occur through this transport mechanism.");
	}

	for (size_t i = 0; i < n_fluxes(); i++)
		m_vScaleFluxes[i] = scale[i];
}


void IMembraneTransporter::check_and_lock()
{
	// nothing to do if already locked
	if (m_bLocked) return;

	// check that each unknown is either given as a function or constant
	std::vector<size_t> not_ok_ind;
	for (size_t i = 0; i < n_fct; i++)
	{
		if (m_mfInd.find(i) == m_mfInd.end()
			&& m_mConstVal.find(i) == m_mConstVal.end())
		{
			not_ok_ind.push_back(i);
		}
	}

	// throw if not
	if (not_ok_ind.size())
	{
		UG_THROW("Setup for membrane transport mechanism of type \"" << name() << "\" is incorrect:\n"
				 "Any unknown involved must either be given as function in the constructor\n"
				 "or set to a constant value using the set_constant() method." );
	}

	// resize fluxes scaling vector (if necessary)
	m_vScaleFluxes.resize(n_fluxes(), 1.0);

	// lock
	m_bLocked = true;
}

const bool IMembraneTransporter::is_locked() const
{
	return m_bLocked;
}

} // namespace neuro_collection
} // namespace ug
