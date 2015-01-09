/*
 *	Discretization for the IP3R calcium channel in the ER membrane
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__IP3R_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__IP3R_H__

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug{
namespace neuro_collection{


/// Discretization for the IP3R calcium channel in the ER membrane
/**
 * This class implements the MembraneTransport interface to provide flux densities
 * and their derivatives for the De Young & Keizer (1992) model of IP3R channels.
 *
 * Units used in the implementation of this channel:
 * [Ca_cyt]  mM (= mol/m^3)
 * [Ca_er]   mM (= mol/m^3)
 * [IP3]     mM (= mol/m^3)
 *
 * Ca flux   mol/s
 *
 */

class IP3R : public IMembraneTransporter
{
	public:
		enum{_CCYT_=0, _CER_, _IP3_};

	protected:
		const number R;			// universal gas constant
		const number T;			// temperature
		const number F;			// Faraday constant

		const number D1;		// IP3 binding (w/o Ca2+ inhibition)
		const number D2;		// Ca2+ inhibiting binding
		const number D3;		// IP3 binding (w/ Ca2+ inhibition)
		const number D5;		// Ca2+ activating binding
		const number RHO_IP3R;	// average channel density for IP3R channel
		const number MU_IP3R;	// IP3R channel conductance

		const number REF_CA_ER;	// reference endoplasmatic Ca2+ concentration (for conductances)

	public:
		/// constructor
		IP3R(std::vector<std::string> fcts);

		/// destructor
		virtual ~IP3R();

		/// flux output function
		virtual void calc_flux(const std::vector<number>& u, std::vector<number>& flux) const;

		/// flux derivative output function
		virtual void calc_flux_deriv(const std::vector<number>& u, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

		/// return number of unknowns this transport mechanism depends on
		virtual const size_t n_dependencies() const;

		/// return number of fluxes calculated by this machanism
		virtual size_t n_fluxes() const;

		/// from where to where do the fluxes occur
		virtual const std::pair<size_t,size_t> flux_from_to(size_t flux_i) const;

		/// return supplied function names
		virtual const std::string name() const;

		/// Check whether setting the i-th unknown to a constant value of val is allowed.
		/**
		 * UG_THROWs, if not allowed.
		 * @param i		index of the unknown
		 * @param val	constant value to be set
		 */
		virtual void check_constant_allowed(const size_t i, const number val) const;

		/// prints the units this implementation uses
		virtual void print_units() const;
};


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__IP3R_H__

