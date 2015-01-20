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


///@addtogroup plugin_neuro_collection
///@{


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
		const number R;			///< universal gas constant
		const number T;			///< temperature
		const number F;			///< Faraday constant

		const number D1;		///< IP3 binding (w/o Ca2+ inhibition)
		const number D2;		///< Ca2+ inhibiting binding
		const number D3;		///< IP3 binding (w/ Ca2+ inhibition)
		const number D5;		///< Ca2+ activating binding
		const number MU_IP3R;	///< IP3R channel conductance

		const number REF_CA_ER;	///< reference endoplasmic Ca2+ concentration (for conductances)

	public:
		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		IP3R(std::vector<std::string> fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~IP3R();

		/// @copydoc IMembraneTransporter::calc_flux()
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

		/// @copydoc IMembraneTransporter::calc_flux_deriv()
		virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

		/// @copydoc IMembraneTransporter::n_dependencies()
		virtual const size_t n_dependencies() const;

		/// @copydoc IMembraneTransporter::n_fluxes()
		virtual size_t n_fluxes() const;

		/// @copydoc IMembraneTransporter::flux_from_to()
		virtual const std::pair<size_t,size_t> flux_from_to(size_t flux_i) const;

		/// @copydoc IMembraneTransporter::name()
		virtual const std::string name() const;

		/// @copydoc IMembraneTransporter::check_supplied_functions()
		virtual void check_supplied_functions() const;

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const;
};

///@}

} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__IP3R_H__

