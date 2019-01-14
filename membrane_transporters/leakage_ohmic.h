/*
 *	Discretization for an ohmic leakage flux through a membrane
 *
 *  Created on: 2017-11-01
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAKAGE_OHMIC_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAKAGE_OHMIC_H__

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug {
namespace neuro_collection {


///@addtogroup plugin_neuro_collection
///@{


/// Discretization for an ohmic leakage flux through a membrane
/**
 * This class implements an ohmic leakage flux through a membrane,
 * i.e., it assumes a linear dependency of the current from the difference between
 * membrane potential and reversal potential:
 * j = g * (Vm - E_L). The leakage conductance g can be set using the appropriate method.
 *
 * The calculated current is outward.
 *
 * Units used in the implementation of this channel:
 * Potential:    V
 * Conductance:  C/(Vs)
 * Current:      C/s
 */

class OhmicLeakage : public IMembraneTransporter
{
	public:
		enum{_PHII_ = 0, _PHIO_};

	public:
		/// @copydoc IMembraneTransporter::IMembraneTransporter(const std::vector<std::string)
		OhmicLeakage(const std::vector<std::string>& fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		OhmicLeakage(const char* fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~OhmicLeakage();

		/// set leakage conductance
		void set_conductance(number g);

		/// set leakage reversal potential
		void set_reversal_potential(number el);

		/// @copydoc IMembraneTransporter::calc_flux()
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

		/// @copydoc IMembraneTransporter::calc_flux_deriv()
		virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

		/// @copydoc IMembraneTransporter::n_dependencies()
		virtual size_t n_dependencies() const;

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

	protected:
		number m_g;  ///< conductance
		number m_eL; ///< reversal potential
};


/// Discretization for an ohmic leakage flux through a membrane
/**
 * This class implements an ohmic leakage flux through a membrane,
 * i.e., it assumes a linear dependency of the current from the difference between
 * membrane potential and reversal potential:
 * j = g * (Vm - E_L). The leakage conductance g can be set using the appropriate method.
 *
 * The calculated current is outward.
 *
 * Units used in the implementation of this channel:
 * Potential:    V
 * Conductance:  C/(Vs)
 * Current:      C/s
 */

class OhmicLeakageCharges : public IMembraneTransporter
{
	public:
		enum{_RHOI_ = 0, _RHOO_, _PHII_, _PHIO_};

	public:
		/// @copydoc IMembraneTransporter::IMembraneTransporter(const std::vector<std::string)
		OhmicLeakageCharges(const std::vector<std::string>& fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		OhmicLeakageCharges(const char* fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~OhmicLeakageCharges();

		/// set leakage conductance
		void set_conductance(number g);

		/// set leakage reversal potential
		void set_reversal_potential(number el);

		/// @copydoc IMembraneTransporter::calc_flux()
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

		/// @copydoc IMembraneTransporter::calc_flux_deriv()
		virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

		/// @copydoc IMembraneTransporter::n_dependencies()
		virtual size_t n_dependencies() const;

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

	protected:
		number m_g;  ///< conductance
		number m_eL; ///< reversal potential
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__LEAKAGE_OHMIC_H__

