/*
 *	Discretization for a simplicstic NMDA receptor channel.
 *
 *  Created on: 2018-04-19
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS_NMDAR_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS_NMDAR_H__

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug {
namespace neuro_collection {


///@addtogroup plugin_neuro_collection
///@{


/// Discretization for the NMDA receptor channel
/**
 * This class implements the ionic calcium current through an NMDA receptor channel.
 * The current is assumed to be the product of open probability and single-channel current.
 *
 * The probability is given as a time-dependent function
 *   p_o(t) = exp(-(t-t0)/tau)  for t >= t0,
 *   p_o(t) = 0                 otherwise,
 * where t0 is the activation time and tau the decay time constant.
 *
 * The single_channel current is of Goldman-Hodgkin-Katz style, i.e.,
 *   i = p * zF/(RT) * V_m * (c_e - c_i*exp(zF/(RT)*V_m)) / (1.0 - exp(zF/(RT)*V_m))
 * with the permeability p, the ion valency z, Faraday's constant F, the ideal gas constant R,
 * the temperature T, the membrane potential V_m and the extra- and intracellular
 * concentrations c_e, c_i.
 *
 * The channel needs to be given the extracellular and intracellular calcium concentrations
 * as functions (in this order) in the constructor.
 *
 * Units used in the implementation of this channel:
 * [Ca_cyt]  mM (= mol/m^3)
 * [Ca_ext]  mM (= mol/m^3)
 * T         K
 * t0        s
 * tau       s
 * p         m^3/s
 * V_m       V
 *
 * Ca current   mol/s
 */

class NMDAR : public IMembraneTransporter
{
	public:
		enum{_EXT_ = 0, _CYT_};

	public:
		/// @copydoc IMembraneTransporter::IMembraneTransporter(const std::vector<std::string)
		NMDAR(const std::vector<std::string>& fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		NMDAR(const char* fcts);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~NMDAR();

	public:
		/// set activation time
		void set_activation_time(number t0);

		/// set decay time constant
		void set_decay_time(number tau);

		/// set single-channel current
		void set_permeability(number perm);

		/// set membrane potential
		void set_membrane_potential(number vm);

		/// set temperature
		void set_temperature(number temp);

	public:
		/// @copydoc IMembraneTransporter::prep_timestep()
		virtual void prep_timestep(number future_time, const number time, VectorProxyBase* upb);

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

	private:
		const number m_z;
		const number m_F;
		const number m_R;

		number m_T;
		number m_t0;
		number m_tau;
		number m_perm;
		number m_vm;

		number m_time;
};

///@}


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS_NMDAR_H__
