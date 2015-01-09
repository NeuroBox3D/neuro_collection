/*
 *  General channel/pump transport interface class
 *
 *  Created on: 07.01.2015
 *     Authors: mbreit, mstepniewski
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_INTERFACE_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_INTERFACE_H__


#include <common/common.h>
#include <utility>      // std::pair
#include <string>
#include <vector>
#include <map>


namespace ug{
namespace neuro_collection{


///
/**
 * TODO:
 * Each membrane transport mechanism MUST provide information on the units used for
 * each of the unknowns and time and space and so on...!!
 */

class IMembraneTransporter
{
	public:
		/// constructor
		IMembraneTransporter(std::vector<std::string> vFct);

		/// destructor
		virtual ~IMembraneTransporter();

		/// fluxes through this mechanism
		void flux(const std::vector<number>& u, std::vector<number>& flux) const;

		/// derivatives of fluxes through this mechanism
		void flux_deriv(const std::vector<number>& u, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const;

		/// calculation of flux (mechanism-specific)
		virtual void calc_flux(const std::vector<number>& u, std::vector<number>& flux) const = 0;

		/// calculation of flux derivatives (mechanism-specific)
		virtual void calc_flux_deriv(const std::vector<number>& u, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const = 0;

		/// return number of unknowns this transport mechanism depends on
		virtual const size_t n_dependencies() const = 0;

		/// return number of fluxes calculated by this machanism
		virtual size_t n_fluxes() const = 0;

		/// from where to where do the fluxes occur
		virtual const std::pair<size_t,size_t> flux_from_to(size_t flux_i) const = 0;

		/// name of the membrane transport system
		virtual const std::string name() const = 0;

		/// return supplied function names
		const std::vector<std::string>& symb_fcts() const;

		/// Check whether setting the i-th unknown to a constant value of val is allowed.
		/**
		 * UG_THROWs, if not allowed.
		 * @param i		index of the unknown
		 * @param val	constant value to be set
		 */
		virtual void check_constant_allowed(const size_t i, const number val) const;

		/// set a constant value instead of one of the unknowns
		void set_constant(const size_t i, const number val);

		/// returns whether the unknown of an index is set constant; and if so gets this value
		const bool has_constant_value(const size_t i, number& val) const;

		/// returns whether the unknown of an index is set constant
		const bool has_constant_value(const size_t i) const;

		/// returns whether the unknown of an index is given as a function
		const bool allows_flux(const size_t i) const;

		/// prints the units this implementation uses
		virtual void print_units() const;

		/// scaling of the inputs
		void set_scale_inputs(const std::vector<number>& scale);

		/// scaling of the inputs
		void set_scale_fluxes(const std::vector<number>& scale);

		/// check that all values are either given as unknowns or constants
		void check_and_lock();

		/// return whether MembraneTransporter is locked
		const bool is_locked() const;

	private:
		void create_local_vector_with_constants(const std::vector<number>& u, std::vector<number>& u_wc) const;

	private:
		/// local vector of supplied function names
		std::vector<std::string> m_vFct;

		/// indices of unknowns in local vector
		std::map<size_t, size_t> m_mfInd;

		/// constant values map
		std::map<size_t, number> m_mConstVal;

		/// number of functions in total (supplied or constant)
		const size_t n_fct;

		/// scaling factors
		std::vector<number> m_vScaleInputs;
		std::vector<number> m_vScaleFluxes;

		/// lock status
		bool m_bLocked;
};


} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__MEMBRANE_TRANSPORTER_INTERFACE_H__

