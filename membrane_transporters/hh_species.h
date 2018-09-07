/*
 *	Discretization of Hodgkin-Huxley voltage-dependent channels in the plasma membrane
 *
 *  Created on: 2018-09-05
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__HH_SPECIES_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__HH_SPECIES_H__

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug {
namespace neuro_collection {


///@addtogroup plugin_neuro_collection
///@{


/// Discretization of Hodgkin-Huxley voltage-dependent channels in the plasma membrane
/**
 * This class implements the MembraneTransport interface to provide current densities
 * and their derivatives for the Hodgkin-Huxley model of voltage-dependent channels.
 *
 * The gating parameters are kept internally and need to be updated before each step.
 * This means, they are discretized explicitly w.r.t. the membrane potential.
 *
 * Units used in the implementation of this channel:
 * membrane potential     V
 * conductances           S/m^2
 * [current] C/s
 *
 */

template <typename TDomain>
class HHSpecies
: public IMembraneTransporter
{
	public:
		enum{_KI_ = 0, _KO_, _NAI_, _NAO_, _PHII_, _PHIO_};

		static const int dim = TDomain::dim;	 //!< world dimension
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		/// constructor with functions and subsets as vectors of string
		HHSpecies
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			ConstSmartPtr<ISubsetHandler> spSH
		);

		/// constructor with functions and subsets as c-style string
		HHSpecies(const char* fcts, const char* subsets, ConstSmartPtr<ISubsetHandler> spSH);

		/// destructor
		virtual ~HHSpecies();

	public:
		/**
		 * @brief Set HH channel conductances
		 * @param gk   potassium conductance [C/Vs]
		 * @param gna  sodium conductance [C/Vs]
		 */
		void set_conductances(number gk, number gna);

		/**
		 * @brief Set HH channel reversal potentials
		 * @param ek   potassium reversal potential [V]
		 * @param ena  sodium reversal potential [V]
		 */
		void set_reversal_potentials(number ek, number ena);

		/**
		 * @brief Set the reference time used in the simulation.
		 * @todo maybe replace by "time scale"?
		 * @param refTime reference time in units of s
		 */
		void set_reference_time(number refTime);

		/**
		 * @brief Set voltage-explicit exact discretization mode for gates
		 *
		 * The gating parameters m, n and h are governed by a linear
		 * ordinary differential equation
		 *     u' = \frac{u_\infty - u}{\tau_u}
		 * the parameters of which (u_/infty and \tau_u) depend on the
		 * (time-dependent) membrane voltage V_m.
		 * When this dependency is discretized in an explicit manner,
		 * the solution of the ODE can be given in exact terms:
		 *     u(t+\Delta t) = u_\infty - (u_\infty - u(t)) * \exp{(-\frac{\Delta t}{\tau_u})}.
		 * This can be assembled with only changes to the stiffness terms as
		 *     0 = u(t+\Delta t) - u(t) - (u_\infty - u(t)) * (1 - \exp{(-\frac{\Delta t}{\tau_u})}).
		 * Note that it is necessary to scale the new stiffness term with {\Delta t}^{-1},
		 * because of the stiffness scale factor in the instationary case.
		 * This time step size has to be supplied by the user!
		 */
		void use_exact_gating_mode();


	private:
		template <typename TAlgebra, int locDim>
		void update_potential
		(
			ConstSmartPtr<DoFDistribution> dd,
			int si,
			const typename TAlgebra::vector_type& u
		);

        /// initializes the defined channel type
        /**
         * During the initialization, the gating parameters are set to the equilibrium state
         * for the start membrane potential.
        **/
		void init_gating(GridObject* elem);

		/// updates the gating parameters
		/**
		 * This method needs to be called before calc_flux().
		 * @param newTime new point in time
		 */
		void update_gating(GridObject* elem);

		/// updates internal time if necessary
		void update_time(number newTime);

		template <typename TAlgebra, int dim>
		void prep_timestep_with_algebra_type_and_dim
		(
			number future_time,
			const number time,
			const typename TAlgebra::vector_type& u,
			int si
		);

		template <typename TAlgebra>
		void prep_timestep_with_algebra_type
		(
			number future_time,
			const number time,
			const typename TAlgebra::vector_type& u
		);


	// inheritances from IMembraneTransporter
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


	protected:
		number m_gK;    ///< potassium single-channel conductance [C/Vs]
		number m_gNa;   ///< sodium single-channel conductance [C/Vs]
		number m_eK;    ///< potassium reversal potential [V]
		number m_eNa;   ///< sodium reversal potential [V]

		struct GatingInfo
		{
			number vm;
			number n;
			number m;
			number h;
		};
		typedef std::map<GridObject*, GatingInfo> GatingMap;
		GatingMap m_mGating;                        //!< current values for Vm and n, m, h

		ConstSmartPtr<ISubsetHandler> m_spSH;       //!< subset handler
		std::vector<std::string> m_vSubset;         //!< subsets this channel exists on

		number m_refTime;
		number m_time;								//!< current time
		number m_initTime;							//!< time of initialization
		number m_oldTime;							//!< time step before current time

		bool m_bInitiated;							//!< indicates whether channel has been initialized by init()

		bool m_bVoltageExplicitDiscMode;
};

///@}

} // namespace neuro_collection
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__HH_SPECIES_H__

