/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2018-09-05
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__HH_SPECIES_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__HH_SPECIES_H

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

		const number RoverF;

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
		 * @brief Set constant HH channel reversal potentials
		 * @param ek   potassium reversal potential [V]
		 * @param ena  sodium reversal potential [V]
		 */
		void set_reversal_potentials(number ek, number ena);

		/**
		 * @brief Set the temperature
		 * Sets the temperature. Only needed when no constant reversal potential is used.
		 * Default is 310 K.
		 */
		void set_temperature(number t);

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
		virtual void prepare_timestep(number future_time, const number time, VectorProxyBase* upb) override;

		/// @copydoc IMembraneTransporter::calc_flux()
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const override;

		/// @copydoc IMembraneTransporter::calc_flux_deriv()
		virtual void calc_flux_deriv(const std::vector<number>& u, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const override;

		/// @copydoc IMembraneTransporter::n_dependencies()
		virtual size_t n_dependencies() const override;

		/// @copydoc IMembraneTransporter::n_fluxes()
		virtual size_t n_fluxes() const override;

		/// @copydoc IMembraneTransporter::flux_from_to()
		virtual const std::pair<size_t,size_t> flux_from_to(size_t flux_i) const override;

		/// @copydoc IMembraneTransporter::name()
		virtual const std::string name() const override;

		/// @copydoc IMembraneTransporter::check_supplied_functions()
		virtual void check_supplied_functions() const override;

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const override;


	protected:
		number m_gK;    ///< potassium single-channel conductance [C/Vs]
		number m_gNa;   ///< sodium single-channel conductance [C/Vs]

		bool m_bConstNernstPotentials;
		number m_eK;    ///< potassium reversal potential [V]
		number m_eNa;   ///< sodium reversal potential [V]

		number m_T;

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

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__HH_SPECIES_H

