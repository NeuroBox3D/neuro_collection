/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-10-30
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__HH_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__HH_H

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
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
 * Current densities go directly to the inner and outer potentials, respectively, instead
 * of to inner and outer charge densities as in HHCharges or HHSpecies.
 *
 * Units used in the implementation of this channel:
 * membrane potential     V
 * conductances           S/m^2
 * [current] C/s
 *
 */

template <typename TDomain>
class HH
: public IMembraneTransporter,
  public IElemDisc<TDomain>
{
	public:
		enum{_PHII_ = 0, _PHIO_, _N_, _M_, _H_};

		static const int dim = TDomain::dim;	 //!< world dimension
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		/// constructor with functions and subsets as vectors of string
		HH(const std::vector<std::string>& fcts, const std::vector<std::string>& subsets);

		/// constructor with functions and subsets as c-style string
		HH(const char* fcts, const char* subsets);

		/// destructor
		virtual ~HH();

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
		void use_exact_gating_mode(number timeStep);

		void use_gating_explicit_current_mode();

	// inheritances from IMembraneTransporter
	public:
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

	// inheritances from IElemDisc
	public:
		/// type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		/// returns if hanging nodes are used
		virtual bool use_hanging() const;

	// assembling functions
	protected:
		///	prepares the loop over all elements (of a type and subset)
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	prepares the element for assembling
		template<typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

		/// finishes the loop over all elements
		template<typename TElem, typename TFVGeom>
		void fsh_elem_loop();

		///	assembles the local stiffness matrix
		template<typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the local mass matrix
		template<typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the stiffness part of the local defect
		template<typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the mass part of the local defect
		template<typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the local right hand side
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	protected:
		void register_all_fv1_funcs();

		struct RegisterFV1
		{
			RegisterFV1(HH<TDomain>* pThis) : m_pThis(pThis){}
			HH<TDomain>* m_pThis;
			template< typename TElem > void operator()(TElem&)
			{
				if (m_pThis->m_bNonRegularGrid)
					m_pThis->register_fv1_func<TElem, HFV1ManifoldGeometry<TElem, dim> >();
				else
					m_pThis->register_fv1_func<TElem, FV1ManifoldGeometry<TElem, dim> >();
			}
		};

		template <typename TElem, typename TFVGeom>
		void register_fv1_func();

	protected:
		number m_gK;    ///< potassium single-channel conductance [C/Vs]
		number m_gNa;   ///< sodium single-channel conductance [C/Vs]
		number m_eK;    ///< potassium reversal potential [V]
		number m_eNa;   ///< sodium reversal potential [V]

		number m_refTime;

		bool m_bVoltageExplicitDiscMode;
		bool m_bGatingExplicitCurrentMode;
		number m_VEDMdt;

	protected:
		bool m_bNonRegularGrid;
		bool m_bCurrElemIsHSlave;
};

///@}

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__HH_H

