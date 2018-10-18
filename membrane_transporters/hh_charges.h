/*
 *	Discretization of Hodgkin-Huxley voltage-dependent channels in the plasma membrane
 *
 *  Created on: 2017-10-30
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__HH_CHARGES_H__
#define __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__HH_CHARGES_H__

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
 * Units used in the implementation of this channel:
 * membrane potential     V
 * conductances           S/m^2
 * [current] C/s
 *
 */

template <typename TDomain>
class HHCharges
: public IMembraneTransporter,
  public IElemDisc<TDomain>
{
	public:
		enum{_RHOI_ = 0, _RHOO_, _PHII_, _PHIO_, _N_, _M_, _H_};

		static const int dim = TDomain::dim;	 //!< world dimension
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	public:
		/// constructor with functions and subsets as vectors of string
		HHCharges(const std::vector<std::string>& fcts, const std::vector<std::string>& subsets);

		/// constructor with functions and subsets as c-style string
		HHCharges(const char* fcts, const char* subsets);

		/// destructor
		virtual ~HHCharges();

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
			RegisterFV1(HHCharges<TDomain>* pThis) : m_pThis(pThis){}
			HHCharges<TDomain>* m_pThis;
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

#endif // __UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__HH_CHARGES_H__

