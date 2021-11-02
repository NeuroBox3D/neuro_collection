/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-08-16
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_LINEARIZED_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_LINEARIZED_H

#include "membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"


namespace ug {
namespace neuro_collection {


///@addtogroup plugin_neuro_collection
///@{


/// Linearized discretization for the RyR calcium channel in the ER membrane
/**
 * This class implements the IMembraneTransporter to provide Ca2+ flux densities
 * and their derivatives for the Keizer & Levine (1996) RyR model.
 * Additionally, it implements the IElemDisc interface to provide element-local
 * assemblings for the (ordinary) differential equations governing the
 * channel state probabilites.
 *
 * The discretization eliminates non-linearities by suitably choosing where to evaluate
 * functions at the new and where with the old point in time.
 *
 * Units used in the implementation of this channel:
 * [Ca_cyt]  mM (= mol/m^3)
 * [Ca_er]   mM (= mol/m^3)
 *
 * Ca flux   mol/s
 *
 */

template<typename TDomain>
class RyRLinearized
: public IMembraneTransporter,
  public IElemDisc<TDomain>
{
	public:
		enum
		{
			_CCYT_ = 0,
			_CER_  = 1,
			_O2_   = 2,
			_C1_   = 3,
			_C2_   = 4
		};

		static const int dim = TDomain::dim;	//!< world dimension

		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;

	protected:
		const number R;			// universal gas constant
		const number T;			// temperature
		const number F;			// Faraday constant

		const number KAplus;	// C1 --> O1
		const number KBplus;	// O1 --> O2
		const number KCplus;	// O1 --> C2
		const number KAminus;	// C1 <-- O1
		const number KBminus;	// O1 <-- O2
		const number KCminus;	// O1 <-- C2
		const number MU_RYR;	// RyR channel conductance

		const number REF_CA_ER;		// reference endoplasmic Ca2+ concentration (for conductances)

	public:
		/// constructor with vectors
		RyRLinearized(const std::vector<std::string>& fcts, const std::vector<std::string>& subsets);

		/// constructor with c-style strings
		RyRLinearized(const char* fcts, const char* subsets);

		/// destructor
		virtual ~RyRLinearized();

	// inheritances from IMembraneTransporter
	public:
		/// @copydoc IMembraneTransporter::calc_flux()
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const override
		{return calc_flux(u, u, e, flux);}
		virtual void calc_flux(const std::vector<number>& u, const std::vector<number>& uOld, GridObject* e, std::vector<number>& flux) const override;

		/// @copydoc IMembraneTransporter::calc_flux_deriv()
		virtual void calc_flux_deriv(const std::vector<number>& u, const std::vector<number>& uOld, GridObject* e, std::vector<std::vector<std::pair<size_t, number> > >& flux_derivs) const override;

		/// @copydoc IMembraneTransporter::n_dependencies()
		virtual size_t n_dependencies() const override;

		/// @copydoc IMembraneTransporter::n_fluxes()
		virtual size_t n_fluxes() const override;

		/// @copydoc IMembraneTransporter::flux_from_to()
		virtual const std::pair<size_t,size_t> flux_from_to(size_t flux_i) const override;

		/// @copydoc IMembraneTransporter::name()
		virtual const std::string name() const override;

		/// @copydoc IMembraneTransporter::needsPreviousSolution()
		virtual bool needsPreviousSolution() const override;

		/// @copydoc IMembraneTransporter::check_supplied_functions()
		virtual void check_supplied_functions() const override;

		/// @copydoc IMembraneTransporter::print_units()
		virtual void print_units() const override;

	public:
		/// init gating variables to equilibrium
		template <typename TGridFunction> // this is supposed to be some algebra_type::vector_type
		void calculate_steady_state(SmartPtr<TGridFunction> u) const;

	// inheritances from IElemDisc
	public:
		/// type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid) override;

		/// returns if hanging nodes are used
		virtual bool use_hanging() const override;

	// assembling functions
	protected:
		/// preapre a time step
		template <typename TAlgebra>
		void prep_timestep(number future_time, number time, VectorProxyBase* upb);

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
		template <typename TElem, typename TFVGeom>
		void register_fv1_func();

		void register_all_fv1_funcs();

		struct RegisterFV1
		{
			RegisterFV1(RyRLinearized<TDomain>* pThis) : m_pThis(pThis){}
			RyRLinearized<TDomain>* m_pThis;
			template< typename TElem > void operator()(TElem&)
			{
				if (m_pThis->m_bNonRegularGrid)
					m_pThis->register_fv1_func<TElem, HFV1ManifoldGeometry<TElem, dim> >();
				else
					m_pThis->register_fv1_func<TElem, FV1ManifoldGeometry<TElem, dim> >();
			}
		};


		typedef RyRLinearized<TDomain> this_type;
		template <typename List>
		struct RegisterPrepTimestep
		{
			RegisterPrepTimestep(this_type* p)
			{
				static const bool isEmpty = boost::mpl::empty<List>::value;
				(typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type(p));
			}
			struct RegEnd
			{
				RegEnd(this_type*) {}
			};
			struct RegNext
			{
				RegNext(this_type* p)
				{
					typedef typename boost::mpl::front<List>::type AlgebraType;
					typedef typename boost::mpl::pop_front<List>::type NextList;

					size_t aid = bridge::AlgebraTypeIDProvider::instance().id<AlgebraType>();
					p->set_prep_timestep_fct(aid, &this_type::template prep_timestep<AlgebraType>);
					(RegisterPrepTimestep<NextList>(p));
				}
			};
		};

	protected:
		bool m_bNonRegularGrid;
		bool m_bCurrElemIsHSlave;

		LocalVector m_locUOld;
		SmartPtr<VectorProxyBase> m_spOldSolutionProxy;
};

///@}


} // namespace neuro_collection
} // namespace ug

#include "ryr_linearized_impl.h"

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_LINEARIZED_H

