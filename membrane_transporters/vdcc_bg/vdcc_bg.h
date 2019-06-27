/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2013-02-05
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_H

#include "../membrane_transporter_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"  // for FV1ManifoldGeometry
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"  // for HFV1ManifoldGeometry
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"  // for IElemDisc

namespace ug {
namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{

/// Class interface for Borg Graham type VGCCs of the plasma membrane.
/** This class is an interface for the Borg-Graham-type voltage-gated calcium channels
 *	in the plasma membrane (see chapter 8 of "Interpretations of data and mechanisms
 *	for hippocampal pyramidal cell models", Borg-Graham (1998) in Cerebral Cortex,
 *	Vol. 13: Models of Cortical Circuits.
 *
 *	The unknowns of the discretization (the so-called "gates" of the channel) are
 *	not added to the system unknowns, but held separately as vertex attachments, since
 *	they are governed by ODEs and are therefore merely updated by an implicit time
 *	schema once per time step. This will have to be done by the user (in the lua script)
 *	via the ITimeDiscretization method prepare_step_elem() before the execution of any
 *	time step assembling.
 *
 *	The class implements simple versions of N-, L- and T type channels but may very
 *	well be generalized to more complex models thereof. The type of channel can be set
 *	before calling the init() method, N-type channel parameters are taken as default.
 *
 *	This class does not handle the procuration of the values for the membrane
 *	potential. This must be dealt with in a specialization of this class.
 *
 *	Any class specializing this interface _must_ implement the virtual method:
 *	- void update_potential(number newTime);
 *	they _can_ re-implement the virtual method
 *	- void init(number time).
 *
 *
 *	The units required for this discretization are:
 * 		V_m	: mV		membrane voltage (VDCC_BG_VM2UG & VDCC_BG_VM2UG_NEURON)
 * 		V_m	: V			membrane voltage (VDCC_BG & VDCC_BG_UserData)
 * 		t	: s			time
 *  	f	: mol*s^-1	ionic flux
 *
 *  Remarks:
 *  	- Internally, all membrane potentials are attached to elements in [V] for
 *  	  the element discretization!
 *
 *  	- VDCC_BG & VDCC_BG_UserData use [ms] and [mV] in gating and flux calculations!
 *  	  This is due to the use of the original gating parameter sets by Borg-Graham
 *  	  with tau_0 in [ms] and V_12 in [mV].
 *  	  Note that the update potential method already takes care of this, when calling
 *  	  the corresponding gating & flux calculation methods.
 *
**/

template<typename TDomain>
class VDCC_BG
: public IMembraneTransporter,
  public IElemDisc<TDomain>
{
	public:
		/// channel types N, L and T
		enum {BG_Ntype, BG_Ltype, BG_Ttype};
		enum{_CCYT_ = 0, _CEXT_, _M_, _H_};

		static const int dim = TDomain::dim;	//!< world dimension

		// some type definitions
		typedef VDCC_BG<TDomain> this_type;
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		typedef typename elem_t::side side_t;


	protected:
		const number R;			///< universal gas constant
		const number T;			///< temperature
		const number F;			///< Faraday constant

		/// holds the parameters of a channel type
		struct GatingParams
		{
			GatingParams(number _z, number _v, number _t) : z(_z), V_12(_v), tau_0(_t){};
			number z;
			number V_12;  // in mV
			number tau_0; // in ms
		};

	public:
	// inherited from IMembraneTransport
		/**
		 * @brief Constructor for the Borg-Graham type channel interface
		 *
		 * This constructor not only needs information on the functions involved, but also on the subsets
		 * and the approximation space, as it needs to create side attachments for each of the sides in
		 * the subsets involved.
		 *
		 * @param fcts		functions vector
		 * @param subsets	subsets vector
		 * @param approx	underlying approximation space
		 */
		VDCC_BG
		(
			const std::vector<std::string>& fcts,
			const std::vector<std::string>& subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		/**
		 * @brief Constructor for the Borg-Graham type channel interface
		 *
		 * This constructor not only needs information on the functions involved, but also on the subsets
		 * and the approximation space, as it needs to create side attachments for each of the sides in
		 * the subsets involved.
		 *
		 * @param fcts		functions as comma-separated c-string
		 * @param subsets	subsets as comma-separated c-string
		 * @param approx	underlying approximation space
		 */
		VDCC_BG
		(
			const char* fcts,
			const char* subsets,
			SmartPtr<ApproximationSpace<TDomain> > approx
		);

		/// @copydoc IMembraneTransporter::IMembraneTransporter()
		virtual ~VDCC_BG();

		/// @copydoc IMembraneTransporter::prepare_timestep()
		virtual void prepare_timestep
		(
			number future_time, const number time, VectorProxyBase* upb
		);

		/// @copydoc IMembraneTransporter::calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const
		virtual void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

/*
		/// @copydoc IMembraneTransporter::calc_flux(const std::vector<number>&, number&,  flux) const
		virtual number calc_flux(const std::vector<number>& u, size_t index);
*/
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

	// own methods
		/// sets the channel type
		template<int TType> void set_channel_type();

        /**
         * @brief Sets the permeability of this channel
         * @param perm    permeability values
         */
		void set_permeability(const number perm);

        /// initializes the defined channel type
        /** During the initialization, the necessary attachments are attached to the vertices
         *  and their values calculated by the equilibrium state for the start membrane potential.
        **/
		virtual void init(number time);

		/// updates the potential values in the corresponding attachments to new time.
		/**
		 * This method needs to be called before update_gating() if potential is non-constant.
		 * @param newTime new point in time
		 */
		virtual void update_potential(side_t* elem) = 0;

		/// updates internal time if necessary
		virtual void update_time(number newTime);


		// inheritances from IElemDisc
	public:
		/// type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		/// returns if hanging nodes are used
		virtual bool use_hanging() const;
#if 0
		/// @copydoc IElemDisc::approximation_space_changed()
		virtual void approximation_space_changed();
#endif

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
			RegisterFV1(VDCC_BG<TDomain>* pThis) : m_pThis(pThis){}
			VDCC_BG<TDomain>* m_pThis;
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
		bool m_bNonRegularGrid;
		bool m_bCurrElemIsHSlave;


	protected:
		/// calculates the equilibrium state of a gating "particle"
		/** The calculation is done with respect to the given gating parameters set (which represents
		 *	one gating "particle") and the given membrane potential (to be specified in [mV]!).
		**/
		number calc_gating_start(const GatingParams& gp, number Vm) const;

	public:
		/// init gating variables to equilibrium
		template <typename TGridFunction>
		void calculate_steady_state(SmartPtr<TGridFunction> u, number vm) const;

	private:
		void after_construction();

	protected:
		/// whether this channel has an inactivating gate
		bool has_hGate() const {return this->m_channelType == BG_Ntype || this->m_channelType == BG_Ttype;}

	protected:
		SmartPtr<TDomain> m_dom;							//!< underlying domain
		SmartPtr<Grid> m_mg;								//!< underlying multigrid
		SmartPtr<DoFDistribution> m_dd;						//!< underlying surface dof distribution
		ConstSmartPtr<MGSubsetHandler> m_sh;				//!< underlying subset handler
		typename TDomain::position_accessor_type& m_aaPos;	//!< underlying position accessor

		std::vector<std::string> m_vSubset;					//!< subsets this channel exists on
		size_t m_localIndicesOffset;

		ADouble m_Vm;								//!< membrane voltage (in Volt)
		Grid::AttachmentAccessor<side_t, ADouble> m_aaVm;		//!< accessor for membrane potential

		GatingParams m_gpMGate;						//!< gating parameter set for activating gate
		GatingParams m_gpHGate;						//!< gating parameter set for inactivating gate

		number m_time;								//!< current time
		number m_initTime;							//!< time of initialization
		number m_oldTime;							//!< time step before current time

		number m_perm;								//!< channel permeability [m^3/s] (= diff coeff * cross section / membrane thickness)
		int m_mp, m_hp;								//!< powers for gating parameters

		int m_channelType;							//!< channel type

		bool m_initiated;							//!< indicates whether channel has been initialized by init()
};


///@}


} // namespace neuro_collection
} // namespace ug

#include "vdcc_bg_impl.h"

#endif // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__VDCC_BG__VDCC_BG_H
