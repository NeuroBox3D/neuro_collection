/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2017-08-15
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORT_1D_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORT_1D_H


//#include "bindings/lua/lua_user_data.h"
#include "common/util/smart_pointer.h"
#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"
#include "membrane_transporters/membrane_transporter_interface.h"
#include "../cable_neuron/util/diam_attachment_handler.h"	// attachment handling for diameter attachment


namespace ug {
namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{


/// 1d finite volume element discretization for the inner BndCond on a two-sided membrane
/**
 * This class implements the InnerBoundary interface to provide element local
 * assemblings for the unknown-dependent flux over a membrane, where the flowing
 * unknowns are present on both sides of the membrane.
 * This is a special implementation for a 1d "cable" representation of a perfectly
 * rotationally symmetric membrane.
 * The radius at which the membrane is located must be provided, either as an explicit
 * constant or as a diameter attachment in the grid as in cable_neuron applications.
 */
template<typename TDomain>
class MembraneTransport1d
: public IElemDisc<TDomain>
{
	public:
		///	world dimension
		static const int dim = IElemDisc<TDomain>::dim;

		typedef InnerBoundaryFluxCond FluxCond;
		typedef InnerBoundaryFluxDerivCond FluxDerivCond;

	public:
		/// constructor with c-string
		MembraneTransport1d(const char* subsets, SmartPtr<IMembraneTransporter> mt);

		/// constructor with vector
		MembraneTransport1d(const std::vector<std::string>& subsets, SmartPtr<IMembraneTransporter> mt);

		/// destructor
		virtual ~MembraneTransport1d();

		/// adding density information for pumps/channels in membrane
		void set_density_function(SmartPtr<CplUserData<number,dim> > densityFct);

		/// adding density information for pumps/channels in membrane
		void set_density_function(const number dens);

		/// adding density information for pumps/channels in membrane
		void set_density_function(const char* name);

		/// set radius of plasma membrane
		void set_radius(number r);

		/// set plasma membrane radius fraction at which membrane (ERM or PM) is located
		void set_radius_factor(number r);

		/// the flux function
		/**	This is the actual flux function defining the flux density over the boundary
		 *	depending on the unknowns on the boundary;
		 */
		bool fluxDensityFct
		(
			const std::vector<LocalVector::value_type>& u,
			GridObject* e,
			const MathVector<dim>& cc,
			int si,
			FluxCond& fc
		);

		/**	This is the flux derivative function defining the flux density derivatives over the boundary
		 *	depending on the unknowns on the boundary;
		 */
		bool fluxDensityDerivFct
		(
			const std::vector<LocalVector::value_type>& u,
			GridObject* e,
			const MathVector<dim>& cc,
			int si,
			FluxDerivCond& fdc
		);


	public:	// inherited from IElemDisc
		/// @copydoc IElemDisc::approximation_space_changed()
		virtual void approximation_space_changed();

		///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		///	returns if hanging nodes are used
		virtual bool use_hanging() const;

		/// @copydoc IElemDisc<TDomain>::prepare_timestep()
		//template <typename TAlgebra>
		void prep_timestep(number future_time, number time, VectorProxyBase* upb);

		///	prepares the loop over all elements
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	prepares the element for assembling
		template<typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

		///	finishes the loop over all elements
		template<typename TElem, typename TFVGeom>
		void fsh_elem_loop();

		///	assembles the local stiffness matrix using a finite volume scheme
		template<typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

		///	assembles the local mass matrix using a finite volume scheme
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

	private:
		template <typename List>
		struct RegisterPrepTimestepFct
		{
			RegisterPrepTimestepFct(MembraneTransport1d* p)
			{
				static const bool isEmpty = boost::mpl::empty<List>::value;
				(typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type (p));
			}

			struct RegEnd
			{
				RegEnd(MembraneTransport1d*) {}
			};

			struct RegNext
			{
				RegNext(MembraneTransport1d* p)
				{
					typedef typename boost::mpl::front<List>::type AlgebraType;
					typedef typename boost::mpl::pop_front<List>::type NextList;

					size_t aid = bridge::AlgebraTypeIDProvider::instance().id<AlgebraType>();

					// TODO: should be (with a TAlgebra-templated version of this method):
					// p->set_prep_timestep_fct(aid, &MembraneTransport1d::prep_timestep<AlgebraType>);
					p->set_prep_timestep_fct(aid, &MembraneTransport1d::prep_timestep);

					(RegisterPrepTimestepFct<NextList> (p));
				}
			};
		};

		void register_assembling_funcs();

	protected:
		number m_radiusFactor;
		number m_constRadius;
		bool m_bConstRadiusSet;
		ANumber m_aDiameter;									 ///< diameter attachment
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaDiameter;  ///< diameter attachment accessor
		cable_neuron::DiamAttachmentHandler m_dah;							 ///< handler for multigrid usage of diameter attachment

		SmartPtr<CplUserData<number,dim> > m_spDensityFct;
		SmartPtr<IMembraneTransporter> m_spMembraneTransporter;

	private:
		int m_currSI;
};

///@}

} // end namespace neuro_collection
} // end namespace ug


#endif  // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORT_1D_H
