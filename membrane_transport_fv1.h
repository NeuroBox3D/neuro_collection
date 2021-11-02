/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2011-12-20
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORT_FV1_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORT_FV1_H


#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"
#include "common/util/smart_pointer.h"
#include "membrane_transporters/membrane_transporter_interface.h"


namespace ug {
namespace neuro_collection {

// forward declaration of IMembraneTransporter
class IMembraneTransporter;

///@addtogroup plugin_neuro_collection
///@{

/// Finite Volume element discretization for the inner BndCond on a two-sided membrane
/**
 * This class implements the InnerBoundary interface to provide element local
 * assemblings for the unknown-dependent Neumann flux over a membrane, where the flowing
 * unknowns are present on both sides of the membrane.
 *
 * It can be used for the discretization of all kinds of trans-membrane transport mechanisms
 * (channels or pumps), the dynamics of which are defined by objects of the interface class
 * IMembraneTransporter.
 * Such an object can be assigned to an object of this class using the method
 * set_membrane_transporter(). The density of the corresponding channels or pumps needs
 * to be set using set_density_function().
 */
template <typename TDomain>
class MembraneTransportFV1
: public FV1InnerBoundaryElemDisc<MembraneTransportFV1<TDomain>, TDomain>
{
	protected:
		const number R;			// universal gas constant
		const number T;			// temperature
		const number F;			// Faraday constant

		typedef MembraneTransportFV1<TDomain> this_type;
		typedef FV1InnerBoundaryElemDisc<this_type, TDomain> base_type;
		typedef typename base_type::FluxCond FluxCond;
		typedef typename base_type::FluxDerivCond FluxDerivCond;

	public:
		typedef TDomain domain_type;

	///	world dimension
		static const int dim = TDomain::dim;

	public:
	/// constructor with c-string
		MembraneTransportFV1(const char* subsets, SmartPtr<IMembraneTransporter> mt);

	/// constructor with vector
		MembraneTransportFV1(const std::vector<std::string>& subsets, SmartPtr<IMembraneTransporter> mt);

	/// destructor
		virtual ~MembraneTransportFV1();

	public:
	/// adding density information for pumps/channels in membrane
		void set_density_function(SmartPtr<CplUserData<number,dim> > densityFct);

	/// adding density information for pumps/channels in membrane
		void set_density_function(const number dens);

	/// adding density information for pumps/channels in membrane
		void set_density_function(const char* name);

	/// set transport mechanism
		void set_membrane_transporter(SmartPtr<IMembraneTransporter> mt);

		// functions needed by FV1InnerBoundaryElemDisc

		/**	This is the actual flux function defining the flux density over the boundary
		 *	depending on the unknowns on the boundary;
		 */
		bool fluxDensityFct
		(
			const std::vector<LocalVector::value_type>& u,
			GridObject* e,
			const MathVector<dim>& coords,
			int si,
			FluxCond& fc
		);
		bool fluxDensityFct
		(
			const std::vector<LocalVector::value_type>& u,
			const std::vector<LocalVector::value_type>& uOld,
			GridObject* e,
			const MathVector<dim>& coords,
			int si,
			FluxCond& fc
		);
		template <typename CallToMemTransporter>
		bool fluxDensityFctImpl
		(
			const CallToMemTransporter& call,
			GridObject* e,
			const MathVector<dim>& coords,
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
			const MathVector<dim>& coords,
			int si,
			FluxDerivCond& fdc
		);
		bool fluxDensityDerivFct
		(
			const std::vector<LocalVector::value_type>& u,
			const std::vector<LocalVector::value_type>& uOld,
			GridObject* e,
			const MathVector<dim>& coords,
			int si,
			FluxDerivCond& fdc
		);
		template <typename CallToMemTransporter>
		bool fluxDensityDerivFctImpl
		(
			const CallToMemTransporter& call,
			GridObject* e,
			const MathVector<dim>& coords,
			int si,
			FluxDerivCond& fdc
		);

	/// @copydoc IElemDisc<TDomain>::prepare_setting()
		void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid) override;

	/// @copydoc IElemDisc<TDomain>::prep_timestep()
		template <typename TAlgebra>
		void prep_timestep(number future_time, number time, VectorProxyBase* upb);


	protected:
		SmartPtr<CplUserData<number,dim> > m_spDensityFct;
		SmartPtr<IMembraneTransporter> m_spMembraneTransporter;

	private:
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


		void register_all_fv1_funcs();
};

///@}

} // end namespace neuro_collection
} // end namespace ug


#endif  // UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORT_FV1_H
