/*
 *  membrane_transport_fv1.h
 *
 *  Created on: 20.12.2011
 *      Author: mbreit
 */

#ifndef __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__TWO_SIDED_MEMBRANE_TRANSPORT_FV1__
#define __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__TWO_SIDED_MEMBRANE_TRANSPORT_FV1__


#include "lib_disc/spatial_disc/elem_disc/inner_boundary/inner_boundary.h"
#include "bindings/lua/lua_user_data.h"
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
 */
template<typename TDomain>
class MembraneTransportFV1
: public FV1InnerBoundaryElemDisc<TDomain>
{
	protected:
		const number R;			// universal gas constant
		const number T;			// temperature
		const number F;			// Faraday constant

		typedef MembraneTransportFV1<TDomain> this_type;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxCond FluxCond;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxDerivCond FluxDerivCond;

	public:
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

	/// @copydoc FV1InnerBoundary<TDomain>::fluxDensityFct()
		virtual bool fluxDensityFct
		(
			const std::vector<LocalVector::value_type>& u,
			GridObject* e,
			const MathVector<dim>& coords,
			int si,
			FluxCond& fc
		);

	/// @copydoc FV1InnerBoundary<TDomain>::fluxDensityDerivFct()
		bool fluxDensityDerivFct
		(
			const std::vector<LocalVector::value_type>& u,
			GridObject* e,
			const MathVector<dim>& coords,
			int si,
			FluxDerivCond& fdc
		);

	/// @copydoc IElemDisc<TDomain>::prepare_setting()
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	/// @copydoc IElemDisc<TDomain>::prepare_timestep()
		void prep_timestep(number future_time, number time, VectorProxyBase* upb);

	protected:
		SmartPtr<CplUserData<number,dim> > m_spDensityFct;
		SmartPtr<IMembraneTransporter> m_spMembraneTransporter;

	private:
		template <typename List>
		struct Register
		{
			Register(this_type* p)
			{
				static const bool isEmpty = boost::mpl::empty<List>::value;
				(typename boost::mpl::if_c<isEmpty, RegEnd, RegNext>::type (p));
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
					p->set_prep_timestep_fct(aid, &this_type::prep_timestep);

					(Register<NextList> (p));
				}
			};
		};


		void register_all_fv1_funcs();

	private:
		bool m_bNonRegularGrid;
};




/// 1d Finite Volume element discretization for the inner BndCond on a two-sided membrane
/**
 * This class implements the InnerBoundary interface to provide element local
 * assemblings for the unknown-dependent flux over a membrane, where the flowing
 * unknowns are present on both sides of the membrane.
 * This is a special implementation for a 1d "cable" representation of a perfectly
 * rotationally symmetric membrane. The radius at which the membrane is located must be provided.
 */
template<typename TDomain>
class MembraneTransport1d
: public IElemDisc<TDomain>
{
	public:
		///	world dimension
		static const int dim = IElemDisc<TDomain>::dim;

		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxCond FluxCond;
		typedef typename FV1InnerBoundaryElemDisc<TDomain>::FluxDerivCond FluxDerivCond;

	public:
		/// constructor with c-string
		MembraneTransport1d(const char* subsets, SmartPtr<IMembraneTransporter> mt);

		/// constructor with vector
		MembraneTransport1d(const std::vector<std::string>& subsets, SmartPtr<IMembraneTransporter> mt);

		/// adding density information for pumps/channels in membrane
		void set_density_function(SmartPtr<CplUserData<number,dim> > densityFct);

		/// adding density information for pumps/channels in membrane
		void set_density_function(const number dens);

		/// adding density information for pumps/channels in membrane
		void set_density_function(const char* name);

		/// set radius at which membrane is located
		void set_radius(number r);

		/// the flux function
		/**	This is the actual flux function defining the flux density over the boundary
		 *	depending on the unknowns on the boundary;
		 *	shall be defined in a specialized class that is derived from FV1InnerBoundaryElemDisc.
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
		 *	shall be defined in a specialized class that is derived from FV1InnerBoundaryElemDisc.
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
		number m_radius;
		SmartPtr<CplUserData<number,dim> > m_spDensityFct;
		SmartPtr<IMembraneTransporter> m_spMembraneTransporter;

	private:
		int m_currSI;
};


///@}

} // end namespace neuro_collection
} // end namespace ug


#endif //__H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__TWO_SIDED_MEMBRANE_TRANSPORT_FV1__
