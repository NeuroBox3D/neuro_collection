/*
 * dependent_neumann_boundary_fv1.h
 * 
 *  Created on: 18.01.2013
 *      Author: markusbreit
 */

#ifndef __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__DEPENDENT_NEUMANN_BOUNDARY_FV1__
#define __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__DEPENDENT_NEUMANN_BOUNDARY_FV1__

#include <boost/function.hpp>
#include <vector>
#include <string>

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"



namespace ug
{
namespace neuro_collection
{

/// Neumann boundary condition that depends on the unknowns (on the bnd)
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the unknown-dependent Neumann flux over an outer boundary.
 * The equation of this flux should be given in a concretization of this class.
 *
 * This is the one-sided equivalent of inner_boundary_fv1.h.
 *
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */

template<typename TDomain>
class DependentNeumannBoundaryFV1
	: public IElemDisc<TDomain>
{
	public:
	/// struct that holds information about the flux densities and from where to where the flux occurs
		struct NFluxCond
		{
			// vector of fluxFctValues: flux[i] holds the flux density of unknown i
			std::vector<number> flux;
			std::vector<std::size_t> to;
		};

	/// struct that holds information about the derivatives of the flux densities
	/// and from where to where the flux occurs
		struct NFluxDerivCond
		{
			// vector of fluxFctDerivValues: fluxDeriv[i][j] holds the derivative
			// of the flux density of unknown i wrt unknown j
			std::vector<std::vector<number> > fluxDeriv;
			std::vector<std::size_t> to;
		};

	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	own type
		typedef DependentNeumannBoundaryFV1<TDomain> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	/// error estimator type
		typedef MultipleSideAndElemErrEstData<TDomain> err_est_type;

	public:
	
    /// Constructor
		DependentNeumannBoundaryFV1(const char* functions, const char* subsets)
        	: IElemDisc<TDomain>(functions, subsets), m_bNonRegularGrid(false)
        {
        	register_all_fv1_funcs();
        }

	/// destructor
		virtual ~DependentNeumannBoundaryFV1() {};

	
	public:	// inherited from IElemDisc
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns if hanging nodes are used
		virtual bool use_hanging() const;

	private:
	
	/// the flux function
	/**	This is the actual flux function defining the flux density over the boundary
	 *	depending on the unknowns on the boundary;
	 *	shall be defined in a specialized class that is derived from FV1MyNeumannBoundaryElemDisc.
	 */
		virtual bool fluxDensityFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& cc,
									int si, NFluxCond& fc) = 0;

	/**	This is the flux derivative function defining the flux density derivatives over the boundary
	 *	depending on the unknowns on the boundary;
	 *	shall be defined in a specialized class that is derived from FV1MyNeumannBoundaryElemDisc.
	 */
		virtual bool fluxDensityDerivFct(const std::vector<LocalVector::value_type>& u, const MathVector<dim>& cc,
										 int si, NFluxDerivCond& fdc) = 0;
	
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
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
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	prepares the loop over all elements of one type for the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling the error estimator
		template <typename TElem, typename TFVGeom>
		void prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	computes the error estimator contribution for one element
		template <typename TElem, typename TFVGeom>
		void compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	summarizes the contributions of the error estimator in one element
		template <typename TElem, typename TFVGeom>
		number get_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	postprocesses the loop over all elements of one type in the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void fsh_err_est_elem_loop();


	protected:
		// position access
		std::vector<Vertex*> m_vVertices;
		Vertex* m_currVertex;

	private:
		void register_all_fv1_funcs();

		struct RegisterFV1
		{
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
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

		/// struct holding values of shape functions in IPs
		struct ShapeValues
		{
			public:
				void resize(std::size_t nSip, std::size_t _nSh)
				{
					nSh = _nSh;
					sideVals.resize(nSip);
					for (std::size_t i = 0; i < nSip; i++) sideVals[i].resize(nSh);
				}
				number& shapeAtSideIP(std::size_t sh, std::size_t ip) {return sideVals[ip][sh];}
				number* shapesAtSideIP(std::size_t ip) {return &sideVals[ip][0];}
				std::size_t num_sh() {return nSh;}
			private:
				std::size_t nSh;
				std::vector<std::vector<number> > sideVals;
		} m_shapeValues;

	private:
		bool m_bNonRegularGrid;
};

} // end namespace neuro_collection
} // end namspace ug

#endif /*__H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__DEPENDENT_NEUMANN_BOUNDARY_FV1__*/
