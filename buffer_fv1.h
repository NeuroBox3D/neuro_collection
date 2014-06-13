/*
 * buffer.h
 *
 * Finite Volume Element discretization for an arbitrary number of chemical
 * (buffering) two-agents reactions.
 *
 * \tparam	TDomain		Domain
 *
 *  Created on: 07.11.2012
 *      Author: mbreit
 */

#ifndef __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__BUFFER_FV1__
#define __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__BUFFER_FV1__

#include <boost/function.hpp>
#include <vector>
#include <string>

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"


/// discretization for a buffering equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for a buffering process involving one buffer and one "buffee".
 * The equations have the form
 * \f[
 * 	\partial_t c + k_b \cdot c \cdot  b - k_u \cdot \left( b_{tot} - b \right) = 0
 * 	\partial_t b + k_b \cdot c \cdot  b - k_u \cdot \left( b_{tot} - b \right) = 0
 * \f]
 * where
 * <ul>
 * <li>	\f$ c \f$ is the concentration of the buffered substance,
 * <li>	\f$ b \f$ is the concentration of the buffering substance (unbound buffer),
 * <li>	\f$ b_{tot} \f$ is the total concentration of buffer (bound or unbound),
 * <li>	\f$ k_b \f$ is the association rate for the buffering reaction,
 * <li>	\f$ k_u \f$ is the dissociation rate for the buffering reaction.
 * </ul>
 *
 * NOTE: This discretization must always be used in conjunction with another
 * time-dependent process such as diffusion which carries out the discretization
 * of the time derivative. NO MASS ASSEMBLINGS ARE PERFORMED IN THIS DISCRETIZATION.
 * This has the advantage that the discretization can simply be added to a domain
 * disc already containing a diffusion disc e.g.
 * \tparam	TDomain		Domain
 */



namespace ug
{
namespace neuro_collection
{

/// struct that holds information about the unknowns involved in a reaction
/// as well as the kinetics constants for their reaction
template<int dim> struct ReactionInfo
{
	ReactionInfo(std::size_t _b, std::size_t _bd,
				 SmartPtr<CplUserData<number, dim> > _tb,
			     SmartPtr<CplUserData<number, dim> > _kb,
			     SmartPtr<CplUserData<number, dim> > _ku)
		: buffer(_b), buffered(_bd)
	{
		tot_buffer.set_data(_tb);
		k_bind.set_data(_kb);
		k_unbind.set_data(_ku);
	};

	std::size_t buffer;		// index of buffer
	std::size_t buffered;	// index of buffered agent
	DataImport<number, dim> tot_buffer; // data import for total buffer concentration
	DataImport<number, dim> k_bind;		// binding constant
	DataImport<number, dim> k_unbind;	// unbinding constant
};


template<typename TDomain>
class BufferFV1
: public IElemDisc<TDomain>
{
	private:
		///	base class type
		typedef IElemDisc<TDomain> base_type;

		///	own type
		typedef BufferFV1<TDomain> this_type;

	public:
		///	domain type
		typedef typename base_type::domain_type domain_type;

		///	world dimension
		static const int dim = base_type::dim;

		///	position type
		typedef typename base_type::position_type position_type;

		/// error estimator type
		typedef MultipleSideAndElemErrEstData<TDomain> err_est_type;

	public:

		/// constructor
        BufferFV1(const char* subsets)
        	: IElemDisc<TDomain>(NULL, subsets), m_bNonRegularGrid(false)
        {
        	register_all_fv1_funcs();
        	m_reactions.clear();
        }

        /// add a reaction (with DataImports)
        void add_reaction(const char* fct1, const char* fct2,
        			 SmartPtr<CplUserData<number, dim> > tbc,
        			 SmartPtr<CplUserData<number, dim> > k1,
        			 SmartPtr<CplUserData<number, dim> > k2);

        /// add a reaction (with constants)
        void add_reaction(const char* fct1, const char* fct2, number tbc, number k1, number k2);


	private:
		/// reactions information
		std::vector<ReactionInfo<dim> > m_reactions;


    public:	// inherited from IElemDisc

        ///	type of trial space for each function used
		void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

		///	returns if hanging nodes are used
		bool use_hanging() const;


	private:

		///	prepares the loop over all elements
		/**
		 * This method prepares the loop over all elements. It resizes the position
		 * array for the corner coordinates and schedules the local ip positions
		 * at the data imports.
		 */
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

		///	finishes the loop over all elements
		template<typename TElem, typename TFVGeom>
		void fsh_elem_loop();

		///	prepares the element for assembling
		/**
		 * This methods prepares an element for the assembling. The positions of
		 * the element corners are read and the Finite Volume Geometry is updated.
		 * The global ip positions are scheduled at the data imports.
		 */
		template<typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

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

		///	postprocesses the loop over all elements of one type in the computation of the error estimator
		template <typename TElem, typename TFVGeom>
		void fsh_err_est_elem_loop();


	private:
		void register_all_fv1_funcs();

		struct RegisterFV1
		{
				RegisterFV1(this_type* pThis) : m_pThis(pThis) {};
				this_type* m_pThis;
				template<typename TElem > void operator()(TElem&)
				{
					if (m_pThis->m_bNonRegularGrid)
						m_pThis->register_fv1_func<TElem, HFV1Geometry<TElem, dim> >();
					else
						m_pThis->register_fv1_func<TElem, FV1Geometry<TElem, dim> >();
				}
		};

		template <typename TElem, typename TFVGeom> void register_fv1_func();

		/// struct holding values of shape functions in IPs
		struct ShapeValues
		{
			public:
				void resize(std::size_t nEip, std::size_t _nSh)
				{
					nSh = _nSh;
					elemVals.resize(nEip);
					for (std::size_t i = 0; i < nEip; i++) elemVals[i].resize(nSh);
				}
				number& shapeAtElemIP(std::size_t sh, std::size_t ip) {return elemVals[ip][sh];}
				number* shapesAtElemIP(std::size_t ip) {return &elemVals[ip][0];}
				std::size_t num_sh() {return nSh;}
			private:
				std::size_t nSh;
				std::vector<std::vector<number> > elemVals;
		} m_shapeValues;

	private:
		bool m_bNonRegularGrid;
};

} // namespace neuro_collection
} // namespace ug

#include "buffer_fv1_impl.h"

#endif /* __H__UG__PLUGINS__EXPERIMENTAL__NEURO_COLLECTION__BUFFER_FV1__ */
