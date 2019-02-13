/*
 *  ryr_discrete.h
 *
 *  Created on: 2019-02-06
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_DISCRETE_H
#define UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_DISCRETE_H


#include "common/util/smart_pointer.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"  // for IDomainConstraint

#include <string>
#include <vector>


namespace ug {
namespace neuro_collection {


///@addtogroup plugin_neuro_collection
///@{

/// Finite Volume element discretization for discrete RyR channels
/**
 * This class implements the InnerBoundary interface to provide element local
 * assemblings for the unknown-dependent Neumann flux over a membrane, where the flowing
 * unknowns are present on both sides of the membrane.
 */
template <typename TDomain, typename TAlgebra>
class RyRDiscrete
: public IDomainConstraint<TDomain, TAlgebra>
{
	public:
		///	world dimension
		static const int dim = TDomain::dim;

		///	matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

		///	vector type
		typedef typename TAlgebra::vector_type vector_type;

		/// value type
		typedef typename matrix_type::value_type value_type;

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

		enum
		{
			_CCYT_ = 0,
			_CER_  = 1,
			_O2_   = 2,
			_C1_   = 3,
			_C2_   = 4
		};

	public:
	/// constructor with c-string
		RyRDiscrete(const char* functions, const char* subsets);

	/// constructor with vector
		RyRDiscrete(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

	/// destructor
		virtual ~RyRDiscrete();


	// inherited from IConstraint
	public:
		///	adapts defect to enforce constraints
		virtual void adjust_defect
		(
			vector_type& d,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
			const std::vector<number>* vScaleMass = NULL,
			const std::vector<number>* vScaleStiff = NULL
		);

		///	adapts jacobian to enforce constraints
		virtual void adjust_jacobian
		(
			matrix_type& J,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
			const number s_a0 = 1.0
		);


		///	sets the constraints in a solution vector
		virtual void adjust_solution
		(
			vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		);

		///	adapts matrix and rhs (linear case) to enforce constraints
		virtual void adjust_linear
		(
			matrix_type& mat,
			vector_type& rhs,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		);

		///	adapts a rhs to enforce constraints
		virtual void adjust_rhs
		(
			vector_type& rhs,
			const vector_type& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		);

		///	returns the type of constraints
		virtual int type() const;


	// inherited from IDomainConstraint
	public:
		virtual void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	public:
		/// init gating variables to equilibrium
		void calculate_steady_state(SmartPtr<vector_type> u) const;

		/// set open probability below which the channel is supposed to be certainly closed
		void set_cutoff_open_probability(number cutoffProb);

	protected:
		std::vector<std::string> m_vSubsetNames;
		std::vector<int> m_vSI;

		std::vector<std::string> m_vFunctionNames;
		std::vector<size_t> m_vFctMap;

		number m_cutoffOpenProb;
};

///@}

} // end namespace neuro_collection
} // end namespace ug


#endif //UG__PLUGINS__NEURO_COLLECTION__MEMBRANE_TRANSPORTERS__RYR_DISCRETE_H
