/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2019-02-06
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

/**
 *  @brief Discretization for discretely distributed ryanodine receptor channels
 *
 *  The discrete RyR channels are represented as vertices in the grid.
 *  To assemble suitable defect and Jacobian terms there, the usual IMembraneTransporter
 *  interface cannot be used (cannot assemble over 0d elements).
 *  So here we choose a constraint to do the discretization. We write directly to the
 *  defect and matrix and need to take care of proper handling of the mass scales.
 *
 *  In order to really close a channel (non-trivial equilibrium state otherwise), we
 *  supply a method set_cutoff_open_probability that set open probability below which
 *  the channel is supposed to be certainly closed.
 *
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

		enum
		{
			_CCYT_ = 0,
			_CER_  = 1,
			_O2_   = 2,
			_C1_   = 3,
			_C2_   = 4
		};

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
		/// calculate current with given unknowns (method is needed in maxRyRFluxDensity)
		void calc_flux(const std::vector<number>& u, GridObject* e, std::vector<number>& flux) const;

		/// scaling factor of a single input (method is needed in maxRyRFluxDensity)
		number scale_input(const size_t i) const;

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
