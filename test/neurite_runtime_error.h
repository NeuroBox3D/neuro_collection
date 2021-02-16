 /* Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2020-08-26
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_ERROR_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_ERROR_H

#include "common/math/ugmath_types.h"
#include "types.h"
#include <stdexcept>

namespace ug {
	namespace neuro_collection {
		/*!
		 * \brief A list of possible error codes associated with the RTEs
		 * DO NOT USE these codes in production - see please comments below.
		 *
		 * Note that these error codes should not go towards the public API
		 * as they are used only internally for generation of the meshing
		 * statistics in the grid generation subroutines in neuro_collection.
		 * By intention the error codes are not tied to the custom runtime
		 * error classes, as the error codes are used only auxiliary in
		 * offline post-processing of meshes to generate desired statistics.
		 */
		enum NeuriteErrorCode
		{
			NEURITE_RUNTIME_ERROR_CODE_SUCCESS, // 0
			NEURITE_RUNTIME_ERROR_CODE_OTHER, // 1
			NEURITE_RUNTIME_ERROR_CODE_REGULARIZATION_INCOMPLETE, // 2
			NEURITE_RUNTIME_ERROR_CODE_INVALID_BRANCHES, // 3
			NEURITE_RUNTIME_ERROR_CODE_CONTAINS_CYCLES, // 4
			NEURITE_RUNTIME_ERROR_CODE_CYLINDER_CYLINDER_OVERLAP, // 5
			NEURITE_RUNTIME_ERROR_CODE_SOMA_CONNECTION_OVERLAP, // 6
			NEURITE_RUNTIME_ERROR_CODE_TETRAHEDRALIZE_FAILURE, // 7
			NEURITE_RUNTIME_ERROR_CODE_BP_ITERATION_FAILURE, // 8
			NEURITE_RUNTIME_ERROR_CODE_NO_PERMISSIBLE_RENDER_VECTOR_FOUND, // 9
			NEURITE_RUNTIME_ERROR_CODE_HIGH_DIAMETER_VARIABILITY, // 10
			NEURITE_RUNTIME_ERROR_CODE_BRANCHING_POINT_CLUSTERING, // 11
			NEURITE_RUNTIME_ERROR_CODE_SMALL_OR_NEGATIVE_RADIUS // 12
		};

		/*!
		 * \brief common base for all neurite errors
		 */
		struct NeuriteRuntimeError : public std::runtime_error {
			 NeuriteRuntimeError() : std::runtime_error("NeuriteError") {}
			 NeuriteRuntimeError(const std::string& msg) : std::runtime_error(msg) {}
			 virtual ~NeuriteRuntimeError() throw() {}
		};


		/// RegularizationIncomplete
		struct RegularizationIncomplete : public NeuriteRuntimeError {
				RegularizationIncomplete(const std::string& msg)
					: NeuriteRuntimeError(msg) {}
		};

		/// InvalidBranches
		struct InvalidBranches : public NeuriteRuntimeError {
			virtual const char* what() const throw() {
				return "Invalid number of neurite branches n > 3 in 1D geometry.";
			}
		};

		/// ContainsCycle
		struct ContainsCycles : public NeuriteRuntimeError {
			virtual const char* what() const throw() {
				return "1D Geometry contains at least one cycle.";
			}
		};

		/// CylinderCylinderOverlap
		 struct CylinderCylinderOverlap: public NeuriteRuntimeError {
			virtual const char* what() const throw() {
				return "Neurites overlap in 3D geoemtry";
			}
		};

		/// SomaConnectionOverlap
		struct SomaConnectionOverlap : public NeuriteRuntimeError {
			virtual const char* what() const throw() {
				return "Soma connections overlap in 3D geometry";
			}
		};

		/// TetrahedralizeFailure
		struct TetrahedralizeFailure : public NeuriteRuntimeError {
			virtual const char* what() const throw() {
				return "Tetrahedralize failed for 3D geometry";
			}
		};

		/// NoPermissibleRenderVector
		struct NoPermissibleRenderVector : public NeuriteRuntimeError {
			NoPermissibleRenderVector(const std::string& msg)
								: NeuriteRuntimeError(msg) {}
		};

		/// HighDiameterVariability
		struct HighDiameterVariability : public NeuriteRuntimeError {
			HighDiameterVariability(const std::string& msg)
	 							 : NeuriteRuntimeError(msg) {}
		};

		/// BranchingPointClustering
		struct BranchingPointClustering : public NeuriteRuntimeError {
			BranchingPointClustering(const std::string& msg) 
								 : NeuriteRuntimeError(msg) {}
		}; 

		/// SmallOrNegativeRadius
		struct SmallOrNegativeRadius : public NeuriteRuntimeError {
            SmallOrNegativeRadius(const std::string& msg) 
								: NeuriteRuntimeError(msg) {}
		};
	}
}

#endif //  UG__PLUGINS__NEURO_COLLECTION__TEST__NEURITE_ERROR_H
