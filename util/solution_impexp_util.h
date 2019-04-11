/*
 * solution_impexp_util.h
 *
 *  Created on: 31.07.2012
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__SOLUTION_IMPEXP_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__SOLUTION_IMPEXP_UTIL_H

#include "common/types.h"  // for number
#include "common/util/smart_pointer.h"  // for SmartPtr
#include "lib_disc/function_spaces/approximation_space.h"

// configuration file for compile options
#include "nc_config.h"


namespace ug {
namespace neuro_collection {


/// @addtogroup plugin_neuro_collection
/// @{


/**
 * \brief outputs solution for specified functions on specified subsets to file
 *
 *	The solutions (separately for each function) are written to (a) continuing
 *	file(s) specified by the user. A new filename is created for each time step.
 *	At the moment, this is only functional, if UG is not used in parallel mode
 *	- for my lack of parallelizing knowledge.
 *	The function could later serve for saving results of a simulation that can
 *	later be used to "feed" the unknown functions with start values in a
 *	consecutive simulation.
 *
 * \param solution		the vector containing the unknowns of the problem
 * \param approx		the underlying approximation space
 * \param time			the simulation time which the averaged value is taken at
 * 						and shall be accorded to in the output file
 * \param subsetNames	contains the names of the subsets that the averaging
 * 						is to be performed on, separated by commas
 * \param functionNames	contains the names of the functions that the averaging
 * 						is to be performed for, separated by commas;
 * 						for each function, a separate file with the function name
 * 						as a suffix will be created
 * \param outFileName	the name of the output file(s), i.e. their prefix
 *
 * \warning This function is very old and will probably not work properly.
 */
template <typename TGridFunction>
void exportSolution
(
	SmartPtr<TGridFunction> solution,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName
);


#ifdef NC_WITH_MPM
/**
 * \brief Inputs solution for a specified function on specified subsets from file
 * 		  and sets them as a constraint to the solution.
 *
 *	The solutions (separately for each function) are read from a file specified
 *	by the user.
 *	At the moment, this is only functional, if UG is not used in parallel mode.
 *
 * \param solution		the vector the solution is to be written to
 * \param subsetNames	subsets the solution is to be specified on
 * \param functionName	function the solution is to be specified for
 * \param inFileName	the name of the input file
 */
template <typename TGridFunction>
void importSolution
(
	SmartPtr<TGridFunction> solution,
	const char* subsetNames,
	const char* functionName,
	const char* inFileName
);
#endif

/// @}

}  // namspace neuro_collection
}  // namespace ug

#include "solution_impexp_util_impl.h"

#endif  // UG__PLUGINS__NEURO_COLLECTION__UTIL__SOLUTION_IMPEXP_UTIL_H
