/*
 * measurement.h
 *
 *  Created on: 31.05.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__MEASUREMENT_H
#define UG__PLUGINS__NEURO_COLLECTION__MEASUREMENT_H

#include "common/util/smart_pointer.h"
#include "lib_disc/function_spaces/approximation_space.h"


namespace ug {
namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{



/**
 * \brief outputs the volumes of chosen subsets
 *
 *	The result is written to the command line interface.
 *
 * \param approx		the underlying approximation space
 * \param subsetNames	contains the names of the subsets that the volume measuring
 * 						is to be performed for, separated by commas
 *
 * @todo	In the parallel case:
 * 			For lower-dimensional objects it is possible that some are accounted for
 * 			more than once! --> master/slave layout !?
 */
template <typename TDomain>
void computeVolume
(
	ConstSmartPtr<ApproximationSpace<TDomain> > approx,
	const char* subsetNames
);


/**
 * \brief outputs average values of specified unknowns on specified subsets.
 *
 *	The result is written to (a) continuing file(s) specified by the user.
 *
 * \param solution		the vector containing the unknowns of the problem
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
 * @todo	In the parallel case:
 * 			For lower-dimensional objects it is possible that some are accounted for
 * 			more than once! --> master/slave layout !?
 */
template <typename TGridFunction>
void takeMeasurement
(
	SmartPtr<TGridFunction> solution,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName
);

///@}

} // namespace ug
} // namespace neuro_collection

#include "measurement_impl.h"

#endif // UG__PLUGINS__NEURO_COLLECTION__MEASUREMENT_H
