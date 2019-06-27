/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-05-31
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__UTIL__MEASUREMENT_H
#define UG__PLUGINS__NEURO_COLLECTION__UTIL__MEASUREMENT_H

#include "common/util/smart_pointer.h"
#include "lib_disc/function_spaces/approximation_space.h"


namespace ug {
namespace neuro_collection {

///@addtogroup plugin_neuro_collection
///@{


/**
 * \brief calculates the volume of a subset
 *
 * \param approx   the underlying approximation space
 * \param subset   contains the subset index that the volume measurement
 *                 is to be performed for
 *
 * @todo In the parallel case:
 *       For lower-dimensional objects it is possible that some are accounted for
 *       more than once! --> master/slave layout !?
 */
template <typename TDomain>
number computeVolume(ConstSmartPtr<ApproximationSpace<TDomain> > approx, const int subset);


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
number takeMeasurement
(
	SmartPtr<TGridFunction> solution,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName
);

template <typename TGridFunction>
number takeMeasurement
(
	SmartPtr<TGridFunction> solution,
	const number time,
	const char* subsetNames,
	const char* functionNames,
	const char* outFileName,
	const char* outFileExt
);

///@}

} // namespace ug
} // namespace neuro_collection

#include "measurement_impl.h"

#endif // UG__PLUGINS__NEURO_COLLECTION__UTIL__MEASUREMENT_H
