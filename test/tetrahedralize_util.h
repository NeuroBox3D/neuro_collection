/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-07-19
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__TETRAHEDRALIZE_UTIL_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__TETRAHEDRALIZE_UTIL_H

#include <lib_grid/algorithms/geom_obj_util/geom_obj_util.h>
#include <lib_grid/algorithms/remove_duplicates_util.h>

namespace ug {
	namespace neuro_collection {
		/*!
		* \brief Fills a closed surface-grid selection with tetrahedrons
		* You may specify a quality parameter. If this parameter is <= 0, no
		* inner vertices will be created. The default value for quality is 5.
		* The quality resembles the minimal valid dihedral angle.
		* The algorithm should always terminate for this quality. If you choose
		* a lower quality parameter (careful with quality < 1), the algotithm may
		* not terminate. The lower the quality parameter (but > 0), the
		* better the tetrahedron quality. Using Tetgen by Hang Si.
		*
		* \param[in] sel
		* \param[in,out] grid
		* \param[in,out] sh
		* \param[in] quality        number specifiying grid quality between 0 and 18
		* \param[in] preserveBnds   bool to specify if outer boundaries shall be preserved
		* \param[in] preserveAll    bool to specify if outer and inner boundaries shall be preserved
		* \param[in] aPos
		* \param[in] verbosity	  number between 0 and 3 to specify level of verbosity
		*/
		bool Tetrahedralize
		(
			Selector& sel,
			Grid& grid,
			ISubsetHandler* sh,
		    number quality = 5,
		    bool preserveBnds = false,
		    bool preserveAll = false,
		    APosition& aPos = aPosition,
		    int verbosity = 0
		);
	}
}

#endif // UG__PLUIGNS___NEURO_COLLECTION_TEST__TETRAHEDRALIZE_UTIL_H
