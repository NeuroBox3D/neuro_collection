/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Markus Breit
 * Creation date: 2016-12-27
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

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H

#include "common/math/ugmath_types.h"
#include "common/util/smart_pointer.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"
#include <vector>
#include <string>
#include <list>
#include <utility>
#include "types.h"

namespace ug {
namespace neuro_collection {

/*!
 * \brief import a SWC file
 * \param[in] fileName
 * \param[out] vPointsOut
 * \param[in] correct
 * \param[in] scale
 */
void import_swc
(
    const std::string& fileName,
    std::vector<SWCPoint>& vPointsOut,
    bool correct,
    number scale
);

/*!
 * \brief smoothing
 * \param[inout] vPointsInOut
 * \param[in] n
 * \param[in] h
 * \param[in] gamma
 * Note:s Needs debugging.
 */
void smoothing
(
	std::vector<SWCPoint>& vPointsInOut,
	size_t n,
	number h,
	number gamma
);

/*!
 * \brief collapse short edges
 * \param[in] g
 * \param[in] sh
 * Note: Needs debugging
 */
void collapse_short_edges
(
	Grid& g,
	SubsetHandler& sh
);

/*!
 * \brief converts pointlist to a neuritelist
 * \param[in] vPoints
 * \param[in] vPosOut
 * \param[in] vRadOut
 * \param[in] vBPInfoOut
 * \param[in] vRootNeuriteIndsOut
 */
void convert_pointlist_to_neuritelist
(
    const std::vector<SWCPoint>& vPoints,
    std::vector<std::vector<vector3> >& vPosOut,
    std::vector<std::vector<number> >& vRadOut,
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > >& vBPInfoOut,
    std::vector<size_t>& vRootNeuriteIndsOut
);

/*!
 * \brief test smoothing
 * \param[in] fileName
 * \param[in] n
 * \param[in] h
 * \param[in] gamma
 */
void test_smoothing
(
	const std::string& fileName,
	size_t n,
	number h,
	number gamma
);

/*!
 * \brief test import swc
 * \param[in] fileName
 * \param[in] anisotropy
 * \param[in] numRefs
 */
void test_import_swc
(
	const std::string& fileName,
	number anisotropy = 2.0,
	size_t numRefs = 0
);

/*!
 * \brief test import swc with er
 * \param[in] fileName
 * \param[in] correct
 * \param[in] erScaleFactor
 * \param[in] anisotropy
 * \param[in] numRefs
 */
void test_import_swc_with_er
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number erScaleFactor,
	number anisotropy = 2.0,
	size_t numRefs = 0
);

/*!
 * \brief test import swc surf
 * \param[in] fileName
 */
void test_import_swc_surf
(
	const std::string& fileName
);

/*!
 * \brief test import swc 1d
 * \param[in] fileName
 * \param[in] anisotropy
 * \param[in] numRefs
 * \param[in] scale
 */
void test_import_swc_1d
(
	const std::string& fileName,
	number anisotropy = 2.0,
	size_t numRefs = 0,
	number scale = 1e-6
);

/*!
 * \brief shrink quadrilateral copy
 * \param[in] vVrt
 * \param[out] outvVrt
 * \param[in] oldVertices
 * \param[out] outvEdge
 * \param[in] g
 * \param[in] aaPos
 * \param[in] percentage
 * \param[in] createFaces
 * \param[in] outSel
 * \param[in] currentDir
 */
void shrink_quadrilateral_copy
(
	const std::vector<Vertex*>& vVrt,
	std::vector<Vertex*>& outvVrt,
	const std::vector<Vertex*>& oldVertices,
	std::vector<Edge*>& outvEdge,
	Grid& g,
	Grid::VertexAttachmentAccessor<APosition>& aaPos,
	number percentage,
	bool createFaces,
	ISelector* outSel,
	ug::vector3* currentDir
);

/*!
 * \brief test smoothing
 * \param[in] fileName
 * \param[in] n
 * \param[in] h
 * \param[in] gamma
 * \param[in] scale
 */
void test_smoothing
(
	const std::string& fileName,
	size_t n,
	number h,
	number gamma,
	number scale
);

/*!
 * \brief test import swc
 * \param[in] fileName
 * \param[in] correct
 */
void test_import_swc
(
	const std::string& fileName,
	bool correct
);

/*!
 * \brief test import swc scale
 * \param[in] fileName
 * \param[in] correct
 * \param[in] scale
 */
void test_import_swc_scale
(
	const std::string& fileName,
	bool correct,
	number scale
);

/*!
 * \brief test import swc general
 * Builds geometry with or without ER and connects to soma
 * \param[in] fileName
 * \param[in] correct
 * \param[in] erScaleFactor
 * \param[in] withER
 * \param[in] anisotropy
 * \param[in] numRefs
 */
void test_import_swc_general
(
	const std::string& fileName,
	bool correct,
	number erScaleFactor,
	bool withER,
	number anisotropy = 2.0,
	size_t numRefs = 1,
	bool regularize = false
);

/*!
 * \brief test import swc general smooth
 * \param[in] fileName
 * \param[in] correct
 * \param[in] shrinkPercentage
 * \param[in] withER
 */
void test_import_swc_general_smooth
(
	const std::string& fileName,
	bool correct,
	number shrinkPercentage,
	bool withER
);

/*!
 * \brief test neurite projector with four section tube
 */
void test_neurite_projector_with_four_section_tube();

/*!
 * \brief test neurite projector with four section tube and branch point
 */
void test_neurite_projector_with_four_section_tube_and_branch_point();

/*!
 * \brief test shrink geometry
 * \param[in] percentage
 */
void test_shrink_geom
(
	number percentage
);

/*!
 * \brief test shrink geom copy
 * \param[in] percentage
 */
void test_shrink_geom_copy
(
	number percentage
);

/*!
 * \brief test split geom
 * \param[in] percentage
 */
void test_split_geom
(
	number percentage
);

/*!
 * \brief apply neurite projector
 * \param[inout] mg
 * \param[in] neuriteProj
 */
void apply_neurite_projector
(
	MultiGrid& mg,
	SmartPtr<NeuriteProjector> neuriteProj
);

/*!
 * \brief converts swc to ugx
 * \param[in] fileName
 */
void test_convert_swc_to_ugx
(
	const std::string& fileName
);

/*!
 * \brief new strategy
 */
void test_import_swc_general_var(
	const std::string& fileName,
	bool correct,
	number erScaleFactor,
	bool withER,
	number anisotropy,
	size_t numRefs,
	bool regularize,
	number blowUpFactor,
	bool forVR
);

/*!
 * \brief vr strategy
 */
void test_import_swc_general_var_for_vr(
	const std::string& fileName,
	bool correct,
	number erScaleFactor,
	bool withER,
	number anisotropy,
	size_t numRefs,
	bool regularize,
	number blowUpFactor
);

/*!
 * \brief refine swc grid
 * \param[in] inFileName
 * \param[in] outFileName
 */
void refine_swc_grid(
	const std::string& inFileName,
	const std::string& outFileName
);

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H
