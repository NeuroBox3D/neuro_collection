/*!
 * \file test_neurite_proj.h
 * How to run: ugshell -call "test_import_swc_general(\"smith.swc\", false, 0.5, true, 0.5, 0)"
 *
 *  Created on: 27.12.2016
 *      Author: mbreit
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

namespace ug {
namespace neuro_collection {

/// Edge comparable
struct EdgeLengthCompare
{
	bool operator()(const std::pair<Edge*, number> e1, const std::pair<Edge*, number> e2)
	{return e1.second > e2.second;}
};

/// SWC types
enum swc_type
{
    SWC_UNDF = 0,
    SWC_SOMA = 1,
    SWC_AXON = 2,
    SWC_DEND = 3,
    SWC_APIC = 4
};

/// SWC point
struct SWCPoint
{
    vector3 coords;
    number radius;
    swc_type type;
    std::vector<size_t> conns;
};

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
	size_t numRefs = 1
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

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H
