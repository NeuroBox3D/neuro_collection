/*!
 * test_neurite_proj.h
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
 * TODO: Needs debugging
 * \param[inout] vPointsInOut
 * \param[in] n
 * \param[in] h
 * \param[in] gamma
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
 * TODO: Needs debugging
 */
void collapse_short_edges
(
	Grid& g,
	SubsetHandler& sh
);

/*!
 * \brief converts pointlist to a neuritelist
 * TODO: Document parameters
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
 * TODO: Document parameters
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
 * TODO: Document parameters
 */
void test_import_swc
(
	const std::string& fileName,
	number anisotropy = 2.0,
	size_t numRefs = 0
);

/*!
 * \brief test import swc with er
 * TODO: Document parameters
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
 * TODO: Document parameters
 */
void test_import_swc_surf
(
	const std::string& fileName
);

/*!
 * \brief test import swc 1d
 * TODO: Document parameters
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
 * TODO: Document parameters
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
 * TODO: Document parameters
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
 * TODO: Document parameters
 */
void test_import_swc
(
	const std::string& fileName,
	bool correct
);

/*!
 * \brief test import swc scale
 * TODO: Document parameters
 */
void test_import_swc_scale
(
	const std::string& fileName,
	bool correct,
	number scale
);

/*!
 * \brief test import swc general
 * TODO: Document parameters
 */
void test_import_swc_general
(
	const std::string& fileName,
	bool correct,
	number shrinkPercentage,
	bool withER
);

/*!
 * \brief test import swc general smooth
 * TODO: Document parameters
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
 * TODO: Document parameters
 */
void test_neurite_projector_with_four_section_tube();

/*!
 * \brief test neurite projector with four section tube and branch point
 * TODO: Document parameters
 */
void test_neurite_projector_with_four_section_tube_and_branch_point();

/*!
 * \brief test shrink geometry
 * TODO: Document parameters
 */
void test_shrink_geom
(
	number percentage
);

/*!
 * \brief test shrink geom copy
 * TODO: Document parameters
 */
void test_shrink_geom_copy
(
	number percentage
);

/*!
 * \brief test split geom
 * TODO: Document parameters
 */
void test_split_geom
(
	number percentage
);

/*!
 * \brief apply neurite projector
 * TODO: Document parameters
 */
void apply_neurite_projector
(
	MultiGrid& mg,
	SmartPtr<NeuriteProjector> neuriteProj
);

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H
