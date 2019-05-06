/*
 * test_neurite_proj.h
 *
 *  Created on: 27.12.2016
 *      Author: mbreit
 */

#ifndef UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H
#define UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H

#include "common/math/ugmath_types.h"   // vector3
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

/// TODO: Move these to a types class / file
enum swc_type
{
    SWC_UNDF = 0,
    SWC_SOMA = 1,
    SWC_AXON = 2,
    SWC_DEND = 3,
    SWC_APIC = 4
};

struct SWCPoint
{
    vector3 coords;
    number radius;
    swc_type type;
    std::vector<size_t> conns;
};

/*!
 * \brief import an SWC file
 */
void import_swc
(
    const std::string& fileName,
    std::vector<SWCPoint>& vPointsOut,
    bool correct,
    number scale
);

// TODO: both need debugging!
void smoothing(std::vector<SWCPoint>& vPointsInOut, size_t n, number h, number gamma);
void collapse_short_edges(Grid& g, SubsetHandler& sh);

void convert_pointlist_to_neuritelist
(
    const std::vector<SWCPoint>& vPoints,
    std::vector<std::vector<vector3> >& vPosOut,
    std::vector<std::vector<number> >& vRadOut,
    std::vector<std::vector<std::pair<size_t, std::vector<size_t> > > >& vBPInfoOut,
    std::vector<size_t>& vRootNeuriteIndsOut
);

void test_smoothing(const std::string& fileName, size_t n, number h, number gamma);
void test_import_swc(const std::string& fileName, number anisotropy = 2.0, size_t numRefs = 0);
void test_import_swc_with_er
(
	const std::string& fileNameIn,
	const std::string& fileNameOut,
	number erScaleFactor,
	number anisotropy = 2.0,
	size_t numRefs = 0
);
void test_import_swc_surf(const std::string& fileName);
void test_import_swc_1d(const std::string& fileName, number anisotropy = 2.0, size_t numRefs = 0, number scale = 1e-6);
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

void test_smoothing(const std::string& fileName, size_t n, number h, number gamma, number scale);
void test_import_swc(const std::string& fileName, bool correct);
void test_import_swc_scale(const std::string& fileName, bool correct, number scale);
void test_import_swc_general(const std::string& fileName, bool correct, number shrinkPercentage, bool withER);
void test_import_swc_general_smooth(const std::string& fileName, bool correct, number shrinkPercentage, bool withER);
void test_neurite_projector_with_four_section_tube();
void test_neurite_projector_with_four_section_tube_and_branch_point();
void test_shrink_geom(number percentage);
void test_shrink_geom_copy(number percentage);
void test_split_geom(number percentage);

void apply_neurite_projector(MultiGrid& mg, SmartPtr<NeuriteProjector> neuriteProj);

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H
