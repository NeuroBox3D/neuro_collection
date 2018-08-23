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

void test_smoothing(const std::string& fileName, size_t n, number h, number gamma, number scale);
void test_import_swc(const std::string& fileName, bool correct, number scaleER, bool withER);
void test_import_swc_scale(const std::string& fileName, bool correct, number scale);
void test_import_swc_general(const std::string& fileName, bool correct, number shrinkPercentage, bool withER);
void test_neurite_projector_with_four_section_tube();
void test_neurite_projector_with_four_section_tube_and_branch_point();
void test_shrink_geom(number percentage);

void apply_neurite_projector(MultiGrid& mg, SmartPtr<NeuriteProjector> neuriteProj);

void test_cylinder_volume_projector();

} // namespace neuro_collection
} // namespace ug

#endif // UG__PLUGINS__NEURO_COLLECTION__TEST__TEST_NEURITE_PROJ_H
