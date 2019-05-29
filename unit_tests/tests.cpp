/*!
 * \file tests.cpp
 * \brief neuro collection unit tests
 *
 *  Created on: May 28, 2019
 *      Author: stephan
 */

#define BOOST_TEST_MODULE __CPP__UNIT_TESTS__UG__NEURO_COLLECTION__TESTS__

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <common/math/ugmath.h>
#include <lib_grid/grid/grid.h>

#include "../test/neurite_util.h"

using namespace ug;
using namespace ug::neuro_collection;

////////////////////////////////////////////////////////////////////////
/// neuro_collection/test tests
////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(test);

////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(CreatePyramid) {
	// create vertices
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	ug::Vertex *p1, *p2, *p3, *p4;
	p1 = *g.create<RegularVertex>(); p2 = *g.create<RegularVertex>();
	p3 = *g.create<RegularVertex>(); p4 = *g.create<RegularVertex>();
	BOOST_REQUIRE_MESSAGE(p1, "Creating first vertex.");
	BOOST_REQUIRE_MESSAGE(p2, "Creating second vertex.");
	BOOST_REQUIRE_MESSAGE(p3, "Creating third vertex.");
	BOOST_REQUIRE_MESSAGE(p4, "Creating forth vertex.");

	// assign coordinates
	aaPos[p1] = ug::vector3(0, 0, 0); aaPos[p2] = ug::vector3(0, 1, 0);
	aaPos[p3] = ug::vector3(1, 1, 0); aaPos[p4] = ug::vector3(1, 0, 0);

	// create elements and checks
	const Quadrilateral* const quad = *g.create<Quadrilateral>(QuadrilateralDescriptor(p1, p2, p3, p4));
	BOOST_REQUIRE_MESSAGE(quad, "Creating quadrilateral out of supplied vertices.");
	Pyramid* p = create_pyramid(g, quad, aaPos);
	BOOST_REQUIRE_MESSAGE(p, "Creating pyramid out of supplied vertices and created top vertex.");
	vector3 top = aaPos[p->vertex(4)]; //!< last vertex of pyramid is top vertex always
	BOOST_CHECK_EQUAL(top.x(), 0.5);
	BOOST_CHECK_EQUAL(top.y(), 0.5);
	BOOST_CHECK_EQUAL(top.z(), -1.0);
}
////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END();
////////////////////////////////////////////////////////////////////////
/// neuro_collection/test tests
////////////////////////////////////////////////////////////////////////
