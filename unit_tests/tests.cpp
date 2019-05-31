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
#include "fixtures.cpp"

using namespace ug;
using namespace ug::neuro_collection;

////////////////////////////////////////////////////////////////////////
/// neuro_collection/test tests
////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(test);

////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(CreatePyramid, FixtureGrid) {
	Pyramid* p = create_pyramid(g, quad, aaPos);
	BOOST_REQUIRE_MESSAGE(p, "Creating pyramid out of supplied vertices "
			"and created top vertex.");
	vector3 top = aaPos[p->vertex(4)]; //!< last vertex is top vertex
	BOOST_CHECK_EQUAL(top.x(), 0.5);
	BOOST_CHECK_EQUAL(top.y(), 0.5);
	BOOST_CHECK_EQUAL(top.z(), -1.0);
}
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(FindQuadrilateralConstrained, FixtureGrid) {
	// find auadrilateral and checks
	Grid::traits<Quadrilateral>::secure_container quadCont;
	find_quadrilaterals_constrained(g, aaSurfParams, quadCont);
	BOOST_REQUIRE_MESSAGE(quadCont.size() == 1, "Finding quadrilateral "
			"out of supplied vertices and constrained axial and radial "
			"parameters.");

}
////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END();
////////////////////////////////////////////////////////////////////////
/// neuro_collection/test tests
////////////////////////////////////////////////////////////////////////
