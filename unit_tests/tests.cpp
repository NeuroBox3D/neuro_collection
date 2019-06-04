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
BOOST_FIXTURE_TEST_CASE(FindQuadrilateralConstrained, FixtureGrid) {
	// find auadrilateral and checks
	Grid::traits<Quadrilateral>::secure_container quadCont;
	find_quadrilaterals_constrained(g, aaSurfParams, quadCont);
	BOOST_REQUIRE_MESSAGE(quadCont.size() == 1, "Finding quadrilateral "
			"out of supplied vertices and constrained axial and radial "
			"parameters.");
}

////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(ExtractSubGrid) {
	Grid g;
	SubsetHandler sh(g);
	g.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	sh.set_default_subset_index(0);
	// create vertices in subset 0
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

	/// create edges in subsets 1, 2, 3 and 4 of square
	ug::Edge *e1, *e2, *e3, *e4;
	sh.assign_subset(*g.create<RegularEdge>(EdgeDescriptor(p1, p2)), 1);
	sh.assign_subset(*g.create<RegularEdge>(EdgeDescriptor(p2, p3)), 2);
	sh.assign_subset(*g.create<RegularEdge>(EdgeDescriptor(p3, p4)), 3);
	sh.assign_subset(*g.create<RegularEdge>(EdgeDescriptor(p4, p1)), 4);

	std::vector<size_t> vSi; vSi.push_back(0);
	Grid gridOut;
	SubsetHandler destSh(gridOut);
	split_grid_based_on_subset_indices(g, sh, gridOut, destSh, aaPos, vSi);
	BOOST_REQUIRE_MESSAGE(gridOut.num_vertices() == 4, "Requiring four vertices.");
	BOOST_REQUIRE_MESSAGE(g.num_edges() == 4, "Requiring four edges.");
}

////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(ExtendERintoSoma, FixtureGrid) {
	// calculate new vertices manually in normal direction
	ug::vector3 newVertex1, newVertex2, newVertex3, newVertex4;
	ug::vector3 normalDir(0, 0, 1);
	VecScaleAdd(newVertex1, 1.0, aaPos[p1], 1.0, normalDir);
	VecScaleAdd(newVertex2, 1.0, aaPos[p2], 1.0, normalDir);
	VecScaleAdd(newVertex3, 1.0, aaPos[p3], 1.0, normalDir);
	VecScaleAdd(newVertex4, 1.0, aaPos[p4], 1.0, normalDir);

	// extrude with method and check for equality
	std::vector<ug::Vertex*> outVerts;
	extend_ER_within(g, sh, aaPos, aaSurfParams, 0, 1, 1.0, outVerts);
	BOOST_REQUIRE_MESSAGE(outVerts.size() == 4, "Requiring four extruded vertices");
	BOOST_REQUIRE_SMALL(VecDistance(aaPos[outVerts[0]], newVertex1), SMALL);
	BOOST_REQUIRE_SMALL(VecDistance(aaPos[outVerts[1]], newVertex2), SMALL);
	BOOST_REQUIRE_SMALL(VecDistance(aaPos[outVerts[2]], newVertex3), SMALL);
	BOOST_REQUIRE_SMALL(VecDistance(aaPos[outVerts[3]], newVertex4), SMALL);
}

BOOST_AUTO_TEST_SUITE_END();
////////////////////////////////////////////////////////////////////////
/// neuro_collection/test tests
////////////////////////////////////////////////////////////////////////
