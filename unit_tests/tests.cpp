/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-05-28
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

#define BOOST_TEST_MODULE __CPP__UNIT_TESTS__UG__NEURO_COLLECTION__TESTS__

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <common/math/ugmath.h>
#include <lib_grid/grid/grid.h>

#include "../test/neurite_util.h"
#include "fixtures.cpp"

using namespace ug;
using namespace ug::neuro_collection;

////////////////////////////////////////////////////////////////////////////////
/// neuro_collection/test tests
////////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(test);
////////////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(CreatePyramid, FixtureOneGrid) {
	Pyramid* p = create_pyramid(g, quad, aaPos, 1.0, NULL);
	BOOST_REQUIRE_MESSAGE(p, "Creating pyramid out of supplied vertices "
			"and created top vertex.");
	vector3 top = aaPos[p->vertex(4)]; //!< last vertex is top vertex
	BOOST_CHECK_EQUAL(top.x(), 0.5);
	BOOST_CHECK_EQUAL(top.y(), 0.5);
	BOOST_CHECK_EQUAL(top.z(), -1.0);
}

////////////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(FindQuadrilateralConstrained, FixtureOneGrid) {
	// find auadrilateral and checks
	Grid::traits<Quadrilateral>::secure_container quadCont;
	find_quadrilaterals_constrained(g, aaSurfParams, quadCont);
	BOOST_REQUIRE_MESSAGE(quadCont.size() == 1, "Finding quadrilateral "
			"out of supplied vertices and constrained axial and radial "
			"parameters.");
}

////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(ExtendERintoSoma, FixtureOneGrid) {
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

////////////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(MergeTwoGrids, FixtureTwoGrid) {
	// check input for consistency
	BOOST_REQUIRE_MESSAGE(sh1.num<Vertex>(0) == 4, "Requiring 4 vertices in subset 0 in first grid");
	BOOST_REQUIRE_MESSAGE(sh1.num<Quadrilateral>() == 1, "Requiring 1 quadrilateral in subset 0 in first grid");
	BOOST_REQUIRE_MESSAGE(sh2.num<Vertex>(1) == 4, "Requiring 4 vertices in subset 1 of second grid");
	BOOST_REQUIRE_MESSAGE(sh2.num<Quadrilateral>() == 1, "Requiring 1 quadrilateral in subset 1 in first grid");

	// merge grids and check output for basic consistency
    typedef NeuriteProjector::SurfaceParams NPSP;
    UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
            "GlobalAttachment 'npSurfParams' not declared.");
    Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	MergeFirstGrids<Attachment<NPSP> >(g1, g2, sh1, sh2, aSP, true);
	BOOST_REQUIRE_MESSAGE(g1.num_vertices() == 8, "Requiring 8 vertices in merged grid (first grid)");
	BOOST_REQUIRE_MESSAGE(g1.num_edges() == 8, "Requiring 8 edges in merged grid (first grid)");

	// additional checks
	BOOST_CHECK_MESSAGE(sh1.num<Vertex>(0) == 4, "Checking 4 vertices in subset 0 of merged grid");
	BOOST_CHECK_MESSAGE(sh1.num<Vertex>(1) == 4, "Checking 4 vertices in subset 1 of merged grid");
	BOOST_CHECK_MESSAGE(sh1.num<Quadrilateral>() == 2, "Checking 2 quadrilaterals in merged grid");
}

////////////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(CopyTwoGrids, FixtureOneGrid) {
	Grid g2;
	SubsetHandler sh2(g2);
	SubsetHandler sh1(g);
	AInt aInt;
    Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
    UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
            "GlobalAttachment 'npSurfParams' not declared.");
    Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	g2.attach_to_vertices(aPosition);
    CopyGrid<APosition>(g, g2, sh1, sh2, aPosition);
	g2.attach_to_vertices(aSP);
	aaSurfParams.access(g, aSP);
	// require that copied attachment values from grid g1 (axial and radial) matches values in copied grid g2
	// Grid g1 has value 0.5 for axial parameter and value 0.1 for radial parameter
    for (VertexIterator iter = g2.vertices_begin(); iter != g2.vertices_end(); ++iter) {
    	BOOST_REQUIRE_MESSAGE(aaSurfParams[*iter].axial == 0.0, "Requiring axial 0.0 in copied grid");
    	BOOST_REQUIRE_MESSAGE(aaSurfParams[*iter].radial == 1.0, "Requiring radial 1.0 in copied grid");
    }
}

////////////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(SelectElementsInUnitSphere, FixtureSphere) {
	Selector sel(g);
	SelectElementsInSphere<ug::Vertex>(g, sel, ug::vector3(0, 0, 0), 1, aaPos);
	InvertSelection(sel);
	BOOST_REQUIRE_MESSAGE(g.num<Vertex>() == 12, "Requiring 12 vertices in original grid");
	InvertSelection(sel);
	EraseSelectedObjects(sel);
	BOOST_REQUIRE_MESSAGE(g.num_vertices() == 0, "Requiring empty grid (all vertices erased)");
	BOOST_REQUIRE_MESSAGE(g.num_edges() == 0, "Requiring empty grid (all vertices erased thus no further edge elements)");
	BOOST_REQUIRE_MESSAGE(g.num_faces() == 0, "Requiring empty grid (all vertices erased thus no further face elements)");
	BOOST_REQUIRE_MESSAGE(g.num_volumes() == 0, "Requiring empty grid (all vertices erased thus no further volume elements)");
}

////////////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(SelectElementsInUnitSphereByAxialPosition, FixtureSphereAxial) {
	Selector sel(g);
	BOOST_REQUIRE_MESSAGE(g.num<Face>() == 20, "Requiring 20 faces in total in original grid");

	/// positive test (should select all elements of unit sphere since they have aaSurfParams.axial = -1.0)
	SelectElementsByAxialPosition<Face>(g, sel, 0.0, aaPos, aaSurfParams);
	BOOST_REQUIRE_MESSAGE(sel.num<Face>() == 20, "Requiring 20 faces in total in selection.");
	CloseSelection(sel);
	BOOST_REQUIRE_MESSAGE(sel.num<Vertex>() == 12, "Requiring 12 vertices in total in selection.");

	/// negative test
	sel.clear();
	SelectElementsByAxialPosition<Face>(g, sel, -2.0, aaPos, aaSurfParams);
	BOOST_REQUIRE_MESSAGE(sel.num<Face>() == 0, "Requiring 0 faces in total in selection.");
	CloseSelection(sel);
	BOOST_REQUIRE_MESSAGE(sel.num<Vertex>() == 0, "Requiring 0 vertices in total in selection.");
}

////////////////////////////////////////////////////////////////////////////////
BOOST_FIXTURE_TEST_CASE(DeleteInnerEdgesFromQuadrilateralWithInnerEdges, FixtureQuadrilateralWithInnerEdges) {
	BOOST_REQUIRE_MESSAGE(g.num<Edge>() == 6, "Requiring a quadrilateral with two additional edges from p1->p3 and p2->p4.");
	BOOST_REQUIRE_MESSAGE(g.num<Quadrilateral>() == 1, "Requiring one quadrilateral.");
	DeleteInnerEdgesFromQuadrilaterals(g, sh, 0);
	BOOST_REQUIRE_MESSAGE(g.num<Edge>() == 4, "Requiring four edges from one quadrilateral.");
	BOOST_REQUIRE_MESSAGE(g.num<Quadrilateral>() == 1, "Requiring one quadrilateral.");
}

BOOST_FIXTURE_TEST_CASE(DeleteInnerEdgesFromQuadrilateralWithoutInnerEdges, FixtureQuadrilateralWithoutInnerEdges) {
	BOOST_REQUIRE_MESSAGE(g.num<Edge>() == 4, "Requiring a quadrilateral with no additional inner edges.");
	BOOST_REQUIRE_MESSAGE(g.num<Quadrilateral>() == 1, "Requiring one quadrilateral.");
	DeleteInnerEdgesFromQuadrilaterals(g, sh, 0);
	BOOST_REQUIRE_MESSAGE(g.num<Edge>() == 4, "Requiring four edges from one quadrilateral.");
	BOOST_REQUIRE_MESSAGE(g.num<Quadrilateral>() == 1, "Requiring one quadrilateral.");
}

BOOST_AUTO_TEST_SUITE_END();
////////////////////////////////////////////////////////////////////////////////
/// neuro_collection/test tests
////////////////////////////////////////////////////////////////////////////////
