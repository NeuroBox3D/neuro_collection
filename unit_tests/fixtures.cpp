/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-03-31
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

#include <boost/test/included/unit_test.hpp>
#include <common/math/ugmath.h>
#include <lib_grid/grid/grid.h>
#include <lib_grid/global_attachments.h>
#include <lib_grid/algorithms/grid_generation/icosahedron.h>

using namespace ug;

/*!
 * \brief empty grid fixture
 */
struct FixtureEmptyGrid {
	Grid g;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	FixtureEmptyGrid() {
		g.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g, aPosition);
		sh = SubsetHandler(g);
		sh.set_default_subset_index(0);
		BOOST_REQUIRE_MESSAGE(g.num<Vertex>() == 0, "Grid not empty.");
	}
};

/*!
 * \brief one grid fixture
 */
struct FixtureOneGrid {
	typedef NeuriteProjector::SurfaceParams NPSP;
	Grid g;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
    Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	const Quadrilateral* quad;
	number axial;
	size_t si;
	number radial;
	ug::Vertex *p1, *p2, *p3, *p4;

	/*!
	 * \brief set up a quadrilateral with optional axial and radial parameters
	 */
	FixtureOneGrid() : axial(0.0), si (0), radial(1.0) {
		// create vertices
		g.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g, aPosition);
		sh = SubsetHandler(g);
		sh.set_default_subset_index(si);
		p1 = *g.create<RegularVertex>(); p2 = *g.create<RegularVertex>();
		p3 = *g.create<RegularVertex>(); p4 = *g.create<RegularVertex>();
		BOOST_REQUIRE_MESSAGE(p1, "Creating first vertex.");
		BOOST_REQUIRE_MESSAGE(p2, "Creating second vertex.");
		BOOST_REQUIRE_MESSAGE(p3, "Creating third vertex.");
		BOOST_REQUIRE_MESSAGE(p4, "Creating forth vertex.");

		// assign coordinates
		aaPos[p1] = ug::vector3(0, 0, 0); aaPos[p2] = ug::vector3(0, 1, 0);
		aaPos[p3] = ug::vector3(1, 1, 0); aaPos[p4] = ug::vector3(1, 0, 0);

		// create quadrilateral
		quad = *g.create<Quadrilateral>(QuadrilateralDescriptor(p1, p2, p3, p4));
		BOOST_REQUIRE_MESSAGE(quad, "Creating quadrilateral out of supplied vertices.");

	    // assign axial and radial arameter
		GlobalAttachments::declare_attachment<Attachment<NPSP> >("npSurfParams", true);
	    UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
	            "GlobalAttachment 'npSurfParams' not declared.");
	    Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	    if (!g.has_vertex_attachment(aSP))
	        g.attach_to_vertices(aSP);
	    aaSurfParams.access(g, aSP);
	    aaSurfParams[p1].axial = axial; aaSurfParams[p2].axial = axial;
	    aaSurfParams[p3].axial = axial; aaSurfParams[p4].axial = axial;

	    BOOST_REQUIRE_MESSAGE(aaSurfParams[p1].axial == axial, "Setting axial parameter for first vertex");
	    BOOST_REQUIRE_MESSAGE(aaSurfParams[p2].axial == axial, "Setting axial parameter for second vertex");
	    BOOST_REQUIRE_MESSAGE(aaSurfParams[p3].axial == axial, "Setting axial parameter for third vertex");
	    BOOST_REQUIRE_MESSAGE(aaSurfParams[p4].axial == axial, "Setting axial parameter for forth vertex");

	    aaSurfParams[p1].radial = radial; aaSurfParams[p2].radial = radial;
	    aaSurfParams[p3].radial = radial; aaSurfParams[p4].radial = radial;

	    BOOST_REQUIRE_MESSAGE(aaSurfParams[p1].radial == radial, "Setting radial parameter for first vertex");
	    BOOST_REQUIRE_MESSAGE(aaSurfParams[p2].radial == radial, "Setting radial parameter for second vertex");
	    BOOST_REQUIRE_MESSAGE(aaSurfParams[p3].radial == radial, "Setting radial parameter for third vertex");
	    BOOST_REQUIRE_MESSAGE(aaSurfParams[p4].radial == radial, "Setting radial parameter for forth vertex");
	}

	/*!
	 * \brief undeclare the npSurfParam attachment
	 */
	~FixtureOneGrid() {
		GlobalAttachments::undeclare_attachment<Attachment<NPSP> >("npSurfParams");
	}
};

/*!
 * \brief two grid fixture
 */
struct FixtureTwoGrid {
	Grid g1, g2;
	SubsetHandler sh1, sh2;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	/*!
	 * \brief creates two grids with non-overlapping quadrilaterals in distinct subsets
	 */
	FixtureTwoGrid() {
		// create vertices for first grid
		g1.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g1, aPosition);
		sh1 = SubsetHandler(g1);
		sh1.set_default_subset_index(0);

		GlobalAttachments::declare_attachment<Attachment<NeuriteProjector::SurfaceParams> >("npSurfParams", true);

		ug::Vertex *p1, *p2, *p3, *p4;
		p1 = *g1.create<RegularVertex>(); p2 = *g1.create<RegularVertex>();
		p3 = *g1.create<RegularVertex>(); p4 = *g1.create<RegularVertex>();
		BOOST_REQUIRE_MESSAGE(p1, "Creating first vertex.");
		BOOST_REQUIRE_MESSAGE(p2, "Creating second vertex.");
		BOOST_REQUIRE_MESSAGE(p3, "Creating third vertex.");
		BOOST_REQUIRE_MESSAGE(p4, "Creating forth vertex.");


		// assign coordinates
		aaPos[p1] = ug::vector3(2, 2, 0); aaPos[p2] = ug::vector3(2, 3, 0);
		aaPos[p3] = ug::vector3(3, 3, 0); aaPos[p4] = ug::vector3(3, 2, 0);
		Quadrilateral* quad = *g1.create<Quadrilateral>(QuadrilateralDescriptor(p1, p2, p3, p4));
		BOOST_REQUIRE_MESSAGE(quad, "Creating quadrilateral out of supplied vertices.");
		sh1.assign_grid(g1);
		sh1.assign_subset(p1, 0);
		sh1.assign_subset(p2, 0);
		sh1.assign_subset(p3, 0);
		sh1.assign_subset(p4, 0);
		sh1.assign_subset(quad, 0);

		// create vertices
		g2.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g2, aPosition);
		sh2 = SubsetHandler(g2);
		sh2.set_default_subset_index(1);
		ug::Vertex *p12, *p22, *p32, *p42;
		p12 = *g2.create<RegularVertex>(); p22 = *g2.create<RegularVertex>();
		p32 = *g2.create<RegularVertex>(); p42 = *g2.create<RegularVertex>();
		BOOST_REQUIRE_MESSAGE(p12, "Creating first vertex.");
		BOOST_REQUIRE_MESSAGE(p22, "Creating second vertex.");
		BOOST_REQUIRE_MESSAGE(p32, "Creating third vertex.");
		BOOST_REQUIRE_MESSAGE(p42, "Creating forth vertex.");

		// assign coordinates
		aaPos[p12] = ug::vector3(0, 0, 0); aaPos[p22] = ug::vector3(0, 1, 0);
		aaPos[p32] = ug::vector3(1, 1, 0); aaPos[p42] = ug::vector3(1, 0, 0);
		Quadrilateral* quad2 = *g2.create<Quadrilateral>(QuadrilateralDescriptor(p12, p22, p32, p42));
		BOOST_REQUIRE_MESSAGE(quad2, "Creating quadrilateral out of supplied vertices.");
		sh2.assign_grid(g2);
		sh2.assign_subset(p12, 1);
		sh2.assign_subset(p22, 1);
		sh2.assign_subset(p32, 1);
		sh2.assign_subset(p42, 1);
		sh2.assign_subset(quad2, 1);
	}
	~FixtureTwoGrid() {
		GlobalAttachments::undeclare_attachment<Attachment<NeuriteProjector::SurfaceParams> >("npSurfParams");
	}
};

/*!
 * \brief unit sphere fixture
 */
struct FixtureSphere {
	Grid g;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	ug::vector3 center;
	number radius;
	/*!
	 * \brief creates a sphere
	 */
	FixtureSphere() : center(ug::vector3(0, 0, 0)), radius(1.0) {
		g.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g, aPosition);
		sh = SubsetHandler(g);
		sh.set_default_subset_index(0);
		Selector sel(g);
		GenerateIcosphere(g, center, radius, 0, aPosition, &sel);
	}
};

/*!
 * \brief unit sphere fixture axial
 */
struct FixtureSphereAxial {
	typedef NeuriteProjector::SurfaceParams NPSP;
	Grid g;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	ug::vector3 center;
	number radius;
    Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	/*!
	 * \brief creates a sphere with axial parameters
	 */
	FixtureSphereAxial() : center(ug::vector3(0, 0, 0)), radius(1.0) {
		g.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g, aPosition);
		sh = SubsetHandler(g);
		sh.set_default_subset_index(0);
		Selector sel(g);
		GenerateIcosphere(g, center, radius, 0, aPosition, &sel);
		GlobalAttachments::declare_attachment<Attachment<NPSP> >("npSurfParams", true);
	    UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
	            "GlobalAttachment 'npSurfParams' not declared.");
	    Attachment<NPSP> aSP = GlobalAttachments::attachment<Attachment<NPSP> >("npSurfParams");
	    if (!g.has_vertex_attachment(aSP))
	        g.attach_to_vertices(aSP);
	    aaSurfParams.access(g, aSP);
	    for(VertexIterator iter = g.vertices_begin(); iter != g.vertices_end(); ++iter) {
	    	aaSurfParams[*iter].axial = -1.0;
	    }
	}
};

/*!
 * \brief quadrilateral with inner edges fixture
 */
struct FixtureQuadrilateralWithInnerEdges {
	Grid g;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	/*!
	 * \brief creates a quadrilateral with two inner edges p1->p3 and p2->p4
	 */
	FixtureQuadrilateralWithInnerEdges() {
		g.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g, aPosition);
		sh = SubsetHandler(g);
		sh.set_default_subset_index(0);
		Vertex* p1, *p2, *p3, *p4;
		p1 = *g.create<RegularVertex>(); p2 = *g.create<RegularVertex>();
		p3 = *g.create<RegularVertex>(); p4 = *g.create<RegularVertex>();
		BOOST_REQUIRE_MESSAGE(p1, "Creating first vertex.");
		BOOST_REQUIRE_MESSAGE(p2, "Creating second vertex.");
		BOOST_REQUIRE_MESSAGE(p3, "Creating third vertex.");
		BOOST_REQUIRE_MESSAGE(p4, "Creating forth vertex.");
		aaPos[p1] = ug::vector3(0, 0, 0); aaPos[p2] = ug::vector3(0, 1, 0);
		aaPos[p3] = ug::vector3(1, 1, 0); aaPos[p4] = ug::vector3(1, 0, 0);
		Quadrilateral* quad = *g.create<Quadrilateral>(QuadrilateralDescriptor(p1, p2, p3, p4));
		BOOST_REQUIRE_MESSAGE(quad, "Creating quadrilateral.");
		RegularEdge* e1 = *g.create<RegularEdge>(EdgeDescriptor(p1, p3));
		RegularEdge* e2 = *g.create<RegularEdge>(EdgeDescriptor(p2, p4));
		BOOST_REQUIRE_MESSAGE(e1, "Creating first additional edge.");
		BOOST_REQUIRE_MESSAGE(e2, "Creating second additional edge.");
		sh.assign_grid(g);
		sh.assign_subset(quad, 0);
		sh.assign_subset(e1, 0);
		sh.assign_subset(e2, 0);
	}
};

/*!
 * \brief quadrilateral without inner edges fixture
 */
struct FixtureQuadrilateralWithoutInnerEdges {
	Grid g;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	/*!
	 * \brief creates a quadrilateral with two inner edges p1->p3 and p2->p4
	 */
	FixtureQuadrilateralWithoutInnerEdges() {
		g.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g, aPosition);
		sh = SubsetHandler(g);
		sh.set_default_subset_index(0);
		Vertex* p1, *p2, *p3, *p4;
		p1 = *g.create<RegularVertex>(); p2 = *g.create<RegularVertex>();
		p3 = *g.create<RegularVertex>(); p4 = *g.create<RegularVertex>();
		BOOST_REQUIRE_MESSAGE(p1, "Creating first vertex.");
		BOOST_REQUIRE_MESSAGE(p2, "Creating second vertex.");
		BOOST_REQUIRE_MESSAGE(p3, "Creating third vertex.");
		BOOST_REQUIRE_MESSAGE(p4, "Creating forth vertex.");
		aaPos[p1] = ug::vector3(0, 0, 0); aaPos[p2] = ug::vector3(0, 1, 0);
		aaPos[p3] = ug::vector3(1, 1, 0); aaPos[p4] = ug::vector3(1, 0, 0);
		Quadrilateral* quad = *g.create<Quadrilateral>(QuadrilateralDescriptor(p1, p2, p3, p4));
		BOOST_REQUIRE_MESSAGE(quad, "Creating quadrilateral.");
		sh.assign_grid(g);
		sh.assign_subset(quad, 0);
	}
};
