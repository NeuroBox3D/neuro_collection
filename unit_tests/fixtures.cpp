/*!
 * \file fixtures.cpp
 * \brief neuro collection fixtures
 *
 *  Created on: May 31, 2019
 *      Author: stephan
 */

#include <boost/test/included/unit_test.hpp>
#include <common/math/ugmath.h>
#include <lib_grid/grid/grid.h>
#include <lib_grid/global_attachments.h>

using namespace ug;

/*!
 * \brief grid fixture
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
	 * \param[in] axial
	 */
	FixtureOneGrid(number axial=0.0, size_t si=0, number radial=1.0) : axial(axial), si (si), radial(radial) {
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
};
