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
struct FixtureGrid {
	typedef NeuriteProjector::SurfaceParams NPSP;
	Grid g;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
    Grid::VertexAttachmentAccessor<Attachment<NPSP> > aaSurfParams;
	const Quadrilateral* quad;
	number axial;
	size_t si;
	number radial;

	/*!
	 * \brief set up a quadrilateral with optional axial and radial parameters
	 * \param[in] axial
	 */
	FixtureGrid(number axial=0.0, size_t si=0, number radial=1.0) : axial(axial), si (si), radial(radial) {
		// create vertices
		g.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(g, aPosition);
		sh = SubsetHandler(g);
		sh.set_default_subset_index(si);
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
	 * \brief
	 */
	~FixtureGrid() {
		GlobalAttachments::undeclare_attachment<Attachment<NPSP> >("npSurfParams");
	}
};
