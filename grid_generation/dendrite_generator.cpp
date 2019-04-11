/*
 * dendrite_generator.cpp
 *
 *  Created on: 09.08.2017
 *      Author: mbreit
 */

#include "dendrite_generator.h"

#include <vector>                                                // for vector

#include "common/error.h"                                        // for UG_THROW
#include "common/log.h"                                          // for UG_LOGN
#include "common/math/misc/math_constants.h"                     // for PI
#include "common/math/ugmath_types.h"                            // for vector3
#include "common/util/file_util.h"                               // for FindDirInStandardPaths ...
#include "common/util/string_util.h"                             // for GetFilenameExtension...
#include "lib_grid/common_attachments.h"                         // for aPosition
#include "lib_grid/algorithms/extrusion/extrude.h"               // for Extrude
#include "lib_grid/algorithms/subset_util.h"                     // for AssignSubsetColors, EraseEmptySubsets
#include "lib_grid/file_io/file_io_ugx.h"                        // for GridWriterUGX
#include "lib_grid/grid/grid.h"                                  // for Grid
#include "lib_grid/grid/grid_base_objects.h"                     // for Face, Edge, EdgeDescriptor
#include "lib_grid/grid_objects/grid_objects_0d.h"               // for RegularVertex
#include "lib_grid/tools/subset_handler_grid.h"                  // for SubsetHandler


namespace ug {
namespace neuro_collection {


DendriteGenerator::DendriteGenerator()
: m_dendrite_length(50.0),
  m_dendrite_radius(0.5),
  m_er_radius(0.158),
  m_synapse_width(0.5),
  m_numSegments(100),
  m_bNumSegSet(false),
  m_bBobbelER(false),
  m_numERBlockSegments(1),
  m_numHoleBlockSegments(1),
  CYT_SI(0),
  ER_SI(1),
  PM_SI(2),
  ERM_SI(3),
  RYR_SI(4),
  SYN_SI(5),
  BND_CYT_SI(6),
  BND_ER_SI(7)
{}

void DendriteGenerator::set_dendrite_length(number l)
{
	m_dendrite_length = l;
}

void DendriteGenerator::set_dendrite_radius(number r)
{
	// recover synapse area (mod 2pi)
	number synArea = m_synapse_width * m_dendrite_radius;

	// set new radius
	m_dendrite_radius = r;

	// adapt synapse width
	m_synapse_width = synArea / m_dendrite_radius;
}

void DendriteGenerator::set_er_radius(number r)
{
	m_er_radius = r;
}

void DendriteGenerator::set_synapse_area(number a)
{
	m_synapse_width = a / (2*PI*m_dendrite_radius);
}

void DendriteGenerator::set_num_segments(size_t n)
{
	m_numSegments = 2 * (n/2); // number of segments must be even
	m_bNumSegSet = true;
}

number DendriteGenerator::num_segments() const
{
	return m_numSegments;
}

void DendriteGenerator::set_bobbel_er(size_t erSegLength, size_t holeSegLength)
{
	m_bBobbelER = true;
	m_numERBlockSegments = erSegLength;
	m_numHoleBlockSegments = holeSegLength;


}


void DendriteGenerator::create_dendrite_middle_influx(const std::string& filename)
{
	typedef Grid::VertexAttachmentAccessor<APosition> AAPosition;

	// if #segments has not been explicitly chosen, calculate a sensible number
	if (!m_bNumSegSet)
		m_numSegments = (size_t) floor(m_dendrite_length / m_dendrite_radius);

	// number of segments must be multiple of 2
	// or even of 2*periodLength in case of bobbel ER
	if (m_bBobbelER)
	{
		size_t periodLength = m_numERBlockSegments + m_numHoleBlockSegments;
		m_numSegments = 2*periodLength * m_numSegments / (2 * periodLength);

		m_numSegments = std::max(m_numSegments, 2*periodLength);
	}
	else
	{
		m_numSegments = 2 * m_numSegments / 2;
		m_numSegments = std::max(m_numSegments, (size_t) 2);
	}

	// create grid etc.
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	AAPosition aaPos = AAPosition(g, aPosition);

	// create start vertices and edges at left end
	Vertex* vTmpBndIn = *g.create<RegularVertex>();
	aaPos[vTmpBndIn] = vector3(-0.5*m_dendrite_length, 0, 0);
	int siTmpBndIn = ER_SI;

	Vertex* vTmpERM = *g.create<RegularVertex>();
	aaPos[vTmpERM] = vector3(-0.5*m_dendrite_length, m_er_radius, 0);
	int siTmpERM = ERM_SI;

	Vertex* vTmpPM = *g.create<RegularVertex>();
	aaPos[vTmpPM] = vector3(-0.5*m_dendrite_length, m_dendrite_radius, 0);
	int siTmpPM = PM_SI;

	Edge* eTmpER = *g.create<RegularEdge>(EdgeDescriptor(vTmpBndIn, vTmpERM));
	int siTmpER = BND_ER_SI;

	Edge* eTmpCyt = *g.create<RegularEdge>(EdgeDescriptor(vTmpERM, vTmpPM));
	int siTmpCyt = BND_CYT_SI;

// extrude left half
	std::vector<Vertex*> vrts(3);
	vrts[0] = vTmpBndIn;
	vrts[1] = vTmpERM;
	vrts[2] = vTmpPM;
	std::vector<Edge*> edges(2);
	edges[0] = eTmpER;
	edges[1] = eTmpCyt;

	vector3 extrudeDir(0.0);
	extrudeDir.coord(0) = (m_dendrite_length - m_synapse_width) / m_numSegments;

	if (!m_bBobbelER)
	{
		// prepare subsets for extrusion of ER blocks
		sh.assign_subset(vrts[0], ER_SI);
		sh.assign_subset(vrts[1], ERM_SI);
		sh.assign_subset(vrts[2], PM_SI);
		sh.assign_subset(edges[0], ER_SI);
		sh.assign_subset(edges[1], CYT_SI);

		// extrude
		for (size_t i = 0; i < m_numSegments/2; ++i)
			Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

		// repair original subsets
		sh.assign_subset(vTmpBndIn, siTmpBndIn);
		sh.assign_subset(vTmpERM, siTmpERM);
		sh.assign_subset(vTmpPM, siTmpPM);
		sh.assign_subset(eTmpER, siTmpER);
		sh.assign_subset(eTmpCyt, siTmpCyt);

		// remember correct subsets for current end
		vTmpBndIn = vrts[0];
		vTmpERM = vrts[1];
		vTmpPM = vrts[2];
		eTmpER = edges[0];
		eTmpCyt = edges[1];
		siTmpBndIn = ER_SI;
		siTmpERM = ERM_SI;
		siTmpPM = SYN_SI;
		siTmpER = ER_SI;
		siTmpCyt = CYT_SI;
	}
	else
	{
		size_t numPeriods =  m_numSegments / (2*(m_numERBlockSegments + m_numHoleBlockSegments));
		for (size_t i = 0; i < numPeriods; ++i)
		{
			// prepare subsets for extrusion of ER blocks
			sh.assign_subset(vrts[0], ER_SI);
			sh.assign_subset(vrts[1], ERM_SI);
			sh.assign_subset(vrts[2], PM_SI);
			sh.assign_subset(edges[0], ER_SI);
			sh.assign_subset(edges[1], CYT_SI);

			// extrude a block of ER
			for (size_t j = 0; j < m_numERBlockSegments; ++j)
				Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

			// repair original subsets
			sh.assign_subset(vTmpBndIn, siTmpBndIn);
			sh.assign_subset(vTmpERM, siTmpERM);
			sh.assign_subset(vTmpPM, siTmpPM);
			sh.assign_subset(eTmpER, siTmpER);
			sh.assign_subset(eTmpCyt, siTmpCyt);

			// remember correct subsets for current end
			vTmpBndIn = vrts[0];
			vTmpERM = vrts[1];
			vTmpPM = vrts[2];
			eTmpER = edges[0];
			eTmpCyt = edges[1];
			siTmpBndIn = ERM_SI;
			siTmpERM = ERM_SI;
			siTmpPM = PM_SI;
			siTmpER = ERM_SI;
			siTmpCyt = CYT_SI;
			

			// prepare subsets for extrusion of hole blocks
			sh.assign_subset(vrts[0], CYT_SI);
			sh.assign_subset(vrts[1], CYT_SI);
			sh.assign_subset(vrts[2], PM_SI);
			sh.assign_subset(edges[0], CYT_SI);
			sh.assign_subset(edges[1], CYT_SI);
			
			// extrude a block of hole
			for (size_t j = 0; j < m_numERBlockSegments; ++j)
				Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

			// repair original subsets
			sh.assign_subset(vTmpBndIn, siTmpBndIn);
			sh.assign_subset(vTmpERM, siTmpERM);
			sh.assign_subset(vTmpPM, siTmpPM);
			sh.assign_subset(eTmpER, siTmpER);
			sh.assign_subset(eTmpCyt, siTmpCyt);

			// remember correct subsets for current end
			vTmpBndIn = vrts[0];
			vTmpERM = vrts[1];
			vTmpPM = vrts[2];
			eTmpER = edges[0];
			eTmpCyt = edges[1];
			siTmpBndIn = ERM_SI;
			siTmpERM = ERM_SI;
			siTmpPM = PM_SI;
			siTmpER = ERM_SI;
			siTmpCyt = CYT_SI;
		}

		vTmpBndIn = vrts[0];
		vTmpERM = vrts[1];
		vTmpPM = vrts[2];
		eTmpER = edges[0];
		eTmpCyt = edges[1];
		siTmpBndIn = ERM_SI;
		siTmpERM = ERM_SI;
		siTmpPM = SYN_SI;
		siTmpER = ERM_SI;
		siTmpCyt = CYT_SI;
	}


// make synapse
	size_t nExtr = std::max((size_t) 1, (size_t) floor(m_synapse_width / (m_dendrite_radius - m_er_radius)));
	extrudeDir.coord(0) = m_synapse_width / nExtr;

	// prepare subsets for extrusion of synapse blocks
	sh.assign_subset(vrts[0], ER_SI);
	sh.assign_subset(vrts[1], ERM_SI);
	sh.assign_subset(vrts[2], SYN_SI);
	sh.assign_subset(edges[0], ER_SI);
	sh.assign_subset(edges[1], CYT_SI);

	for (size_t i = 0; i < nExtr; ++i)
		Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

	// repair original subsets
	sh.assign_subset(vTmpBndIn, siTmpBndIn);
	sh.assign_subset(vTmpERM, siTmpERM);
	sh.assign_subset(vTmpPM, siTmpPM);
	sh.assign_subset(eTmpER, siTmpER);
	sh.assign_subset(eTmpCyt, siTmpCyt);

	// remember correct subsets for current end
	vTmpBndIn = vrts[0];
	vTmpERM = vrts[1];
	vTmpPM = vrts[2];
	eTmpER = edges[0];
	eTmpCyt = edges[1];
	if (!m_bBobbelER)
	{
		siTmpBndIn = ER_SI;
		siTmpERM = ERM_SI;
		siTmpPM = SYN_SI;
		siTmpER = ER_SI;
		siTmpCyt = CYT_SI;
	}
	else
	{
		siTmpBndIn = ERM_SI;
		siTmpERM = ERM_SI;
		siTmpPM = SYN_SI;
		siTmpER = ERM_SI;
		siTmpCyt = CYT_SI;
	}


// extrude right half
	extrudeDir.coord(0) = (m_dendrite_length - m_synapse_width) / m_numSegments;

	if (!m_bBobbelER)
	{
		// prepare subsets for extrusion of ER blocks
		sh.assign_subset(vrts[0], ER_SI);
		sh.assign_subset(vrts[1], ERM_SI);
		sh.assign_subset(vrts[2], PM_SI);
		sh.assign_subset(edges[0], ER_SI);
		sh.assign_subset(edges[1], CYT_SI);

		// extrude
		for (size_t i = 0; i < m_numSegments/2; ++i)
			Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

		// repair original subsets
		sh.assign_subset(vTmpBndIn, siTmpBndIn);
		sh.assign_subset(vTmpERM, siTmpERM);
		sh.assign_subset(vTmpPM, siTmpPM);
		sh.assign_subset(eTmpER, siTmpER);
		sh.assign_subset(eTmpCyt, siTmpCyt);
	}
	else
	{
		size_t numPeriods =  m_numSegments / (2*(m_numERBlockSegments + m_numHoleBlockSegments));
		for (size_t i = 0; i < numPeriods; ++i)
		{
			// prepare subsets for extrusion of hole blocks
			sh.assign_subset(vrts[0], CYT_SI);
			sh.assign_subset(vrts[1], CYT_SI);
			sh.assign_subset(vrts[2], PM_SI);
			sh.assign_subset(edges[0], CYT_SI);
			sh.assign_subset(edges[1], CYT_SI);

			// extrude a block of hole
			for (size_t j = 0; j < m_numERBlockSegments; ++j)
				Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

			// repair original subsets
			sh.assign_subset(vTmpBndIn, siTmpBndIn);
			sh.assign_subset(vTmpERM, siTmpERM);
			sh.assign_subset(vTmpPM, siTmpPM);
			sh.assign_subset(eTmpER, siTmpER);
			sh.assign_subset(eTmpCyt, siTmpCyt);

			// remember correct subsets for current end
			vTmpBndIn = vrts[0];
			vTmpERM = vrts[1];
			vTmpPM = vrts[2];
			eTmpER = edges[0];
			eTmpCyt = edges[1];
			siTmpBndIn = ERM_SI;
			siTmpERM = ERM_SI;
			siTmpPM = PM_SI;
			siTmpER = ERM_SI;
			siTmpCyt = CYT_SI;


			// prepare subsets for extrusion of ER blocks
			sh.assign_subset(vrts[0], ER_SI);
			sh.assign_subset(vrts[1], ERM_SI);
			sh.assign_subset(vrts[2], PM_SI);
			sh.assign_subset(edges[0], ER_SI);
			sh.assign_subset(edges[1], CYT_SI);

			// extrude a block of ER
			for (size_t j = 0; j < m_numERBlockSegments; ++j)
				Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

			// repair original subsets
			sh.assign_subset(vTmpBndIn, siTmpBndIn);
			sh.assign_subset(vTmpERM, siTmpERM);
			sh.assign_subset(vTmpPM, siTmpPM);
			sh.assign_subset(eTmpER, siTmpER);
			sh.assign_subset(eTmpCyt, siTmpCyt);

			// remember correct subsets for current end
			vTmpBndIn = vrts[0];
			vTmpERM = vrts[1];
			vTmpPM = vrts[2];
			eTmpER = edges[0];
			eTmpCyt = edges[1];
			siTmpBndIn = ERM_SI;
			siTmpERM = ERM_SI;
			siTmpPM = PM_SI;
			siTmpER = ERM_SI;
			siTmpCyt = CYT_SI;
		}
	}

	// correct subsets at end
	sh.assign_subset(vrts[0], ER_SI);
	sh.assign_subset(vrts[1], ERM_SI);
	sh.assign_subset(vrts[2], PM_SI);
	sh.assign_subset(edges[0], BND_ER_SI);
	sh.assign_subset(edges[1], BND_CYT_SI);


	// subset names
	sh.set_subset_name("cyt", CYT_SI);
	sh.set_subset_name("er", ER_SI);
	sh.set_subset_name("pm", PM_SI);
	sh.set_subset_name("erm", ERM_SI);
	sh.set_subset_name("syn", SYN_SI);
	sh.set_subset_name("bnd_cyt", BND_CYT_SI);
	sh.set_subset_name("bnd_er", BND_ER_SI);
	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);

	// save to .ugx file
	// check that filename ends in ".ugx"
	std::string useFileName = filename;
	if (GetFilenameExtension(filename) != std::string("ugx"))
	{
		UG_LOGN("File name extension needs to be '.ugx' - appending extension.")
		useFileName.append(".ugx");
	}
	std::string filePath = FindDirInStandardPaths(PathFromFilename(useFileName).c_str());
	if (filePath.empty())
		UG_THROW("Directory '" << PathFromFilename(useFileName) << "' could not be located. "
				"The file cannot be written there.");

	std::string fileName = FilenameWithoutPath(useFileName);
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	if (!ugxWriter.write_to_file((filePath+fileName).c_str()))
		UG_THROW("Grid could not be written to file '" << filePath << fileName << "'.");
}



void DendriteGenerator::create_dendrite(const std::string& filename)
{
	typedef Grid::VertexAttachmentAccessor<APosition> AAPosition;

	// if #segments has not been explicitly chosen, calculate a sensible number
	if (!m_bNumSegSet)
		m_numSegments = (size_t) floor(m_dendrite_length / m_dendrite_radius);

	// number of segments must be multiple of 2
	// or even of 2*periodLength in case of bobbel ER
	if (m_bBobbelER)
	{
		size_t periodLength = m_numERBlockSegments + m_numHoleBlockSegments;
		m_numSegments = 2*periodLength * m_numSegments / (2 * periodLength);

		m_numSegments = std::max(m_numSegments, 2*periodLength);
	}
	else
	{
		m_numSegments = 2 * m_numSegments / 2;
		m_numSegments = std::max(m_numSegments, (size_t) 2);
	}

	// create grid etc.
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	AAPosition aaPos = AAPosition(g, aPosition);

	// create start vertices and edges at left end
	Vertex* vTmpBndIn = *g.create<RegularVertex>();
	aaPos[vTmpBndIn] = vector3(-0.5*m_dendrite_length, 0, 0);
	int siTmpBndIn = ER_SI;

	Vertex* vTmpERM = *g.create<RegularVertex>();
	aaPos[vTmpERM] = vector3(-0.5*m_dendrite_length, m_er_radius, 0);
	int siTmpERM = ERM_SI;

	Vertex* vTmpPM = *g.create<RegularVertex>();
	aaPos[vTmpPM] = vector3(-0.5*m_dendrite_length, m_dendrite_radius, 0);
	int siTmpPM = PM_SI;

	Edge* eTmpER = *g.create<RegularEdge>(EdgeDescriptor(vTmpBndIn, vTmpERM));
	int siTmpER = ER_SI;

	Edge* eTmpCyt = *g.create<RegularEdge>(EdgeDescriptor(vTmpERM, vTmpPM));
	int siTmpCyt = SYN_SI;

// extrude left half
	std::vector<Vertex*> vrts(3);
	vrts[0] = vTmpBndIn;
	vrts[1] = vTmpERM;
	vrts[2] = vTmpPM;
	std::vector<Edge*> edges(2);
	edges[0] = eTmpER;
	edges[1] = eTmpCyt;

	vector3 extrudeDir(0.0);
	extrudeDir.coord(0) = m_dendrite_length / m_numSegments;

	if (!m_bBobbelER)
	{
		// prepare subsets for extrusion of ER blocks
		sh.assign_subset(vrts[0], ER_SI);
		sh.assign_subset(vrts[1], ERM_SI);
		sh.assign_subset(vrts[2], PM_SI);
		sh.assign_subset(edges[0], ER_SI);
		sh.assign_subset(edges[1], CYT_SI);

		// extrude
		for (size_t i = 0; i < m_numSegments; ++i)
			Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

		// repair original subsets
		sh.assign_subset(vTmpBndIn, siTmpBndIn);
		sh.assign_subset(vTmpERM, siTmpERM);
		sh.assign_subset(vTmpPM, siTmpPM);
		sh.assign_subset(eTmpER, siTmpER);
		sh.assign_subset(eTmpCyt, siTmpCyt);

		// correct subsets at end
		sh.assign_subset(vrts[0], ER_SI);
		sh.assign_subset(vrts[1], ERM_SI);
		sh.assign_subset(vrts[2], PM_SI);
		sh.assign_subset(edges[0], ER_SI);
		sh.assign_subset(edges[1], BND_CYT_SI);
	}
	else
	{
		size_t numPeriods =  m_numSegments / (m_numERBlockSegments + m_numHoleBlockSegments);
		for (size_t i = 0; i < numPeriods; ++i)
		{
			// prepare subsets for extrusion of ER blocks
			sh.assign_subset(vrts[0], ER_SI);
			sh.assign_subset(vrts[1], ERM_SI);
			sh.assign_subset(vrts[2], PM_SI);
			sh.assign_subset(edges[0], ER_SI);
			sh.assign_subset(edges[1], CYT_SI);

			// extrude a block of ER
			for (size_t j = 0; j < m_numERBlockSegments; ++j)
				Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

			// repair original subsets
			sh.assign_subset(vTmpBndIn, siTmpBndIn);
			sh.assign_subset(vTmpERM, siTmpERM);
			sh.assign_subset(vTmpPM, siTmpPM);
			sh.assign_subset(eTmpER, siTmpER);
			sh.assign_subset(eTmpCyt, siTmpCyt);

			// remember correct subsets for current end
			vTmpBndIn = vrts[0];
			vTmpERM = vrts[1];
			vTmpPM = vrts[2];
			eTmpER = edges[0];
			eTmpCyt = edges[1];
			siTmpBndIn = ERM_SI;
			siTmpERM = ERM_SI;
			siTmpPM = PM_SI;
			siTmpER = ERM_SI;
			siTmpCyt = CYT_SI;


			// prepare subsets for extrusion of hole blocks
			sh.assign_subset(vrts[0], CYT_SI);
			sh.assign_subset(vrts[1], CYT_SI);
			sh.assign_subset(vrts[2], PM_SI);
			sh.assign_subset(edges[0], CYT_SI);
			sh.assign_subset(edges[1], CYT_SI);

			// extrude a block of hole
			for (size_t j = 0; j < m_numERBlockSegments; ++j)
				Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

			// repair original subsets
			sh.assign_subset(vTmpBndIn, siTmpBndIn);
			sh.assign_subset(vTmpERM, siTmpERM);
			sh.assign_subset(vTmpPM, siTmpPM);
			sh.assign_subset(eTmpER, siTmpER);
			sh.assign_subset(eTmpCyt, siTmpCyt);

			// remember correct subsets for current end
			vTmpBndIn = vrts[0];
			vTmpERM = vrts[1];
			vTmpPM = vrts[2];
			eTmpER = edges[0];
			eTmpCyt = edges[1];
			siTmpBndIn = ERM_SI;
			siTmpERM = ERM_SI;
			siTmpPM = PM_SI;
			siTmpER = ERM_SI;
			siTmpCyt = CYT_SI;
		}

		// correct subsets at end
		sh.assign_subset(vrts[0], CYT_SI);
		sh.assign_subset(vrts[1], CYT_SI);
		sh.assign_subset(vrts[2], PM_SI);
		sh.assign_subset(edges[0], CYT_SI);
		sh.assign_subset(edges[1], BND_CYT_SI);
	}


	// subset names
	sh.set_subset_name("cyt", CYT_SI);
	sh.set_subset_name("er", ER_SI);
	sh.set_subset_name("pm", PM_SI);
	sh.set_subset_name("erm", ERM_SI);
	sh.set_subset_name("act", SYN_SI);
	sh.set_subset_name("meas", BND_CYT_SI);
	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);

	// save to .ugx file
	// check that filename ends in ".ugx"
	std::string useFileName = filename;
	if (GetFilenameExtension(filename) != std::string("ugx"))
	{
		UG_LOGN("File name extension needs to be '.ugx' - appending extension.")
		useFileName.append(".ugx");
	}
	std::string filePath = FindDirInStandardPaths(PathFromFilename(useFileName).c_str());
	if (filePath.empty())
		UG_THROW("Directory '" << PathFromFilename(useFileName) << "' could not be located. "
				"The file cannot be written there.");

	std::string fileName = FilenameWithoutPath(useFileName);
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	if (!ugxWriter.write_to_file((filePath+fileName).c_str()))
		UG_THROW("Grid could not be written to file '" << filePath << fileName << "'.");
}



void DendriteGenerator::create_dendrite_1d(const std::string& filename)
{
	typedef Grid::VertexAttachmentAccessor<APosition> AAPosition;

	// if #segments has not been explicitly chosen, calculate a sensible number
	if (!m_bNumSegSet)
		m_numSegments = (size_t) floor(m_dendrite_length / m_dendrite_radius);

	// number of segments must be multiple of 2
	m_numSegments = 2 * m_numSegments / 2;
	m_numSegments = std::max(m_numSegments, (size_t) 2);

	// create grid etc.
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	AAPosition aaPos = AAPosition(g, aPosition);

	// create start vertiex at left end
	Vertex* vTmpBndIn = *g.create<RegularVertex>();
	aaPos[vTmpBndIn] = vector3(-0.5*m_dendrite_length, 0, 0);
	sh.assign_subset(vTmpBndIn, 0);

	// extrude
	std::vector<Vertex*> vrts(1);
	vrts[0] = vTmpBndIn;
	vector3 extrudeDir(0.0);
	extrudeDir.coord(0) = m_dendrite_length / m_numSegments;
	for (size_t i = 0; i < m_numSegments; ++i)
		Extrude(g, &vrts, NULL, NULL, extrudeDir, aaPos, 0, NULL);

	// correct subsets at both ends
	sh.assign_subset(vTmpBndIn, 1);
	sh.assign_subset(vrts[0], 2);

	// subset names
	sh.set_subset_name("dend", 0);
	sh.set_subset_name("act", 1);
	sh.set_subset_name("meas", 2);
	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);

	// save to .ugx file
	// check that filename ends in ".ugx"
	std::string useFileName = filename;
	if (GetFilenameExtension(filename) != std::string("ugx"))
	{
		UG_LOGN("File name extension needs to be '.ugx' - appending extension.")
		useFileName.append(".ugx");
	}
	std::string filePath = FindDirInStandardPaths(PathFromFilename(useFileName).c_str());
	if (filePath.empty())
		UG_THROW("Directory '" << PathFromFilename(useFileName) << "' could not be located. "
				"The file cannot be written there.");

	std::string fileName = FilenameWithoutPath(useFileName);
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	if (!ugxWriter.write_to_file((filePath+fileName).c_str()))
		UG_THROW("Grid could not be written to file '" << filePath << fileName << "'.");
}


void DendriteGenerator::create_dendrite_discreteRyR(const std::string& filename, number channelDistance)
{
	typedef Grid::VertexAttachmentAccessor<APosition> AAPosition;

	const size_t nSegMin = std::ceil(m_dendrite_length / channelDistance);
	if (m_numSegments < nSegMin)
		m_numSegments = nSegMin;

	const size_t ryrElemDist = round(m_numSegments / nSegMin);
	const number segLength = channelDistance / ryrElemDist;
	m_numSegments = ryrElemDist * nSegMin;

	// create grid etc.
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	AAPosition aaPos = AAPosition(g, aPosition);

	// create start vertices and edges at left end
	Vertex* vTmpBndIn = *g.create<RegularVertex>();
	aaPos[vTmpBndIn] = vector3(0, 0, 0);

	Vertex* vTmpERM = *g.create<RegularVertex>();
	aaPos[vTmpERM] = vector3(0, m_er_radius, 0);

	Vertex* vTmpPM = *g.create<RegularVertex>();
	aaPos[vTmpPM] = vector3(0, m_dendrite_radius, 0);

	Edge* eTmpER = *g.create<RegularEdge>(EdgeDescriptor(vTmpBndIn, vTmpERM));
	Edge* eTmpCyt = *g.create<RegularEdge>(EdgeDescriptor(vTmpERM, vTmpPM));

	std::vector<Vertex*> vrts(3);
	vrts[0] = vTmpBndIn;
	vrts[1] = vTmpERM;
	vrts[2] = vTmpPM;
	std::vector<Edge*> edges(2);
	edges[0] = eTmpER;
	edges[1] = eTmpCyt;

	vector3 extrudeDir(0.0);
	extrudeDir.coord(0) = segLength;

	// prepare subsets for extrusion of ER blocks
	sh.assign_subset(vrts[0], ER_SI);
	sh.assign_subset(vrts[1], ERM_SI);
	sh.assign_subset(vrts[2], PM_SI);
	sh.assign_subset(edges[0], ER_SI);
	sh.assign_subset(edges[1], CYT_SI);

	for (size_t n = 0; n < nSegMin; ++n)
	{
		// extrude
		for (size_t i = 0; i < ryrElemDist; ++i)
			Extrude(g, &vrts, &edges, NULL, extrudeDir, aaPos, EO_CREATE_FACES, NULL);

		// set RyR subset
		sh.assign_subset(vTmpERM, RYR_SI);
		vTmpERM = vrts[1];
	}

	// correct subsets at ends
	sh.assign_subset(eTmpCyt, SYN_SI);

	sh.assign_subset(vrts[0], ER_SI);
	sh.assign_subset(vrts[1], RYR_SI);
	sh.assign_subset(vrts[2], PM_SI);
	sh.assign_subset(edges[0], ER_SI);
	sh.assign_subset(edges[1], BND_CYT_SI);

	// subset names
	sh.set_subset_name("cyt", CYT_SI);
	sh.set_subset_name("er", ER_SI);
	sh.set_subset_name("pm", PM_SI);
	sh.set_subset_name("erm", ERM_SI);
	sh.set_subset_name("ryr", RYR_SI);
	sh.set_subset_name("act", SYN_SI);
	sh.set_subset_name("meas", BND_CYT_SI);
	EraseEmptySubsets(sh);
	AssignSubsetColors(sh);


	// save to .ugx file
	// check that filename ends in ".ugx"
	std::string useFileName = filename;
	if (GetFilenameExtension(filename) != std::string("ugx"))
	{
		UG_LOGN("File name extension needs to be '.ugx' - appending extension.")
		useFileName.append(".ugx");
	}
	std::string filePath = FindDirInStandardPaths(PathFromFilename(useFileName).c_str());
	if (filePath.empty())
		UG_THROW("Directory '" << PathFromFilename(useFileName) << "' could not be located. "
				"The file cannot be written there.");

	std::string fileName = FilenameWithoutPath(useFileName);
	GridWriterUGX ugxWriter;
	ugxWriter.add_grid(g, "defGrid", aPosition);
	ugxWriter.add_subset_handler(sh, "defSH", 0);
	if (!ugxWriter.write_to_file((filePath+fileName).c_str()))
		UG_THROW("Grid could not be written to file '" << filePath << fileName << "'.");
}



} // namespace neuro_collection
} // namespace ug
