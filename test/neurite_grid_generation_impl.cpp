/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2019-05-06
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

#include "neurite_grid_generation.h"
#include "lib_algebra/common/operations_vec.h"
#include "lib_grid/algorithms/extrusion/extrude.h"
#include "lib_grid/grid/neighborhood_util.h"
#include "lib_disc/quadrature/gauss_legendre/gauss_legendre.h"
#include <algorithm>

namespace ug {
	namespace neuro_collection {
	////////////////////////////////////////////////////////////////////////
	/// create_neurite_with_er
	////////////////////////////////////////////////////////////////////////
	void create_neurite_with_er(
		const std::vector<NeuriteProjector::Neurite>& vNeurites,
		const std::vector<std::vector<vector3> >& vPos,
		const std::vector<std::vector<number> >& vR, size_t nid,
		number erScaleFactor, number anisotropy, Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		Grid::VertexAttachmentAccessor<
		Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		SubsetHandler& sh, std::vector<Vertex*>* connectingVrts,
		std::vector<Edge*>* connectingEdges,
		std::vector<Face*>* connectingFaces,
		number initialOffset,
		std::vector<Vertex*>* outVerts,
		std::vector<Vertex*>* outVertsInner,
		std::vector<number>* outRads,
		std::vector<number>* outRadsInner
	) {
	const NeuriteProjector::Neurite& neurite = vNeurites[nid];
	const std::vector<vector3>& pos = vPos[nid];
	const std::vector<number>& r = vR[nid];

	number neurite_length = 0.0;
	for (size_t i = 1; i < pos.size(); ++i)
		neurite_length += VecDistance(pos[i], pos[i - 1]);

	size_t nSec = neurite.vSec.size();

	const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit =
			vBR.begin();
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end =
			vBR.end();

	std::vector<Vertex*> vVrt;
	std::vector<Edge*> vEdge;
	std::vector<Face*> vFace;
	vVrt.resize(16);
	vEdge.resize(24);
	vFace.resize(9);

	vector3 vel;
	const NeuriteProjector::Section& sec = neurite.vSec[0];
	number h = sec.endParam;
	vel[0] = -3.0 * sec.splineParamsX[0] * h * h
			- 2.0 * sec.splineParamsX[1] * h - sec.splineParamsX[2];
	vel[1] = -3.0 * sec.splineParamsY[0] * h * h
			- 2.0 * sec.splineParamsY[1] * h - sec.splineParamsY[2];
	vel[2] = -3.0 * sec.splineParamsZ[0] * h * h
			- 2.0 * sec.splineParamsZ[1] * h - sec.splineParamsZ[2];

	vector3 projRefDir;
	VecNormalize(vel, vel);
	number fac = VecProd(neurite.refDir, vel);
	VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
	VecNormalize(projRefDir, projRefDir);

	vector3 thirdDir;
	VecCross(thirdDir, vel, projRefDir);

	number angleOffset = 0.0;

	number t_start = 0.0;
	number t_end = 0.0;

	if (connectingVrts && connectingEdges && connectingFaces) {
		vVrt = *connectingVrts;
		vEdge = *connectingEdges;
		vFace = *connectingFaces;

		// calculate start angle offset
		vector3 center(0.0);
		for (size_t i = 0; i < 4; ++i)
			VecAdd(center, center, aaPos[(*connectingVrts)[i]]);
		center /= 4;

		vector3 centerToFirst;
		VecSubtract(centerToFirst, aaPos[(*connectingVrts)[0]], center);

		vector2 relCoord;
		VecScaleAdd(centerToFirst, 1.0, centerToFirst,
				-VecProd(centerToFirst, vel), vel);
		relCoord[0] = VecProd(centerToFirst, projRefDir);
		VecScaleAdd(centerToFirst, 1.0, centerToFirst, -relCoord[0],
				projRefDir);
		relCoord[1] = VecProd(centerToFirst, thirdDir);
		VecNormalize(relCoord, relCoord);

		if (fabs(relCoord[0]) < 1e-8)
			angleOffset = relCoord[1] < 0 ? 1.5 * PI : 0.5 * PI;
		else
			angleOffset =
					relCoord[0] < 0 ?
							PI - atan(-relCoord[1] / relCoord[0]) :
							atan(relCoord[1] / relCoord[0]);
		if (angleOffset < 0)
			angleOffset += 2.0 * PI;

		// ignore first branching region (the connecting region)
		++brit;

		// apply initial offset (to prevent first segment being shorter than the others)
		t_end = initialOffset / neurite_length;
		for (size_t i = 0; i < 4; ++i) {
			aaSurfParams[(*connectingVrts)[i]].axial = t_end;
			aaSurfParams[(*connectingVrts)[i]].angular =
					0.5 * PI * i + angleOffset < 2 * PI ?
							0.5 * PI * i + angleOffset :
							0.5 * PI * i + angleOffset - 2 * PI;
			aaSurfParams[(*connectingVrts)[i]].radial = erScaleFactor;
		}
	} else {
		// create first layer of vertices/edges //

		// ER vertices
		for (size_t i = 0; i < 4; ++i) {
			Vertex* v = *g.create<RegularVertex>();
			vVrt[i] = v;
			number angle = 0.5 * PI * i;
			VecScaleAdd(aaPos[v], 1.0, pos[0],
					erScaleFactor * r[0] * cos(angle), projRefDir,
					erScaleFactor * r[0] * sin(angle), thirdDir);

			aaSurfParams[v].neuriteID = nid;
			aaSurfParams[v].axial = 0.0;
			aaSurfParams[v].angular = angle;
			aaSurfParams[v].radial = erScaleFactor;
			sh.assign_subset(v, 3);
			if (outVertsInner) {
				outVertsInner->push_back(v);
			}

			if (outRadsInner) {
				outRadsInner->push_back(erScaleFactor * r[0]);
			}
		}

		for (size_t i = 0; i < 12; ++i) {
			Vertex* v = *g.create<RegularVertex>();
			vVrt[i + 4] = v;
			number angle = PI * ((number) i / 6);
			VecScaleAdd(aaPos[v], 1.0, pos[0], r[0] * cos(angle), projRefDir,
					r[0] * sin(angle), thirdDir);

			aaSurfParams[v].neuriteID = nid;
			aaSurfParams[v].axial = 0.0;
			aaSurfParams[v].angular = angle;
			aaSurfParams[v].radial = 1.0;
			sh.assign_subset(v, 2);
			if (outVerts) {
				outVerts->push_back(v);
			}

			if (outRads) {
				outRads->push_back(r[0]);
			}
		}

		// edges
		for (size_t i = 0; i < 4; ++i) {
			vEdge[i] = *g.create<RegularEdge>(
					EdgeDescriptor(vVrt[i], vVrt[(i + 1) % 4]));
			vEdge[i + 4] = *g.create<RegularEdge>(
					EdgeDescriptor(vVrt[i], vVrt[5 + 3 * i]));
			vEdge[i + 8] = *g.create<RegularEdge>(
					EdgeDescriptor(vVrt[(i + 1) % 4], vVrt[6 + 3 * i]));

			sh.assign_subset(vEdge[i], 3);
			sh.assign_subset(vEdge[i + 4], 0);
		}
		for (size_t i = 0; i < 12; ++i) {
			vEdge[i + 12] = *g.create<RegularEdge>(
					EdgeDescriptor(vVrt[i + 4], vVrt[(i + 1) % 12 + 4]));
			sh.assign_subset(vEdge[i + 12], 2);
		}

		// faces
		vFace[0] = *g.create<Quadrilateral>(
				QuadrilateralDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3]));
		sh.assign_subset(vFace[0], 1);
		for (size_t i = 0; i < 4; ++i) {
			vFace[i + 1] = *g.create<Quadrilateral>(
					QuadrilateralDescriptor(vVrt[i],
							vVrt[(3 * i + 11) % 12 + 4], vVrt[3 * i + 4],
							vVrt[3 * i + 5]));
			vFace[i + 5] = *g.create<Quadrilateral>(
					QuadrilateralDescriptor(vVrt[i], vVrt[3 * i + 5],
							vVrt[3 * i + 6], vVrt[(i + 1) % 4]));
			sh.assign_subset(vFace[i + 1], 0);
			sh.assign_subset(vFace[i + 5], 0);
		}
	}

	// Now create dendrite to the next branching point and iterate this process.
	// We want to create each of the segments with approx. the same aspect ratio.
	// To that end, we first calculate the length of the section to be created (in units of radius)
	// and then divide this number by 2^n where n is the number of anisotropic refinements to
	// be performed to make all segments (more or less) isotropic. The result is the number
	// of segments to be used for the section.
	vector3 lastPos = pos[0];
	size_t curSec = 0;

	while (true) {
		t_start = t_end;

		// if there is another BP, store its extensions here
		number bp_start = 1.0;
		number bp_end = 0.0;

		// initial branch offsets (to prevent the first segment being shorter than the following)
		std::vector<number> branchOffset;
		number surfBPoffset;

		// last section: create until tip
		if (brit == brit_end)
			t_end = 1.0;

		// otherwise: section goes to next branching point
		else {
			// calculate the exact position of the branching point,
			// i.e., the axial position of the intersection of the branching neurite's
			// spline with the surface of the current neurite
			// this is necessary esp. when the branching angle is very small
			const std::vector<uint32_t>& vBranchInd = brit->bp->vNid;
			size_t nBranches = vBranchInd.size();

			branchOffset.resize(brit->bp->vNid.size(), 0.0);
			for (size_t br = 1; br < nBranches; ++br) {
				uint32_t brInd = vBranchInd[br];

				// get position and radius of first point of branch
				const number brRadSeg1 = vR[brInd][0];

				// get position and radius of branching point
				const number bpTPos = brit->t;
				size_t brSec = curSec;
				for (; brSec < nSec; ++brSec) {
					const NeuriteProjector::Section& sec = neurite.vSec[brSec];
					if (bpTPos - sec.endParam < 1e-6 * bpTPos)
						break;
				}
				UG_COND_THROW(brSec == nSec,
						"Could not find section containing branching point "
						"at t = " << bpTPos << ".");
				const number bpRad = vR[nid][brSec + 1];

				// calculate branch and neurite directions
				vector3 branchDir;
				const NeuriteProjector::Section& childSec =
						vNeurites[brInd].vSec[0];
				number te = childSec.endParam;

				const number* s = &childSec.splineParamsX[0];
				number& v0 = branchDir[0];
				v0 = -3.0 * s[0] * te - 2.0 * s[1];
				v0 = v0 * te - s[2];

				s = &childSec.splineParamsY[0];
				number& v1 = branchDir[1];
				v1 = -3.0 * s[0] * te - 2.0 * s[1];
				v1 = v1 * te - s[2];

				s = &childSec.splineParamsZ[0];
				number& v2 = branchDir[2];
				v2 = -3.0 * s[0] * te - 2.0 * s[1];
				v2 = v2 * te - s[2];

				VecNormalize(branchDir, branchDir);

				vector3 neuriteDir;
				const NeuriteProjector::Section& sec = neurite.vSec.at(brSec);
				vel[0] = -sec.splineParamsX[2];
				vel[1] = -sec.splineParamsY[2];
				vel[2] = -sec.splineParamsZ[2];
				number velNorm = sqrt(VecNormSquared(vel));
				VecScale(neuriteDir, vel, 1.0 / velNorm);

				// calculate offset of true branch position, which is r1/sqrt(2)*cot(alpha)
				const number brScProd = VecProd(neuriteDir, branchDir);
				const number sinAlphaInv = 1.0
						/ sqrt(1.0 - brScProd * brScProd);
				surfBPoffset = 0.5 * sqrt(2.0) * bpRad * brScProd * sinAlphaInv;

				// calculate offset of new branch, which is r1/sqrt(2)/sin(alpha)
				branchOffset[br] = 0.5 * sqrt(2.0) * bpRad * sinAlphaInv;

				// calculate true half-length of BP, which is r2/sqrt(2)/sin(alpha)
				number surfBPhalfLength = brRadSeg1 * sinAlphaInv;

				// finally set bp start and end
				bp_start = std::min(bp_start,
						bpTPos - surfBPhalfLength / neurite_length);
				bp_end = std::max(bp_end,
						bpTPos + surfBPhalfLength / neurite_length);
			}

			t_end = bp_start;
		}

		// calculate total length in units of radius
		// = integral from t_start to t_end over: ||v(t)|| / r(t) dt
		number lengthOverRadius = calculate_length_over_radius(t_start, t_end,
				neurite, curSec);

		// to reach the desired anisotropy on the surface in the refinement limit,
		// it has to be multiplied by pi/2 h
		size_t nSeg = (size_t) floor(
				lengthOverRadius / (anisotropy * 0.5 * PI));
		if (!nSeg)
			nSeg = 1;
		number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
		std::vector<number> vSegAxPos(nSeg);
		calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite,
				curSec, segLength);

		// add the branching point to segment list (if present)
		if (brit != brit_end) {
			vSegAxPos.resize(nSeg + 1);
			vSegAxPos[nSeg] = bp_end;
			++nSeg;
		}

		// in case we construct to a BP, find out the branching angle
		// in order to adjust this neurites offset bit by bit
		number addOffset = 0.0;
		size_t child_nid;
		size_t connFaceInd = 0;
		if (brit != brit_end) {
			// find branching child neurite
			SmartPtr<NeuriteProjector::BranchingPoint> bp = brit->bp;
			UG_COND_THROW(bp->vNid.size() > 2,
					"This implementation can only handle branching points with one branching child.");

			if (bp->vNid[0] != nid)
				child_nid = bp->vNid[0];
			else
				child_nid = bp->vNid[1];

			// find out branching child neurite initial direction
			vector3 childDir;
			const NeuriteProjector::Section& childSec =
					vNeurites[child_nid].vSec[0];
			number te = childSec.endParam;

			const number* sp = &childSec.splineParamsX[0];
			number& vc0 = childDir[0];
			vc0 = -3.0 * sp[0] * te - 2.0 * sp[1];
			vc0 = vc0 * te - sp[2];

			sp = &childSec.splineParamsY[0];
			number& vc1 = childDir[1];
			vc1 = -3.0 * sp[0] * te - 2.0 * sp[1];
			vc1 = vc1 * te - sp[2];

			sp = &childSec.splineParamsZ[0];
			number& vc2 = childDir[2];
			vc2 = -3.0 * sp[0] * te - 2.0 * sp[1];
			vc2 = vc2 * te - sp[2];

			// find out neurite direction in next BP
			number bpAxPos = vSegAxPos[nSeg - 1];
			size_t tmpSec = curSec;
			for (; tmpSec < nSec; ++tmpSec) {
				const NeuriteProjector::Section& sec = neurite.vSec[tmpSec];
				if (sec.endParam >= bpAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[tmpSec];
			number monom = sec.endParam - bpAxPos;
			sp = &sec.splineParamsX[0];
			number& v0 = vel[0];
			v0 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v0 = v0 * monom - sp[2];

			sp = &sec.splineParamsY[0];
			number& v1 = vel[1];
			v1 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v1 = v1 * monom - sp[2];

			sp = &sec.splineParamsZ[0];
			number& v2 = vel[2];
			v2 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v2 = v2 * monom - sp[2];

			VecNormalize(vel, vel);

			// calculate offset
			number fac = VecProd(neurite.refDir, vel);
			VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
			VecNormalize(projRefDir, projRefDir);
			VecCross(thirdDir, vel, projRefDir);

			vector2 relCoord;
			VecScaleAppend(childDir, -VecProd(childDir, vel), vel);
			relCoord[0] = VecProd(childDir, projRefDir);
			VecScaleAppend(childDir, -relCoord[0], projRefDir);
			relCoord[1] = VecProd(childDir, thirdDir);
			VecNormalize(relCoord, relCoord);

			number branchOffset = 0.0;
			if (fabs(relCoord[0]) < 1e-8)
				branchOffset = relCoord[1] < 0 ? 1.5 * PI : 0.5 * PI;
			else
				branchOffset =
						relCoord[0] < 0 ?
								PI - atan(-relCoord[1] / relCoord[0]) :
								atan(relCoord[1] / relCoord[0]);

			addOffset = branchOffset - angleOffset;
			connFaceInd = floor(
					std::fmod(addOffset + 4 * PI, 2 * PI) / (PI / 2));
			addOffset = std::fmod(
					addOffset - (connFaceInd * PI / 2 + PI / 4) + 4 * PI,
					2 * PI);
			if (addOffset > PI)
				addOffset -= 2 * PI;
			addOffset /= nSeg - 1;
		}

		// create mesh for segments
		Selector sel(g);
		for (size_t s = 0; s < nSeg; ++s) {
			// get exact position, velocity and radius of segment end
			number segAxPos = vSegAxPos[s];
			for (; curSec < nSec; ++curSec) {
				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				if (sec.endParam >= segAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			vector3 curPos;
			number monom = sec.endParam - segAxPos;
			const number* sp = &sec.splineParamsX[0];
			number& p0 = curPos[0];
			number& v0 = vel[0];
			p0 = sp[0] * monom + sp[1];
			p0 = p0 * monom + sp[2];
			p0 = p0 * monom + sp[3];
			v0 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v0 = v0 * monom - sp[2];

			sp = &sec.splineParamsY[0];
			number& p1 = curPos[1];
			number& v1 = vel[1];
			p1 = sp[0] * monom + sp[1];
			p1 = p1 * monom + sp[2];
			p1 = p1 * monom + sp[3];
			v1 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v1 = v1 * monom - sp[2];

			sp = &sec.splineParamsZ[0];
			number& p2 = curPos[2];
			number& v2 = vel[2];
			p2 = sp[0] * monom + sp[1];
			p2 = p2 * monom + sp[2];
			p2 = p2 * monom + sp[3];
			v2 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v2 = v2 * monom - sp[2];

			sp = &sec.splineParamsR[0];
			number radius;
			radius = sp[0] * monom + sp[1];
			radius = radius * monom + sp[2];
			radius = radius * monom + sp[3];

			VecNormalize(vel, vel);

			// calculate reference dir projected to normal plane of velocity
			number fac = VecProd(neurite.refDir, vel);
			VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
			VecNormalize(projRefDir, projRefDir);

			VecCross(thirdDir, vel, projRefDir);

			// usual segment: extrude
			if (s != nSeg - 1 || brit == brit_end) {
				// apply additional offset
				angleOffset = std::fmod(angleOffset + addOffset + 2 * PI,
						2 * PI);

				// extrude from last pos to new pos
				vector3 extrudeDir;
				VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, lastPos);
				std::vector<Volume*> vVol;
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos,
						EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j) {
					number angle = 0.5 * PI * j + angleOffset;
					if (angle > 2 * PI)
						angle -= 2 * PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, erScaleFactor * radius * cos(angle),
							projRefDir, erScaleFactor * radius * sin(angle),
							thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = erScaleFactor;
				}
				for (size_t j = 0; j < 12; ++j) {
					number angle = PI * ((number) j / 6) + angleOffset;
					if (angle > 2 * PI)
						angle -= 2 * PI;
					Vertex* v = vVrt[j + 4];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius * cos(angle), projRefDir,
							radius * sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// ensure correct volume orientation
				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);
			}

			// BP segment: create BP with tetrahedra/pyramids and create whole branch
			else {
				std::vector<Volume*> vBPVols;
				vBPVols.reserve(27);

				// correct vertex offsets to reflect angle at which child branches
				VecScaleAppend(aaPos[vVrt[(connFaceInd) % 4]],
						erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 1) % 4]],
						erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 2) % 4]],
						-erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 3) % 4]],
						-erScaleFactor * surfBPoffset, vel);
				aaSurfParams[vVrt[(connFaceInd) % 4]].axial += erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 1) % 4]].axial += erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 2) % 4]].axial -= erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 3) % 4]].axial -= erScaleFactor
						* surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd) % 4)]],
						surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]],
						surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 2) % 4)]],
						-surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 3) % 4)]],
						-surfBPoffset, vel);
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd) % 4)]].axial +=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]].axial +=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 3) % 4)]].axial -=
						surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd) % 4)]],
						1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd) % 4)]],
						1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 2) % 4)]],
						-1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 2) % 4)]],
						-1.366 * surfBPoffset, vel);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd) % 4)]].axial += 1.366
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd) % 4)]].axial += 1.366
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						1.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						1.366 * surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 1) % 4)]],
						0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 1) % 4)]],
						-0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 3) % 4)]],
						-0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 3) % 4)]],
						0.366 * surfBPoffset, vel);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 1) % 4)]].axial +=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 1) % 4)]].axial -=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 3) % 4)]].axial -=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 3) % 4)]].axial +=
						0.366 * surfBPoffset / neurite_length;

				// add branch neurite ID to its initial vertices
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);

				// prepare branch vertices
				std::vector<Vertex*> vBranchVrts(16);
				vBranchVrts[4] = vVrt[4 + 3 * ((connFaceInd + 1) % 4)];
				vBranchVrts[13] = vVrt[4 + 3 * ((connFaceInd) % 4)];
				vBranchVrts[14] = vVrt[5 + 3 * ((connFaceInd) % 4)];
				vBranchVrts[15] = vVrt[6 + 3 * ((connFaceInd) % 4)];

				// extrude to first third of BP //
				segAxPos = 0.5 * (1.0 + erScaleFactor) * vSegAxPos[s - 1]
				                                                   + 0.5 * (1.0 - erScaleFactor) * vSegAxPos[s];
				vector3 firstPos;
				VecScaleAdd(firstPos, 0.5 * (1.0 + erScaleFactor), lastPos,
						0.5 * (1.0 - erScaleFactor), curPos);
				vector3 extrudeDir;
				VecScaleAdd(extrudeDir, 1.0, firstPos, -1.0, lastPos);
				std::vector<Volume*> vVol;
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos,
						EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);
				for (size_t j = 0; j < 9; ++j)
					vBPVols.push_back(vVol[j]);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j) {
					number angle = 0.5 * PI * j + angleOffset;
					if (angle > 2 * PI)
						angle -= 2 * PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, erScaleFactor * radius * cos(angle),
							projRefDir, erScaleFactor * radius * sin(angle),
							thirdDir);
					VecAdd(aaPos[v], firstPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = erScaleFactor;
				}
				for (size_t j = 0; j < 12; ++j) {
					number angle = PI * ((number) j / 6) + angleOffset;
					if (angle > 2 * PI)
						angle -= 2 * PI;
					Vertex* v = vVrt[j + 4];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius * cos(angle), projRefDir,
							radius * sin(angle), thirdDir);
					VecAdd(aaPos[v], firstPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// correct vertex offsets to reflect angle at which child branches
				VecScaleAppend(aaPos[vVrt[(connFaceInd) % 4]],
						erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 1) % 4]],
						erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 2) % 4]],
						-erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 3) % 4]],
						-erScaleFactor * surfBPoffset, vel);
				aaSurfParams[vVrt[(connFaceInd) % 4]].axial += erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 1) % 4]].axial += erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 2) % 4]].axial -= erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 3) % 4]].axial -= erScaleFactor
						* surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd) % 4)]],
						surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]],
						surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 2) % 4)]],
						-surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 3) % 4)]],
						-surfBPoffset, vel);
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd) % 4)]].axial +=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]].axial +=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 3) % 4)]].axial -=
						surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd) % 4)]],
						1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd) % 4)]],
						1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 2) % 4)]],
						-1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 2) % 4)]],
						-1.366 * surfBPoffset, vel);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd) % 4)]].axial += 1.366
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd) % 4)]].axial += 1.366
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						1.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						1.366 * surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 1) % 4)]],
						0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 1) % 4)]],
						-0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 3) % 4)]],
						-0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 3) % 4)]],
						0.366 * surfBPoffset, vel);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 1) % 4)]].axial +=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 1) % 4)]].axial -=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 3) % 4)]].axial -=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 3) % 4)]].axial +=
						0.366 * surfBPoffset / neurite_length;

				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);

				// add branch neurite ID to its initial vertices
				aaSurfParams[vVrt[(connFaceInd) % 4]].neuriteID += ((brit
						- vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[(connFaceInd + 1) % 4]].neuriteID += ((brit
						- vBR.begin()) << 20) + (1 << 28);

				aaSurfParams[vVrt[4 + 3 * ((connFaceInd) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd) % 4)]].neuriteID =
						child_nid;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd) % 4)]].neuriteID =
						child_nid;

				// prepare branch vertices
				vBranchVrts[0] = vVrt[6 + 3 * ((connFaceInd) % 4)];
				vBranchVrts[3] = vVrt[5 + 3 * ((connFaceInd) % 4)];
				vBranchVrts[5] = vVrt[4 + 3 * ((connFaceInd + 1) % 4)];
				vBranchVrts[12] = vVrt[4 + 3 * ((connFaceInd) % 4)];

				// extrude to second third of BP //
				segAxPos = 0.5 * (1.0 - erScaleFactor) * vSegAxPos[s - 1]
				                                                   + 0.5 * (1.0 + erScaleFactor) * vSegAxPos[s];
				vector3 secondPos;
				VecScaleAdd(secondPos, 0.5 * (1.0 - erScaleFactor), lastPos,
						0.5 * (1.0 + erScaleFactor), curPos);
				VecScaleAdd(extrudeDir, 1.0, secondPos, -1.0, firstPos);
				vVol.clear();
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos,
						EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);
				for (size_t j = 0; j < 9; ++j)
					vBPVols.push_back(vVol[j]);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j) {
					number angle = 0.5 * PI * j + angleOffset;
					if (angle > 2 * PI)
						angle -= 2 * PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, erScaleFactor * radius * cos(angle),
							projRefDir, erScaleFactor * radius * sin(angle),
							thirdDir);
					VecAdd(aaPos[v], secondPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = erScaleFactor;
				}
				for (size_t j = 0; j < 12; ++j) {
					number angle = PI * ((number) j / 6) + angleOffset;
					if (angle > 2 * PI)
						angle -= 2 * PI;
					Vertex* v = vVrt[j + 4];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius * cos(angle), projRefDir,
							radius * sin(angle), thirdDir);
					VecAdd(aaPos[v], secondPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// correct vertex offsets to reflect angle at which child branches
				VecScaleAppend(aaPos[vVrt[(connFaceInd) % 4]],
						erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 1) % 4]],
						erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 2) % 4]],
						-erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 3) % 4]],
						-erScaleFactor * surfBPoffset, vel);
				aaSurfParams[vVrt[(connFaceInd) % 4]].axial += erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 1) % 4]].axial += erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 2) % 4]].axial -= erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 3) % 4]].axial -= erScaleFactor
						* surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd) % 4)]],
						surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]],
						surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 2) % 4)]],
						-surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 3) % 4)]],
						-surfBPoffset, vel);
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd) % 4)]].axial +=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]].axial +=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 3) % 4)]].axial -=
						surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd) % 4)]],
						1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd) % 4)]],
						1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 2) % 4)]],
						-1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 2) % 4)]],
						-1.366 * surfBPoffset, vel);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd) % 4)]].axial += 1.366
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd) % 4)]].axial += 1.366
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						1.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						1.366 * surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 1) % 4)]],
						0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 1) % 4)]],
						-0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 3) % 4)]],
						-0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 3) % 4)]],
						0.366 * surfBPoffset, vel);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 1) % 4)]].axial +=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 1) % 4)]].axial -=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 3) % 4)]].axial -=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 3) % 4)]].axial +=
						0.366 * surfBPoffset / neurite_length;

				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);

				// add branch neurite ID to its initial vertices
				aaSurfParams[vVrt[(connFaceInd) % 4]].neuriteID += ((brit
						- vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[(connFaceInd + 1) % 4]].neuriteID += ((brit
						- vBR.begin()) << 20) + (1 << 28);

				aaSurfParams[vVrt[4 + 3 * ((connFaceInd) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd) % 4)]].neuriteID =
						child_nid;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd) % 4)]].neuriteID =
						child_nid;

				// prepare branch vertices
				vBranchVrts[1] = vVrt[6 + 3 * ((connFaceInd) % 4)];
				vBranchVrts[2] = vVrt[5 + 3 * ((connFaceInd) % 4)];
				vBranchVrts[6] = vVrt[4 + 3 * ((connFaceInd + 1) % 4)];
				vBranchVrts[11] = vVrt[4 + 3 * ((connFaceInd) % 4)];

				// extrude to end of BP //
				segAxPos = vSegAxPos[s];
				VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, secondPos);
				vVol.clear();
				Extrude(g, &vVrt, &vEdge, &vFace, extrudeDir, aaPos,
						EO_CREATE_FACES | EO_CREATE_VOLUMES, &vVol);
				for (size_t j = 0; j < 9; ++j)
					vBPVols.push_back(vVol[j]);

				// set new positions and param attachments
				for (size_t j = 0; j < 4; ++j) {
					number angle = 0.5 * PI * j + angleOffset;
					if (angle > 2 * PI)
						angle -= 2 * PI;
					Vertex* v = vVrt[j];
					vector3 radialVec;
					VecScaleAdd(radialVec, erScaleFactor * radius * cos(angle),
							projRefDir, erScaleFactor * radius * sin(angle),
							thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = erScaleFactor;
				}
				for (size_t j = 0; j < 12; ++j) {
					number angle = PI * ((number) j / 6) + angleOffset;
					if (angle > 2 * PI)
						angle -= 2 * PI;
					Vertex* v = vVrt[j + 4];
					vector3 radialVec;
					VecScaleAdd(radialVec, radius * cos(angle), projRefDir,
							radius * sin(angle), thirdDir);
					VecAdd(aaPos[v], curPos, radialVec);

					aaSurfParams[v].neuriteID = nid;
					aaSurfParams[v].axial = segAxPos;
					aaSurfParams[v].angular = angle;
					aaSurfParams[v].radial = 1.0;
				}

				// correct vertex offsets to reflect angle at which child branches
				VecScaleAppend(aaPos[vVrt[(connFaceInd) % 4]],
						erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 1) % 4]],
						erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 2) % 4]],
						-erScaleFactor * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[(connFaceInd + 3) % 4]],
						-erScaleFactor * surfBPoffset, vel);
				aaSurfParams[vVrt[(connFaceInd) % 4]].axial += erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 1) % 4]].axial += erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 2) % 4]].axial -= erScaleFactor
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[(connFaceInd + 3) % 4]].axial -= erScaleFactor
						* surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd) % 4)]],
						surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]],
						surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 2) % 4)]],
						-surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[4 + 3 * ((connFaceInd + 3) % 4)]],
						-surfBPoffset, vel);
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd) % 4)]].axial +=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]].axial +=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						surfBPoffset / neurite_length;
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 3) % 4)]].axial -=
						surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd) % 4)]],
						1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd) % 4)]],
						1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 2) % 4)]],
						-1.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 2) % 4)]],
						-1.366 * surfBPoffset, vel);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd) % 4)]].axial += 1.366
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd) % 4)]].axial += 1.366
						* surfBPoffset / neurite_length;
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						1.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 2) % 4)]].axial -=
						1.366 * surfBPoffset / neurite_length;

				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 1) % 4)]],
						0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 1) % 4)]],
						-0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[5 + 3 * ((connFaceInd + 3) % 4)]],
						-0.366 * surfBPoffset, vel);
				VecScaleAppend(aaPos[vVrt[6 + 3 * ((connFaceInd + 3) % 4)]],
						0.366 * surfBPoffset, vel);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 1) % 4)]].axial +=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 1) % 4)]].axial -=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd + 3) % 4)]].axial -=
						0.366 * surfBPoffset / neurite_length;
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd + 3) % 4)]].axial +=
						0.366 * surfBPoffset / neurite_length;

				FixOrientation(g, vVol.begin(), vVol.end(), aaPos);

				// add branch neurite ID to its initial vertices
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[4 + 3 * ((connFaceInd + 1) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[5 + 3 * ((connFaceInd) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);
				aaSurfParams[vVrt[6 + 3 * ((connFaceInd) % 4)]].neuriteID +=
						((brit - vBR.begin()) << 20) + (1 << 28);

				// prepare branch vertices
				vBranchVrts[7] = vVrt[4 + 3 * ((connFaceInd + 1) % 4)];
				vBranchVrts[8] = vVrt[6 + 3 * ((connFaceInd) % 4)];
				vBranchVrts[9] = vVrt[5 + 3 * ((connFaceInd) % 4)];
				vBranchVrts[10] = vVrt[4 + 3 * ((connFaceInd) % 4)];

				// prepare connectingEdges, -Faces for branch
				std::vector<EdgeDescriptor> vED(24);
				for (size_t i = 0; i < 4; ++i) {
					vED[i] = EdgeDescriptor(vBranchVrts[i],
							vBranchVrts[(i + 1) % 4]);
					vED[i + 4] = EdgeDescriptor(vBranchVrts[i],
							vBranchVrts[5 + 3 * i]);
					vED[i + 8] = EdgeDescriptor(vBranchVrts[(i + 1) % 4],
							vBranchVrts[6 + 3 * i]);
				}
				for (size_t i = 0; i < 12; ++i)
					vED[i + 12] = EdgeDescriptor(vBranchVrts[i + 4],
							vBranchVrts[(i + 1) % 12 + 4]);

				std::vector<FaceDescriptor> vFD(9);
				vFD[0] = FaceDescriptor(vBranchVrts[0], vBranchVrts[1],
						vBranchVrts[2], vBranchVrts[3]);
				for (size_t i = 0; i < 4; ++i) {
					vFD[i + 1] = FaceDescriptor(vBranchVrts[i],
							vBranchVrts[(3 * i + 11) % 12 + 4],
							vBranchVrts[3 * i + 4], vBranchVrts[3 * i + 5]);
					vFD[i + 5] = FaceDescriptor(vBranchVrts[i],
							vBranchVrts[3 * i + 5], vBranchVrts[3 * i + 6],
							vBranchVrts[(i + 1) % 4]);
				}

				typedef Grid::traits<Face>::secure_container faceCont;
				std::vector<Face*> vBranchFaces(9);
				for (size_t j = 0; j < 9; ++j) {
					const FaceDescriptor& qDesc = vFD[j];

					faceCont fl;
					bool found = false;
					for (size_t k = 0; k < 27; ++k) {
						g.associated_elements(fl, vBPVols[k]);
						const size_t flSz = fl.size();
						for (size_t f = 0; f < flSz; ++f) {
							if (CompareVertices(fl[f], &qDesc)) {
								vBranchFaces[j] = fl[f];
								found = true;
								break;
							}
						}
						if (found)
							break;
					}
					UG_COND_THROW(!found,
							"Connecting face " << j << " not found.")
				}

				typedef Grid::traits<Edge>::secure_container edgeCont;
				std::vector<Edge*> vBranchEdges(24);
				for (size_t j = 0; j < 24; ++j) {
					const EdgeDescriptor& eDesc = vED[j];

					edgeCont el;
					bool found = false;
					for (size_t k = 0; k < 9; ++k) {
						g.associated_elements(el, vBranchFaces[k]);
						const size_t elSz = el.size();
						for (size_t e = 0; e < elSz; ++e) {
							if (CompareVertices(el[e], &eDesc)) {
								vBranchEdges[j] = el[e];
								found = true;
								break;
							}
						}
						if (found)
							break;
					}
					UG_COND_THROW(!found,
							"Connecting edge " << j << " not found.")
				}

				// correct connecting volume subset indices
				{
					Volume* connVol = vBPVols[connFaceInd + 14];
					sh.assign_subset(connVol, 1);

					typedef Grid::traits<Face>::secure_container faceCont;
					faceCont fl;
					g.associated_elements(fl, connVol);
					const size_t flSz = fl.size();
					for (size_t f = 0; f < flSz; ++f) {
						Face* sideFace = fl[f];
						Volume* opp = GetConnectedNeighbor(g, sideFace,
								connVol);
						if (!opp || sh.get_subset_index(opp) == 1)
							sh.assign_subset(sideFace, 1);
						else {
							sh.assign_subset(sideFace, 3);

							typedef Grid::traits<Edge>::secure_container edgeCont;
							edgeCont el;
							g.associated_elements(el, sideFace);
							const size_t elSz = el.size();
							for (size_t e = 0; e < elSz; ++e) {
								Edge* sideEdge = el[e];
								sh.assign_subset(sideEdge, 3);
								sh.assign_subset(sideEdge->vertex(0), 3);
								sh.assign_subset(sideEdge->vertex(1), 3);
							}
						}
					}
				}
				for (size_t j = 1; j < 9; ++j) {
					Face* connFace = vBranchFaces[j];
					sh.assign_subset(connFace, 0);
				}
				for (size_t j = 4; j < 12; ++j) {
					Edge* connEdge = vBranchEdges[j];
					sh.assign_subset(connEdge, 0);
				}

				// recursively build branch
				create_neurite_with_er(vNeurites, vPos, vR, child_nid,
						erScaleFactor, anisotropy, g, aaPos, aaSurfParams, sh,
						&vBranchVrts, &vBranchEdges, &vBranchFaces,
						branchOffset[1], NULL, NULL, NULL, NULL);
			}

			lastPos = curPos;
		}

		// update t_end and curSec
		if (brit != brit_end)
			t_end = bp_end;

		for (; curSec < nSec; ++curSec) {
			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			if (sec.endParam >= t_end)
				break;
		}

		// check whether tip has been reached
		if (brit == brit_end)
			break;
		else
			++brit;
		}
	}

	////////////////////////////////////////////////////////////////////////
	/// create_neurite_surf
	////////////////////////////////////////////////////////////////////////
	void create_neurite_surf(
		const std::vector<NeuriteProjector::Neurite>& vNeurites,
		const std::vector<std::vector<vector3> >& vPos,
		const std::vector<std::vector<number> >& vR, size_t nid,
		number anisotropy, Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		Grid::VertexAttachmentAccessor<
		Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		std::vector<Vertex*>* connectingVrts,
		std::vector<Edge*>* connectingEdges,
		number initialOffset) {
	const NeuriteProjector::Neurite& neurite = vNeurites[nid];
	const std::vector<vector3>& pos = vPos[nid];
	const std::vector<number>& r = vR[nid];

	number neurite_length = 0.0;
	for (size_t i = 1; i < pos.size(); ++i)
		neurite_length += VecDistance(pos[i], pos[i - 1]);

	size_t nSec = neurite.vSec.size();

	const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit =
			vBR.begin();
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end =
			vBR.end();

	std::vector<Vertex*> vVrt;
	std::vector<Edge*> vEdge;
	vVrt.resize(4);
	vEdge.resize(4);

	vector3 vel;
	const NeuriteProjector::Section& sec = neurite.vSec[0];
	number h = sec.endParam;
	vel[0] = -3.0 * sec.splineParamsX[0] * h * h
			- 2.0 * sec.splineParamsX[1] * h - sec.splineParamsX[2];
	vel[1] = -3.0 * sec.splineParamsY[0] * h * h
			- 2.0 * sec.splineParamsY[1] * h - sec.splineParamsY[2];
	vel[2] = -3.0 * sec.splineParamsZ[0] * h * h
			- 2.0 * sec.splineParamsZ[1] * h - sec.splineParamsZ[2];

	vector3 projRefDir;
	VecNormalize(vel, vel);
	number fac = VecProd(neurite.refDir, vel);
	VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
	VecNormalize(projRefDir, projRefDir);

	vector3 thirdDir;
	VecCross(thirdDir, vel, projRefDir);

	number angleOffset = 0.0;

	number t_start = 0.0;
	number t_end = 0.0;

	if (connectingVrts && connectingEdges) {
		vVrt = *connectingVrts;
		vEdge = *connectingEdges;

		// calculate start angle offset
		vector3 center(0.0);
		for (size_t i = 0; i < 4; ++i)
			VecAdd(center, center, aaPos[(*connectingVrts)[i]]);
		center /= 4;

		vector3 centerToFirst;
		VecSubtract(centerToFirst, aaPos[(*connectingVrts)[0]], center);

		vector2 relCoord;
		VecScaleAdd(centerToFirst, 1.0, centerToFirst,
				-VecProd(centerToFirst, vel), vel);
		relCoord[0] = VecProd(centerToFirst, projRefDir);
		VecScaleAdd(centerToFirst, 1.0, centerToFirst, -relCoord[0],
				projRefDir);
		relCoord[1] = VecProd(centerToFirst, thirdDir);
		VecNormalize(relCoord, relCoord);

		if (fabs(relCoord[0]) < 1e-8)
			angleOffset = relCoord[1] < 0 ? 1.5 * PI : 0.5 * PI;
		else
			angleOffset =
					relCoord[0] < 0 ?
							PI - atan(-relCoord[1] / relCoord[0]) :
							atan(relCoord[1] / relCoord[0]);
		if (angleOffset < 0)
			angleOffset += 2.0 * PI;

		// ignore first branching region (the connecting region)
		++brit;

		// apply initial offset (to prevent first segment being shorter than the others)
		t_end = initialOffset / neurite_length;
	} else {
		// create first layer of vertices/edges
		for (size_t i = 0; i < 4; ++i) {
			Vertex* v = *g.create<RegularVertex>();
			vVrt[i] = v;
			number angle = 0.5 * PI * i;
			VecScaleAdd(aaPos[v], 1.0, pos[0], r[0] * cos(angle), projRefDir,
					r[0] * sin(angle), thirdDir);

			aaSurfParams[v].neuriteID = nid;
			aaSurfParams[v].axial = 0.0;
			aaSurfParams[v].angular = angle;
			aaSurfParams[v].radial = 1.0;
		}
		for (size_t i = 0; i < 4; ++i)
			vEdge[i] = *g.create<RegularEdge>(
					EdgeDescriptor(vVrt[i], vVrt[(i + 1) % 4]));
	}

	// Now create dendrite to the next branching point and iterate this process.
	// We want to create each of the segments with approx. the same aspect ratio.
	// To that end, we first calculate the length of the section to be created (in units of radius)
	// and then divide this number by 2^n where n is the number of anisotropic refinements to
	// be performed to make all segments (more or less) isotropic. The result is the number
	// of segments to be used for the section.
	vector3 lastPos = pos[0];
	size_t curSec = 0;

	while (true) {
		t_start = t_end;

		// if there is another BP, store its extensions here
		number bp_start = 1.0;
		number bp_end = 0.0;

		// initial branch offsets (to prevent the first segment being shorter than the following)
		std::vector<number> branchOffset;

		// last section: create until tip
		if (brit == brit_end)
			t_end = 1.0;

		// otherwise: section goes to next branching point
		else {
			// calculate the exact position of the branching point,
			// i.e., the axial position of the intersection of the branching neurite's
			// spline with the surface of the current neurite
			// this is necessary esp. when the branching angle is very small
			const std::vector<uint32_t>& vBranchInd = brit->bp->vNid;
			size_t nBranches = vBranchInd.size();

			branchOffset.resize(brit->bp->vNid.size(), 0.0);
			for (size_t br = 1; br < nBranches; ++br) {
				uint32_t brInd = vBranchInd[br];

				// get position and radius of first point of branch
				const number brRadSeg1 = vR[brInd][0];

				// get position and radius of branching point
				const number bpTPos = brit->t;
				size_t brSec = curSec;
				for (; brSec < nSec; ++brSec) {
					const NeuriteProjector::Section& sec = neurite.vSec[brSec];
					if (bpTPos - sec.endParam < 1e-6 * bpTPos)
						break;
				}
				UG_COND_THROW(brSec == nSec,
						"Could not find section containing branching point "
						"at t = " << bpTPos << ".");
				const number bpRad = vR[nid][brSec + 1];

				// calculate branch and neurite directions
				vector3 branchDir;
				const NeuriteProjector::Section& childSec =
						vNeurites[brInd].vSec[0];
				number te = childSec.endParam;

				const number* s = &childSec.splineParamsX[0];
				number& v0 = branchDir[0];
				v0 = -3.0 * s[0] * te - 2.0 * s[1];
				v0 = v0 * te - s[2];

				s = &childSec.splineParamsY[0];
				number& v1 = branchDir[1];
				v1 = -3.0 * s[0] * te - 2.0 * s[1];
				v1 = v1 * te - s[2];

				s = &childSec.splineParamsZ[0];
				number& v2 = branchDir[2];
				v2 = -3.0 * s[0] * te - 2.0 * s[1];
				v2 = v2 * te - s[2];

				VecNormalize(branchDir, branchDir);

				vector3 neuriteDir;
				const NeuriteProjector::Section& sec = neurite.vSec.at(brSec);
				vel[0] = -sec.splineParamsX[2];
				vel[1] = -sec.splineParamsY[2];
				vel[2] = -sec.splineParamsZ[2];
				number velNorm = sqrt(VecNormSquared(vel));
				VecScale(neuriteDir, vel, 1.0 / velNorm);

				// calculate offset of true branch position, which is sqrt(2)/2*r1*cot(alpha)
				number brScProd = VecProd(neuriteDir, branchDir);
				number surfBPoffset = 0.5 * sqrt(2.0) * bpRad * brScProd
						/ sqrt(1.0 - brScProd * brScProd);

				// calculate offset of new branch, which is sqrt(2)/2*r1/sin(alpha)
				branchOffset[br] = 0.5 * sqrt(2.0) * bpRad
						/ sqrt(1.0 - brScProd * brScProd);

				// calculate true length of BP, which is sqrt(2)*r2/sin(alpha)
				number surfBPhalfLength = 0.5 * sqrt(2.0) * brRadSeg1
						/ sqrt(1.0 - brScProd * brScProd);

				// finally set bp start and end
				bp_start = std::min(bp_start,
						bpTPos
						+ (surfBPoffset - surfBPhalfLength)
						/ neurite_length);
				bp_end = std::max(bp_end,
						bpTPos
						+ (surfBPoffset + surfBPhalfLength)
						/ neurite_length);
			}

			t_end = bp_start;
		}

		// calculate total length in units of radius
		// = integral from t_start to t_end over: ||v(t)|| / r(t) dt
		number lengthOverRadius = calculate_length_over_radius(t_start, t_end,
				neurite, curSec);

		// to reach the desired anisotropy on the surface in the refinement limit,
		// it has to be multiplied by pi/2 h
		size_t nSeg = (size_t) floor(
				lengthOverRadius / (anisotropy * 0.5 * PI));
		if (!nSeg)
			nSeg = 1;
		number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
		std::vector<number> vSegAxPos(nSeg);
		calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite,
				curSec, segLength);

		// add the branching point to segment list (if present)
		if (brit != brit_end) {
			vSegAxPos.resize(nSeg + 1);
			vSegAxPos[nSeg] = bp_end;
			++nSeg;
		}

		// create mesh for segments
		Selector sel(g);
		for (size_t s = 0; s < nSeg; ++s) {
			// get exact position, velocity and radius of segment end
			number segAxPos = vSegAxPos[s];
			for (; curSec < nSec; ++curSec) {
				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				if (sec.endParam >= segAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			vector3 curPos;
			number monom = sec.endParam - segAxPos;
			const number* sp = &sec.splineParamsX[0];
			number& p0 = curPos[0];
			number& v0 = vel[0];
			p0 = sp[0] * monom + sp[1];
			p0 = p0 * monom + sp[2];
			p0 = p0 * monom + sp[3];
			v0 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v0 = v0 * monom - sp[2];

			sp = &sec.splineParamsY[0];
			number& p1 = curPos[1];
			number& v1 = vel[1];
			p1 = sp[0] * monom + sp[1];
			p1 = p1 * monom + sp[2];
			p1 = p1 * monom + sp[3];
			v1 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v1 = v1 * monom - sp[2];

			sp = &sec.splineParamsZ[0];
			number& p2 = curPos[2];
			number& v2 = vel[2];
			p2 = sp[0] * monom + sp[1];
			p2 = p2 * monom + sp[2];
			p2 = p2 * monom + sp[3];
			v2 = -3.0 * sp[0] * monom - 2.0 * sp[1];
			v2 = v2 * monom - sp[2];

			sp = &sec.splineParamsR[0];
			number radius;
			radius = sp[0] * monom + sp[1];
			radius = radius * monom + sp[2];
			radius = radius * monom + sp[3];

			VecNormalize(vel, vel);

			// calculate reference dir projected to normal plane of velocity
			number fac = VecProd(neurite.refDir, vel);
			VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
			VecNormalize(projRefDir, projRefDir);

			VecCross(thirdDir, vel, projRefDir);

			// extrude from last pos to new pos
			if (s == nSeg - 1 && brit != brit_end)
				sel.enable_autoselection(true); // for last segment (BP), select new elems

			vector3 extrudeDir;
			VecScaleAdd(extrudeDir, 1.0, curPos, -1.0, lastPos);
			Extrude(g, &vVrt, &vEdge, NULL, extrudeDir, aaPos, EO_CREATE_FACES,
					NULL);

			sel.enable_autoselection(false);

			// set new positions and param attachments; also ensure correct face orientation
			for (size_t j = 0; j < 4; ++j) {
				number angle = 0.5 * PI * j + angleOffset;
				if (angle > 2 * PI)
					angle -= 2 * PI;
				Vertex* v = vVrt[j];
				vector3 radialVec;
				VecScaleAdd(radialVec, radius * cos(angle), projRefDir,
						radius * sin(angle), thirdDir);
				VecAdd(aaPos[v], curPos, radialVec);

				aaSurfParams[v].neuriteID = nid;
				aaSurfParams[v].axial = segAxPos;
				aaSurfParams[v].angular = angle;
				aaSurfParams[v].radial = 1.0;

				Grid::traits<Face>::secure_container faceCont;
				g.associated_elements(faceCont, vEdge[j]); // faceCont must contain exactly one face
				vector3 normal;
				CalculateNormal(normal, faceCont[0], aaPos);
				if (VecProd(normal, radialVec) < 0)
					g.flip_orientation(faceCont[0]);
			}

			lastPos = curPos;
		}

		// connect branching neurites if present
		if (brit != brit_end) {
			// find branching child neurite
			SmartPtr<NeuriteProjector::BranchingPoint> bp = brit->bp;
			UG_COND_THROW(bp->vNid.size() > 2,
					"This implementation can only handle branching points with one branching child.");

			size_t child_nid;
			if (bp->vNid[0] != nid)
				child_nid = bp->vNid[0];
			else
				child_nid = bp->vNid[1];

			// find out branching child neurite initial direction
			vector3 childDir;
			const NeuriteProjector::Section& childSec =
					vNeurites[child_nid].vSec[0];
			number te = childSec.endParam;

			const number* s = &childSec.splineParamsX[0];
			number& v0 = childDir[0];
			v0 = -3.0 * s[0] * te - 2.0 * s[1];
			v0 = v0 * te - s[2];

			s = &childSec.splineParamsY[0];
			number& v1 = childDir[1];
			v1 = -3.0 * s[0] * te - 2.0 * s[1];
			v1 = v1 * te - s[2];

			s = &childSec.splineParamsZ[0];
			number& v2 = childDir[2];
			v2 = -3.0 * s[0] * te - 2.0 * s[1];
			v2 = v2 * te - s[2];

			// now choose best side of hexahedron to connect to
			vector3 normal;
			Face* best = NULL;
			number bestProd = 0.0;
			Selector::traits<Face>::iterator fit = sel.faces_begin();
			Selector::traits<Face>::iterator fit_end = sel.faces_end();
			for (; fit != fit_end; ++fit) {
				CalculateNormal(normal, *fit, aaPos);
				number prod = VecProd(normal, childDir);
				if (prod > bestProd) {
					best = *fit;
					bestProd = prod;
				}
			}
			UG_COND_THROW(!best,
					"None of the branching point faces pointed in a suitable direction.")
			sel.deselect(sel.faces_begin(), sel.faces_end());

			// add branch neurite ID to its initial vertices
			for (size_t j = 0; j < 4; ++j) {
				Vertex* brVrt = best->vertex(j);
				aaSurfParams[brVrt].neuriteID += (brit - vBR.begin()) << 20; // add branching region index
				aaSurfParams[brVrt].neuriteID += 1 << 28; // add child ID (always 0, since there can only be one child here)
			}

			// remove face and call creation of child neurite (recursion)
			std::vector<Vertex*> vrts(4);
			UG_COND_THROW(best->num_vertices() != 4,
					"Hexaeder face does not have 4 vertices!");
			for (size_t j = 0; j < 4; ++j)
				vrts[j] = best->vertex(j);
			std::vector<Edge*> edges(4);

			Grid::traits<Edge>::secure_container edgeCont;
			g.associated_elements(edgeCont, best);
			size_t esz = edgeCont.size();

			for (size_t j = 0; j < 4; ++j) {
				Vertex* first = vrts[j];
				Vertex* second = vrts[(j + 1) % 4];

				size_t k = 0;
				for (; k < esz; ++k) {
					if ((edgeCont[k]->vertex(0) == first
							&& edgeCont[k]->vertex(1) == second)
							|| (edgeCont[k]->vertex(0) == second
									&& edgeCont[k]->vertex(1) == first)) {
						edges[j] = edgeCont[k];
						break;
					}
				}
				UG_COND_THROW(k == esz,
						"Connecting edges for child neurite could not be determined.");
			}

			g.erase(best);

			// TODO: create prism to connect to in case the branching angle is small or big
			create_neurite_surf(vNeurites, vPos, vR, child_nid, anisotropy, g,
					aaPos, aaSurfParams, &vrts, &edges, branchOffset[1]);
		}

		// update t_end and curSec
		if (brit != brit_end)
			t_end = bp_end;

		for (; curSec < nSec; ++curSec) {
			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			if (sec.endParam >= t_end)
				break;
		}

		// check whether tip has been reached
		if (brit == brit_end)
			break;
		else
			++brit;
	}

	// close the tip of the neurite
	const NeuriteProjector::Section& lastSec = neurite.vSec[nSec - 1];
	vel = vector3(-lastSec.splineParamsX[2], -lastSec.splineParamsY[2],
			-lastSec.splineParamsZ[2]);
	number radius = lastSec.splineParamsR[3];
	VecScale(vel, vel, radius / sqrt(VecProd(vel, vel)));
	Extrude(g, &vVrt, &vEdge, NULL, vel, aaPos, EO_CREATE_FACES, NULL);
	vector3 center = CalculateBarycenter(vVrt.begin(), vVrt.end(), aaPos);
	MergeMultipleVertices(g, vVrt.begin(), vVrt.end());

	Vertex* v = *vVrt.begin();
	aaPos[v] = center;

	aaSurfParams[v].neuriteID = nid;
	aaSurfParams[v].axial = 2.0;
	aaSurfParams[v].angular = 0.0;
	aaSurfParams[v].radial = 1.0;
	}

	////////////////////////////////////////////////////////////////////////
	/// create_neurite_1d
	////////////////////////////////////////////////////////////////////////
	void create_neurite_1d(
		const std::vector<NeuriteProjector::Neurite>& vNeurites,
		const std::vector<std::vector<vector3> >& vPos,
		const std::vector<std::vector<number> >& vR, size_t nid,
		number anisotropy, Grid& g,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		Grid::VertexAttachmentAccessor<
		Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
		Grid::VertexAttachmentAccessor<Attachment<number> >& aaDiam,
		Vertex* connectingVrt
	) {
	const NeuriteProjector::Neurite& neurite = vNeurites[nid];
	const std::vector<vector3>& pos = vPos[nid];
	const std::vector<number>& r = vR[nid];

	number neurite_length = 0.0;
	for (size_t i = 1; i < pos.size(); ++i)
		neurite_length += VecDistance(pos[i], pos[i - 1]);

	size_t nSec = neurite.vSec.size();

	const std::vector<NeuriteProjector::BranchingRegion>& vBR = neurite.vBR;
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit =
			vBR.begin();
	std::vector<NeuriteProjector::BranchingRegion>::const_iterator brit_end =
			vBR.end();

	if (connectingVrt) {
		// ignore first branching region (the connecting region)
		++brit;
	} else {
		// create first vertex
		connectingVrt = *g.create<RegularVertex>();
		aaPos[connectingVrt] = pos[0];
		aaDiam[connectingVrt] = r[0];

		aaSurfParams[connectingVrt].neuriteID = nid;
		aaSurfParams[connectingVrt].axial = 0.0;
		aaSurfParams[connectingVrt].angular = 0.0;
		aaSurfParams[connectingVrt].radial = 0.0;
	}

	// Now create dendrite to the next branching point and iterate this process.
	// We want to create each of the segments with approx. the same aspect ratio.
	// To that end, we first calculate the length of the section to be created (in units of radius)
	// and then divide this number by 2^n where n is the number of anisotropic refinements to
	// be performed to make all segments (more or less) isotropic. The result is the number
	// of segments to be used for the section.
	number t_start = 0.0;
	number t_end = 0.0;
	size_t curSec = 0;

	while (true) {
		t_start = t_end;

		// last section: create until tip
		if (brit == brit_end)
			t_end = 1.0;

		// otherwise: section goes to next branching point
		else
			t_end = brit->t;

		// calculate total length in units of radius
		// = integral from t_start to t_end over: ||v(t)|| / r(t) dt
		number lengthOverRadius = calculate_length_over_radius(t_start, t_end,
				neurite, curSec);

		// to reach the desired anisotropy on the surface in the refinement limit,
		// it has to be multiplied by pi/2 h
		size_t nSeg = (size_t) floor(
				lengthOverRadius / (anisotropy * 0.5 * PI));
		if (nSeg < 1 || lengthOverRadius < 0)
			nSeg = 1;
		number segLength = lengthOverRadius / nSeg;	// segments are between 8 and 16 radii long
		std::vector<number> vSegAxPos(nSeg);
		calculate_segment_axial_positions(vSegAxPos, t_start, t_end, neurite,
				curSec, segLength);

		// create mesh for segments
		Selector sel(g);
		for (size_t s = 0; s < nSeg; ++s) {
			// get exact position and radius of segment end
			number segAxPos = vSegAxPos[s];
			for (; curSec < nSec; ++curSec) {
				const NeuriteProjector::Section& sec = neurite.vSec[curSec];
				if (sec.endParam >= segAxPos)
					break;
			}

			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			vector3 curPos;
			number monom = sec.endParam - segAxPos;
			const number* sp = &sec.splineParamsX[0];
			number& p0 = curPos[0];
			p0 = sp[0] * monom + sp[1];
			p0 = p0 * monom + sp[2];
			p0 = p0 * monom + sp[3];

			sp = &sec.splineParamsY[0];
			number& p1 = curPos[1];
			p1 = sp[0] * monom + sp[1];
			p1 = p1 * monom + sp[2];
			p1 = p1 * monom + sp[3];

			sp = &sec.splineParamsZ[0];
			number& p2 = curPos[2];
			p2 = sp[0] * monom + sp[1];
			p2 = p2 * monom + sp[2];
			p2 = p2 * monom + sp[3];

			number curRad;
			sp = &sec.splineParamsR[0];
			curRad = sp[0] * monom + sp[1];
			curRad = curRad * monom + sp[2];
			curRad = curRad * monom + sp[3];

			// create new vertex and connect with edge
			Vertex* newVrt = *g.create<RegularVertex>();
			*g.create<RegularEdge>(EdgeDescriptor(connectingVrt, newVrt));

			// position and diameter
			aaPos[newVrt] = curPos;
			aaDiam[newVrt] = 2 * curRad;

			// set new param attachments
			aaSurfParams[newVrt].neuriteID = nid;
			aaSurfParams[newVrt].axial = segAxPos;
			aaSurfParams[newVrt].angular = 0.0;
			aaSurfParams[newVrt].radial = 0.0;

			// update
			connectingVrt = newVrt;
		}

		// connect branching neurites if present
		if (brit != brit_end) {
			// find branching child neurite
			SmartPtr<NeuriteProjector::BranchingPoint> bp = brit->bp;
			UG_COND_THROW(bp->vNid.size() > 2,
					"This implementation can only handle branching points with one branching child.");

			size_t child_nid;
			if (bp->vNid[0] != nid)
				child_nid = bp->vNid[0];
			else
				child_nid = bp->vNid[1];

			// add branch neurite ID to its initial vertex
			aaSurfParams[connectingVrt].neuriteID += (brit - vBR.begin()) << 20; // add branching region index
			aaSurfParams[connectingVrt].neuriteID += 1 << 28; // add child ID (always 0, since there can only be one child here)

			create_neurite_1d(vNeurites, vPos, vR, child_nid, anisotropy, g,
					aaPos, aaSurfParams, aaDiam, connectingVrt);
		}

		// update curSec
		for (; curSec < nSec; ++curSec) {
			const NeuriteProjector::Section& sec = neurite.vSec[curSec];
			if (sec.endParam >= t_end)
				break;
		}

		// check whether tip has been reached
		if (brit == brit_end)
			break;
		else
			++brit;
	}
}

////////////////////////////////////////////////////////////////////////
/// calculate_length_over_radius
////////////////////////////////////////////////////////////////////////
number calculate_length_over_radius
(
	number t_start,
	number t_end,
	const NeuriteProjector::Neurite& neurite,
	size_t startSec
)
{
	GaussLegendre gl(5);
	size_t nPts = gl.size();

	std::vector<NeuriteProjector::Section>::const_iterator sec_it = neurite.vSec.begin() + startSec;
	std::vector<NeuriteProjector::Section>::const_iterator sec_end = neurite.vSec.end();

	while (sec_it->endParam < t_start && sec_it->endParam < 1.0)
		++sec_it;

	// check that startSec was correct
	number sec_tstart = startSec > 0 ? (sec_it - 1)->endParam : 0.0;
	number sec_tend = sec_it->endParam;
	UG_COND_THROW(sec_tend < t_start || sec_tstart > t_start,
		"Wrong section iterator given to calculate_length_over_radius().\n"
		"Section goes from " << (sec_it-1)->endParam << " to " << sec_it->endParam
		<< ", but t_start is " << t_start << ".");

                number integral = 0.0;
                while (sec_it != sec_end)
                {
                        // integrate from t_start to min{t_end, sec_tend}
                        const NeuriteProjector::Section& sec = *sec_it;
                        sec_tstart = std::max(t_start, sec_it != neurite.vSec.begin() ? (sec_it - 1)->endParam : 0.0);
                        sec_tend = std::min(t_end, sec.endParam);
                        number dt = sec_tend - sec_tstart;
                        number sec_integral = 0.0;
                        for (size_t i = 0; i < nPts; ++i)
                        {
                                number t = sec.endParam - (sec_tstart + dt*gl.point(i)[0]);

                                vector3 vel;
                                const number* s = &sec.splineParamsX[0];
                                number& v0 = vel[0];
                                v0 = -3.0*s[0]*t - 2.0*s[1];
                                v0 = v0*t - s[2];

                                s = &sec.splineParamsY[0];
                                number& v1 = vel[1];
                                v1 = -3.0*s[0]*t - 2.0*s[1];
                                v1 = v1*t - s[2];

                                s = &sec.splineParamsZ[0];
                                number& v2 = vel[2];
                                v2 = -3.0*s[0]*t - 2.0*s[1];
                                v2 = v2*t - s[2];

                                s = &sec.splineParamsR[0];
                                number r = s[0]*t + s[1];
                                r = r*t + s[2];
                                r = r*t + s[3];

                                UG_COND_THROW(r*r <= VecNormSquared(vel)*1e-12, "r = " << r << " at t = " << t << "!");

                                sec_integral += gl.weight(i) * sqrt(VecNormSquared(vel)) / r;
                        }

                        integral += dt * sec_integral;


                        // update lower bound and iterator
                        t_start = sec_tend;
                        if (t_start >= t_end) break;
                        ++sec_it;
                }

                return integral;
                }

                ////////////////////////////////////////////////////////////////////////
                /// calculate_segment_axial_positions
                ////////////////////////////////////////////////////////////////////////
                void calculate_segment_axial_positions
                (
                        std::vector<number>& segAxPosOut,
                        number t_start,
                        number t_end,
                        const NeuriteProjector::Neurite& neurite,
                        size_t startSec,
                        number segLength
                )
                {
                const size_t nSeg = segAxPosOut.size();

                GaussLegendre gl(5);
                size_t nPts = gl.size();

                std::vector<NeuriteProjector::Section>::const_iterator sec_it = neurite.vSec.begin() + startSec;
                std::vector<NeuriteProjector::Section>::const_iterator sec_end = neurite.vSec.end();

            	while (sec_it->endParam < t_start && sec_it->endParam < 1.0)
            		++sec_it;

                // check that startSec was correct
                number sec_tstart = startSec > 0 ? (sec_it - 1)->endParam : 0.0;
                number sec_tend = sec_it->endParam;
                UG_COND_THROW(sec_tend < t_start || sec_tstart > t_start,
					"Wrong section iterator given to calculate_segment_axial_positions()."
					"Section goes from " << (sec_it-1)->endParam << " to " << sec_it->endParam
					<< ", but t_start is " << t_start << ".");

                number integral = 0.0;
                size_t seg = 0;
                while (sec_it != sec_end)
                {
                        // integrate from t_start to min{t_end, sec_tend}
                        const NeuriteProjector::Section& sec = *sec_it;
                        sec_tstart = std::max(t_start, sec_it != neurite.vSec.begin() ? (sec_it - 1)->endParam : 0.0);
                        sec_tend = std::min(t_end, sec.endParam);
                        number dt = sec_tend - sec_tstart;
                        number sec_integral = 0.0;
                        for (size_t i = 0; i < nPts; ++i)
                        {
                                number t = sec.endParam - (sec_tstart + dt*gl.point(i)[0]);

                                vector3 vel;
                                const number* s = &sec.splineParamsX[0];
                                number& v0 = vel[0];
                                v0 = -3.0*s[0]*t - 2.0*s[1];
                                v0 = v0*t - s[2];

                                s = &sec.splineParamsY[0];
                                number& v1 = vel[1];
                                v1 = -3.0*s[0]*t - 2.0*s[1];
                                v1 = v1*t - s[2];

                                s = &sec.splineParamsZ[0];
                                number& v2 = vel[2];
                                v2 = -3.0*s[0]*t - 2.0*s[1];
                                v2 = v2*t - s[2];

                                s = &sec.splineParamsR[0];
                                number r = s[0]*t + s[1];
                                r = r*t + s[2];
                                r = r*t + s[3];

                                UG_COND_THROW(r*r <= VecNormSquared(vel)*1e-12, "r = " << r << " at t = " << t << "!");

                                sec_integral += gl.weight(i) * sqrt(VecNormSquared(vel)) / r;
                        }
                        integral += dt * sec_integral;

                        // calculate exact position by linear interpolation, whenever integral has surpassed it
                        while (integral >= (seg+1)*segLength)
                        {
                                number lastIntegral = integral - dt * sec_integral;
                                segAxPosOut[seg] = t_start + ((seg+1)*segLength - lastIntegral) / sec_integral;
                                ++seg;
                        }

                        // update lower bound and iterator
                        t_start = sec_tend;
                        if (t_start >= t_end) break;
                        ++sec_it;
                }

                // rounding errors may make this necessary
                if (seg+1 == nSeg && (nSeg*segLength - integral)/integral < 1e-6)
                {
                        segAxPosOut[nSeg-1] = t_end;
                        ++seg;
                }

                UG_ASSERT(seg == nSeg, "seg = " << seg << " != " << nSeg << " = nSeg");
                }

                ////////////////////////////////////////////////////////////////////////
                /// create_neurite_with_er
                ////////////////////////////////////////////////////////////////////////
                void create_neurite_with_er
                (
                        const std::vector<NeuriteProjector::Neurite>& vNeurites,
                        const std::vector<std::vector<vector3> >& vPos,
                        const std::vector<std::vector<number> >& vR,
                        size_t nid,
                        number erScaleFactor,
                        number anisotropy,
                        Grid& g,
                        Grid::VertexAttachmentAccessor<APosition>& aaPos,
                        Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
                        SubsetHandler& sh,
                        std::vector<Vertex*>* outVerts,
                        std::vector<Vertex*>* outVertsInner,
                        std::vector<number>* outRads,
                        std::vector<number>* outRadsInner
                )
                {
                        create_neurite_with_er(vNeurites, vPos, vR, nid, erScaleFactor, anisotropy, g, aaPos, aaSurfParams, sh, NULL, NULL, NULL, 0, outVerts, outVertsInner, outRads, outRadsInner);
                }

                ////////////////////////////////////////////////////////////////////////
                /// create_neurite_root_vertices
                ////////////////////////////////////////////////////////////////////////
                void create_neurite_root_vertices
                (
                        const std::vector<NeuriteProjector::Neurite>& vNeurites,
                        const std::vector<std::vector<vector3> >& vPos,
                        const std::vector<std::vector<number> >& vR,
                        size_t nid,
                        Grid& g,
                        SubsetHandler& sh,
                        number erScaleFactor,
                        Grid::VertexAttachmentAccessor<APosition>& aaPos,
                        std::vector<Vertex*>* outVerts,
                        std::vector<Vertex*>* outVertsInner,
                        std::vector<number>* outRads,
                        std::vector<number>* outRadsInner,
                        bool withER
                )
                {
                        const NeuriteProjector::Neurite& neurite = vNeurites[nid];
                        const std::vector<vector3>& pos = vPos[nid];
                        const std::vector<number>& r = vR[nid];

                        vector3 vel;
                        const NeuriteProjector::Section& sec = neurite.vSec[0];
                        number h = sec.endParam;
                        vel[0] = -3.0 * sec.splineParamsX[0] * h * h
                                        - 2.0 * sec.splineParamsX[1] * h - sec.splineParamsX[2];
                        vel[1] = -3.0 * sec.splineParamsY[0] * h * h
                                        - 2.0 * sec.splineParamsY[1] * h - sec.splineParamsY[2];
                        vel[2] = -3.0 * sec.splineParamsZ[0] * h * h
                                        - 2.0 * sec.splineParamsZ[1] * h - sec.splineParamsZ[2];

                        vector3 projRefDir;
                        VecNormalize(vel, vel);
                        number fac = VecProd(neurite.refDir, vel);
                        VecScaleAdd(projRefDir, 1.0, neurite.refDir, -fac, vel);
                        VecNormalize(projRefDir, projRefDir);

                        vector3 thirdDir;
                        VecCross(thirdDir, vel, projRefDir);

                        std::vector<Vertex*> vVrt;
                        std::vector<Edge*> vEdge;
                        std::vector<Face*> vFace;

                        if (withER) {
                                vVrt.resize(16);
                                vEdge.resize(24);
                                vFace.resize(9);
                        } else {
                                vVrt.resize(12);
                                vEdge.resize(12);
                                vFace.resize(8);
                        }

                        if (withER) {
                                // ER vertices and radii
                                for (size_t i = 0; i < 4; ++i) {
                                        Vertex* v = *g.create<RegularVertex>();
                                        vVrt[i] = v;
                                        number angle = 0.5 * PI * i;
                                        VecScaleAdd(aaPos[v], 1.0, pos[0],
                                                        erScaleFactor * r[0] * cos(angle), projRefDir,
                                                        erScaleFactor * r[0] * sin(angle), thirdDir);

                                        if (outVertsInner) {
                                                outVertsInner->push_back(v);
                                        }

                                        if (outRadsInner) {
                                                outRadsInner->push_back(erScaleFactor * r[0]);
                                        }
                                        sh.assign_subset(v, 3);
                                }
                        }

                        // PM vertices and radii
                        for (size_t i = 0; i < 12; ++i) {
                                Vertex* v = *g.create<RegularVertex>();
                                vVrt[i + 4] = v;
                                number angle = PI * ((number) i / 6);
                                VecScaleAdd(aaPos[v], 1.0, pos[0], r[0] * cos(angle), projRefDir,
                                                r[0] * sin(angle), thirdDir);

                                if (outVerts) {
                                        outVerts->push_back(v);
                                }

                                if (outRads) {
                                        outRads->push_back(r[0]);
                                }
                                sh.assign_subset(v, 2);
                        }

                        // edges
                        if (withER) {
                                for (size_t i = 0; i < 4; ++i) {
                                        vEdge[i] = *g.create<RegularEdge>(
                                                        EdgeDescriptor(vVrt[i], vVrt[(i + 1) % 4]));
                                        vEdge[i + 4] = *g.create<RegularEdge>(
                                                        EdgeDescriptor(vVrt[i], vVrt[5 + 3 * i]));
                                        vEdge[i + 8] = *g.create<RegularEdge>(
                                                        EdgeDescriptor(vVrt[(i + 1) % 4], vVrt[6 + 3 * i]));

                                        sh.assign_subset(vEdge[i], 3);
                                        sh.assign_subset(vEdge[i + 4], 0);
                                }
                        }

                        for (size_t i = 0; i < 12; ++i) {
                                vEdge[i + 12] = *g.create<RegularEdge>(
                                                EdgeDescriptor(vVrt[i + 4], vVrt[(i + 1) % 12 + 4]));
                                sh.assign_subset(vEdge[i + 12], 2);
                        }

                        // faces
                        if (withER) {
                                vFace[0] = *g.create<Quadrilateral>(
                                QuadrilateralDescriptor(vVrt[0], vVrt[1], vVrt[2], vVrt[3]));
                                sh.assign_subset(vFace[0], 1);
                        }

                        for (size_t i = 0; i < 4; ++i) {
                                vFace[i + 1] = *g.create<Quadrilateral>(
                                        QuadrilateralDescriptor(vVrt[i],
                                                vVrt[(3 * i + 11) % 12 + 4], vVrt[3 * i + 4],
                                                vVrt[3 * i + 5]));
                                vFace[i + 5] = *g.create<Quadrilateral>(
                                        QuadrilateralDescriptor(vVrt[i], vVrt[3 * i + 5],
                                                vVrt[3 * i + 6], vVrt[(i + 1) % 4]));
                                sh.assign_subset(vFace[i + 1], 0);
                                sh.assign_subset(vFace[i + 5], 0);
                        }
                }

                ////////////////////////////////////////////////////////////////////////
                /// constrained_smoothing TODO: test
                ////////////////////////////////////////////////////////////////////////
                void constrained_smoothing
                (
                        std::vector<SWCPoint>& vPointsIn,
                        size_t n,
                        number h,
                        number gamma,
                        number minAngle,
                        number maxRadiusRatio
                )
                {
                        // store first vertices after branch which should not be smoothed
                        std::vector<size_t> firstVerticesAfterBranch;

                        // find neurite root vertices
                        const size_t nP = vPointsIn.size();
                        std::vector<size_t> rootVrts;
                        std::vector<bool> treated(nP, false);
                        for (size_t i = 0; i < nP; ++i)
                        {
                                if (treated[i]) continue;
                                treated[i] = true;
                                if (vPointsIn[i].type != SWC_SOMA) continue;

                                // here, we have a soma point;
                                // find first non-soma point in all directions
                                std::queue<size_t> q;
                                q.push(i);

                                while (!q.empty())
                                {
                                        const size_t ind = q.front();
                                        const SWCPoint& pt = vPointsIn[ind];
                                        q.pop();

                                        if (pt.type == SWC_SOMA)
                                        {
                                                const size_t nConn = pt.conns.size();
                                                for (size_t j = 0; j < nConn; ++j)
                                                        if (!treated[pt.conns[j]])
                                                                q.push(pt.conns[j]);
                                        }
                                        else
                                                rootVrts.push_back(ind);

                                        treated[ind] = true;
                                }
                        }

                        // starting at root vertices, smooth the entire tree(s),
                        // but leave out soma vertices as well as branching points
                        std::vector<vector3> newPos(nP);
                        for (size_t i = 0; i < n; ++i)
                        {
                                treated.clear();
                                treated.resize(nP, false);

                                std::stack<size_t> stack;
                                for (size_t rv = 0; rv < rootVrts.size(); ++rv)
                                        stack.push(rootVrts[rv]);

                                while (!stack.empty())
                                {
                                        size_t ind = stack.top();
                                        stack.pop();
                                        const SWCPoint& pt = vPointsIn[ind];
                                        const vector3& x = pt.coords;
                                        vector3& x_new = newPos[ind];

                                        UG_COND_THROW(treated[ind], "Circle detected in supposedly tree-shaped neuron!\n"
                                                "Position: " << vPointsIn[ind].coords);
                                        treated[ind] = true;

                                        // somata are not smoothed and not iterated over
                                        if (pt.type == SWC_SOMA)
                                        {
                                                x_new = x;
                                                continue;
                                        }

                                        // mark branching points as treated
                                        std::vector<size_t> conns = pt.conns;
                                        for (size_t c = 0; c < conns.size(); ++c) {
                                                if (!treated[conns[c]]) {
                                                        stack.push(conns[c]);
                                                }
                                        }

                                        // actually treat branching points
                                        if (conns.size() != 2) {
                                                std::vector<number> radii;
                                                for (size_t i = 0; i < conns.size(); ++i) {
                                                   radii.push_back(vPointsIn[conns[i]].radius);
                                                }

                                                // largest radius of all connecting points originating from
                                                // current point x (pt.conns.coords) deemed as root branch
                                                // continuation of current neurite the algorithm is tracing
                                                int max = std::distance(radii.begin(), std::max_element(radii.begin(), radii.end()));
                                                UG_COND_THROW(max < maxRadiusRatio, "Diameter of neurite does not satify main branch criterion.");

                                                // branch with min angle is assumed to be root neurite branch continuation
                                                number min = PI/2;
                                                std::pair<size_t, size_t> minPair;
                                                for (size_t i = 0; conns.size(); ++i) {
                                                        for (size_t j = 0; conns.size(); ++j) {
                                                                ug::vector3 v1;
                                                                ug::vector3 v2;
                                                                VecSubtract(v1, vPointsIn[conns[i]].coords, x);
                                                                VecSubtract(v2, vPointsIn[conns[j]].coords, x);
                                                                VecNormalize(v1, v1);
                                                                VecNormalize(v2, v2);
                                                                number angle = acos(VecDot(v1, v2));
                                                                if (angle < min) {
                                                                        min = angle;
                                                                        minPair = std::pair<size_t, size_t>(i, j);
                                                                }
                                                        }
                                                }

                                                /// min angle branch which matches the largest diameter deemed as root branch continuation
                                                if ((minPair.second == (size_t) max) || (minPair.first == (size_t) max)) {
                                                        conns.erase(conns.begin() + max);
                                                }

                                                // remember first vertex after or before branch of non root branch of current neurite
                                                for (size_t i = 0; i < conns.size(); i++) {
                                                        firstVerticesAfterBranch.push_back(conns[i]);
                                                }

                                                continue;
				}

				// ignore non-root branch first vertices for smoothing
				if(std::find(firstVerticesAfterBranch.begin(), firstVerticesAfterBranch.end(),
						pt.conns[0]) != firstVerticesAfterBranch.end()) {
					continue;
				}

				if(std::find(firstVerticesAfterBranch.begin(), firstVerticesAfterBranch.end(),
						pt.conns[1]) != firstVerticesAfterBranch.end()) {
					continue;
				}

				// here we have a non-branching, non-end, non-soma point: smooth
				const vector3& x1 = vPointsIn[pt.conns[0]].coords;
				const vector3& x2 = vPointsIn[pt.conns[1]].coords;

				number d1 = VecDistanceSq(x1, x);
				number d2 = VecDistanceSq(x2, x);
				number w1 = std::exp(-d1/(h*h));
				number w2 = std::exp(-d2/(h*h));

				// only really smooth if both adjacent edges are short
				number w = std::min(w1, w2);

				// correction
				vector3 corr;
				VecScaleAdd(corr, w, x1, -2*w, x, w, x2);
				VecScale(corr, corr, 1.0 / (1.0 + 2*w));

				// take only the part orthogonal to x1 - x2,
				// we do not want to shift x towards the nearer neighbor
				VecSubtract(x_new, x1, x2); // using x_new as intermediate variable
				number normSq = VecNormSquared(x_new);
				VecScaleAdd(corr, 1.0, corr, - VecProd(corr, x_new) / normSq, x_new);
				VecScaleAdd(x_new, 1.0, x, gamma, corr);
			}

			// assign new positions
			for (size_t p = 0; p < nP; ++p)
				if (treated[p]) { // soma points may not have been treated
					vPointsIn[p].coords = newPos[p];
				}
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// regularize_bps
		////////////////////////////////////////////////////////////////////////
		void regularize_bps
		(
			std::vector<SWCPoint>& vPoints,
			bool orthogonalize
		) {
			size_t nPts = vPoints.size();
			std::vector<bool> ptProcessed(nPts, false);
			size_t nProcessed = 0;

			while (nProcessed != nPts)
			{
				// find first soma's root point in geometry and save its index as i
				size_t i = 0;
				for (; i < nPts; ++i) {
					if (vPoints[i].type == SWC_SOMA && !ptProcessed[i])
						break;
				}

				// collect neurite root points
				std::vector<std::pair<size_t, size_t> > rootPts;
				std::queue<std::pair<size_t, size_t> > soma_queue;
				soma_queue.push(std::make_pair((size_t)-1, i));
				while (!soma_queue.empty())
				{
					size_t pind = soma_queue.front().first;
					size_t ind = soma_queue.front().second;
					soma_queue.pop();

					const SWCPoint& pt = vPoints[ind];

					if (pt.type == SWC_SOMA)
					{
						ptProcessed[ind] = true;
						++nProcessed;

						size_t nConn = pt.conns.size();
						for (size_t j = 0; j < nConn; ++j) {
							if (pt.conns[j] != pind) {
								soma_queue.push(std::make_pair(ind, pt.conns[j]));
							}
						}
					}
					else
					{
						rootPts.push_back(std::make_pair(pind, ind));
					}
				}

				std::stack<std::pair<size_t, size_t> > processing_stack;
				for (size_t i = 0; i < rootPts.size(); ++i)
					processing_stack.push(rootPts[i]);

				while (!processing_stack.empty())
				{
					size_t pind = processing_stack.top().first;
					size_t ind = processing_stack.top().second;
					processing_stack.pop();

					SWCPoint& pt = vPoints[ind];
					UG_COND_THROW(pt.type == SWC_SOMA, "Detected neuron with more than one soma.");
					size_t nConn = pt.conns.size();
					++nProcessed;

					UG_COND_THROW(nConn > 3, "nConn > 3 not supported yet.");

					// branching point
					if (nConn > 2)
					{
						// branch with minimal angle will continue current branch
						vector3 parentDir;
						VecSubtract(parentDir, pt.coords, vPoints[pind].coords);
						VecNormalize(parentDir, parentDir);

						size_t parentToBeDiscarded = 0;
						size_t minAngleInd = 0;
						number minAngle = std::numeric_limits<number>::infinity();

						for (size_t i = 0; i < nConn; ++i)
						{
							if (pt.conns[i] == pind)
							{
								parentToBeDiscarded = i;
								continue;
							}

							vector3 dir;
							VecSubtract(dir, vPoints[pt.conns[i]].coords, pt.coords);
							VecNormalize(dir, dir);

							number angle = acos(VecProd(dir, parentDir));
							if (angle < minAngle)
							{
								minAngle = angle;
								minAngleInd = i;
							}
						}

						/// project only branching point onto branching edge
						// current pt (branching point) = B
						// parentToBeDiscarded = P
						// minAngleInd = Q
						// pt.conns[i] = element of vector R
						// add to vPoints
						const vector3& B = pt.coords;
						vector3 Bprime;
						const vector3& Q = vPoints[pt.conns[minAngleInd]].coords;
						const vector3& P = vPoints[pt.conns[parentToBeDiscarded]].coords;
						ProjectPointToLine(Bprime, B, P, Q);
						pt.coords = Bprime;

						number QtoPt = VecDistance(Q, pt.coords) / VecDistance(Q, P);
						number PtoPt = VecDistance(P, pt.coords) / VecDistance(Q, P);
						pt.radius = (vPoints[pt.conns[minAngleInd]].radius * QtoPt)
								+ (vPoints[pt.conns[parentToBeDiscarded]].radius * PtoPt);

						/// "orthogonalize" branching points
						if (orthogonalize) {
#if 0
							SWCPoint pointA;
							vector3 dir;
							VecCross(dir, Bprime, Q); // Bprime is projected point onto line
							VecNormalize(dir, dir);

							vector3 A, vector, axis;
							VecAdd(A, Bprime, dir); // new point A orthogonal to Bprime and edge PQ

							VecSubtract(axis, A, Bprime);
							VecSubtract(vector, Q, Bprime);

							/// rotated vectors
							vector3 xPrime;
							vector3 xPrime2;
							rotate_vector_around_axis(vector, axis, Bprime, xPrime, xPrime2, 90);

							Grid g;
							SubsetHandler sh;
							Grid::VertexAttachmentAccessor<APosition> aaPos;
							g.attach_to_vertices(aPosition);
							aaPos = Grid::VertexAttachmentAccessor<APosition>(g, aPosition);
							sh = SubsetHandler(g);
							sh.set_default_subset_index(0);

							Vertex* p1, *p2, *p3, *p4;
							p1 = *g.create<RegularVertex>(); p2 = *g.create<RegularVertex>();
							p3 = *g.create<RegularVertex>(); p4 = *g.create<RegularVertex>();
							aaPos[p1] = Bprime;
							aaPos[p2] = xPrime;
							aaPos[p3] = A;
							aaPos[p4] = Q;

							ug::RegularEdge* e3 = *g.create<RegularEdge>(EdgeDescriptor(p1, p2));
							ug::RegularEdge* e4 = *g.create<RegularEdge>(EdgeDescriptor(p1, p3));
							ug::RegularEdge* e5 = *g.create<RegularEdge>(EdgeDescriptor(p1, p4));

							sh.assign_grid(g);
							AssignSubsetColors(sh);
							SaveGridToFile(g, sh, "test_rotation_with_real_grid.ugx");

							size_t continuation;
							for (size_t j = 0; j < nConn; ++j)
							{
								if (j == parentToBeDiscarded || j == minAngleInd)
								{
									continue;
								}
								continuation=j;
							}

							pointA.coords = xPrime;
#endif

#if 0
							vector3 xPrime;
							vector3 xPrime2;
							number minDist = std::numeric_limits<number>::infinity();
							vector3 minPt;
							for (int i = 0; i < 180; i++) {
								rotate_vector_around_axis(axis, dir, Bprime, xPrime, xPrime2, i);
								number dist1 = VecDistance(vPoints[pt.conns[continuation]].coords, xPrime);
								number dist2 = VecDistance(vPoints[pt.conns[continuation]].coords, xPrime2);
								if (dist1 < dist2) {
									if (dist1 < minDist) {
										minPt = xPrime;
										minDist = dist1;
									}
									if (dist2 < minDist) {
										minPt = xPrime2;
										minDist = dist2;
									}
								}
								/*
								 minDist = std::min(dist1, dist2);
								 minPt = (minDist == dist1 ? xPrime : xPrime2);

								 */
							}

							pointA.coords = minPt;
#endif


							#if 1
							SWCPoint pointA;
							vector3 dir, dirAlt;
							VecCross(dir, P, Q);
							VecNormalize(dir, dir);


							size_t continuation;
							for (size_t j = 0; j < nConn; ++j)
							{
								if (j == parentToBeDiscarded || j == minAngleInd)
								{
									continue;
								}
								continuation=j;
							}

							/// Could scale with distance from original point B to continuation point of branch
							number scaleFactor = VecDistance(vPoints[pt.conns[continuation]].coords, B);
							VecScale(dir, dir, scaleFactor);
							///VecScale(dir, dir, 1.5); // was: 1
							/// VecDistance(vPoints[pt.conns[continuation]].coords, B);
							/*UG_COND_THROW(VecLength(dir) < 1-SMALL, "VecLength(dir) < 1."
									" Minimum length of dir is one however actual"
									" length is " << VecLength(dir));
									*/

							vector3 A, Aalt;
							VecAdd(A, Bprime, dir);
							dirAlt = -dir;
							VecAdd(Aalt, Bprime, dirAlt);

							/// Quick fix for 3 way branches - TODO: rotate around axis to find best direction
							if (VecDistance(vPoints[pt.conns[continuation]].coords, A) <
								VecDistance(vPoints[pt.conns[continuation]].coords, Aalt)) {
								pointA.coords = A;
							} else {
								pointA.coords = Aalt;
							}
							#endif

							pointA.radius = pt.radius; // vPoints[pt.conns[continuation]].radius;
							pointA.type = vPoints[pt.conns[minAngleInd]].type;
							size_t newIndex = vPoints.size();
							std::vector<size_t> Bs;
							size_t Rindex;

							for (size_t j = 0; j < nConn; ++j)
							{
								if (j == parentToBeDiscarded || j == minAngleInd)
								{
									continue;
								}

								/// Removes current index of B (ind) from connected
								/// branching child neurite neighbor (which is neither
								/// the point after branching point which continues the
								/// main branch or the point before the branching point
								/// called parent). The connected branching child neurite
								/// is connected to the new intermediate point pointA
								/// with index newIndex instead and not to B anymore
								vPoints[pt.conns[j]].conns.push_back(newIndex);
								vPoints[pt.conns[j]].conns.erase(std::remove(vPoints[pt.conns[j]].conns.begin(), vPoints[pt.conns[j]].conns.end(), ind), vPoints[pt.conns[j]].conns.end());
								pointA.conns.push_back(pt.conns[j]);
								nPts++;
								Rindex = j;
							}

							/// old branching point B (now with new position) is
							/// still connected to main branch continuation with
							/// index minAngleInd and is still connected to the
							/// parent point with index parentToBeDiscarded. in
							/// addition the old branching point is connected to
							/// the new intermediate point but not to the other
							/// point after the branching point which starts a
							/// new child neurite branch with index pt.conns[j].
							std::vector<size_t> Rs;
							Rs.resize(nConn);
							Rs[minAngleInd] = pt.conns[minAngleInd];
							Rs[parentToBeDiscarded] = pt.conns[parentToBeDiscarded];
							Rs[Rindex] = newIndex;
							pt.conns = Rs;

							/// pointA is the new intermediate point connected to
							/// old branching point B and to the start of a new
							/// branching neurite child point after the original
							/// branching point B. The intermediate point is only
							/// connected to these two points: old branching point
							/// B with index ind and to the branching child neighbor
							/// with index pt.conns[j] where j is not minAngleIndex
							/// or parentIndexToDiscard which is the main branch
							/// continuation and the parent point respectively
							pointA.conns.push_back(ind);
							vPoints.push_back(pointA);
							ptProcessed.push_back(false);
						}

						for (size_t j = 0; j < nConn; ++j)
						{
							if (j == parentToBeDiscarded || j == minAngleInd)
							{
								continue;
							}

							// push new neurite starting point index to stack
							// connected to new intermediate point with
							// newIndex = pt.conns[j]
							processing_stack.push(std::make_pair(ind, pt.conns[j]));

						}

						// push next index of the current neurite to stack
						processing_stack.push(std::make_pair(ind, pt.conns[minAngleInd]));
					}

					// root point
					else if (nConn == 1)
					{
						// nothing to do
					}
					// normal point
					else
					{
						for (size_t i = 0; i < nConn; ++i)
						{
							if (pt.conns[i] != pind)
							{
								processing_stack.push(std::make_pair(ind, pt.conns[i]));
							}
						}
					}
				}
			}
		}
	}
}
