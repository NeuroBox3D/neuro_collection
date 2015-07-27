/*
 * bouton_generator.cpp
 *
 *  Created on: 03.09.2013
 *      Author: Martin Stepniewski
 */

#include "bouton_generator.h"


namespace ug{
namespace neuro_collection{


////////////////////////////////////////////////////////////////////////////////////////////
//	DistributeNPointsOnASphere
////////////////////////////////////////////////////////////////////////////////////////////
void GetNEvenlyDistributedSphereCoords(vector<vector3>& coords, int N, double radius)
{
/*
 * Check, if N defines a platonic solid
 */

/*
//	Tetrahedron
	if(N == 4)
	{
		const number A = sqrt(2);
		const number B = sqrt(6);
		const number vrts[4][3] = {	{ 0,        0      ,  1},
									{ 2.0/3*A,  0      , -1.0/3},
									{-1.0/3*A,  1.0/3*B, -1.0/3},
									{-1.0/3*A, -1.0/3*B, -1.0/3}};

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;
			coords.push_back(tmpCoords);
		}
	}


//	Octahedron
	else if(N == 6)
	{
		const number vrts[6][3] = {	{1,0,0},
									{-1,0,0},
									{0,1,0},
									{0,0,1},
									{0,-1,0},
									{0,0,-1}
								  };

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;
			coords.push_back(tmpCoords);
		}
	}


//	Cube
	else if(N == 8)
	{
		const number c = 1.0/sqrt(3);
		const number vrts[8][3] = {	{c,c,c},
									{-c,c,c},
									{c,-c,c},
									{c,c,-c},
									{-c,-c,-c},
									{-c,-c,c},
									{-c,c,-c},
									{c,-c,-c}
								  };

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;
			coords.push_back(tmpCoords);
		}
	}


//	Icosahedron
	else if(N == 12)
	{
		const number A = 0.85065080835204;
		const number B = 0.525731112119134;
		const number vrts[12][3] = {	{-B, A, 0},
										{0, B, A},
										{B, A, 0},
										{0, B, -A},
										{-A, 0, B},
										{A, 0, B},
										{A, 0, -B},
										{-A, 0, -B},
										{-B, -A, 0},
										{0, -B, A},
										{B, -A, 0},
										{0, -B, -A}
									};

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;
			coords.push_back(tmpCoords);
		}
	}


//	Dodecahedron
	else if(N == 20)
	{
		const number c = 1.0/sqrt(3);
		const number x = (sqrt(5)+1)/(2*sqrt(3));
		const number y = (sqrt(5)-1)/(2*sqrt(3));
		const number vrts[20][3] = {	{c,c,c}, {c,-c,c}, {c,c,-c}, {c,-c,-c},
										{-c,c,c}, {-c,-c,c}, {-c,c,-c}, {-c,-c,-c},
										{x,y,0}, {-x,y,0}, {x,-y,0}, {-x,-y,0},
										{0,x,y}, {0,-x,y}, {0,x,-y}, {0,-x,-y},
										{y,0,x}, {-y,0,x}, {y,0,-x}, {-y,0,-x}};

		vector3 tmpCoords;
		for(int i = 0; i < N; i++)
		{
			tmpCoords = vector3(vrts[i][0], vrts[i][1], vrts[i][2]);
			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;
			coords.push_back(tmpCoords);
		}
	}
*/

/*
 * otherwise use heuristic
 */

	//else
	//{
		/***********************
		 * Golden spiral section
		 * (2nd best algorithm)
		 **********************/
		/*
		int n = N;
		double r, phi;
		vector3 tmpCoords;

		for(int i = 0; i < n; i++)
		{
			tmpCoords.y() 	= i * (double)2/n - 1 + (double)1/n;
			r				= sqrt(1.0 - tmpCoords.y() * tmpCoords.y());
			phi 			= i * M_PI * (3.0 - sqrt(5));
			tmpCoords.x		= r * cos(phi);
			tmpCoords.z() 	= r * sin(phi);

			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;

			coords.push_back(tmpCoords);
		}
		*/


		/****************************************
		 * "Distributing many points on a sphere"
		 * by E.B. Saff and A.B.J. Kuijlaars
		 ***************************************/

		/*
		int n = N;
		double theta[n+1], phi[n+1], h;
		vector3 tmpCoords;

		for(int i = 1; i <= n; i++)
		{
			h = -1 + 2 * (double)(i-1)/(n-1);
			theta[i] = acos(h);

			if(i == 1 || i == n)
				phi[i] = 0;
			else
				phi[i] = fmod(phi[i-1] + 3.6/sqrt(n*(1.0-h*h)), (2*M_PI));

			tmpCoords.x() = sin(theta[i])*cos(phi[i]);
			tmpCoords.y() = sin(phi[i])*sin(theta[i]);
			tmpCoords.z() = cos(theta[i]);

			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;

			coords.push_back(tmpCoords);
		}
		*/


		/****************************************
		 * "Distributing many points on a sphere"
		 * by E.B. Saff and A.B.J. Kuijlaars
		 *
		 * Modification #1
		 ***************************************/

		/*
		int n = N;
		vector3 tmpCoords;
		double p = 0.5;
		double a = 1 - 2 * p / (n-3);
		double b = p * (double)(n+1)/(n-1);
		double theta[n+1], phi[n+1], h[n+1], r[n+1];

		for(int i = 1; i <= n; i++)
		{
			if(i == 1)
			{
				r[1] 		= 0.0;
				theta[1] 	= M_PI;
				phi[1] 		= 0.0;
			}

			else if (i == n)
			{
				theta[n] 	= 0.0;
				phi[n] 		= 0.0;
			}

			else
			{
				double k 	= a * i + b;
				h[i] 		= -1 + 2 * (k-1)/(n-1);
				r[i] 		= sqrt(1 - h[i]*h[i]);
				theta[i] 	= acos(h[i]);
				phi[i] 		= fmod(phi[i-1] + 3.6/sqrt(n)*(double)2/(r[i-1]+r[i]), (2*M_PI));
			}

			tmpCoords.x() = sin(theta[i])*cos(phi[i]);
			tmpCoords.y() = sin(phi[i])*sin(theta[i]);
			tmpCoords.z() = cos(theta[i]);

			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;

			coords.push_back(tmpCoords);
		}
		*/


		/****************************************
		 * "Distributing many points on a sphere"
		 * by E.B. Saff and A.B.J. Kuijlaars
		 *
		 * Modification #2
		 * (best algorithm)
		 ***************************************/

		int n = N;
		vector3 tmpCoords;
		const int size = n+1;
		vector<double> theta(size);
		vector<double> phi(size);
		vector<double> h(size);

		for(int i = 1; i <= n; i++)
		{
			h[i] 		= -1 + (double)(2*i-1)/n;
			theta[i] 	= acos(h[i]);
			phi[i] 		= sqrt(n*PI)*theta[i];

			tmpCoords.x() = sin(theta[i])*cos(phi[i]);
			tmpCoords.y() = sin(phi[i])*sin(theta[i]);
			tmpCoords.z() = cos(theta[i]);

			tmpCoords.x() *= radius;
			tmpCoords.y() *= radius;
			tmpCoords.z() *= radius;

			coords.push_back(tmpCoords);
		}
	//}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	BuildBouton
////////////////////////////////////////////////////////////////////////////////////////////
void BuildBouton(	bool bExtSpace, number radius, int numRefinements, int numReleaseSites,
					number TbarHeight,
					number TbarLegRadius,
					number TbarTopRadius,
					number TbarTopHeight,
					string fileName)
{
//	Initial grid management setup
	Grid grid;
	AInt aInt;
	grid.attach_to_vertices(aInt);
	grid.attach_to_vertices(aPosition);
	grid.attach_to_vertices(aNormal);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	Grid::VertexAttachmentAccessor<ANormal> aaNorm(grid, aNormal);

	SubsetHandler sh(grid);

	Selector sel(grid);
	Selector tmpSel(grid);


//	Final subset management

	// Basic subsets:
	const int si_bouton_bnd 	= 0;
	const int si_mature_AZ 		= 1;
	const int si_immature_AZ 	= 2;
	const int si_CChCl 			= 3;
	const int si_Tbars_bnd 		= 4;
	const int si_mit_bnd 		= 5;
	const int si_cyt_int 		= 6;
	const int si_locks			= 7;
	const int si_Tbars_int 		= 8;
	const int si_mit_int		= 9;
	const int si_ext_int		=10;
	const int si_ext_bnd		=11;

	// Probe subsets:
	if(numReleaseSites < 2 || numReleaseSites%2 != 0)
		UG_THROW("ERROR in bouton_generator: numReleaseSites has to be EVEN and >= 2.");

	vector<int> si_probe(numReleaseSites/2);
	vector<int> si_recycling(numReleaseSites/2);
	int probeSubsetsBegin = 10;

	if(bExtSpace)
		probeSubsetsBegin = 12;

	for(int i = 0; i < numReleaseSites/2; ++i)
	{
		si_probe[i] = i + probeSubsetsBegin;
	}

	// Recycling subsets:
	for(int i = 0; i < numReleaseSites/2; ++i)
	{
		si_recycling[i] = i + probeSubsetsBegin + numReleaseSites/2;
	}


////
//	Generate raw icosphere
////
	vector3 center(0.0, 0.0, 0.0);
	sh.set_default_subset_index(si_bouton_bnd);
	GenerateIcosphere(grid, center, radius, numRefinements, aPosition);


////
//	Create mitochondrium
////
	number mit_radius = 0.25;
	sh.set_default_subset_index(si_mit_bnd);
	GenerateIcosphere(grid, center, mit_radius, 1, aPosition);


//	Get #numReleaseSites evenly distributed sphere coordinates
	vector<vector3> coords;
	GetNEvenlyDistributedSphereCoords(coords, numReleaseSites, radius);

/*
 * 	Rotate evenly distributed sphere coordinates around x axes by "a" degrees
 * 	so that every t-bar is of hexagonal shape
 */
	double a = 10.0;
	for(size_t i = 0; i < coords.size(); ++i)
	{
		double y, z;
		y = coords[i].y();
		z = coords[i].z();
		coords[i].y() = y * cos(a) - z * sin(a);
		coords[i].z() = y * sin(a) + z * cos(a);
	}


//	Testwise creation of evenly distributed vertices
	/*
	Vertex* vrts[numReleaseSites];
	for(size_t i = 0; i < numReleaseSites; ++i)
	{
		vrts[i] = *grid.create<RegularVertex>();
		aaPos[vrts[i]] = coords[i];
	}
	*/


//	Find corresponding vertices on the icosphere and assign evenly distributed vertices for release sites
	number minDist, tmpMinDist;

	for(VertexIterator vIter = grid.vertices_begin(); vIter != grid.vertices_end(); ++vIter)
	{
		Vertex* vrt = *vIter;
		sel.select(vrt);
	}

	if(grid.num<RegularVertex>() > 0)
	{
		for(size_t i = 0; i < coords.size(); ++i)
		{
			bool gotOne = false;
			Vertex* tmpVrt;

			for(VertexIterator vIter = grid.vertices_begin(); vIter != grid.vertices_end(); ++vIter)
			{
				Vertex* vrt = *vIter;
				tmpMinDist = VecDistance(aaPos[vrt], coords[i]);

				if(((!gotOne) || (tmpMinDist < minDist)) && sel.is_selected(vrt))
				{
					minDist = tmpMinDist;
					tmpVrt = vrt;
					gotOne = true;
				}
			}

			// throw if !gotOne? -- otherwise the following instructions could be
			// executed with tmpVrt uninit'ed or with a value from the previous iteration!
			//if (!gotOne) UG_THROW("No vertex found.");

			sel.deselect(tmpVrt);

		//	spiral distribution of mature and immature release sites
			if(i % 2 == 0)
			{
				sh.assign_subset(tmpVrt, si_Tbars_bnd);
			}
			else
			{
				sh.assign_subset(tmpVrt, si_immature_AZ);
			}
		}
	}


//	Calculate vertex normals for t-bar extrusion
	CalculateVertexNormals(grid, aaPos, aaNorm);


//	store evenly distributed vertices of mature release sites to T-bars_bnd subset
	vector<Vertex*> vrts;
	for(VertexIterator vIter = sh.begin<Vertex>(si_Tbars_bnd); vIter != sh.end<Vertex>(si_Tbars_bnd); ++vIter)
	{
		Vertex* vrt = *vIter;
		vrts.push_back(vrt);
	}


////
//	Create mature release sites
////
	vector3 vExtrDir;
	vector3 vUnitExtrDir;
	vector<Edge*> vExtrusionEdges;
	vector<Edge*> vTmpExtrusionEdges;

	for(size_t i = 0; i < vrts.size(); ++i)
	{
		tmpSel.clear();
		vExtrusionEdges.clear();
		vTmpExtrusionEdges.clear();

	//	Calculate extrusion norm for probe subset
		vector3 negNorm;
		Vertex* vrt = vrts[i];
		VecScale(negNorm, aaNorm[vrt], -1.0);

	//	Hack especially for 1b boutons
		if(numReleaseSites == 20)
		{
			QualityGridGeneration(grid, sh.begin<Face>(si_bouton_bnd), sh.end<Face>(si_bouton_bnd), aaPos, 30.0);
			AdaptSurfaceGridToCylinder(sel, grid, vrt, negNorm, 0.25, 0.05);
			Refine(grid, sel);
			sel.clear();
		}

		AdaptSurfaceGridToCylinder(sel, grid, vrt, negNorm, 0.25, 0.05);

		for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, si_mature_AZ);
			for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f); eIter != grid.associated_edges_end(f); ++eIter)
			{
				Edge* e = *eIter;
				sh.assign_subset(e, si_mature_AZ);
			}

			for(size_t i = 0; i < f->num_vertices(); ++i)
			{
				Vertex* vrt = f->vertex(i);
				sh.assign_subset(vrt, si_mature_AZ);
			}
		}

	//	Select mature_AZ boundary for extrusion of probe subset
		SelectAreaBoundary(tmpSel, sel.begin<Face>(), sel.end<Face>());
		//	No refinement here in favor of a coarser grid
		//Refine(grid, sel);
		for(EdgeIterator eIter = tmpSel.begin<Edge>(); eIter != tmpSel.end<Edge>(); ++eIter)
		{
			Edge* e = *eIter;
			vExtrusionEdges.push_back(e);
			vTmpExtrusionEdges.push_back(e);

			sh.assign_subset(e, si_probe[i]);
			sh.assign_subset(e->vertex(0), si_probe[i]);
			sh.assign_subset(e->vertex(1), si_probe[i]);
		}

		vUnitExtrDir = aaNorm[vrt];
		VecScale(vExtrDir, vUnitExtrDir, -1.0 * (TbarHeight+TbarTopHeight));
		//	No double extrusion here in favor of a coarser grid
		// 	Info: factor 0.5 because of following double extrusion!
		//VecScale(vExtrDir, vUnitExtrDir, -1.0 * (TbarHeight+TbarTopHeight)*0.5);
		//Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);
		Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);

	//	Reassign root extrusion edges to mature_AZ subset
		for(size_t i = 0; i < vTmpExtrusionEdges.size(); ++i)
		{
			Edge* e = vTmpExtrusionEdges[i];

			sh.assign_subset(e, si_mature_AZ);
			sh.assign_subset(e->vertex(0), si_mature_AZ);
			sh.assign_subset(e->vertex(1), si_mature_AZ);
		}

	//	Triangulate top of probe subset
		tmpSel.clear();
		for(size_t j = 0; j < vExtrusionEdges.size(); ++j)
		{
			Edge* e = vExtrusionEdges[j];
			tmpSel.select(e);
		}
		sh.set_default_subset_index(si_probe[i]);
		TriangleFill_SweepLine(grid, tmpSel.begin<Edge>(), tmpSel.end<Edge>(), aPosition, aInt, &sh, si_probe[i]);

		//	No refinement here in favor of a coarser grid
		//Refine(grid, sel);
	}

//	Optimize probe subset triangulation
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		sel.clear();
		for(FaceIterator fIter = sh.begin<Face>(si_probe[i]); fIter != sh.end<Face>(si_probe[i]); ++fIter)
		{
			Face* f = *fIter;
			sel.select(f);
		}
		sh.set_default_subset_index(si_probe[i]);
		Refine(grid, sel);
		QualityGridGeneration(grid, sh.begin<Face>(si_probe[i]), sh.end<Face>(si_probe[i]), aaPos, 30.0);
	}


////
//	Create CChCl
////
	sel.clear();
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		vector3 negNorm;
		Vertex* vrt = vrts[i];
		VecScale(negNorm, aaNorm[vrt], -1.0);
		AdaptSurfaceGridToCylinder(sel, grid, vrt, negNorm, 0.1, 0.02);

		for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, si_CChCl);
			for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f); eIter != grid.associated_edges_end(f); ++eIter)
			{
				Edge* e = *eIter;
				sh.assign_subset(e, si_CChCl);
			}

			for(size_t i = 0; i < f->num_vertices(); ++i)
			{
				Vertex* vrt = f->vertex(i);
				sh.assign_subset(vrt, si_CChCl);
			}
		}

		Refine(grid, sel);
	}


////
//	Create T-bars_bnd
////
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		Vertex* vrt = vrts[i];
		BuildTbar(grid, sh, vrt, aaPos, aaNorm, si_Tbars_bnd, TbarHeight*0.25, TbarLegRadius, TbarTopRadius, TbarTopHeight);
	}

//	Reassign CChCl boundary edges to CChCl subset
	sel.clear();
	SelectAreaBoundary(sel, sh.begin<Face>(si_CChCl), sh.end<Face>(si_CChCl));
	for(EdgeIterator eIter = sel.begin<Edge>(); eIter != sel.end<Edge>(); ++eIter)
	{
		Edge* e = *eIter;
		sh.assign_subset(e, si_CChCl);
		sh.assign_subset(e->vertex(0), si_CChCl);
		sh.assign_subset(e->vertex(1), si_CChCl);
	}

//	Reassign T-bars_bottom
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(vrts[i]); eIter != grid.associated_edges_end(vrts[i]); ++eIter)
		{
			Edge* e = *eIter;
			sh.assign_subset(e, si_locks);
		}

		for(Grid::AssociatedFaceIterator fIter = grid.associated_faces_begin(vrts[i]); fIter != grid.associated_faces_end(vrts[i]); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, si_locks);
		}

		sh.assign_subset(vrts[i], si_locks);
	}


////
//	Create immature_AZ
////
	vector<Vertex*> immature_vrts;
	for(VertexIterator vIter = sh.begin<Vertex>(si_immature_AZ); vIter != sh.end<Vertex>(si_immature_AZ); ++vIter)
	{
		Vertex* vrt = *vIter;
		immature_vrts.push_back(vrt);
	}

	for(size_t i = 0; i < immature_vrts.size(); ++i)
	{
	//	Calculate extrusion norm of immatureAZ subset for procedure AdaptSurfaceGridToCylinder
		vector3 negNorm;
		Vertex* vrt = immature_vrts[i];
		VecScale(negNorm, aaNorm[vrt], -1.0);

	//	Hack especially for 1b boutons
		if(numReleaseSites == 20)
		{
			QualityGridGeneration(grid, sh.begin<Face>(si_bouton_bnd), sh.end<Face>(si_bouton_bnd), aaPos, 30.0);
			AdaptSurfaceGridToCylinder(sel, grid, vrt, negNorm, 0.25, 0.05);
			Refine(grid, sel);
			sel.clear();
		}

		AdaptSurfaceGridToCylinder(sel, grid, vrt, negNorm, 0.25, 0.05);

		for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
		{
			Face* f = *fIter;
			sh.assign_subset(f, si_immature_AZ);

			for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f); eIter != grid.associated_edges_end(f); ++eIter)
			{
				Edge* e = *eIter;
				sh.assign_subset(e, si_immature_AZ);
			}

			for(size_t i = 0; i < f->num_vertices(); ++i)
			{
				Vertex* vrt = f->vertex(i);
				sh.assign_subset(vrt, si_immature_AZ);
			}

		}

	//	WARNING: use of refine procedure before SelectAreaBoundary doesn't result in the expected way
		//Refine(grid, sel);

	//	Select immature_AZ boundary for extrusion of recycling subset
		tmpSel.clear();
		vExtrusionEdges.clear();
		vTmpExtrusionEdges.clear();
		SelectAreaBoundary(tmpSel, sel.begin<Face>(), sel.end<Face>());

		//	Abstain from refining here in favour of a coarser coarse grid (#1)
		//Refine(grid, sel);

		for(EdgeIterator eIter = tmpSel.begin<Edge>(); eIter != tmpSel.end<Edge>(); ++eIter)
		{
			Edge* e = *eIter;
			vExtrusionEdges.push_back(e);
			vTmpExtrusionEdges.push_back(e);

			sh.assign_subset(e, si_recycling[i]);
			sh.assign_subset(e->vertex(0), si_recycling[i]);
			sh.assign_subset(e->vertex(1), si_recycling[i]);
		}

		vUnitExtrDir = aaNorm[vrt];
		VecScale(vExtrDir, vUnitExtrDir, -1.0 * (TbarHeight+TbarTopHeight));
		//	No double extrusion here in favor of a coarser grid
		// 	Info: factor 0.5 because of following double extrusion!
		//VecScale(vExtrDir, vUnitExtrDir, -1.0 * (TbarHeight+TbarTopHeight));
		//Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);
		Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);

	//	Reassign root extrusion edges to mature_AZ subset
		for(size_t i = 0; i < vTmpExtrusionEdges.size(); ++i)
		{
			Edge* e = vTmpExtrusionEdges[i];

			sh.assign_subset(e, si_immature_AZ);
			sh.assign_subset(e->vertex(0), si_immature_AZ);
			sh.assign_subset(e->vertex(1), si_immature_AZ);
		}

	//	Triangulate top of probe subset
		tmpSel.clear();
		for(size_t i = 0; i < vExtrusionEdges.size(); ++i)
		{
			Edge* e = vExtrusionEdges[i];
			tmpSel.select(e);
		}
		sh.set_default_subset_index(si_recycling[i]);
		TriangleFill_SweepLine(grid, tmpSel.begin<Edge>(), tmpSel.end<Edge>(), aPosition, aInt, &sh, si_recycling[i]);
	}

//	Optimize recycling subset triangulation
	for(size_t i = 0; i < immature_vrts.size(); ++i)
	{
		sel.clear();
		for(FaceIterator fIter = sh.begin<Face>(si_recycling[i]); fIter != sh.end<Face>(si_recycling[i]); ++fIter)
		{
			Face* f = *fIter;
			sel.select(f);
		}
		sh.set_default_subset_index(si_recycling[i]);

		//	Abstain from refining here in favour of a coarser coarse grid (#1)
		Refine(grid, sel);

		QualityGridGeneration(grid, sh.begin<Face>(si_recycling[i]), sh.end<Face>(si_recycling[i]), aaPos, 30.0);
	}


////
//	Optimize triangulation
////
	QualityGridGeneration(grid, sh.begin<Face>(si_immature_AZ), sh.end<Face>(si_immature_AZ), aaPos, 30.0);
	QualityGridGeneration(grid, sh.begin<Face>(si_mature_AZ), sh.end<Face>(si_mature_AZ), aaPos, 30.0);
	QualityGridGeneration(grid, sh.begin<Face>(si_bouton_bnd), sh.end<Face>(si_bouton_bnd), aaPos, 25.0);


////
//	Create extracellular space
////
	if(bExtSpace)
	{
		sh.set_default_subset_index(si_ext_bnd);
		GenerateIcosphere(grid, center, 1.2*radius, numRefinements, aPosition);
	}


////
//	Volume grid generation
////
	sh.set_default_subset_index(-1);
	Tetrahedralize(grid, sh, 18.0, false, false);
	//SeparateSubsetsByLowerDimSubsets<Volume>(grid, sh, true);


////
//	New volume subset separation with SelectRegion using IsNotInSubset callback
////

//	si_mit_int
	SelectRegion<Volume>(sel, center, aaPos, IsNotInSubset(sh, -1));
	for(VolumeIterator vIter = sel.begin<Volume>(); vIter != sel.end<Volume>(); ++vIter)
	{
		Volume* v = *vIter;
		sh.assign_subset(v, si_mit_int);
	}

//	si_cyt_int
	number x = (radius + mit_radius) / 2;
	vector3 coord_in_cyt_int(x, 0, 0);
	SelectRegion<Volume>(sel, coord_in_cyt_int, aaPos, IsNotInSubset(sh, -1));
	for(VolumeIterator vIter = sel.begin<Volume>(); vIter != sel.end<Volume>(); ++vIter)
	{
		Volume* v = *vIter;
		sh.assign_subset(v, si_cyt_int);
	}

//	si_Tbars_int
	vector3 coord;
	for(size_t i = 0; i < vrts.size(); i++)
	{
		sel.clear();
		VecScale(coord, aaNorm[vrts[i]], -0.5*TbarHeight);
		VecAdd(coord, aaPos[vrts[i]], coord);
		SelectRegion<Volume>(sel, coord, aaPos, IsNotInSubset(sh, -1));

		for(VolumeIterator vIter = sel.begin<Volume>(); vIter != sel.end<Volume>(); ++vIter)
		{
			Volume* v = *vIter;
			sh.assign_subset(v, si_Tbars_int);
		}
	}

//	si_ext_int
	if(bExtSpace)
	{
		x = (radius + 1.2*radius) / 2;
		vector3 coord_in_ext_int(x, 0, 0);
		SelectRegion<Volume>(sel, coord_in_ext_int, aaPos, IsNotInSubset(sh, -1));
		for(VolumeIterator vIter = sel.begin<Volume>(); vIter != sel.end<Volume>(); ++vIter)
		{
			Volume* v = *vIter;
			sh.assign_subset(v, si_ext_int);
		}
	}

//	si_recycling
	for(size_t i = 0; i < immature_vrts.size(); ++i)
	{
		sel.clear();
		VecScale(coord, aaNorm[immature_vrts[i]], -0.5*TbarHeight);
		VecAdd(coord, aaPos[immature_vrts[i]], coord);
		SelectRegion<Volume>(sel, coord, aaPos, IsNotInSubset(sh, -1));

		for(VolumeIterator vIter = sel.begin<Volume>(); vIter != sel.end<Volume>(); ++vIter)
		{
			Volume* v = *vIter;
			sh.assign_subset(v, si_recycling[i]);
		}
	}

//	si_probe
	for(size_t i = 0; i < vrts.size(); ++i)
	{
		sel.clear();
		VecScale(coord, aaNorm[vrts[i]], -1.0*(0.1+TbarHeight+TbarTopHeight)*0.5);
		VecAdd(coord, aaPos[vrts[i]], coord);
		SelectRegion<Volume>(sel, coord, aaPos, IsNotInSubset(sh, -1));

		for(VolumeIterator vIter = sel.begin<Volume>(); vIter != sel.end<Volume>(); ++vIter)
		{
			Volume* v = *vIter;
			sh.assign_subset(v, si_probe[i]);
		}
	}


//	Final subset management (Instead of AdjustSubsetsForSimulation we use: )
	//AdjustSubsetsForSimulation(sh, true);

	//	INFO:
	//	faces subsets indices have priority over volume subset indices.
	//	Thus copy faces subset indices to sides first and perform the copy-mechanism
	//	on the whole grid afterwards
	//	(note that we only copy subset-indices to unassigned sides)

	sh.assign_subset(grid.vertices_begin(), grid.vertices_end(), -1);
	sh.assign_subset(grid.edges_begin(), grid.edges_end(), -1);

	for(int i = 0; i < sh.num_subsets(); ++i){
		CopySubsetIndicesToSides(sh, sh.begin<Face>(i), sh.end<Face>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Edge>(i), sh.end<Edge>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Vertex>(i), sh.end<Vertex>(i), true);
	}

	for(int i = 0; i < sh.num_subsets(); ++i){
		CopySubsetIndicesToSides(sh, sh.begin<Volume>(i), sh.end<Volume>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Face>(i), sh.end<Face>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Edge>(i), sh.end<Edge>(i), true);
		CopySubsetIndicesToSides(sh, sh.begin<Vertex>(i), sh.end<Vertex>(i), true);
	}

	AssignSubsetColors(sh);


//	Reassign mature_AZ boundary edges to mature_AZ subset
	sel.clear();
	SelectAreaBoundary(sel, sh.begin<Face>(si_mature_AZ), sh.end<Face>(si_mature_AZ));
	for(EdgeIterator eIter = sel.begin<Edge>(); eIter != sel.end<Edge>(); ++eIter)
	{
		Edge* e = *eIter;
		sh.assign_subset(e, si_mature_AZ);
		sh.assign_subset(e->vertex(0), si_mature_AZ);
		sh.assign_subset(e->vertex(1), si_mature_AZ);
	}

//	Reassign immature_AZ boundary edges to immature_AZ subset
	sel.clear();
	SelectAreaBoundary(sel, sh.begin<Face>(si_immature_AZ), sh.end<Face>(si_immature_AZ));
	for(EdgeIterator eIter = sel.begin<Edge>(); eIter != sel.end<Edge>(); ++eIter)
	{
		Edge* e = *eIter;
		sh.assign_subset(e, si_immature_AZ);
		sh.assign_subset(e->vertex(0), si_immature_AZ);
		sh.assign_subset(e->vertex(1), si_immature_AZ);
	}

//	Reassign CChCl boundary edges to CChCl subset
	sel.clear();
	SelectAreaBoundary(sel, sh.begin<Face>(si_CChCl), sh.end<Face>(si_CChCl));
	for(EdgeIterator eIter = sel.begin<Edge>(); eIter != sel.end<Edge>(); ++eIter)
	{
		Edge* e = *eIter;
		sh.assign_subset(e, si_CChCl);
		sh.assign_subset(e->vertex(0), si_CChCl);
		sh.assign_subset(e->vertex(1), si_CChCl);
	}


//	Name subsets
	sh.set_subset_name("bouton_bnd", 	si_bouton_bnd);
	sh.set_subset_name("mature_AZ",	 	si_mature_AZ);
	sh.set_subset_name("immature_AZ", 	si_immature_AZ);
	sh.set_subset_name("CChCl", 		si_CChCl);
	sh.set_subset_name("T-bars_bnd", 	si_Tbars_bnd);
	sh.set_subset_name("mit_bnd", 		si_mit_bnd);
	sh.set_subset_name("cyt",	 		si_cyt_int);
	sh.set_subset_name("T-bars_int", 	si_Tbars_int);
	sh.set_subset_name("mit_int", 		si_mit_int);
	sh.set_subset_name("locks", 		si_locks);

	if(bExtSpace)
	{
		sh.set_subset_name("ext_bnd", 		si_ext_bnd);
		sh.set_subset_name("ext_int", 		si_ext_int);
	}

	stringstream probeStream;
	stringstream recyclingStream;

	for(int i = 0; i < numReleaseSites/2; ++i)
	{
		probeStream << "probe_" << i;
		recyclingStream << "recycling_" << i;

		sh.set_subset_name(probeStream.str().c_str(), 		si_probe[i]);
		sh.set_subset_name(recyclingStream.str().c_str(),	si_recycling[i]);

	//	Clear stringstreams
		probeStream.str(std::string());
		recyclingStream.str(std::string());
	}

//	For now: erase Tbar and mitochondrial volume subsets
	grid.erase(sh.begin<Vertex>(si_mit_int), sh.end<Vertex>(si_mit_int));
	grid.erase(sh.begin<Edge>(si_mit_int), sh.end<Edge>(si_mit_int));
	grid.erase(sh.begin<Face>(si_mit_int), sh.end<Face>(si_mit_int));
	grid.erase(sh.begin<Volume>(si_mit_int), sh.end<Volume>(si_mit_int));

	grid.erase(sh.begin<Vertex>(si_Tbars_int), sh.end<Vertex>(si_Tbars_int));
	grid.erase(sh.begin<Edge>(si_Tbars_int), sh.end<Edge>(si_Tbars_int));
	grid.erase(sh.begin<Face>(si_Tbars_int), sh.end<Face>(si_Tbars_int));
	grid.erase(sh.begin<Volume>(si_Tbars_int), sh.end<Volume>(si_Tbars_int));

//	Erase subset no. 8 two times
	sh.erase_subset(si_Tbars_int);
	sh.erase_subset(si_Tbars_int);


//	Vertex* vrt;
//	Grid::edge_traits::secure_container edges;
//	grid.associated_elements(edges, vrt);

//	Write file
	stringstream ss;
	//ss << fileName << "dnmj_bouton" << "_" << numReleaseSites << "AZ" << ".ugx";
	ss << fileName;
	string outfile = ss.str();

	SaveGridToFile(grid, sh, outfile.c_str());
}


////////////////////////////////////////////////////////////////////////////////////////////
//	BuildTbar
////////////////////////////////////////////////////////////////////////////////////////////
void BuildTbar(	Grid& grid, SubsetHandler& sh_orig, Vertex* vrt,
				Grid::VertexAttachmentAccessor<APosition>& aaPos,
				Grid::VertexAttachmentAccessor<ANormal>& aaNorm,
				int si,
				number TbarHeight,
				number TbarLegRadius,
				number TbarTopRadius,
				number TbarTopHeight)
{
	SubsetHandler sh(grid);
	AInt aInt;
	grid.attach_to_vertices(aInt);


//	Temporal subset index specs
	const int si_Tbars_bnd 		= 0;
	const int si_Tbars_post 	= 1;
	const int si_Tbars_bottom 	= 2;
	const int si_Tbars_sides 	= 3;
	const int si_Tbars_tabletop	= 4;


//	Instantiate temporal subset handler sh
	for(VertexIterator vIter = sh_orig.begin<Vertex>(si); vIter != sh_orig.end<Vertex>(si); ++vIter)
	{
		Vertex* v = *vIter;
		sh.assign_subset(v, 0);
	}

	for(EdgeIterator eIter = sh_orig.begin<Edge>(si); eIter != sh_orig.end<Edge>(si); ++eIter)
	{
		Edge* e = *eIter;
		sh.assign_subset(e, 0);
	}

	for(FaceIterator fIter = sh_orig.begin<Face>(si); fIter != sh_orig.end<Face>(si); ++fIter)
	{
		Face* f = *fIter;
		sh.assign_subset(f, 0);
	}


//	setup extrusion tools
	vector<Edge*> vExtrusionEdges;
	vector<Edge*> vTmpExtrusionEdges;
	vector3 vZero(0.0,0.0,0.0);

	vector3 vExtrDir;
	vector3 vUnitExtrDir;

	Selector sel(grid);
	Selector tmpSel(grid);

	sel.clear();
	tmpSel.clear();
	vExtrusionEdges.clear();
	vTmpExtrusionEdges.clear();


//	create cylinder surface around specified vertex
	vUnitExtrDir = aaNorm[vrt];
	VecScale(vExtrDir, vUnitExtrDir, -1.0 * TbarHeight);
	AdaptSurfaceGridToCylinder(sel, grid, vrt, vExtrDir, TbarLegRadius, 0.01);


//	assign closure of cylinder surface to subset si_Tbars_post
	for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
	{
		Face* f = *fIter;
		sh.assign_subset(f, si_Tbars_post);

		for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f); eIter != grid.associated_edges_end(f); ++eIter)
		{
			Edge* e = *eIter;
			sh.assign_subset(e, si_Tbars_post);
		}

		for(size_t i = 0; i < f->num_vertices(); ++i)
		{
			Vertex* vrt = f->vertex(i);
			sh.assign_subset(vrt, si_Tbars_post);
		}
	}


//	select si_Tbars_post boundary edges and extrude them and assign extrusion to subset si_Tbars_bnd
	SelectAreaBoundary(tmpSel, sel.begin<Face>(), sel.end<Face>());

	for(EdgeIterator eIter = tmpSel.begin<Edge>(); eIter != tmpSel.end<Edge>(); ++eIter)
	{
		Edge* e = *eIter;
		vExtrusionEdges.push_back(e);
		vTmpExtrusionEdges.push_back(e);

		sh.assign_subset(e, si_Tbars_bnd);
		sh.assign_subset(e->vertex(0), si_Tbars_bnd);
		sh.assign_subset(e->vertex(1), si_Tbars_bnd);
	}

	Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);
	Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);
	Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);
	Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);


//	assign upper ring of T-bar posts to subset si_Tbars_bottom
	vTmpExtrusionEdges.clear();
	for(size_t i = 0; i < vExtrusionEdges.size(); i++)
	{
		Edge* e = vExtrusionEdges[i];
		sh.assign_subset(e, si_Tbars_bottom);

		vTmpExtrusionEdges.push_back(vExtrusionEdges[i]);
	}


//	extrude and scale table top horizontally
	tmpSel.clear();
	Extrude(grid, NULL, &vExtrusionEdges, NULL, vZero, EO_CREATE_FACES);

	for(size_t i = 0; i < vExtrusionEdges.size(); i++)
	{
		Edge* e = vExtrusionEdges[i];
		tmpSel.select(e->vertex(0));
		tmpSel.select(e->vertex(1));
		tmpSel.select(e);

		sh.assign_subset(e, si_Tbars_bottom);
		sh.assign_subset(e->vertex(0), si_Tbars_bottom);
		sh.assign_subset(e->vertex(1), si_Tbars_bottom);
	}


//	seperate post from table bottom (reassign upper ring of T-bar table post to subset si_Tbars_bnd)
	for(size_t i = 0; i < vTmpExtrusionEdges.size(); i++)
	{
		Edge* e = vTmpExtrusionEdges[i];
		sh.assign_subset(e, si_Tbars_bnd);
		sh.assign_subset(e->vertex(0), si_Tbars_bnd);
		sh.assign_subset(e->vertex(1), si_Tbars_bnd);
	}

	vector3 baryCenter = CalculateBarycenter(tmpSel.begin<Vertex>(), tmpSel.end<Vertex>(), aaPos);

	for(VertexIterator vIter = tmpSel.begin<Vertex>(); vIter != tmpSel.end<Vertex>(); ++vIter)
	{
		Vertex* vrt = *vIter;
		VecSubtract(vExtrDir, aaPos[vrt], baryCenter);
		VecNormalize(vExtrDir, vExtrDir);
		VecScale(vExtrDir, vExtrDir, TbarTopRadius);
		VecAdd(aaPos[vrt], aaPos[vrt], vExtrDir);
	}


//	select edges to extrude vertically and assign table top side edges to subset si_Tbars_sides
	tmpSel.clear();
	vExtrusionEdges.clear();
	SelectAreaBoundary(tmpSel, sh.begin<Face>(si_Tbars_bottom), sh.end<Face>(si_Tbars_bottom));
	for(EdgeIterator eIter = tmpSel.begin<Edge>(); eIter != tmpSel.end<Edge>(); ++eIter)
	{
		Edge* e = *eIter;
		if(NumAssociatedFaces(grid, e) == 1)
		{
			vExtrusionEdges.push_back(e);
			sh.assign_subset(e, si_Tbars_sides);
			sh.assign_subset(e->vertex(0), si_Tbars_sides);
			sh.assign_subset(e->vertex(1), si_Tbars_sides);
		}
	}


//	extrude table top vertically and assign upper edge ring to subset si_Tbars_tabletop
	VecScale(vExtrDir, vUnitExtrDir, -1.0 * TbarTopHeight);
	Extrude(grid, NULL, &vExtrusionEdges, NULL, vExtrDir, EO_CREATE_FACES);

	tmpSel.clear();
	for(size_t i = 0; i < vExtrusionEdges.size(); i++)
	{
		Edge* e = vExtrusionEdges[i];
		tmpSel.select(e);
		sh.assign_subset(e, si_Tbars_tabletop);
	}

	sh.set_default_subset_index(si_Tbars_tabletop);
	TriangleFill_SweepLine(grid, tmpSel.begin<Edge>(), tmpSel.end<Edge>(), aPosition, aInt, &sh, si_Tbars_tabletop);

	sh.set_default_subset_index(-1);


//	optimize table top
	tmpSel.clear();
	for(FaceIterator fIter = sh.begin<Face>(si_Tbars_tabletop); fIter != sh.end<Face>(si_Tbars_tabletop); ++fIter)
	{
		Face* f = *fIter;
		tmpSel.select(f);
	}
	sh.set_default_subset_index(si_Tbars_tabletop);
	Refine(grid, tmpSel);
	QualityGridGeneration(grid, sh.begin<Face>(si_Tbars_tabletop), sh.end<Face>(si_Tbars_tabletop), aaPos, 30.0);


//	optimize table bottom
	tmpSel.clear();
	Triangulate(grid, sh.begin<Quadrilateral>(si_Tbars_bottom), sh.end<Quadrilateral>(si_Tbars_bottom));

	for(FaceIterator fIter = sh.begin<Face>(si_Tbars_bottom); fIter != sh.end<Face>(si_Tbars_bottom); ++fIter)
	{
		Face* f = *fIter;
		tmpSel.select(f);
	}
	sh.set_default_subset_index(6);
	QualityGridGeneration(grid, sh.begin<Face>(si_Tbars_bottom), sh.end<Face>(si_Tbars_bottom), aaPos, 20.0);


//	optimize table sides
	sel.clear();
	sel.enable_autoselection(true);
	tmpSel.clear();

	ABool aBoolMarked = false;
	grid.attach_to_faces(aBoolMarked);
	Grid::FaceAttachmentAccessor<ABool> aaBoolMarked(grid, aBoolMarked);

	int counter = 1;

	sh.set_default_subset_index(si_Tbars_sides);
	Triangulate(grid, sh.begin<Quadrilateral>(si_Tbars_sides), sh.end<Quadrilateral>(si_Tbars_sides));

	for(FaceIterator fIter = sh.begin<Face>(si_Tbars_sides); fIter != sh.end<Face>(si_Tbars_sides); ++fIter)
	{
		Face* f = *fIter;
		sel.select(f);
	}

	for(FaceIterator fIter = sel.begin<Face>(); fIter != sel.end<Face>(); ++fIter)
	{
		Face* f = *fIter;
		if(aaBoolMarked[f] == false)
		{
			tmpSel.select(f);
			SelectLinkedFlatFaces(tmpSel, 1.0);

			for(FaceIterator fMarkedIter = tmpSel.begin<Face>(); fMarkedIter != tmpSel.end<Face>(); ++fMarkedIter)
			{
				aaBoolMarked[*fMarkedIter] = true;
				sh.assign_subset(*fMarkedIter, counter + 4);
			}

			tmpSel.clear();
			counter++;
		}
	}

	for(int i = 0; i < counter; i++)
		QualityGridGeneration(grid, sh.begin<Face>(i + 4), sh.end<Face>(i + 4), aaPos, 30.0);


//	optimize table top and bottom again
	tmpSel.clear();
	SelectAreaBoundary(tmpSel, sh.begin<Face>(si_Tbars_tabletop), sh.end<Face>(si_Tbars_tabletop));
	SelectAreaBoundary(tmpSel, sh.begin<Face>(si_Tbars_bottom), sh.end<Face>(si_Tbars_bottom));
	QualityGridGeneration(grid, sh.begin<Face>(si_Tbars_tabletop), sh.end<Face>(si_Tbars_tabletop), aaPos, 25.0, IsSelected(tmpSel));
	QualityGridGeneration(grid, sh.begin<Face>(si_Tbars_bottom), sh.end<Face>(si_Tbars_bottom), aaPos, 30.0, IsSelected(tmpSel));


//	force triangulation of table post
	Triangulate(grid, sh.begin<Quadrilateral>(si_Tbars_bnd), sh.end<Quadrilateral>(si_Tbars_bnd));


//	join temporal t-bar subsets created during construction
	while(sh.num_subsets() > 1)
	{
		sh.join_subsets(0, 0, 1, true);
	}


//	transfer new t-bar elements to original subset handler sh_orig
	for(VertexIterator vIter = sh.begin<Vertex>(si_Tbars_bnd); vIter != sh.end<Vertex>(si_Tbars_bnd); ++vIter)
	{
		Vertex* v = *vIter;
		sh_orig.assign_subset(v, si);
	}

	for(EdgeIterator eIter = sh.begin<Edge>(si_Tbars_bnd); eIter != sh.end<Edge>(si_Tbars_bnd); ++eIter)
	{
		Edge* e = *eIter;
		sh_orig.assign_subset(e, si);
	}

	for(FaceIterator fIter = sh.begin<Face>(si_Tbars_bnd); fIter != sh.end<Face>(si_Tbars_bnd); ++fIter)
	{
		Face* f = *fIter;
		sh_orig.assign_subset(f, si);
	}
}









////////////////////////////////////////////////////////////////////////////////////////////
//	SaveSelectionStatesToFile
////////////////////////////////////////////////////////////////////////////////////////////
void SaveSelectionStatesToFile(Grid& mg, Selector& msel, const char* filename)
{
//	create a subset handler which holds different subsets for the different selection states
	//MultiGrid& mg = *msel.multi_grid();
	SubsetHandler sh(mg);


	for(Selector::traits<Volume>::iterator iter = msel.begin<Volume>();
		iter != msel.end<Volume>(); ++iter)
	{
		sh.assign_subset(*iter, msel.get_selection_status(*iter));
	}

	for(Selector::traits<Face>::iterator iter = msel.begin<Face>();
		iter != msel.end<Face>(); ++iter)
	{
		sh.assign_subset(*iter, msel.get_selection_status(*iter));
	}

	for(Selector::traits<Edge>::iterator iter = msel.begin<Edge>();
		iter != msel.end<Edge>(); ++iter)
	{
		sh.assign_subset(*iter, msel.get_selection_status(*iter));
	}

	for(Selector::traits<Vertex>::iterator iter = msel.begin<Vertex>();
		iter != msel.end<Vertex>(); ++iter)
	{
		sh.assign_subset(*iter, msel.get_selection_status(*iter));
	}


	const char* subsetNames[] = {"one", "two"};

	for(int i = 0; i < 2; ++i)
		sh.subset_info(i).name = subsetNames[i];

	AssignSubsetColors(sh);
	EraseEmptySubsets(sh);
	//SaveGridHierarchyTransformed(mg, sh, filename, LG_DISTRIBUTION_Z_OUTPUT_TRANSFORM);

	SaveGridToFile(mg, sh, filename);
}




} // neuro_collection
} // namespace ug




