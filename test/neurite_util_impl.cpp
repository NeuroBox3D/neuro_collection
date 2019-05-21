/*!
 * \file neurite_util_impl.cpp
 *
 *  Created on: Apr 22, 2019
 *      Author: stephan
 */

#include "neurite_util.h"
#include "common/log.h"
#include "common/error.h"
#include "bridge/domain_bridges/selection_bridge.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h"
#include "lib_grid/refinement/regular_refinement.h"
#include <cmath>
#include <boost/lexical_cast.hpp>

namespace ug {
	namespace neuro_collection {
		////////////////////////////////////////////////////////////////////////
		/// calculate_angle
		////////////////////////////////////////////////////////////////////////
		float calculate_angle
		(
			const vector3& pos,
			const vector3& origin,
			const vector3& point,
			const vector3& n
		) {
			vector3 v1, v2;
			VecSubtract(v1, point, pos);
			VecSubtract(v2, origin, pos);
			number dot = VecDot(v1, v2);
			vector3 cross, normal;
			VecCross(cross, v1, v2);
			VecNormalize(normal, n);
			number det = VecDot(normal, cross);
			return rad_to_deg(std::atan2(det, dot));
		}

		////////////////////////////////////////////////////////////////////////
		/// calculate_angle
		////////////////////////////////////////////////////////////////////////
		float calculate_angle
		(
			const vector3& pos,
			const vector3& origin,
			const vector3& point
		) {
			vector3 v1, v2;
			VecSubtract(v1, point, pos);
			VecSubtract(v2, origin, pos);
			number dot = VecDot(v1, v2);
			number v1_len = VecLengthSq(v1);
			number v2_len = VecLengthSq(v2);
			number x = (dot / (std::sqrt(v1_len*v2_len)));
			UG_LOGN("acos(" << x << ")");
			return rad_to_deg(acos(x));
		}

		////////////////////////////////////////////////////////////////////////
		/// deg_to_full_range
		////////////////////////////////////////////////////////////////////////
		float deg_to_full_range
		(
			float angle
		) {
			return fmod(angle+360, 360);
		}

		////////////////////////////////////////////////////////////////////////
		/// calculate_angles
		////////////////////////////////////////////////////////////////////////
		void calculate_angles
		(
			const vector3& originCenter,
			const std::vector<ug::vector3>& points,
			std::vector<number>& angles,
			size_t refIndex
		) {
			/// reference point
			ug::vector3	ref = points[(refIndex > points.size() ? 0 : refIndex)];

			/// calculate each angle
			for (std::vector<ug::vector3>::const_iterator it = points.begin(); it != points.end(); ++it) {
				angles.push_back(calculate_angle(originCenter, ref, *it));
			}
		}


		////////////////////////////////////////////////////////////////////////
		/// calculate_angles
		////////////////////////////////////////////////////////////////////////
		void calculate_angles
		(
			const vector3& originCenter,
			const std::vector<ug::vector3>& points,
			std::vector<number>& angles,
			const std::vector<ug::vector3>& normals,
			const ug::vector3& ref
		) {
			size_t i = 0;
			/// calculate each angle
			for (std::vector<ug::vector3>::const_iterator it = points.begin(); it != points.end(); ++it) {
				angles.push_back(calculate_angle(originCenter, ref, *it, normals[i]));
				i++;
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// calculate_angles
		////////////////////////////////////////////////////////////////////////
		void calculate_angles
		(
			const vector3& originCenter,
			const std::vector<ug::vector3>& points,
			std::vector<number>& angles,
			const std::vector<ug::vector3>& normals,
			size_t refIndex
		) {
			/// reference point
			ug::vector3	ref = points[(refIndex > points.size() ? 0 : refIndex)];

			size_t i = 0;
			/// calculate each angle
			for (std::vector<ug::vector3>::const_iterator it = points.begin(); it != points.end(); ++it) {
				angles.push_back(calculate_angle(originCenter, ref, *it, normals[i]));
				i++;
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// sortIndices
		////////////////////////////////////////////////////////////////////////
		void sortIndices
		(
			const std::vector<number>& values,
			std::vector<size_t>& indices
		) {
			std::vector<std::pair<size_t, size_t> > a;
			for (size_t i = 0 ; i < values.size() ; i++) {
				a.push_back(std::make_pair(values[i], i));
			}

			std::sort(a.begin(), a.end());
			for (size_t i = 0; i < a.size(); i++) {
				indices.push_back(a[i].second);
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// sortbysec
		////////////////////////////////////////////////////////////////////////
		bool sortbysec
		(
			const std::pair<ug::Vertex*, number> &a,
	        const std::pair<ug::Vertex*, number> &b
        ) {
			return (a.second < b.second);
		}

		////////////////////////////////////////////////////////////////////////
		/// deg360
		////////////////////////////////////////////////////////////////////////
		number deg360
		(
			number a
		) {
			if (a > 180) {
				a -= 360;
			}

			if (a < -180) {
				a += 360;
			}

			if (a < 0) {
				a = a + 360;
			}

			return a;
		}

		////////////////////////////////////////////////////////////////////////
		/// connect_outer_and_inner_root_neurites_to_outer_soma
		////////////////////////////////////////////////////////////////////////
		void connect_outer_and_inner_root_neurites_to_outer_soma
		(
			size_t somaIndex,
			size_t numQuads,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			std::vector<ug::Vertex*>& rootNeurites, /// Note: is just one single array => need to iterate by stride 4
			std::vector<ug::Vertex*>& rootNeuritesInner /// Note: is just one single array => need to iterate by stride 4
		) {
			Selector::traits<Edge>::iterator eit;
			Selector::traits<Edge>::iterator eit_end;
			Selector sel(g);

			/// TODO Project root neurite vertices to outer soma surface
			/// Find minimal angle vertices... and merge them at outer soma surface
			/// (don't need to introduce factor 1.05 for outer soma, can let at 1.00)
			std::vector<std::vector<ug::Vertex*> > projectedVertices;
			std::vector<std::vector<ug::Vertex*> > projectedVertices2;
			std::vector<std::vector<ug::vector3> > projected;
			std::vector<std::vector<ug::vector3> > allNormals;
			projected.resize(numQuads);
			projectedVertices.resize(numQuads);
			projectedVertices2.resize(numQuads);
			for (size_t i = 1; i < numQuads+1; i++) {
				UG_LOGN("Selecting now subset: " << somaIndex+i);
				sel.clear();
				/// Select outer soma inner quad
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+i, true);
				ug::Vertex* v0 = *(sel.vertices_begin());
				UG_LOGN("First vertex of subset: " << aaPos[v0]);
				SelectSubsetElements<Edge>(sel, sh, somaIndex+i, true);
				eit = sel.edges_begin();
				eit_end = sel.edges_end();
				std::vector<std::pair<ug::Vertex*, ug::Vertex*> > es;
				size_t count = 0;

				UG_LOGN("Trying to find edges...");
				for (; eit != eit_end; ++eit) {
					Edge* e = *eit;
					UG_LOGN("trying e->vertex(0)");
					std::pair<ug::Vertex*, ug::Vertex*> p;
					if (e->vertex(0) == v0) {
						p.first = v0;
						p.second = e->vertex(1);
						count++;
						es.push_back(p);
					}

					UG_LOGN("trying e->vertex(1)");
					if (e->vertex(1) == v0) {
						p.first = e->vertex(0);
						p.second = v0;
						es.push_back(p);
						count++;
					}
				}
				UG_LOGN("Found edges");
				UG_COND_THROW(count != 2, "Number of edges has to be two!");

				ug::vector3 v1, v2;
				VecSubtract(v1, aaPos[es[0].first], aaPos[es[0].second]);
				UG_LOGN("Subtracted v1");
				VecSubtract(v2, aaPos[es[1].first], aaPos[es[1].second]);
				UG_LOGN("Subtracted v2");
				ug::vector3 normal;
				VecCross(normal, v1, v2);
				VecNormalize(normal, normal);
				sel.clear();
				UG_LOGN("calculated cross product");

				size_t numVerts = 4;
				ug::vector3 vProjected;
				std::vector<ug::vector3> normals;
				for (size_t l = 0; l < numVerts; l++) {
					ug::Vertex* vit = rootNeurites[(i-1)*numVerts+l];
					ug::vector3 v;
					VecSubtract(v, aaPos[vit], v1);
					number dot = VecDot(v, normal);
					vProjected = aaPos[vit];
					VecScaleAdd(vProjected, 1.0, vProjected, dot, normal);
					ug::Vertex* projVert = *g.create<ug::RegularVertex>();
					const ug::vector3& n = -normal;
					ProjectPointToPlane(vProjected, aaPos[vit], aaPos[es[0].first], n);
					aaPos[projVert] = vProjected;
					projected[i-1].push_back(vProjected);
					projectedVertices[i-1].push_back(vit); /// save original vertex from which we projected (root neurite)
					projectedVertices2[i-1].push_back(projVert); /// the actual projected vertex on the soma surface
					normals.push_back(normal);
				}
				allNormals.push_back(normals);
				UG_LOGN("First projection!");
			}

			/// Find the corresponding pairs of projected and unprojected respectively original quad vertex and unprojected vertex
			std::vector<std::vector<number> > allAngles;
			std::vector<std::vector<number> > allAnglesInner;

			Selector::traits<Vertex>::iterator vit;
			Selector::traits<Vertex>::iterator vit_end;
			size_t j = 1;
			/// TODO: Projected vertices have to be centered around the outer soma's quads (ER and PM)
			for (std::vector<std::vector<ug::vector3> >::const_iterator it = projected.begin(); it != projected.end(); ++it) {
				ug::vector3 centerOut;
				std::vector<number> angles;
				CalculateCenter(centerOut, &(*it)[0], it->size());
				calculate_angles(centerOut, *it, angles, allNormals[j-1], (*it)[0]);
				allAngles.push_back(angles);

				sel.clear();
				std::vector<number> angles2;
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+j, true);
				vit = sel.vertices_begin();
				vit_end = sel.vertices_end();
				std::vector<ug::vector3> verts;
				for (; vit != vit_end; ++vit) {
					verts.push_back(aaPos[*vit]);
				}
				calculate_angles(centerOut, verts, angles2, allNormals[j-1], (*it)[0]);
				allAnglesInner.push_back(angles2);
				j++;
			}
			UG_LOGN("Found angles");

			for (size_t i = 0; i < 1; i++) {
				for (size_t j = 0; j < allAngles[i].size(); j++) {
					UG_LOGN("(old angles): " << allAngles[i][j]); /// allAngles are projected inner vertices (from outer vertices)
				}
			}

			for (size_t i = 0; i < 1; i++) {
				for (size_t j = 0; j < allAnglesInner[i].size(); j++) {
					UG_LOGN("(old angles inner): " << allAnglesInner[i][j]); /// allAnglesInner are unprojected inner vertices
				}
			}

			for (size_t i = 0; i < allAngles.size(); i++) {
				for (size_t j = 0; j < allAngles[i].size(); j++) {
					number a = allAngles[i][j];
					number b = allAnglesInner[i][j];
					/// Convert from -180, 180 to 0,360 interval -> then sort array based on this
					if (a > 180) {
						a -= 360;
					}
					if (a < -180) {
						a += 360;
					}

					if (a < 0) {
						a = a + 360;
					}

					if (b > 180) {
						b -= 360;
					}
					if (b < -180) {
						b += 360;
					}

					if (b < 0) {
						b = b + 360;
					}

					allAngles[i][j] = a;
					allAnglesInner[i][j] = b;
				}
			}

			for (size_t i = 0; i < 1; i++) {
				for (size_t j = 0; j < allAngles[i].size(); j++) {
					UG_LOGN("(new angles): " << allAngles[i][j]); /// projected inner vertices (from outer vertices of root neurite)
				}
			}

			for (size_t i = 0; i < 1; i++) {
				for (size_t j = 0; j < allAnglesInner[i].size(); j++) {
					UG_LOGN("(new angles inner): " << allAnglesInner[i][j]); /// unprojected inner vertices
				}
			}

			/// get mapping of outer vertices to inner (Unprojected) vertices
			std::map<Vertex*, Vertex*> innerToOuter;
			std::map<Vertex*, Vertex*> innerToOuter2;
			for (size_t i = 0 ; i < projected.size(); i++) {
				/// calculate angles for inner original soma vertices and projected vertices (projectedVertices is ug::Vertex*)
				std::vector<std::pair<ug::Vertex*, number> > anglesOfProjectedInnerVertices;
				std::vector<std::pair<ug::Vertex*, number> > anglesOfOrginalSomaInnerVertices;

				/// TODO: Convert angles from -180, 180 to 0, 360 degree interval with deg360 -> check if this method works.
				/// these are the inner projected vertices
				ug::vector3 centerOut;
				std::vector<number> angles;
				CalculateCenter(centerOut, &projected[i][0], 4);
				calculate_angles(centerOut, projected[i], angles, allNormals[i], projected[i][0]);
				for (size_t n = 0; n < 4; n++) {
					anglesOfProjectedInnerVertices.push_back(std::make_pair(projectedVertices[i][n], deg360(angles[n])));
				}

				/// these are the inner original vertices
				sel.clear();
				std::vector<number> angles2;
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+i+1, true);
				vit = sel.vertices_begin();
				vit_end = sel.vertices_end();
				std::vector<ug::vector3> verts;
				std::vector<ug::Vertex*> vertsVtx;
				for (; vit != vit_end; ++vit) {
					verts.push_back(aaPos[*vit]);
					vertsVtx.push_back(*vit);
				}
				calculate_angles(centerOut, verts, angles2, allNormals[i], projected[i][0]);
				for (size_t n = 0; n < 4; n++) {
					anglesOfOrginalSomaInnerVertices.push_back(std::make_pair(vertsVtx[n], deg360(angles2[n])));
				}

				/// Store: A mapping of the projected vertices on soma surface to unprojected vertices (root neurite vertices)
				std::map<ug::Vertex*, ug::Vertex*> projectedToUnprojected;
				std::map<ug::Vertex*, ug::Vertex*> projectedToUnprojectedInner;
				for (size_t l = 0; l < 4; l++) {
					projectedToUnprojected[projectedVertices[i][l]] = rootNeurites[(i)*4+l];
					projectedToUnprojectedInner[projectedVertices[i][l]] = rootNeuritesInner[(i)*4+l];
					UG_COND_THROW(!rootNeuritesInner[(i)*4+l], "Root neurite inner vertex not available");
				}

				std::sort(anglesOfOrginalSomaInnerVertices.begin(), anglesOfOrginalSomaInnerVertices.end(), sortbysec);
				std::sort(anglesOfProjectedInnerVertices.begin(), anglesOfProjectedInnerVertices.end(), sortbysec);

				for (size_t l = 0; l < anglesOfOrginalSomaInnerVertices.size(); l++) {
					Vertex* p1 = anglesOfOrginalSomaInnerVertices[l].first; // original soma vertex
					Vertex* temp = anglesOfProjectedInnerVertices[l].first; // projected vertex on somasurface
					Vertex* p2 = projectedToUnprojected[temp]; // look up the unprojected vertex (root neurite vertex)

					UG_COND_THROW(!p1, "P1 not found!");
					UG_COND_THROW(!p2, "P2 not found!");
					innerToOuter[p1] = p2;

					/// find soma inner quads to root neurites inner quads for the first neurite i
					sel.clear();
					SelectSubsetElements<ug::Vertex>(sel, sh, somaIndex+1+numQuads+i, true);
					vit = sel.begin<Vertex>();
					vit_end = sel.end<Vertex>();
					std::vector<Vertex*> array;
					for (; vit != vit_end; ++vit) { array.push_back(*vit); }
					int index = FindClosestVertexInArray(array, aaPos[p1], aaPos, 10); /// closest soma vertex of inner quad to outer soma quad: this is safe if inner and outer quads have same number of vertices

					UG_COND_THROW(!array[index], "Not found!");
					UG_COND_THROW(!projectedToUnprojectedInner[p2], "Not found!");
					innerToOuter2[array[index]] = projectedToUnprojectedInner[projectedToUnprojected[temp]]; /// soma inner quad vertex => rootneurite inner vertex
				}
			}

			/// connect outer quad with outer neurite
			for (std::map<ug::Vertex*, ug::Vertex*>::iterator it = innerToOuter.begin(); it != innerToOuter.end(); ++it) {
				/// outer quads to outer neurite
				ug::Vertex* p1 = it->first; /// soma vertex
				ug::Vertex* p2 = it->second; /// neurite vertex

				/// Merge
				// MergeVertices(g, p1, p2);
				// aaPos[p2] = aaPos[p1]; /// p1 = p2 => merge at neurite, p2 = p1 => merge at soma surface

				/// TODO: Beides mal ist die zuordnung der inneren vertices nicht gesichert
				/// wenn man die angles nicht für beide quads inner und außen mit gleicher referenz berechnet,
				/// bzgl. EINES mittelpunktes und EINER normalen, dann muss es gecentered werden um diesen
				/// mittelpunkt: einfacher oben die zuordnung treffen und diese methode nur 1 mal nutzen hier bzw.
				/// diese methode oben erweitern!
				/// Just edge for debuggging purposes
				ug::Edge* e1 = *g.create<RegularEdge>(EdgeDescriptor(p1, p2));
				UG_COND_THROW(!e1, "Edge (connecting inner quads with inner neurite) was not created");
			}

			for (std::map<ug::Vertex*, ug::Vertex*>::iterator it = innerToOuter2.begin(); it != innerToOuter2.end(); ++it) {
				/// inner quads to inner neurite
				ug::Vertex* p1 = it->first; /// soma vertex
				ug::Vertex* p2 = it->second; /// neurite vertex

				ug::Edge* e1 = *g.create<RegularEdge>(EdgeDescriptor(p1, p2));
				UG_COND_THROW(!e1, "Edge (conencting outer quads with outer neurite) was not created");
			}

			UG_LOGN("Inner done!");
			UG_LOGN("Next merge these vertices from above!");
		}


   		////////////////////////////////////////////////////////////////////////
		/// connect_inner_neurites_to_inner_soma
		////////////////////////////////////////////////////////////////////////
		void connect_inner_neurites_to_inner_soma
		(
			size_t somaIndex, /// inner soma index: beginning of inner sphere's quads (ER) is somaIndex+1, outer sphere's quads (ER) is somaIndex-numQuads-1
			size_t numQuads, /// number of total surface quads or neurites to connect to
		    Grid& g,
		    Grid::VertexAttachmentAccessor<APosition>& aaPos,
		    SubsetHandler& sh
		) {
			Selector::traits<Vertex>::iterator vit;
			Selector::traits<Vertex>::iterator vit_end;
			Selector::traits<Edge>::iterator eit;
			Selector::traits<Edge>::iterator eit_end;
			Selector sel(g);

			/// the projected vertices coordinates as a vector3 and the corresponding grid vertices
			std::vector<std::vector<ug::vector3> > projected;
			std::vector<std::vector<ug::Vertex*> > projectedVertices;
			std::vector<ug::vector3> normals;
			projected.resize(numQuads);
			projectedVertices.resize(numQuads);
			SaveGridToFile(g, "before_projections_inner.ugx");

			/// find all edges for each inner sphere's surface quad - take two starting at the same vertex to get two edges for normal calculation
			for (size_t i = 1; i < numQuads+1; i++) {
				UG_LOGN("Selecting now subset: " << somaIndex+i);
				sel.clear();
				/// Select inner soma quad
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+i, true);
				ug::Vertex* v0 = *(sel.vertices_begin());
				UG_LOGN("First vertex of subset: " << aaPos[v0]);
				SelectSubsetElements<Edge>(sel, sh, somaIndex+i, true);
				eit = sel.edges_begin();
				eit_end = sel.edges_end();
				std::vector<std::pair<ug::Vertex*, ug::Vertex*> > es;
				size_t count = 0;

				UG_LOGN("Trying to find edges...");
				for (; eit != eit_end; ++eit) {
					Edge* e = *eit;
					UG_LOGN("trying e->vertex(0)");
					std::pair<ug::Vertex*, ug::Vertex*> p;
					if (e->vertex(0) == v0) {
						p.first = v0;
						p.second = e->vertex(1);
						count++;
						es.push_back(p);
					}

					UG_LOGN("trying e->vertex(1)");
					if (e->vertex(1) == v0) {
						p.first = e->vertex(0);
						p.second = v0;
						es.push_back(p);
						count++;
					}
				}

				/// Two edges found starting in same vertex
				UG_LOGN("Found edges");
				UG_COND_THROW(count != 2, "Number of edges has to be two to calculate a normal");

				/// Now calculate the normal for this plane / inner sphere's quad
				ug::vector3 v1, v2;
				VecSubtract(v1, aaPos[es[0].first], aaPos[es[0].second]);
				UG_LOGN("Subtracted v1");
				VecSubtract(v2, aaPos[es[1].first], aaPos[es[1].second]);
				UG_LOGN("Subtracted v2");
				ug::vector3 normal;
				VecCross(normal, v1, v2);
				VecNormalize(normal, normal);
				sel.clear();
				UG_LOGN("calculated cross product");

				/// Project each of the outer sphere's quad (ER) to the plane described by the normal of the inner sphere's quad
				ug::vector3 vProjected;
				SelectSubsetElements<Vertex>(sel, sh, somaIndex-numQuads+i-1, true);
				vit = sel.vertices_begin();
				vit_end = sel.vertices_end();
				for (; vit != vit_end; ++vit) {
					ug::vector3 v;
					VecSubtract(v, aaPos[*vit], v1);
					number dot = VecDot(v, normal);
					vProjected = aaPos[*vit];
					VecScaleAdd(vProjected, 1.0, vProjected, dot, normal);
					ug::Vertex* projVert = *g.create<ug::RegularVertex>();
					const ug::vector3& n = -normal;
					ProjectPointToPlane(vProjected, aaPos[*vit], aaPos[es[0].first], n);
					aaPos[projVert] = vProjected;
					projected[i-1].push_back(vProjected);
					projectedVertices[i-1].push_back(*vit); /// save original vertex from which we projected
					g.erase(projVert); /// delete projected vertex (only used during debugging)
				}
				normals.push_back(normal);
				UG_LOGN("First projection!");
			}

			/// Note: Could also calculate an averaged plane, e.g. calculate two plane normals for each quad, average them
			/// Note: Should get normal not from the two points of each inner quad but define the normal to be the edge through the center of the inner sphere's quad (ER) and outer sphere's quad (ER) part

			/// Note: Strategie für innere Verbindungen
			/// 1. Berechne normale (Axis) durch den Mittelpunkt der beiden Oberflächenquads (inner soma und äußeres soma) für den ER Teil
			/// 2. Projiziere äußere soma quad vertices des ER auf inneres Soma quad ebene definiert durch normal und mittelpunkt
			/// 3. Projiziere innere soma quad vertices des ER auf inneres soma quad ebene definiert durch normal und mittelpunkt
			/// 4. Zentriere projizierte äußere soma quad vertices des er um den mittelpunkt des quads auf innerem soma
			/// 5. Berechne winkel und finde paare mit minimalem winkel offset (oder closest distance geht evt. auch)
			/// 6. Move die alten inneren soma quad vertices auf die position der projizierten vertices von äußerem soma quad
			/// 7. merke die paare von projizierten vertices, starte mit jewiels 2 paaren (2 projizierte vertices auf innerem soma , 2 bewegte vertices auf dem inneren soma)
			/// 8. builde faces für 4 vertices immer -> hexaeder. ==> Verbindungshexaeder muss dann gespeichert werden für refinement! ebenso soma radien für beide branching point anschlüsse!!!

			/// Note: Strategie für äußere Verbindungen
			/// Projiziere ebenso auf ebene, jetzt aber auf den äußeren quad vertices des ERs (äußeres Soma)
			/// Projiziere Neuritenstartknoten ebenso darauf. Finde kleinsten Winkel, dann merge diese Vertices!

			/// Calculate all angles with respect to a reference point for the projected outer sphere's quad (ER)
			/// vertices to the inner sphere's quad and the inner sphere's vertices
			/// First move projected vertices to center of inner soma quad -> then do angle calculation accordingly or distance calculation
			size_t j = 1;
			for (std::vector<std::vector<ug::vector3> >::iterator it = projected.begin(); it != projected.end(); ++it) {
				sel.clear();
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+j, true);
				vit = sel.vertices_begin();
				vit_end = sel.vertices_end();
				std::vector<ug::vector3> verts;
				for (; vit != vit_end; ++vit) {
					verts.push_back(aaPos[*vit]);
				}
				ug::vector3 centerOut;
				ug::vector3 centerOut2;
				ug::vector3 dir;
				CalculateCenter(centerOut, &(*it)[0], it->size()); // center of projected inner soma surface quad vertices
				CalculateCenter(centerOut2, &(verts[0]), verts.size()); /// center of unprojected inner soma surface quad
				VecSubtract(dir, centerOut, centerOut2);
				UG_LOGN("Centerout: " << centerOut);
				UG_LOGN("Centerout2: " << centerOut2);
				for (std::vector<ug::vector3>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
					/// center each projected vertex to center of inner soma quad
					VecSubtract(*it2, *it2, dir);
				}
			}
			SaveGridToFile(g, sh, "projection_after_centered.ugx");

			std::vector<std::vector<number> > allAngles;
			std::vector<std::vector<number> > allAnglesInner;
			j = 1;
			for (std::vector<std::vector<ug::vector3> >::const_iterator it = projected.begin(); it != projected.end(); ++it) {
				ug::vector3 centerOut;
				std::vector<number> angles;
				CalculateCenter(centerOut, &(*it)[0], it->size());
				calculate_angles(centerOut, *it, angles, normals, (*it)[0]);
				allAngles.push_back(angles);

				sel.clear();
				std::vector<number> angles2;
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+j, true);
				vit = sel.vertices_begin();
				vit_end = sel.vertices_end();
				std::vector<ug::vector3> verts;
				for (; vit != vit_end; ++vit) {
					verts.push_back(aaPos[*vit]);
				}
				calculate_angles(centerOut, verts, angles2, normals, (*it)[0]);
				allAnglesInner.push_back(angles2);
				j++;
			}

			for (std::vector<std::vector<number> >::const_iterator it = allAngles.begin(); it != allAngles.end(); ++it) {
				UG_LOGN("Quad angles projected from outer to inner...");
				for (std::vector<number>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
					UG_LOGN(*it2);
				}
				UG_LOGN("---");
			}

			for (std::vector<std::vector<number> >::const_iterator it = allAnglesInner.begin(); it != allAnglesInner.end(); ++it) {
				UG_LOGN("Quad angles inner...");
				for (std::vector<number>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
					UG_LOGN(*it2);
				}
				UG_LOGN("---");
			}

			/// TODO: Angle calculation here (Conversion from 0, 180 to 0, 360 interval) is wrong: use from connect_outer_* method to convert
			/// Remember pairs: For each projected vertices to the inner sphere's quad
			/// plane a corresponding vertices we projected from the outer sphere's quad
			/// exist these have to be connected by edges / faces to create a hexaeder
			std::vector<std::vector<std::pair<size_t, size_t > > > pairs;
			for (size_t k = 0; k < numQuads; k++) {
				std::vector<std::pair<size_t, size_t> > pair;
				for (size_t i = 0; i < allAngles[k].size(); i++) {
					number dist = std::numeric_limits<number>::infinity();
					size_t smallest = 0;
					for (size_t j = 0 ; j < allAnglesInner[k].size(); j++) {
						number altDist = allAnglesInner[k][j] - allAngles[k][i];
						altDist += (altDist>180) ? -360 : (altDist<-180) ? 360 : 0;
						altDist = std::abs(altDist);
						if (altDist < dist) {
							dist = altDist;
							smallest = j;
						}
					}
					pair.push_back(std::make_pair<size_t, size_t>(i, smallest));
				}
				pairs.push_back(pair);
			}

			for (std::vector<std::vector<std::pair<size_t, size_t > > >::const_iterator it = pairs.begin(); it != pairs.end(); ++it) {
				UG_LOGN("***");
				for (std::vector<std::pair<size_t, size_t> >::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
					UG_LOGN("Pair " << it2->first << " -> " << it2->second);
				}
				UG_LOGN("***");
			}

			/// find closest vertex instead of minimum angle difference: this should be safe for the inner sphere and outer sphere ER part connection. Note reference point has to be the same to make sense.
			j = 1;
			std::vector<std::vector<std::pair<ug::vector3, ug::vector3> > > myPairs;
			std::map<Vertex*, Vertex*> myPairs2;
			for (size_t i = 0; i < projected.size(); i++) { /// each outer quad projected vertices (now on inner soma)
				sel.clear();
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+j, true);
				vit = sel.vertices_begin();
				vit_end = sel.vertices_end();
				std::vector<ug::vector3> vertsInner; /// each inner quad vertices
				std::vector<ug::Vertex*> vertsInnerVtx; /// each inner quad vertices

				for (; vit != vit_end; ++vit) {
					vertsInner.push_back(aaPos[*vit]);
					vertsInnerVtx.push_back(*vit);
				}

				std::vector<std::pair<ug::vector3, ug::vector3> > pair;
				std::vector<std::pair<ug::Vertex*, ug::Vertex*> > pair2;
				for (size_t k = 0; k < projected[i].size(); k++) {
					number dist = std::numeric_limits<number>::infinity();
					ug::vector3 p1;
					ug::vector3 p2;
					size_t p1i;
					size_t p2i;
					for (size_t l = 0; l < vertsInner.size(); l++) {
						number altDist = VecDistance(projected[i][k], vertsInner[l]);
						if (altDist < dist) {
							dist = altDist;
							p1 = projected[i][k];
							p2 = vertsInner[l];
							p1i = k;
							p2i = l;
						}
					}
					pair.push_back(std::make_pair(p1, p2));
					myPairs2[vertsInnerVtx[p2i]] = projectedVertices[i][p1i];
				}
				myPairs.push_back(pair);
				j++;
			}

			UG_LOGN("vector pairs...");
			for (std::vector<std::vector<std::pair<ug::vector3, ug::vector3 > > >::const_iterator it = myPairs.begin(); it != myPairs.end(); ++it) {
				UG_LOGN("***");
				for (std::vector<std::pair<ug::vector3, ug::vector3> >::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
					UG_LOGN("Pair " << it2->first << " -> " << it2->second);
				}
				UG_LOGN("***");
			}

			/// Note: Could also project the inner soma's quad vertices onto the plane defined by two of the inner soma's quad vertices
			/// (Since two points where taken from each inner soma's quad not all points lie in the plane defined by the normal and the mentioned points)
			/// Iterate over all neurite connections (numQuads) and get the vertices of the
			/// inner sphere's quad edge each and find the corresponding unprojected (outer sphere's quad vertices) and form a face
			/// It is also possible to do the same procedure with the sorted angle differences above to create these faces if angles are correct
			SaveGridToFile(g, sh, "before_projections_inner_connections.ugx");
			for (size_t i = 1; i < numQuads+1; i++) {
				sel.clear();
				UG_LOGN("Selecting now subset: " << somaIndex+i);
				/// Select inner soma quad
				SelectSubsetElements<Edge>(sel, sh, somaIndex+i, true);
				eit = sel.edges_begin();
				eit_end = sel.edges_end();
				for (; eit != eit_end; ++eit) {
					UG_LOGN("Edge!");
					Edge* e = *eit;
					ug::Vertex* p1 = e->vertex(0);
					ug::Vertex* p2 = e->vertex(1);
					ug::Vertex* p3 = myPairs2[e->vertex(0)];
					ug::Vertex* p4 = myPairs2[e->vertex(1)];
					UG_COND_THROW( ! ((p1 != p2) && (p3 != p4)), "Non-unique vertices provided to create quadrilateral.");
					ug::Face* f = *g.create<Quadrilateral>(QuadrilateralDescriptor(p1, p3, p4, p2));
					UG_COND_THROW(!f, "Quadrilateral for connecting inner soma sphere (ER) with inner neurite conneting to outer sphere (PM)");
				}
			}

			SaveGridToFile(g, sh, "after_projections_inner.ugx");

			/// TODO: Verdrehung kann beseitigt werden wenn man die Knoten des inneren Soma Oberflächenquads
			/// auf die Ebene projiziert welche durch das Oberflächenquads (innere) des äußeren Somas projiziert
			/// werden, evt. muss dann nach der Projektion auf das Zentrum des äußeren Somaoberflächen Quads
			/// zentriert werden, oder auf den Mittelpunkt des äußeren Quads des äußeren Somas
		}

		////////////////////////////////////////////////////////////////////////
		/// find_quad_verts_on_soma
		////////////////////////////////////////////////////////////////////////
		std::vector<ug::vector3> find_quad_verts_on_soma
		(
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			std::vector<ug::Vertex*> verticesOld,
			std::vector<std::vector<number> > outRads,
			size_t si,
	   	    SubsetHandler& sh,
	   	    number rimSnapThresholdFactor,
	   	    size_t numQuads,
	   	    size_t numVerts
		) {
			std::vector<ug::vector3> centers;
			UG_LOGN("3. AdaptSurfaceGridToCylinder")
			Selector sel(g);
			for (size_t i = 0; i < verticesOld.size(); i++) {
				sel.clear();
				ug::vector3 normal;
				CalculateVertexNormal(normal, g, verticesOld[i], aaPos);
				number radius = outRads[i][0];
				AdaptSurfaceGridToCylinder(sel, g, verticesOld[i], normal, radius, 1.0*rimSnapThresholdFactor, aPosition);
				UG_LOGN("Adaption done");

				sel.clear();
				sel.select(verticesOld[i]);
				ExtendSelection(sel, 1, true);
				CloseSelection(sel);
				sel.deselect(verticesOld[i]);
				g.erase(verticesOld[i]);

				UG_LOGN("num edges: " << sel.num<Edge>());
				size_t numEdges = sel.num<Edge>();
				size_t j = 0;
				while (numEdges > numVerts) {
					SubsetHandler::traits<Edge>::iterator eit = sel.begin<Edge>();
					SubsetHandler::traits<Edge>::iterator end = sel.end<Edge>();
					number bestLength = -1;
					Edge* eBest = NULL;
					for (; eit != end; ++eit) {
						const Edge* ee = *eit;
						Vertex* const* verts = ee->vertices();
						if (bestLength == -1) {
							bestLength = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
							eBest = *eit;
						} else {
							number length = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
							if (length < bestLength) {
								eBest = *eit;
								bestLength = length;
							}
						}
					}
					CollapseEdge(g, eBest, eBest->vertex(0));
					numEdges--;
					j++;
				}
				UG_LOGN("Collapsing done");

				std::vector<ug::vector3> vertices;
				Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
				Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
				UG_LOGN("Pushing vertices");
				for (; vit != vit_end; ++vit) {
					vertices.push_back(aaPos[*vit]);
				}

				UG_LOGN("Number of vertices: " << vertices.size());
				ug::vector3 centerOut;
				CalculateCenter(centerOut, &vertices[0], sel.num<ug::Vertex>());
				UG_LOGN("centerOut: " << centerOut);
				centers.push_back(centerOut);
				}
			return centers;
		}

		////////////////////////////////////////////////////////////////////////
		/// connect_neurites_with_soma
		////////////////////////////////////////////////////////////////////////
		void connect_neurites_with_soma
		(
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			std::vector<Vertex*> outVerts,
			std::vector<Vertex*> outVertsInner,
			std::vector<number> outRads,
			std::vector<Vertex*>& smallerQuadVerts,
			int si,
			SubsetHandler& sh,
			const std::string& fileName,
			number rimSnapThresholdFactor,
			std::vector<std::pair<size_t, std::pair<ug::vector3, ug::vector3> > >& axisVectors,
			std::vector<NeuriteProjector::Neurite>& vNeurites,
			std::vector<std::vector<ug::Vertex*> >& connectingVertices,
			std::vector<std::vector<ug::Vertex*> >& connectingVerticesInner,
			std::vector<std::vector<ug::Edge*> >& connectingEdges,
			std::vector<std::vector<ug::Edge*> >& connectingEdgesInner,
			bool createInner,
			number alpha,
			int numIterations,
			number resolveThreshold,
			number scale,
			size_t numVerts,
			size_t numQuads
		) {
			UG_LOGN("1. Find the vertices representing dendrite connection to soma.");
			/// 1. Finde die 4 Vertices die den Dendritenanschluss darstellen zum Soma
			std::vector<std::vector<ug::vector3> > quads;
			std::vector<number> quadsRadii;
			numQuads = 1;

			for (size_t i = 0; i < numQuads; i++) {
				std::vector<ug::vector3> temp;
				for (size_t j = 0; j < numVerts; j++) {
					temp.push_back(aaPos[outVerts[(i*12)+j]]);
				}
				UG_LOGN("push a quad!");
				quads.push_back(temp);
			}

			UG_LOGN("2. Calculate center of each quad, find next surface vertex on soma.")
			/// 2. Berechne den Schwerpunkt jedes Quads und finde den nächstgelegenen
			///    Vertex des Oberflächengitters vom Soma
			std::vector<ug::vector3> centerOuts;
			std::vector<ug::vector3> centerOuts2;
			std::vector<Vertex*> bestVertices;
			for (size_t i = 0; i < numQuads; i++) {
				const ug::vector3* pointSet = &(quads[i][0]);
				ug::vector3 centerOut;
				CalculateCenter(centerOut, pointSet, numVerts);
				centerOuts.push_back(centerOut);
				Selector sel(g);
				SelectSubsetElements<Vertex>(sel, sh, si, true);
				Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
				Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
				number best = std::numeric_limits<number>::max();
				ug::Vertex* best_vertex = NULL;
				for (; vit != vit_end; ++vit) {
					number dist = VecDistance(aaPos[*vit], centerOut);
					if (dist < best) {
						best = dist;
						best_vertex = *vit;
					}
				}
				UG_COND_THROW(!best_vertex, "No best vertex found for quad >>" << i << "<<.");
				bestVertices.push_back(best_vertex);
			}

			AssignSubsetColors(sh);
			std::stringstream ss;
			ss << fileName << "_best_vertices.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			UG_LOGN("3. AdaptSurfaceGridToCylinder")
			UG_LOGN("Best vertices size: " << bestVertices.size());
			/// 3. Für jeden Vertex v führe AdaptSurfaceGridToCylinder mit Radius entsprechend
			///    dem anzuschließenden Dendritenende aus. Dadurch entsteht auf der Icosphere
			///    um jedes v ein trianguliertes 6- bzw. 5-Eck.
			Selector sel(g);
			for (size_t i = 0; i < bestVertices.size(); i++) {
				sel.clear();
				ug::vector3 normal;
				CalculateVertexNormal(normal, g, bestVertices[i], aaPos);
				number radius = outRads[i];
				AdaptSurfaceGridToCylinder(sel, g, bestVertices[i], normal, radius, 1.0*rimSnapThresholdFactor, aPosition);
			}

			AssignSubsetColors(sh);
			ss << fileName << "_before_deleting_center_vertices.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			UG_LOGN("5. MergeVertices")
			/// 5. Wandle die stückweise linearen Ringe um die Anschlusslöcher per
			///    MergeVertices zu Vierecken um.
			sel.clear();
			for (std::vector<Vertex*>::iterator it = bestVertices.begin(); it != bestVertices.end(); ++it) {
				sel.select(*it);
				ExtendSelection(sel, 1, true);
				CloseSelection(sel);
				AssignSelectionToSubset(sel, sh, sh.num_subsets()+1);
				sel.clear();
			}

			AssignSubsetColors(sh);
			ss << fileName << "_before_getting_neighborhoods.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			UG_LOGN("4. Remove each vertex. Creates holes in soma")
			/// 4. Lösche jedes v, sodass im Soma Anschlusslöcher für die Dendriten entstehen.
			sel.clear();
			for (std::vector<Vertex*>::iterator it = bestVertices.begin(); it != bestVertices.end(); ++it) {
				sel.select(*it);
			}

			size_t numSubsets = sh.num_subsets();
			AssignSelectionToSubset(sel, sh, numSubsets);
			EraseElements<Vertex>(g, sh.begin<Vertex>(numSubsets), sh.end<Vertex>(numSubsets));
			EraseEmptySubsets(sh);
			AssignSubsetColors(sh);
			ss << fileName << "_after_deleting_center_vertices.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			/// Collapse now edges and take smallest edges first
			int beginningOfQuads = si+1; // subset index where quads are stored in
			for (size_t i = 0; i < numQuads; i++) {
				int si = beginningOfQuads+i;
				size_t numEdges = sh.num<Edge>(si);
				size_t j = 0;
				while (numEdges > numVerts) {
					SubsetHandler::traits<Edge>::iterator eit = sh.begin<Edge>(si);
					SubsetHandler::traits<Edge>::iterator end = sh.end<Edge>(si);
					number bestLength = -1;
					Edge* eBest = NULL;
					for (; eit != end; ++eit) {
						const Edge* ee = *eit;
						Vertex* const* verts = ee->vertices();
						if (bestLength == -1) {
							bestLength = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
							eBest = *eit;
						} else {
							number length = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
							if (length < bestLength) {
								eBest = *eit;
								bestLength = length;
							}
						}
					}
					CollapseEdge(g, eBest, eBest->vertex(0));
					numEdges--;
					j++;
					std::stringstream ss;
					ss << fileName << "_after_collapse_number_" << j << "_for_quad_" << i << ".ugx";
					SaveGridToFile(g, sh, ss.str().c_str());
				}

				SubsetHandler::traits<Vertex>::iterator vit = sh.begin<Vertex>(si);
				SubsetHandler::traits<Vertex>::iterator vend = sh.end<Vertex>(si);

				/// Outer Soma is assumed to start at -5 * outRads[i] of corresponding neurite and inner soma is assumed to start at -1
				for (; vit != vend; ++vit) {
					// aaSurfParams[*vit].soma = true;  // soma is never used
					if (createInner)
						aaSurfParams[*vit].axial = -5 * outRads[i];
					else
						aaSurfParams[*vit].axial = -1;
				}
			}

			if (createInner) {
				/// Shrink each quad on the outer soma surface
				for (size_t i = 0; i < numQuads; i++) {
					sel.clear();
					int si = beginningOfQuads+i;

					SelectSubsetElements<Vertex>(sel, sh, si, true);
					UG_LOGN("Selecting subset index (si): " << si);
					UG_LOGN("Number of vertices: " << sel.num<Vertex>(si));
					std::vector<Vertex*> vrts;
					sel.clear();
					UG_LOGN("verts size: " << vrts.size());

					std::vector<Edge*> edges;
					SelectSubsetElements<Edge>(sel, sh, si, true);
					edges.assign(sel.edges_begin(), sel.edges_end());
					sel.clear();
					UG_LOGN("edges size: " << edges.size());

					vrts.push_back(edges[0]->vertex(0));
					vrts.push_back(edges[0]->vertex(1));
					Vertex* prevVertex = edges[0]->vertex(1);
					std::vector<size_t> indices;
					edges.erase(edges.begin());
					UG_LOGN("number of edges: " << edges.size());

					while (!edges.empty()) {
						UG_LOGN("Still running: edges.size(): " << edges.size());
						for (size_t i = 0; i < edges.size(); i++) {
							Edge* nextEdge = edges[i];
							if (nextEdge->vertex(0) == prevVertex) {
								UG_LOGN("push first if");
								vrts.push_back(nextEdge->vertex(1));
								prevVertex = nextEdge->vertex(1);
								edges.erase(edges.begin()+i);
								break;
							}
							UG_LOGN("in between");
							if (nextEdge->vertex(1) == prevVertex) {
								UG_LOGN("push second if")
		            			vrts.push_back(nextEdge->vertex(0));
								prevVertex = nextEdge->vertex(0);
								edges.erase(edges.begin()+i);
								break;
							}
						}
					}
					vrts.erase(vrts.end()-1);

					std::vector<ug::Edge*> vEdgeOut;
					std::vector<ug::Vertex*> vVrtOut;

					Selector selToAssign(g);
					UG_LOGN("Vrts.size(): " << vrts.size());

					/// UG_COND_THROW(vrts.size() != 4, "Non-quadrilateral encountered. Cannot shrink a non-quadrilateral!");
					/// shrink_quadrilateral_copy(vrts, vVrtOut, vVrtOut, vEdgeOut, g, aaPos, -scale, false, &selToAssign, NULL);
					UG_COND_WARNING(vrts.size() == 0, "Polygon to shrink has no vertices.");
					shrink_polygon_copy(vrts, vVrtOut, vVrtOut, vEdgeOut, g, aaPos, -scale, false, &selToAssign);
					for (std::vector<Vertex*>::const_iterator it = vVrtOut.begin(); it != vVrtOut.end(); ++it) {
						smallerQuadVerts.push_back(*it);
					}
					AssignSelectionToSubset(selToAssign, sh, si+numQuads);
					sel.clear();
					selToAssign.clear();

					SubsetHandler::traits<Vertex>::iterator vit = sh.begin<Vertex>(si);
					SubsetHandler::traits<Vertex>::iterator vend = sh.end<Vertex>(si);

					/// Inner soma is assumed to start at -5 * outRads[i] of corresponding neurite too
					for (; vit != vend; ++vit) {
						// aaSurfParams[*vit].soma = true;  // soma is never used
						aaSurfParams[*vit].axial =  -5 * outRads[i];
						//vNeurites[i].somaStart = 5 * outRads[i];  // somaStart is never used
					}
				}
			}

			/// Collapses only inner polygon to a quadrilateral - outer polygon needs to have 12 vertices
			if (createInner) {
				/// Collapse now edges and take smallest edges first
				int beginningOfQuads = si+1+numQuads; // subset index where inner quads are stored in
				for (size_t i = 0; i < numQuads; i++) {
					int si = beginningOfQuads+i;
					size_t numEdges = sh.num<Edge>(si);
					size_t j = 0;
					while (numEdges > 4) { // while more than 4 edges available
						SubsetHandler::traits<Edge>::iterator eit = sh.begin<Edge>(si);
						SubsetHandler::traits<Edge>::iterator end = sh.end<Edge>(si);
						number bestLength = -1;
						Edge* eBest = NULL;
						for (; eit != end; ++eit) {
							const Edge* ee = *eit;
							Vertex* const* verts = ee->vertices();
							if (bestLength == -1) {
								bestLength = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
								eBest = *eit;
							} else {
								number length = VecDistance(aaPos[verts[0]], aaPos[verts[1]]);
								if (length < bestLength) {
									eBest = *eit;
									bestLength = length;
								}
							}
						}
						CollapseEdge(g, eBest, eBest->vertex(0));
						numEdges--;
						j++;
						std::stringstream ss;
						ss << fileName << "_after_collapse_number_" << j << "_for_quad_" << i << ".ugx";
						SaveGridToFile(g, sh, ss.str().c_str());
					}
				}

				/// refine outer polygon once to get 12 vertices
				beginningOfQuads = si+1; // subset index where outer quads are stored
				for (size_t i = 0; i < numQuads; i++) {
					int si = beginningOfQuads+i;
					sel.clear();
					SelectSubsetElements<Vertex>(sel, sh, si, true);
					SelectSubsetElements<Edge>(sel, sh, si, true);
					Refine(g, sel, NULL, false);
				}
			}

			EraseEmptySubsets(sh);
			AssignSubsetColors(sh);

			ss << fileName << "_after_merging_cylinder_vertices.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			UG_LOGN("8. TangentialSmooth");
			/// Note: TangentialSmooth -> alpha has to be corrected for different geometries.
			/// TangentialSmooth(g, g.vertices_begin(), g.vertices_end(), aaPos, alpha, numIterations);

			ss << fileName << "_after_merging_cylinder_vertices_and_tangential_smooth.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			UG_LOGN("6. Extrude rings along normal")
			/// 6. Extrudiere die Ringe entlang ihrer Normalen mit Höhe 0 (Extrude mit
			///    aktivierter create faces Option).
			sel.clear();
			std::vector<std::vector<Vertex*> > somaVerts;
			std::vector<std::vector<Vertex*> > allVerts;

			std::vector<std::vector<Vertex*> > somaVertsInner;
			std::vector<std::vector<Vertex*> > allVertsInner;


			for (size_t i = 0; i < numQuads; i++) {
				int si = beginningOfQuads+i;
				ug::vector3 normal;
				CalculateVertexNormal(normal, g, *sh.begin<Vertex>(si), aaPos);
				UG_LOGN("normal (outer): " << normal);
				ug::vector3 axisVector;
				CalculateCenter(sh.begin<Vertex>(si), sh.end<Vertex>(si), aaPos);
				/// indicate soma posiiton
				/*for (SubsetHandler::traits<Vertex>::iterator it = sh.begin<Vertex>(si); it != sh.end<Vertex>(si); ++it) {
					UG_LOGN("setting axial to -1!");
					aaSurfParams[*it].axial = -1;
				}*/

				VecAdd(axisVector, axisVector, normal);
				std::vector<Edge*> edges;
				std::vector<Vertex*> vertices;

				for (SubsetHandler::traits<Edge>::const_iterator it = sh.begin<Edge>(si); it != sh.end<Edge>(si); ++it) {
					edges.push_back(*it);
				}

				for (SubsetHandler::traits<Vertex>::const_iterator it = sh.begin<Vertex>(si); it != sh.end<Vertex>(si); ++it) {
					vertices.push_back(*it);
				}

				if (createInner) {
					for (SubsetHandler::traits<Edge>::const_iterator it = sh.begin<Edge>(si+numQuads); it != sh.end<Edge>(si+numQuads); ++it) {
						edges.push_back(*it);
					}
					for (SubsetHandler::traits<Vertex>::const_iterator it = sh.begin<Vertex>(si+numQuads); it != sh.end<Vertex>(si+numQuads); ++it) {
						vertices.push_back(*it);
					}
				}

				SelectSubsetElements<Vertex>(sel, sh, si, true);
				std::vector<Vertex*> temp;
				temp.assign(sel.vertices_begin(), sel.vertices_end());
				somaVerts.push_back(temp);
				sel.clear();
				temp.clear();

				if (createInner) {
					SelectSubsetElements<Vertex>(sel, sh, si+numQuads, true);
					temp.assign(sel.vertices_begin(), sel.vertices_end());
					somaVertsInner.push_back(temp);
					sel.clear();
					temp.clear();
				}

				SelectSubsetElements<Vertex>(sel, sh, si, true);

				if (createInner) {
					SelectSubsetElements<Vertex>(sel, sh, si+numQuads, true);
				}

				/// Extrude(g, &vertices, &edges, NULL, normal, aaPos, EO_CREATE_FACES, NULL); (not needed anymore)

				/// store found soma vertices as connectingvertices for initial neurite vertices (used in create_neurite_general)
				for (size_t j= 0; j < vertices.size(); j++) {
					if (sh.get_subset_index(vertices[j]) == si) {
						connectingVertices[i].push_back(vertices[j]);
					}
					if (sh.get_subset_index(vertices[j]) == si+(int)numQuads) {
						connectingVerticesInner[i].push_back(vertices[j]);
					}
				}

				/// store found soma edges as connectingedges for initial neurite vertices (used in create_neurite_general)
				for (size_t j= 0; j < edges.size(); j++) {
					if (sh.get_subset_index(edges[j]) == si) {
						connectingEdges[i].push_back(edges[j]);
					}
					if (sh.get_subset_index(edges[j]) == si+(int)numQuads) {
						connectingEdgesInner[i].push_back(edges[j]);
					}
				}

				ug::vector3 centerOut2 = CalculateCenter(vertices.begin(), vertices.end(), aaPos);
				centerOuts2.push_back(centerOut2);
				sel.clear();

				/// indicate start of neurite with axial 0 and soma false explicitly: neurite start is 0 for outer soma and for inner soma it is -1 + radOut[i] * 5;
				for (std::vector<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
					aaSurfParams[*it].axial = -1 + 5 * outRads[i];
					// vNeurites[i].somaStart = 5 * outRads[i];  // somaStart is never used
					// aaSurfParams[*it].soma = false;  // soma is never used
				}

				SelectSubsetElements<Vertex>(sel, sh, si, true);
				std::vector<Vertex*> temp2;
				temp2.assign(sel.vertices_begin(), sel.vertices_end());
				allVerts.push_back(temp2);
				sel.clear();
				temp2.clear();

				if (createInner) {
					SelectSubsetElements<Vertex>(sel, sh, si+(int)numQuads, true);
					temp2.assign(sel.vertices_begin(), sel.vertices_end());
					allVertsInner.push_back(temp2);
					sel.clear();
					temp2.clear();
				}

				ug::vector3 cylinderCenter = CalculateCenter(sh.begin<Vertex>(si), sh.end<Vertex>(si), aaPos);
				axisVectors.push_back(make_pair(si, make_pair(axisVector, cylinderCenter)));
			}

			ss << fileName << "_after_extruding_cylinders_before_removing_common_vertices.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			/// delete common vertices, thus keep only newly extruded vertices (used in next step 7.)
			for (size_t i = 0; i < numQuads; i++) {
				size_t numSomaVerts = somaVerts[i].size();
				for (size_t j = 0; j < numSomaVerts; j++) {
					allVerts[i].erase(std::remove(allVerts[i].begin(), allVerts[i].end(), somaVerts[i][j]), allVerts[i].end());
				}
			}

			if (createInner) {
				/// delete common vertices, thus keep only newly extruded vertices (used in next step 7.)
				for (size_t i = 0; i < numQuads; i++) {
					size_t numSomaVerts = somaVertsInner[i].size();
					for (size_t j = 0; j < numSomaVerts; j++) {
						allVertsInner[i].erase(std::remove(allVertsInner[i].begin(), allVertsInner[i].end(), somaVertsInner[i][j]), allVertsInner[i].end());
					}
				}
			}

			ss << fileName << "_after_extruding_cylinders.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			UG_LOGN("7. Calculate convex hull and connect (Not needed anymore)")
			/// 7. Vereine per MergeVertices die Vertices der in 6. extrudierten Ringe jeweils
			///    mit den zu ihnen nächstgelegenen Vertices des entsprechenden Dendritenendes.
			si = beginningOfQuads;
			sel.clear();

			EraseEmptySubsets(sh);
			AssignSubsetColors(sh);
			ss << fileName << "_after_extruding_cylinders_and_merging.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
			ss.str(""); ss.clear();

			UG_LOGN("9. Resolve potentially generated intersection(s)")
			ResolveTriangleIntersections(g, g.begin<ug::Triangle>(), g.end<ug::Triangle>(), resolveThreshold, aPosition);

			ss << fileName << "_final.ugx";
			SaveGridToFile(g, sh, ss.str().c_str());
		}

		////////////////////////////////////////////////////////////////////////
		/// shrink_quadrilateral_copy
		////////////////////////////////////////////////////////////////////////
		void shrink_quadrilateral_copy
		(
			const std::vector<Vertex*>& vVrt,
			std::vector<Vertex*>& outvVrt,
			const std::vector<Vertex*>& oldVertices,
			std::vector<Edge*>& outvEdge,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			bool createFaces,
			ISelector* outSel,
			ug::vector3* currentDir
		) {
			Selector sel(g);
	    	for (size_t i = 0; i < 4; ++i) {
	    		sel.select(vVrt[i]);
			}

	    	ug::vector3 center;
	    	center = CalculateBarycenter(sel.vertices_begin(), sel.vertices_end(), aaPos);
	    	sel.clear();
	    	for (size_t i = 0; i < 4; ++i)
	    	{
	    	ug::Vertex* v = *g.create<RegularVertex>();
	       	   aaPos[v] = aaPos[vVrt[i]];
	       	   ug::vector3 dir;
	       	   VecSubtract(dir, aaPos[vVrt[i]], center);

	       	   /// Note: Shrink towards predscribed direction (currentDir) - otherwise shrink towards barycenter
	       	   if (currentDir) {
	       		   UG_THROW("Only shrinkage towards barycenter currently supported. Provided currentDir: " << currentDir);
	    	   	   /// 1. Get Edge e starting from i % 4 to i+1 % 4
	    	   	   /// 2. Check if e is parallel or anti-parallel to currentDir
	    	   	   /// 3. If true then calculate dir1, dir2 from vertex i, i+1 to center
	    	   	   /// 4. Project dir1 to edge e if e was parallel to currentDir
	    	   	   ///    otherwise project dir1 to edge -e if e was antiparallel to currentDir
	    	   	   /// 5. Project dir2 to edge -e if e was parallel to currentDir
	    	   	   ///    otherwise project dir2 to edge e if e was antiparallel to currentDir
	       	   }

	       	   UG_LOGN("dir:" << dir)
	       	   VecScaleAdd(aaPos[v], 1.0, aaPos[v], percentage, dir);

	       	   if (percentage > 1) {
	       		   UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
	       	   }
	       	   outvVrt.push_back(v);
	       	   if (outSel) outSel->select<ug::Vertex>(v);
	    	}

	    	/// create new edges in new (small) quad
	    	for (size_t i = 0; i < 4; ++i) {
	    	ug::Edge* e = *g.create<RegularEdge>(EdgeDescriptor(outvVrt[i], outvVrt[(i+1)%4]));
	       	   outvEdge.push_back(e);
	       	   if (outSel) outSel->select<ug::Edge>(e);
	    	}

	    	if (createFaces) {
	    		/// create new faces
	    		for (size_t i = 0; i < 4; ++i) {
	    			g.create<Quadrilateral>(QuadrilateralDescriptor(outvVrt[i], outvVrt[(i+1)%4], oldVertices[(i+1)%4], oldVertices[i]));
	    			/// Note: Do we need to flip faces here to wrt radial vector?
	    		}
	    	}
		}

		////////////////////////////////////////////////////////////////////////
		/// shrink_polygon_copy
		////////////////////////////////////////////////////////////////////////
		void shrink_polygon_copy
		(
			const std::vector<Vertex*>& vVrt,
			std::vector<Vertex*>& outvVrt,
			const std::vector<Vertex*>& oldVertices,
			std::vector<Edge*>& outvEdge,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			bool createFaces,
			ISelector* outSel
		) {
			Selector sel(g);
	    	for (size_t i = 0; i < vVrt.size(); ++i) {
	    		sel.select(vVrt[i]);
			}

	    	ug::vector3 center;
	    	center = CalculateBarycenter(sel.vertices_begin(), sel.vertices_end(), aaPos);
	    	sel.clear();
	    	for (size_t i = 0; i < vVrt.size(); ++i)
	    	{
	    	   ug::Vertex* v = *g.create<RegularVertex>();
	       	   aaPos[v] = aaPos[vVrt[i]];
	       	   ug::vector3 dir;
	       	   VecSubtract(dir, aaPos[vVrt[i]], center);

	       	   UG_LOGN("dir:" << dir)
	       	   VecScaleAdd(aaPos[v], 1.0, aaPos[v], percentage, dir);

	       	   if (percentage > 1) {
	       		   UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
	       	   }
	       	   outvVrt.push_back(v);
	       	   if (outSel) outSel->select<ug::Vertex>(v);
		    	}

			   	/// create new edges in new (small) quad
			   	for (size_t i = 0; i < vVrt.size(); ++i) {
			   		ug::Edge* e = *g.create<RegularEdge>(EdgeDescriptor(outvVrt[i], outvVrt[(i+1)%vVrt.size()]));
			    	   outvEdge.push_back(e);
			     	   if (outSel) outSel->select<ug::Edge>(e);
			   	}

			   	if (createFaces) {
			   		/// create new faces
			   		for (size_t i = 0; i < vVrt.size(); ++i) {
			   			g.create<Quadrilateral>(QuadrilateralDescriptor(outvVrt[i], outvVrt[(i+1)%4], oldVertices[(i+1)%4], oldVertices[i]));
			   			/// Note: Do we need to flip faces here to wrt radial vector?
		   		}
		   	}
		}

		////////////////////////////////////////////////////////////////////////
		/// shrink_quadrilateral
		////////////////////////////////////////////////////////////////////////
		void shrink_quadrilateral
		(
			std::vector<Vertex*> vVrt,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage
		) {
			Selector sel(g);
			vVrt.resize(4);
			for (size_t i = 0; i < 4; ++i) {
				sel.select(vVrt[i]);
			}

			ug::vector3 center;
			/// Note: Using barycenter, could also use vertex or area centroid for this:
			/// https://en.wikipedia.org/wiki/Quadrilateral#Remarkable_points_and_lines_in_a_convex_quadrilateral
			center = CalculateBarycenter(sel.vertices_begin(), sel.vertices_end(), aaPos);
			sel.clear();

			for (size_t i = 0; i < 4; ++i)
			{
				ug::vector3 dir;
				VecSubtract(dir, aaPos[vVrt[i]], center);
				UG_LOGN("dir:" << dir)
				VecScaleAdd(aaPos[vVrt[i]], 1.0, aaPos[vVrt[i]], percentage, dir);
				if (percentage > 1) {
					UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
				}
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// create_soma
		////////////////////////////////////////////////////////////////////////
		void create_soma
		(
			const std::vector<SWCPoint>& somaPts,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			size_t si,
			size_t numRefs
		) {
			UG_COND_WARNING(somaPts.size() != 1, "Currently only one soma point is allowed by this implementation");
			Selector sel(g);
			GenerateIcosphere(g, somaPts.front().coords, somaPts.front().radius, numRefs, aPosition, &sel);
			AssignSelectionToSubset(sel, sh, si);
		}

		////////////////////////////////////////////////////////////////////////
		/// create_soma
		////////////////////////////////////////////////////////////////////////
		void create_soma
		(
			const std::vector<SWCPoint>& somaPts,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos
		) {
			if (somaPts.size() == 1) {
			// create soma as icosahedron
			GenerateIcosahedron(g, somaPts[0].coords, somaPts[0].radius, aPosition);
			} else {
				// Note: Generalize this: Could take recipe from the following link to define
				/// and generate a discrete deformated icosphere respectively icosahedron:
				// http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
				UG_THROW("Currently only one soma point is allowed by this implementation.");
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// split_quadrilateral_along_edges
		////////////////////////////////////////////////////////////////////////
		void split_quadrilateral_along_edges
		(
			std::vector<Vertex*> vVrt,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			ug::vector3 vecDir,
			std::vector<ug::Vertex*>& vertices,
			std::vector<ug::Edge*>& edges,
			bool conservative
		) {
			/// "middle edges"
			std::vector<ug::Vertex*> from;
			std::vector<ug::Vertex*> to;
			Selector sel(g);
			vVrt.resize(4);
			size_t numPar = 0;
			for (size_t i = 0; i < 4; ++i) {
				ug::vector3 diffVec;
				VecSubtract(diffVec, aaPos[vVrt[i]], aaPos[vVrt[(i+1)%4]]);
				VecNormalize(diffVec, diffVec);
				VecNormalize(vecDir, vecDir);
				UG_LOGN("Parallel? " << VecDot(vecDir, diffVec));
				if (abs(VecDot(vecDir, diffVec)) > (1-0.1)) {
					numPar++;
					UG_LOGN("Parallel:" << VecDot(vecDir, diffVec));
					Edge* e = g.get_edge(vVrt[i], vVrt[(i+1)%4]);
					ug::RegularVertex* newVertex = SplitEdge<ug::RegularVertex>(g, e, conservative);
					ug::vector3 dir;
					VecSubtract(dir, aaPos[vVrt[i]], aaPos[vVrt[(i+1)%4]]);
					VecScaleAdd(aaPos[newVertex], 1.0, aaPos[vVrt[i]], percentage, dir);
					e = g.get_edge(newVertex, vVrt[(i+1)%4]);
					ug::RegularVertex* newVertex2 = SplitEdge<ug::RegularVertex>(g, e, conservative);
					VecScaleAdd(aaPos[newVertex2], 1.0, aaPos[vVrt[(i+1)%4]], -percentage, dir);
					from.push_back(newVertex);
					to.push_back(newVertex2);
				}
			}

			/// Note: Verify this is correct - seems okay.
			edges.push_back(g.get_edge(to[0], from[0]));
			ug::RegularEdge* e1 = *g.create<RegularEdge>(EdgeDescriptor(to[0], from[1]));
			edges.push_back(e1);
			edges.push_back(g.get_edge(to[1], from[1]));
			ug::RegularEdge* e2 = *g.create<RegularEdge>(EdgeDescriptor(to[1], from[0]));
			edges.push_back(e2);
			vertices.push_back(from[0]);
			vertices.push_back(to[0]);
			vertices.push_back(from[1]);
			vertices.push_back(to[1]);
			UG_COND_THROW(numPar != 2, "Shrinking of connecting quadrilateral failed!");
		}


		////////////////////////////////////////////////////////////////////////
		/// shrink_quadrilateral_center
		////////////////////////////////////////////////////////////////////////
		void shrink_quadrilateral_center
		(
			std::vector<Vertex*>& vVrt,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number percentage,
			ug::vector3& center
		) {
			for (size_t i = 0; i < 4; ++i)
		    {
		       ug::vector3 dir;
		       VecSubtract(dir, aaPos[vVrt[i]], center);

		       UG_LOGN("dir:" << dir)
		       VecScaleAdd(aaPos[vVrt[i]], 1.0, aaPos[vVrt[i]], percentage, dir);

		       if (percentage > 1) {
		          UG_WARNING("Moving vertex beyond center. Will create degenerated elements." << std::endl);
		       }
		    }
		}

		////////////////////////////////////////////////////////////////////////
		/// reorder_connecting_elements
		////////////////////////////////////////////////////////////////////////
		void reorder_connecting_elements
		(
			std::vector<ug::Vertex*>& v,
			std::vector<ug::Edge*> e
		) {
			std::vector<ug::Vertex*> sorted;
			sorted.push_back(v[0]);

			for (size_t j = 1; j < v.size(); ++j) {
				Vertex* next = v[j];
				for (size_t i = 0; i < e.size(); i++) {
					if ( (e[i]->vertex(0) == next && e[i]->vertex(1) == sorted[j-1]) ||
							(e[i]->vertex(1) == next && e[i]->vertex(0) == sorted[j-1])) {
						sorted.push_back(next);
						break;
					}
				}
			}
			UG_COND_THROW(sorted.size() != 4, "Did not find vertices to sort...");
			v = sorted;
		}

		////////////////////////////////////////////////////////////////////////
		/// correct_edges
		////////////////////////////////////////////////////////////////////////
		void correct_edges
		(
			std::vector<ug::Vertex*>& verts,
			std::vector<ug::Edge*>& edges,
			std::vector<ug::Vertex*>& oldVertsSorted,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number scale
		) {
			sort(verts.begin(), verts.end(), CompareBy<&NeuriteProjector::SurfaceParams::axial>(aaSurfParams));
			oldVertsSorted = verts;

			Edge* e1 = g.get_edge(verts[0], verts[2]);
			if (!e1) e1 = g.get_edge(verts[0], verts[3]);
			Edge* e2 = g.get_edge(verts[1], verts[2]);
			if (!e2) e2 = g.get_edge(verts[1], verts[3]);

			/// "bottom vertices of connecting inner face": e1->vertex(0) - newVertex1 - newVertex2 - e1->vertex(1)
			vector3 dir;
			VecSubtract(dir, aaPos[e1->vertex(1)], aaPos[e1->vertex(0)]);
			ug::RegularVertex* newVertex1 = *g.create<ug::RegularVertex>();
			ug::RegularVertex* newVertex2 = *g.create<ug::RegularVertex>();
			aaPos[newVertex1] = aaPos[e1->vertex(0)];
			aaPos[newVertex2] = aaPos[e1->vertex(1)];
			VecScaleAdd(aaPos[newVertex1], 1.0, aaPos[newVertex1], scale/2.0, dir);
			aaSurfParams[newVertex1] = aaSurfParams[e1->vertex(0)];
			aaSurfParams[newVertex1].axial = aaSurfParams[e1->vertex(0)].axial + scale/2*(aaSurfParams[e1->vertex(1)].axial - aaSurfParams[e1->vertex(0)].axial);
			aaSurfParams[newVertex1].neuriteID = aaSurfParams[e1->vertex(0)].neuriteID;
			//aaSurfParams[newVertex1].scale = aaSurfParams[e1->vertex(0)].scale;  // scale is never used
			VecScaleAdd(aaPos[newVertex2], 1.0, aaPos[newVertex2], -scale/2.0, dir);
			aaSurfParams[newVertex2] = aaSurfParams[e1->vertex(1)];
			aaSurfParams[newVertex2].axial = aaSurfParams[e1->vertex(1)].axial - scale/2*(aaSurfParams[e1->vertex(1)].axial - aaSurfParams[e1->vertex(0)].axial);
			aaSurfParams[newVertex2].neuriteID = aaSurfParams[e1->vertex(1)].neuriteID;
			//aaSurfParams[newVertex2].scale = aaSurfParams[e1->vertex(1)].scale;  // scale is never used

			/// "top vertices of connecting inner face": e2->vertex(0) - newVertex3 - newVertex4 - e2->vertex(1)
			vector3 dir2;
			VecSubtract(dir, aaPos[e2->vertex(1)], aaPos[e2->vertex(0)]);
			VecSubtract(dir2, aaPos[e2->vertex(1)], aaPos[e2->vertex(0)]);
			ug::RegularVertex* newVertex3 = *g.create<ug::RegularVertex>();
			ug::RegularVertex* newVertex4 = *g.create<ug::RegularVertex>();
			aaPos[newVertex3] = aaPos[e2->vertex(0)];
			aaPos[newVertex4] = aaPos[e2->vertex(1)];
			VecScaleAdd(aaPos[newVertex3], 1.0, aaPos[newVertex3], scale/2.0, dir);
			aaSurfParams[newVertex3] = aaSurfParams[e2->vertex(0)];
			aaSurfParams[newVertex3].axial =  aaSurfParams[e2->vertex(0)].axial + scale/2*(aaSurfParams[e2->vertex(1)].axial - aaSurfParams[e2->vertex(0)].axial);
			aaSurfParams[newVertex3].neuriteID = aaSurfParams[e2->vertex(0)].neuriteID;
			//aaSurfParams[newVertex3].scale = aaSurfParams[e2->vertex(0)].scale;  // scale is never used
			VecScaleAdd(aaPos[newVertex4], 1.0, aaPos[newVertex4], -scale/2.0, dir);
			aaSurfParams[newVertex4] = aaSurfParams[e2->vertex(1)];
			aaSurfParams[newVertex4].axial = aaSurfParams[e2->vertex(1)].axial - scale/2*(aaSurfParams[e2->vertex(1)].axial - aaSurfParams[e2->vertex(0)].axial);
			aaSurfParams[newVertex4].neuriteID = aaSurfParams[e2->vertex(1)].neuriteID;
			//aaSurfParams[newVertex4].scale = aaSurfParams[e2->vertex(1)].scale;  // scale is never used

			ug::RegularEdge* e31 = *g.create<RegularEdge>(EdgeDescriptor(newVertex1, newVertex3));
			g.create<Quadrilateral>(QuadrilateralDescriptor(e1->vertex(0), newVertex1, newVertex3, e2->vertex(0)));
			ug::RegularEdge* e24 = *g.create<RegularEdge>(EdgeDescriptor(newVertex4, newVertex2));
			g.create<Quadrilateral>(QuadrilateralDescriptor(e1->vertex(1), newVertex2, newVertex4, e2->vertex(1)));
			ug::RegularEdge* e12 =  *g.create<RegularEdge>(EdgeDescriptor(newVertex2, newVertex1));
			ug::RegularEdge* e43 =  *g.create<RegularEdge>(EdgeDescriptor(newVertex3, newVertex4));

			/// Verify edges are quasi parallel (should never happen but you never know)
			VecNormalize(dir, dir);
			VecNormalize(dir2, dir2);
			number dotProd = VecDot(dir, dir2) / (VecLength(dir) * VecLength(dir2));
			UG_COND_THROW( !( fabs(dotProd-1) < SMALL), "Edges need to be quasi parallel during splitting a hexaeder: " << dotProd);

			/// erase old edges
			g.erase(e1);
			g.erase(e2);

			/// set new face vertices for connection
			verts.clear();
		 	verts.push_back(newVertex1);
		 	verts.push_back(newVertex3);
		 	verts.push_back(newVertex4);
		 	verts.push_back(newVertex2);

		 	/// set new edge vertices for connection
		 	edges.clear();
		 	edges.push_back(e31);
		 	edges.push_back(e43);
		 	edges.push_back(e24);
		 	edges.push_back(e12);
		}

		////////////////////////////////////////////////////////////////////////
		/// correct_edges_all
		////////////////////////////////////////////////////////////////////////
		void correct_edges_all
		(
			std::vector<ug::Vertex*>& verts,
			std::vector<ug::Vertex*>& vertsOpp,
			std::vector<ug::Edge*>& edges,
			std::vector<ug::Edge*>& edgesOpp,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number scale
		) {
			/// the connecting vertices are needed later
			UG_LOGN("correcting edges connecting...")
				std::vector<ug::Vertex*> oldVertsSorted;
			correct_edges(verts, edges, oldVertsSorted, aaSurfParams, g, aaPos, scale);
			UG_LOGN("correcting edges opposing...")

			/// backside not needed
			std::vector<ug::Vertex*> oldVertsSortedOpp;
			correct_edges(vertsOpp, edgesOpp, oldVertsSortedOpp, aaSurfParams, g, aaPos, scale);
			g.create<Quadrilateral>(QuadrilateralDescriptor(vertsOpp[0], vertsOpp[1], vertsOpp[2], vertsOpp[3]));

			// connect the sides of splitted edges and fill top center and bottom center holes
			g.create<RegularEdge>(EdgeDescriptor(verts[0], vertsOpp[1]));
			g.create<RegularEdge>(EdgeDescriptor(verts[1], vertsOpp[0]));
			g.create<Quadrilateral>(QuadrilateralDescriptor(verts[0], verts[3], vertsOpp[2], vertsOpp[1]));
			g.create<RegularEdge>(EdgeDescriptor(verts[2], vertsOpp[3]));
			g.create<RegularEdge>(EdgeDescriptor(verts[3], vertsOpp[2]));
			g.create<Quadrilateral>(QuadrilateralDescriptor(verts[1], verts[2], vertsOpp[3], vertsOpp[0]));

			/// fill 4 more holes to the left and right on top and left and right on bottom
			g.create<Quadrilateral>(QuadrilateralDescriptor(verts[0], vertsOpp[1], oldVertsSortedOpp[1], oldVertsSorted[0]));
			g.create<Quadrilateral>(QuadrilateralDescriptor(verts[1], vertsOpp[0], oldVertsSortedOpp[0], oldVertsSorted[1]));
			g.create<Quadrilateral>(QuadrilateralDescriptor(verts[2], vertsOpp[3], oldVertsSortedOpp[3], oldVertsSorted[2]));
			g.create<Quadrilateral>(QuadrilateralDescriptor(verts[3], vertsOpp[2], oldVertsSortedOpp[2], oldVertsSorted[3]));
		}

		////////////////////////////////////////////////////////////////////////
		/// correct_axial_offset
		////////////////////////////////////////////////////////////////////////
		void correct_axial_offset
		(
			std::vector<ug::Vertex*>& verts,
			Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >& aaSurfParams,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			number scale
		) {
			// check for consistency
			UG_COND_THROW(verts.size() != 4, "Exactly 4 vertices are necessary on coarse grid level.");
			// sort to find min and max axial values
			sort(verts.begin(), verts.end(), CompareBy<&NeuriteProjector::SurfaceParams::axial>(aaSurfParams));
			number length = aaSurfParams[verts[2]].axial - aaSurfParams[verts[0]].axial;
			UG_LOGN("length TIMES scale/2: " << length*scale/2)
			// update surface parameters
			aaSurfParams[verts[0]].axial = aaSurfParams[verts[0]].axial + length*scale/2;
			aaSurfParams[verts[1]].axial = aaSurfParams[verts[1]].axial + length*scale/2;
			aaSurfParams[verts[2]].axial = aaSurfParams[verts[2]].axial - length*scale/2;
			aaSurfParams[verts[3]].axial = aaSurfParams[verts[3]].axial - length*scale/2;
		}


		////////////////////////////////////////////////////////////////////////
		/// add_soma_surface_to_swc
		////////////////////////////////////////////////////////////////////////
		void add_soma_surface_to_swc
		(
			const size_t& lines,
			const std::string& fn_precond,
			const std::string& fn_precond_with_soma,
			const std::vector<ug::vector3>& vPointsSomaSurface
		) {
			std::ifstream inFile(fn_precond.c_str());
			UG_COND_THROW(!inFile, "SWC input file '" << fn_precond << "' could not be opened for reading.");
			std::ofstream outFile(fn_precond_with_soma.c_str());
			UG_COND_THROW(!outFile, "SWC output file '" << fn_precond_with_soma << "' could not be opened for reading.");

			size_t lineCnt = 1;
			std::string line;
			int somaIndex;
			std::vector<SWCPoint> swcPoints;
			std::vector<number> rads;
			size_t j = 0;
			while (std::getline(inFile, line)) {
				// trim whitespace
				line = TrimString(line);

				// ignore anything from possible '#' onwards
				size_t nChar = line.size();
				for (size_t i = 0; i < nChar; ++i)
				{
					if (line.at(i) == '#')
					{
						line = line.substr(0, i);
						break;
					}
				}

				// empty lines can be ignored
				if (line.empty()) continue;

				// split the line into tokens
				std::istringstream buf(line);
				std::istream_iterator<std::string> beg(buf), end;
				std::vector<std::string> strs(beg, end);

				// assert number of tokens is correct
				UG_COND_THROW(strs.size() != 7, "Error reading SWC file '" << fn_precond
			          << "': Line " << lineCnt << " does not contain exactly 7 values.");

				// type
				if (boost::lexical_cast<int>(strs[6]) == -1) {
				   somaIndex = boost::lexical_cast<int>(strs[0]);
				    outFile << strs[0] << " " << strs[1] << " " << strs[2] << " "
				        		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
				        		<< strs[6] << std::endl;
				} else {
				 int index = boost::lexical_cast<int>(strs[6]);
			     if (index == somaIndex) {
			        	number rad = boost::lexical_cast<number>(strs[5]);
			        	outFile << lineCnt << " 3 " << vPointsSomaSurface[j].x() << " "
			        			<< vPointsSomaSurface[j].y() << " " << vPointsSomaSurface[j].z()
			        			<< " " << rad << " " << somaIndex << std::endl;
			        	outFile << lineCnt+1 << " " << strs[1] << " " << strs[2] << " "
					   		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
					   		<< lineCnt << std::endl;
			        	lineCnt++;
			        	j++;
			     } else {
			    	 int newIndex = boost::lexical_cast<int>(strs[6])+j;
			    	 outFile << lineCnt << " " << strs[1] << " " << strs[2] << " "
 	 					   		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
   	 					   		<< newIndex << std::endl;
			     }
			   }
			 lineCnt++;
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// get_closest_poins_to_soma
		////////////////////////////////////////////////////////////////////////
		void get_closest_points_to_soma
		(
			const std::string& fn_precond,
			std::vector<ug::vector3>& vPos,
			size_t& lines
		) {
			std::ifstream inFile(fn_precond.c_str());
			UG_COND_THROW(!inFile, "SWC input file '" << fn_precond << "' could not be opened for reading.");

	    	size_t lineCnt = 0;
	    	std::string line;
	    	int somaIndex;
	    	std::vector<SWCPoint> swcPoints;
	    	while (std::getline(inFile, line)) {
	    		lineCnt++;
	    		// trim whitespace
	        	line = TrimString(line);

	        	// ignore anything from possible '#' onwards
	        	size_t nChar = line.size();
	        	for (size_t i = 0; i < nChar; ++i)
	        	{
	        		if (line.at(i) == '#')
	            	{
	            		line = line.substr(0, i);
	                	break;
	            	}
	        	}

	        	// empty lines can be ignored
		    	if (line.empty()) continue;

		    	// split the line into tokens
		    	std::istringstream buf(line);
		    	std::istream_iterator<std::string> beg(buf), end;
		    	std::vector<std::string> strs(beg, end);

		    	// assert number of tokens is correct
		    	UG_COND_THROW(strs.size() != 7, "Error reading SWC file '" << fn_precond
		    		<< "': Line " << lineCnt << " does not contain exactly 7 values.");

		      // type
		      if (boost::lexical_cast<int>(strs[6]) == -1) {
			   	   somaIndex = boost::lexical_cast<int>(strs[0]);
		   	   } else {
		        	if (boost::lexical_cast<int>(strs[6]) == somaIndex) {
		        		number x = boost::lexical_cast<number>(strs[2]);
		        		number y = boost::lexical_cast<number>(strs[3]);
		        		number z = boost::lexical_cast<number>(strs[4]);
		        		vPos.push_back(ug::vector3(x, y, z));
		        	}
		   	   }
	    	}
	    	lines = lineCnt;
		}

		////////////////////////////////////////////////////////////////////////
		/// get_closet_vertices_on_soma
		////////////////////////////////////////////////////////////////////////
		void get_closest_vertices_on_soma
		(
			const std::vector<ug::vector3>& vPos,
			std::vector<ug::Vertex*>& vPointsSomaSurface, Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			size_t si
		) {
			UG_LOGN("Finding now: " << vPos.size());
			for (size_t i = 0; i < vPos.size(); i++) {
				const ug::vector3* pointSet = &vPos[i];
				ug::vector3 centerOut;
				CalculateCenter(centerOut, pointSet, 1);
				Selector sel(g);
				SelectSubsetElements<Vertex>(sel, sh, si, true);
				UG_LOGN("selected vertices: " << sel.num<Vertex>());
				Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
				Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
				number best = -1;
				ug::Vertex* best_vertex = NULL;
				for (; vit != vit_end; ++vit) {
					number dist = VecDistance(aaPos[*vit], centerOut);
					if (best == -1) {
						best = dist;
						best_vertex = *vit;
					} else if (dist < best) {
						best = dist;
						best_vertex = *vit;
					}
				}
				UG_COND_THROW(!best_vertex, "No best vertex found for root neurite >>" << i << "<<.");
				vPointsSomaSurface.push_back(best_vertex);
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// replace_first_root_neurite_vertex_in_swc
		////////////////////////////////////////////////////////////////////////
		void get_closest_points_on_soma
		(
			const std::vector<ug::vector3>& vPos,
			std::vector<ug::vector3>& vPointsSomaSurface, Grid& g,
			Grid::VertexAttachmentAccessor<APosition>& aaPos,
			SubsetHandler& sh,
			size_t si
		) {
			UG_LOGN("finding now: " << vPos.size());
			for (size_t i = 0; i < vPos.size(); i++) {
				const ug::vector3* pointSet = &vPos[i];
				ug::vector3 centerOut;
				CalculateCenter(centerOut, pointSet, 1);
				Selector sel(g);
				SelectSubsetElements<Vertex>(sel, sh, si, true);
				UG_LOGN("selected vertices: " << sel.num<Vertex>());
				Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
				Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
				number best = -1;
				ug::Vertex* best_vertex = NULL;
				for (; vit != vit_end; ++vit) {
					number dist = VecDistance(aaPos[*vit], centerOut);
					if (best == -1) {
						best = dist;
						best_vertex = *vit;
					} else if (dist < best) {
						best = dist;
						best_vertex = *vit;
					}
				}
				UG_COND_THROW(!best_vertex, "No best vertex found for root neurite >>" << i << "<<.");
				vPointsSomaSurface.push_back(aaPos[best_vertex]);
			}
		}

		////////////////////////////////////////////////////////////////////////
		/// replace_first_root_neurite_vertex_in_swc
		////////////////////////////////////////////////////////////////////////
		void replace_first_root_neurite_vertex_in_swc
		(
			const size_t& lines,
			const std::string& fn_precond,
			const std::string& fn_precond_with_soma,
			const std::vector<ug::vector3>& vPointsSomaSurface
		) {
			std::ifstream inFile(fn_precond.c_str());
			UG_COND_THROW(!inFile, "SWC input file '" << fn_precond << "' could not be opened for reading.");
			std::ofstream outFile(fn_precond_with_soma.c_str());
	    	UG_COND_THROW(!outFile, "SWC output file '" << fn_precond_with_soma << "' could not be opened for reading.");

	    	size_t lineCnt = 1;
	    	std::string line;
	    	int somaIndex;
	    	std::vector<SWCPoint> swcPoints;
	    	std::vector<number> rads;
	    	size_t j = 0;
	    	while (std::getline(inFile, line)) {
	    		// trim whitespace
	    		line = TrimString(line);

	    		// ignore anything from possible '#' onwards
	    		size_t nChar = line.size();
	    		for (size_t i = 0; i < nChar; ++i)
	    		{
	    			if (line.at(i) == '#')
	    			{
	    				line = line.substr(0, i);
	    				break;
	    			}
	    		}

	    		// empty lines can be ignored
	    		if (line.empty()) continue;

	    		// split the line into tokens
	    		std::istringstream buf(line);
	    		std::istream_iterator<std::string> beg(buf), end;
	    		std::vector<std::string> strs(beg, end);

	    		// assert number of tokens is correct
	    		UG_COND_THROW(strs.size() != 7, "Error reading SWC file '" << fn_precond
			          << "': Line " << lineCnt << " does not contain exactly 7 values.");

	    		// type
	    		if (boost::lexical_cast<int>(strs[6]) == -1) {
				   somaIndex = boost::lexical_cast<int>(strs[0]);
				    outFile << strs[0] << " " << strs[1] << " " << strs[2] << " "
				        		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
				        		<< strs[6] << std::endl;
	    		} else {
				 int index = boost::lexical_cast<int>(strs[6]);
			     if (index == somaIndex) {
			        	number rad = boost::lexical_cast<number>(strs[5]);
			        	outFile << lineCnt << " 3 " << vPointsSomaSurface[j].x() << " "
			        			<< vPointsSomaSurface[j].y() << " " << vPointsSomaSurface[j].z()
			        			<< " " << rad << " " << somaIndex << std::endl;
			        	///lineCnt++;
			        	j++;
			     } else {
			    	 int newIndex = boost::lexical_cast<int>(strs[6]);
			    	 outFile << lineCnt << " " << strs[1] << " " << strs[2] << " "
 	 					   		<< strs[3] << " " << strs[4] << " " << strs[5] << " "
   	 					   		<< newIndex << std::endl;
			     }
			   }
			 lineCnt++;
	    	}
		}
	}
}
