/*!
 * \file neurite_util.cpp
 *
 *  Created on: Apr 22, 2019
 *      Author: stephan
 */

#include "neurite_util.h"
#include <cmath>
#include "common/log.h"
#include "common/error.h"
#include "../../ugcore/ugbase/bridge/domain_bridges/selection_bridge.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h"

namespace ug {
	namespace neuro_collection {
		////////////////////////////////////////////////////////////////////////
		/// calculate_angle
		////////////////////////////////////////////////////////////////////////
		float calculate_angle(const vector3& pos, const vector3& origin, const vector3& point, const vector3& n) {
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
		float calculate_angle(const vector3& pos, const vector3& origin, const vector3& point) {
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
		float deg_to_full_range(float angle) {
			return fmod(angle+360, 360);
		}

		////////////////////////////////////////////////////////////////////////
		/// calculate_angles
		////////////////////////////////////////////////////////////////////////
		void calculate_angles(const vector3& originCenter, const std::vector<ug::vector3>& points, std::vector<number>& angles, size_t refIndex) {
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
		void calculate_angles(const vector3& originCenter, const std::vector<ug::vector3>& points,  std::vector<number>& angles, std::vector<ug::vector3>& normals, const ug::vector3& ref) {
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
		void calculate_angles(const vector3& originCenter, const std::vector<ug::vector3>& points, std::vector<number>& angles, std::vector<ug::vector3>& normals, size_t refIndex) {
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
		void sortIndices(const std::vector<number>& values, std::vector<size_t>& indices) {
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
		bool sortbysec(const std::pair<ug::Vertex*, number> &a,
	              const std::pair<ug::Vertex*, number> &b)
		{
			return (a.second < b.second);
		}

		////////////////////////////////////////////////////////////////////////
		/// deg360
		////////////////////////////////////////////////////////////////////////
		number deg360(number a) {
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
		)
		{
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
			for (std::vector<std::vector<ug::vector3> >::const_iterator it = projected.begin(); it != projected.end(); ++it) {
				ug::vector3 centerOut;
				std::vector<number> angles;
				CalculateCenter(centerOut, &(*it)[0], it->size());
				/// calculate_angles(centerOut, *it, angles);
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

			/// TODO: Maybe have to center the projected vertices around the soma sphere's large respectively small surface quad center
			for (size_t i = 0; i < allAngles.size(); i++) {
				for (size_t j = 0; j < allAngles[i].size(); j++) {
					number a = allAngles[i][j];
					number b = allAnglesInner[i][j];
					/// Convert from -180,180 to 0,360 -> then sort array based on this
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

			/*
			std::vector<std::vector<size_t> >indices;
			std::vector<std::vector<size_t> >indices2;
			indices.resize(allAngles.size());
			indices2.resize(allAnglesInner.size());
			/// Sort angles to determine the pairs: Pair the smallest angle, then the second smallest angle, and so forth...
			/// This will be used to connect neurite root vertices with the vertices on the outer soma -> same strategy
			/// as before: We remember pairs of unprojected and projected vertices and connect unprojected vertices with
			/// the vertices which correspond to the minimum angle pairs we found below by using the projected vertices
			for (size_t i = 0; i < allAngles.size(); i++) {
				std::vector<number> angleCopy = allAngles[i];
				std::vector<number> angleCopy2 = allAnglesInner[i];
				sortIndices(angleCopy, indices[i]);
				sortIndices(angleCopy2, indices2[i]);
			}
			*/

			/// get mapping of outer vertices to inner (Unprojected) vertices
			std::map<Vertex*, Vertex*> innerToOuter;
			std::map<Vertex*, Vertex*> innerToOuter2;
			for (size_t i = 0 ; i < projected.size(); i++) {
				/// calculate angles for inner original soma vertices and projected vertices (projectedVertices is ug::Vertex*)
				std::vector<std::pair<ug::Vertex*, number> > anglesOfProjectedInnerVertices;
				std::vector<std::pair<ug::Vertex*, number> > anglesOfOrginalSomaInnerVertices;

				/// TODO: convert angles from -180, 180 to 0, 360 degree with deg360 -> check if works.
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
					int index = FindClosestVertexInArray(array, aaPos[p1], aaPos, 10); /// closest soma vertex of inner quad to outer soma quad

					UG_COND_THROW(!array[index], "Not found!");
					UG_COND_THROW(!projectedToUnprojectedInner[p2], "Not found!");
					innerToOuter2[array[index]] = projectedToUnprojectedInner[projectedToUnprojected[temp]]; /// soma inner quad vertex => rootneurite inner vertex
				}
			}

			/*
			UG_LOGN("Sorted indices:")
			for (size_t j = 0; j < 1; j++) {
				for (size_t i = 0; i < indices[j].size(); i++) {
					UG_LOGN("i (indices):" << indices[j][i]);
					UG_LOGN("angles[i] (sorted): " << allAngles[j][indices[j][i]]);
				}
			}

			for (size_t j = 0; j < 1; j++) {
				for (size_t i = 0; i < indices2[j].size(); i++) {
					UG_LOGN("i (indices):" << indices2[j][i]);
					UG_LOGN("angles[i] (sorted): " << allAnglesInner[j][indices2[j][i]]);
				}
			}

			UG_LOGN("sorted angles");
			 */

			/// connect outer quad with outer neurite
			for (std::map<ug::Vertex*, ug::Vertex*>::iterator it = innerToOuter.begin(); it != innerToOuter.end(); ++it) {
				/// outer quads to outer neurite
				ug::Vertex* p1 = it->first; /// soma vertex
				ug::Vertex* p2 = it->second; /// neurite vertex

				/// Merge
				//MergeVertices(g, p1, p2);
				//aaPos[p2] = aaPos[p1]; /// p1 = p2 => merge at neurite, p2 = p1 => merge at soma surface

				/// TODO: beides mal ist die zuordnung der inneren vertices nicht gesicht wenn man die angles nicht für beide quads inner und außen mit gleicher referenz berechnet, bzgl. EINS mittelpunktes und EINER normalen, dann muss es gecentered werden um diesen mittelpunkt: einfacher oben die zuordnung treffen und diese methode nur 1 mal nutzen hier bzw. diese methode oben erweitern!
				/// Just edge for debggging
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
			UG_LOGN("created edges !!! for outer done !!! might be screwed because angles not sorted... looks good")
			UG_LOGN("Next merge these vertices from above!");

		/*
		for (size_t i = 0; i < 1; i++) {
			sel.clear();
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+1+i, true);
			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			std::vector<ug::Vertex*> vertsInnerVtx;

			for (; vit != vit_end; ++vit) {
				vertsInnerVtx.push_back(*vit);
			}

			for (size_t j = 0; j < 2; j++) {
				ug::Vertex* p1 = rootNeurites[(i*4)+indices2[i][j]];
				ug::Vertex* p2 = vertsInnerVtx[indices[i][j]];
				ug::Edge* e = *g.create<RegularEdge>(EdgeDescriptor(p1, p2));
			}
		}
		*/

		/*
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
						pair.push_back(make_pair<size_t, size_t>(indices[k][i], indices2[k][smallest])); /// TODO: use indices from above to map
					}
					pairs.push_back(pair);
				}

			UG_LOGN("found pairs");
			for (size_t i = 0; i < pairs.size(); i++) {
				sel.clear();
				SelectSubsetElements<Vertex>(sel, sh, somaIndex+1+i, true);
				vit = sel.vertices_begin();
				vit_end = sel.vertices_end();
				std::vector<ug::Vertex*> vertsInnerVtx;

				for (; vit != vit_end; ++vit) {
					vertsInnerVtx.push_back(*vit);
				}

				const std::vector<std::pair<size_t, size_t > >& temp = pairs[i]; // first quad
				for (size_t j = 0; j < temp.size(); j++) {
					const std::pair<size_t, size_t >& pair = temp[j]; // get a pair
					ug::Vertex* p1 = rootNeurites[(i*4)+pair.first];
					ug::Vertex* p2 = vertsInnerVtx[pair.second];
					ug::Edge* e = *g.create<RegularEdge>(EdgeDescriptor(p1, p2));
				}
			}
			*/
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
			/// For the test geometry this has to be 10
			/// UG_COND_THROW(somaIndex != 10, "Soma index needs to be 10!");
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
					projectedVertices[i-1].push_back(*vit); /// save original vertex from which we proojected
				}
				normals.push_back(normal);
				UG_LOGN("First projection!");
			}

			/// TODO: Could also calculate an averaged plane, e.g. calculate two plane normals for each quad, average them
			/// TODO: Should get normal not from the two points of each inner quad but define the normal to be the edge through the center of the inner sphere's quad (ER) and outer sphere's quad (ER) part

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

			/// Calculate all angles with respect to a reference point for the projected outer sphere's quad (ER) vertices to the inner sphere's quad and the inner sphere's vertices
			std::vector<std::vector<number> > allAngles;
			std::vector<std::vector<number> > allAnglesInner;
			size_t j = 1;
			for (std::vector<std::vector<ug::vector3> >::const_iterator it = projected.begin(); it != projected.end(); ++it) {
				ug::vector3 centerOut;
				std::vector<number> angles;
				CalculateCenter(centerOut, &(*it)[0], it->size());
				/// calculate_angles(centerOut, *it, angles);
				calculate_angles(centerOut, *it, angles, normals);
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
				calculate_angles(centerOut, verts, angles2, normals);
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

			/// Remember pairs: For each projected vertices to the inner sphere's quad plane a corresponding vertices we projected from the outer sphere's quad exist these have to be connected by edges / faces to create a hexaeder
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
			/// TODO: The angle is not calculate correctly: why? Because two calls to calculate_angles, but the REFERENCE point differs. Has to be precisely the same to make sense.

			/// find closest vertex instead of minimum angle difference: this should be save for the inner sphere and outer sphere ER part connection
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
			/// Nachdem vertices paare gefunden sind, müssen diese gemerkt werden welcher original vertex projiziert wurde.
			/// TODO: dann findet man die edges des inneren somas jeweils und für die edges die korrepsoniderten projizierten vertices => create face!
			/// Note: Winkelstrategie müsste anders behandelt werden wie oben beschrieben

			/// TODO Man kann zusätzlich auch die inneren Vertices auf die Ebene projizieren die zuvor definiert wurde
			/// (denn nicht alle liegen in der Ebene nur die zwei punkte die genommen wurden um die Normale zu berechnen)...
			/// und die projizierten äußeren Soma ER Quad verts auf das Zentrum verschieben siehe unten

			/// Iterate over all neurite connections (numQuads) and get the vertices of the inner sphere's quad edge each and find the corresponding unprojected (outer sphere's quad vertices) and form a face
			/// It is also possible to do the same procedure with the sorted angle differences above to create these faces
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
					ug::Face* f = *g.create<Quadrilateral>(QuadrilateralDescriptor(p1, p3, p4, p2));
					UG_COND_THROW(!f, "Quadrilateral for connecting inner soma sphere (ER) with inner neurite conneting to outer sphere (PM)");
				}
			}

			/// TODO: Vor dem verbinden, verschiebe inner quad vertices auf die
			/// projizierten positionen auf innerem Soma welche das Resultat sind
			/// von der Projektion der äußeren Quad Soma vertices: Sollte nicht nötig
			/// für gutgeartete innere Quads - benötigt falls projiziertes Quad weit
			/// entfernt von dem inneren Quad ist - insebsondere um minimales Angle
			/// Differenz oder minimale Entfernung zu berechnen um Vertices zu verbinden
			SaveGridToFile(g, "after_projections_inner.ugx");


		/// Old strategy with quickhull
		/*
			/// Get vertices attached to inner soma
			SelectSubsetElements<Vertex>(sel, sh, somaIndex+i, true);
			Selector::traits<Vertex>::iterator vit = sel.vertices_begin();
			Selector::traits<Vertex>::iterator vit_end = sel.vertices_end();
			/// select inner vertices connected to outer soma (start vertices of ER)

			vit = sel.vertices_begin();
			vit_end = sel.vertices_end();
			UG_LOGN("selected inner vertices attached to outer soma: " << sel.num());
			for (; vit != vit_end; ++vit) {
				Selector sel2(g);
				sel2.select(*vit);
				ExtendSelection(sel2, 1, true);
				if (sel2.num() == numNeighborHood) {
					temp.push_back(aaPos[*vit]);
					foo2.push_back(aaPos[*vit]);
				}
			}
			sel.clear();

			UG_COND_THROW(temp.size() != 8, "Need 8 vertices for calculating all faces.");
			#ifdef NC_WITH_QHULL
				using ug::neuro_collection::convexhull::gen_face;
				using ug::neuro_collection::convexhull::erase_face;
				gen_face(temp, g, sh, si+i, aaPos);
				erase_face(g, sh, si+i, aaPos, foo);
				erase_face(g, sh, si+i, aaPos, foo2);
			#else
				using ug::neuro_collection::quickhull::gen_face;
				gen_face(temp, temp2, g, sh, si+i, aaPos);
			#endif
		}
		*/
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
	   size_t numQuads
	) {
		size_t numVerts = 4;
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
	   size_t si,
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
	   number scale
) {
	UG_LOGN("1. Find the vertices representing dendrite connection to soma.");
	/// 1. Finde die 4 Vertices die den Dendritenanschluss darstellen zum Soma
	std::vector<std::vector<ug::vector3> > quads;
	std::vector<number> quadsRadii;
	size_t numVerts = 4;
	size_t numQuads = outVerts.size()/numVerts;

	for (size_t i = 0; i < numQuads; i++) {
		std::vector<ug::vector3> temp;
		for (size_t j = 0; j < numVerts; j++) {
			temp.push_back(aaPos[outVerts[(i*4)+j]]);
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
		UG_COND_THROW(!best_vertex, "No best vertex found for quad >>" << i << "<<.");
		bestVertices.push_back(best_vertex);
	}

	UG_LOGN("3. AdaptSurfaceGridToCylinder")
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
	std::stringstream ss;
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
	size_t beginningOfQuads = si+1; // subset index where quads are stored in
	for (size_t i = 0; i < numQuads; i++) {
		size_t si = beginningOfQuads+i;
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
			aaSurfParams[*vit].soma = true;
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
			size_t si = beginningOfQuads+i;

			SelectSubsetElements<Vertex>(sel, sh, si, true);
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
			size_t numIterations = edges.size()-1;
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
			UG_COND_THROW(vrts.size() != 4, "Non-quadrilateral encountered. Cannot shrink a non-quadrilateral!");
			shrink_quadrilateral_copy(vrts, vVrtOut, vVrtOut, vEdgeOut, g, aaPos, -scale, false, &selToAssign, NULL);
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
				aaSurfParams[*vit].soma = true;
				aaSurfParams[*vit].axial =  -5 * outRads[i];
				vNeurites[i].somaStart = 5 * outRads[i];
			}
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
		size_t si = beginningOfQuads+i;
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

		///Extrude(g, &vertices, &edges, NULL, normal, aaPos, EO_CREATE_FACES, NULL);

		/// store found soma vertices as connectingvertices for initial neurite vertices (used in create_neurite_general)
		for (size_t j= 0; j < vertices.size(); j++) {
			if (sh.get_subset_index(vertices[j]) == si) {
				connectingVertices[i].push_back(vertices[j]);
			}
			if (sh.get_subset_index(vertices[j]) == si+numQuads) {
				connectingVerticesInner[i].push_back(vertices[j]);
			}
		}

		/// store found soma edges as connectingedges for initial neurite vertices (used in create_neurite_general)
		for (size_t j= 0; j < edges.size(); j++) {
			if (sh.get_subset_index(edges[j]) == si) {
				connectingEdges[i].push_back(edges[j]);
			}
			if (sh.get_subset_index(edges[j]) == si+numQuads) {
				connectingEdgesInner[i].push_back(edges[j]);
			}
		}

		ug::vector3 centerOut2 = CalculateCenter(vertices.begin(), vertices.end(), aaPos);
		centerOuts2.push_back(centerOut2);
		sel.clear();

		/// indicate start of neurite with axial 0 and soma false explicitly: neurite start is 0 for outer soma and for inner soma it is -1 + radOut[i] * 5;
		for (std::vector<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
			aaSurfParams[*it].axial = -1 + 5 * outRads[i];
			vNeurites[i].somaStart = 5 * outRads[i];
			aaSurfParams[*it].soma = false;
		}

		SelectSubsetElements<Vertex>(sel, sh, si, true);
		std::vector<Vertex*> temp2;
		temp2.assign(sel.vertices_begin(), sel.vertices_end());
		allVerts.push_back(temp2);
		sel.clear();
		temp2.clear();

		if (createInner) {
			SelectSubsetElements<Vertex>(sel, sh, si+numQuads, true);
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


	UG_LOGN("7. Calculate convex hull and connect")
	/// 7. Vereine per MergeVertices die Vertices der in 6. extrudierten Ringe jeweils
	///    mit den zu ihnen nächstgelegenen Vertices des entsprechenden Dendritenendes.
	si = beginningOfQuads;
	sel.clear();
	/*
	for (size_t i = 0; i < numQuads; i++) {
		std::vector<ug::vector3> temp;
		std::vector<ug::vector3> foo;
		std::vector<ug::vector3> foo2;
		std::vector<Vertex*> temp2;
		for (size_t j = 0; j < numVerts; j++) {
			temp.push_back(aaPos[outVerts[i*4+j]]);
			foo.push_back(aaPos[outVerts[i*4+j]]);
		}
		for (size_t j = 0; j < numVerts; j++) {
			temp.push_back(aaPos[allVerts[i][j]]);
			foo2.push_back(aaPos[allVerts[i][j]]);
		}

		UG_COND_THROW(temp.size() != 8, "Need 8 vertices for calculating all faces.");
		#ifdef NC_WITH_QHULL
			using ug::neuro_collection::convexhull::gen_face;
			using ug::neuro_collection::convexhull::erase_face;
			gen_face(temp, g, sh, si+i, aaPos);
			erase_face(g, sh, si+i, aaPos, foo);
			erase_face(g, sh, si+i, aaPos, foo2);
		#else
			using ug::neuro_collection::quickhull::gen_face;
			gen_face(temp, temp2, g, sh, si+i, aaPos);
		#endif
		ug::vector3 center;
		center = CalculateCenter(sh.begin<Vertex>(si+i+1000), sh.end<Vertex>(si+i+1000), aaPos);
		ug::vector3 axis;
		VecSubtract(axis, centerOuts2[i], centerOuts[i]);
		axisVectors.push_back(make_pair(si+i+numQuads*2, make_pair(axis, center)));
		/// numQuads*2 is required: n-inner quads and n-outer quads -> these quads here are the connecting quads
		/// i.e. first come all outer quads, then all inner quads, then the connecting outer quads, then the connecting inner quads
	}
	*/

	/*
	if (!createInner) {
		si = beginningOfQuads+numQuads;
		sel.clear();
		for (size_t i = 0; i < numQuads; i++) {
			UG_LOGN("First quad to connect...: " << i);
			std::vector<ug::vector3> temp;
			std::vector<Vertex*> temp2;
			std::vector<ug::vector3> foo;
			std::vector<ug::vector3> foo2;
			UG_LOGN("Accessing outVertsInner...; " << i);
			for (size_t j = 0; j < numVerts; j++) {
				temp.push_back(aaPos[outVertsInner[i*4+j]]);
				foo.push_back(aaPos[outVertsInner[i*4+j]]);
			}
			UG_LOGN("Accessing allVertsInner...; " << i);
			for (size_t j = 0; j < numVerts; j++) {
				temp.push_back(aaPos[allVertsInner[i][j]]);
				foo2.push_back(aaPos[allVertsInner[i][j]]);
			}

			UG_LOGN("Checking consistency of temp...");
			UG_COND_THROW(temp.size() != 8, "Need 8 vertices for calculating all faces.");
			#ifdef NC_WITH_QHULL
				UG_LOGN("Trying to use convexhull...");
				using ug::neuro_collection::convexhull::gen_face;
				using ug::neuro_collection::convexhull::erase_face;
				gen_face(temp, g, sh, si+i, aaPos);
				erase_face(g, sh, si+i, aaPos, foo);
				erase_face(g, sh, si+i, aaPos, foo2);
			#else
				UG_LOGN("Trying to use quickhull... ")
				using ug::neuro_collection::quickhull::gen_face;
				gen_face(temp, temp2, g, sh, si+i, aaPos);
			#endif
				UG_LOGN("Done with quickhull/convexhull...");
		}
	}
	*/

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
		)
		{
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

	       	   /// TODO: Add the correct shrinkage into a prescribed direction
	       	   if (currentDir) {
	    	   	   /// 1. get Edge e starting from i % 4 to i+1 % 4
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
	    			ug::Face* f = *g.create<Quadrilateral>(QuadrilateralDescriptor(outvVrt[i], outvVrt[(i+1)%4], oldVertices[(i+1)%4], oldVertices[i]));
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
		)
		{
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
		)
		{
			UG_COND_THROW(somaPts.size() != 1, "Currently only one soma point is allowed by this implementation");
			Selector sel(g);
			GenerateIcosphere(g, somaPts.front().coords, somaPts.front().radius, numRefs, aPosition, &sel);
			AssignSelectionToSubset(sel, sh, si);
		}
	}
}