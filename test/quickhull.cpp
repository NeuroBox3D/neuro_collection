/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2018-06-20
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

/*
 * Implementation: Incremental convex hull algorithm
 * Reference: Computational Geometry in C by O'Rourke
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include "common/math/ugmath_types.h"   // vector3
#include "common/util/smart_pointer.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/refinement/projectors/neurite_projector.h"

/// TODO: Delete all

namespace ug {
	namespace neuro_collection {
		using namespace std;
		#define MAXN 1010
		typedef number vtype;

        /* Basic 3D vector implementation */
        struct vec3 {
			vec3() {
              X[0] = X[1] = X[2] = 0;
            }
            vec3(vtype x, vtype y, vtype z) {
               X[0] = x;
               X[1] = y;
               X[2] = z;
            }

            /* 3D cross product */
            vec3 operator*(const vec3& v) const {
                    return vec3(X[1] * v.X[2] - X[2] * v.X[1],
                    X[2] * v.X[0] - X[0] * v.X[2], X[0] * v.X[1] - X[1] * v.X[0]);
            }

            vec3 operator-(const vec3& v) const {
                  return vec3(X[0] - v.X[0], X[1] - v.X[1], X[2] - v.X[2]);
            }

            vec3 operator-() const {
                  return vec3(-X[0], -X[1], -X[2]);
             }

            vtype dot(const vec3& v) const {
                   return X[0] * v.X[0] + X[1] * v.X[1] + X[2] * v.X[2];
            }

            vtype X[3];
            friend std::ostream& operator <<(std::ostream &out, const vec3 &v) {
                   out << v.X[0] << "/" << v.X[1] << "/" << v.X[2];
                   return out;
            }
        };

        namespace quickhull {
        	/* Original points in the input. */
        	vec3 A[MAXN];

        	/* E[i][j] indicates which (up to two) other points combine with the edge i and
        	 * j to make a face in the hull.  Only defined when i < j.
        	 */
        	struct twoset {
                void insert(int x) {
                        (a == -1 ? a : b) = x;
                }
                bool contains(int x) {
                        return a == x || b == x;
                }
                void erase(int x) {
                        (a == x ? a : b) = -1;
                }
                int size() {
                        return (a != -1) + (b != -1);
                }
                int a, b;
        	 } E[MAXN][MAXN];

        	 struct face {
        		 vec3 norm;
                 vtype disc;
                 int I[3];
        	 };

        	 /* Compute the half plane {x : c^T norm < disc}
        	  * defined by the three points A[i], A[j], A[k] where
        	  * A[inside_i] is considered to be on the 'interior' side of the face. */
        	 static face make_qhull_face(int i, int j, int k, int inside_i) {
        		E[i][j].insert(k);
                E[i][k].insert(j);
                E[j][k].insert(i);

                face f;
                f.I[0] = i;
                f.I[1] = j;
                f.I[2] = k;
                f.norm = (A[j] - A[i]) * (A[k] - A[i]);
                f.disc = f.norm.dot(A[i]);
                if (f.norm.dot(A[inside_i]) > f.disc) {
                        f.norm = -f.norm;
                        f.disc = -f.disc;
                }
                return f;
        	 }

        	 static void gen_face
        	 (
        		 const std::vector<ug::vector3>& points,
                 const std::vector<Vertex*>& verts, Grid& g, SubsetHandler& sh,
                 size_t si, Grid::VertexAttachmentAccessor<APosition>& aaPos
             )
        	 {
                 for (size_t i = 0; i < points.size(); i++) {
                        A[i].X[0] = points[i].x();
                        A[i].X[1] = points[i].y();
                        A[i].X[2] = points[i].z();
                }

                /* Initially construct the hull as containing only the first four points. */
                face f;
                vector<face> faces;
                memset(E, -1, sizeof(E));
                for (int i = 0; i < 4; i++)
                        for (int j = i + 1; j < 4; j++)
                                for (int k = j + 1; k < 4; k++) {
                                        faces.push_back(make_qhull_face(i, j, k, 6 - i - j - k));
                                }

                /* Now add a point into the hull one at a time. */
                for (int i = 4; i < points.size(); i++) {
                        /* Find and delete all faces with their outside 'illuminated' by this
                         * point. */
                        for (size_t j = 0; j < faces.size(); j++) {
                                f = faces[j];
                                if (f.norm.dot(A[i]) > f.disc) {
                                        E[f.I[0]][f.I[1]].erase(f.I[2]);
                                        E[f.I[0]][f.I[2]].erase(f.I[1]);
                                        E[f.I[1]][f.I[2]].erase(f.I[0]);
                                        faces[j--] = faces.back();
                                        faces.resize(faces.size() - 1);
                                }
                        }
                        /* Now for any edge still in the hull that is only part of one face
                         * add another face containing the new point and that edge to the hull. */
                        int nfaces = faces.size();
                        for (int j = 0; j < nfaces; j++) {
                                f = faces[j];
                                for (int a = 0; a < 3; a++)
                                        for (int b = a + 1; b < 3; b++) {
                                                int c = 3 - a - b;
                                                if (E[f.I[a]][f.I[b]].size() == 2)
                                                        continue;
                                                faces.push_back(make_qhull_face(f.I[a], f.I[b], i, f.I[c]));
                                 }
                        }
                }

                std::cout << "number of faces found: " << faces.size() << std::endl;
                for (int i = 0; i < faces.size(); i++) {
                        int I[3] = { faces[i].I[0], faces[i].I[1], faces[i].I[2] };
                        std::cout << "Face: " << A[I[0]] << ", " << A[I[1]] << ", " << A[I[2]] << std::endl;
                        Vertex* v1 = verts[I[0]];
                        Vertex* v2 = verts[I[1]];
                        Vertex* v3 = verts[I[2]];
                        Selector sel(g);
                        Triangle* t1 = *g.create<Triangle>(TriangleDescriptor(v1, v2, v3));
                        std::cout << A[I[0]].X[0] << ", " << A[I[0]].X[1] << ", "
                                        << A[I[0]].X[2] << std::endl;
                        Vertex* v1_2 = *g.create<RegularVertex>();
                        aaPos[v1_2] = ug::vector3(A[I[0]].X[0], A[I[0]].X[1], A[I[0]].X[2]);
                        std::cout << ug::vector3(A[I[0]].X[0], A[I[0]].X[1], A[I[0]].X[2])
                        << std::endl;
                        Vertex* v2_2 = *g.create<RegularVertex>();
                        aaPos[v2_2] = ug::vector3(A[I[1]].X[0], A[I[1]].X[1], A[I[1]].X[2]);
                        std::cout << ug::vector3(A[I[1]].X[0], A[I[1]].X[1], A[I[1]].X[2])
                        << std::endl;
                        Vertex* v3_2 = *g.create<RegularVertex>();
                        aaPos[v3_2] = ug::vector3(A[I[2]].X[0], A[I[2]].X[1], A[I[2]].X[2]);
                        std::cout << ug::vector3(A[I[2]].X[0], A[I[2]].X[1], A[I[2]].X[2])
                        << std::endl;
                        Triangle* t2 = *g.create<Triangle>(
                                        TriangleDescriptor(v1_2, v2_2, v3_2));
                        sel.select(t2);
                        sel.select(v1_2);
                        sel.select(v2_2);
                        sel.select(v3_2);
                        AssignSelectionToSubset(sel, sh, si + 1000);
                        sel.clear();
                        sel.select(t1);
                        CloseSelection(sel);
                        AssignSelectionToSubset(sel, sh, si);
                   }
        	 }
        }
    }
 }
