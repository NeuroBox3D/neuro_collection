extern "C" {
	#include <libqhull/qhull_a.h>
}

#include <stdarg.h>
#include <stdio.h>

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

namespace ug {
	namespace neuro_collection {
		namespace convexhull {
			////////////////////////////////////////////////////////////////////
			/// qh_new_qhull
			////////////////////////////////////////////////////////////////////
        	int qh_new_qhull
        	(
        	   int dim,
        	   int numpoints,
        	   coordT *points,
        	   boolT ismalloc,
               char *qhull_cmd,
               FILE *outfile,
               FILE *errfile
            )
        	{
               int exitcode, hulldim;
               boolT new_ismalloc;
               static boolT firstcall = True;
               coordT *new_points;
               if(!errfile){
            	   errfile= stderr;
                }

                if (firstcall) {
                    qh_meminit(errfile);
                    firstcall= False;
                } else {
                    qh_memcheck();
                }

                if (strncmp(qhull_cmd, "qhull ", (size_t)6)) {
                    qh_fprintf(errfile, 6186, "qhull error (qh_new_qhull): start qhull_cmd argument with \"qhull \"\n");
                    return qh_ERRinput;
                }

                qh_initqhull_start(NULL, outfile, errfile);
                if(numpoints==0 && points==NULL){
                    trace1((qh ferr, 1047, "qh_new_qhull: initialize Qhull\n"));
                    return 0;
                }

                trace1((qh ferr, 1044, "qh_new_qhull: build new Qhull for %d %d-d points with %s\n", numpoints, dim, qhull_cmd));
                exitcode = setjmp(qh errexit);
                if (!exitcode) {
                    qh NOerrexit = False;
                    qh_initflags(qhull_cmd);
                    if (qh DELAUNAY)
                      qh PROJECTdelaunay= True;
                    if (qh HALFspace) {
                      /* points is an array of halfspaces,
                         the last coordinate of each halfspace is its offset */
                      hulldim= dim-1;
                      qh_setfeasible(hulldim);
                      new_points= qh_sethalfspace_all(dim, numpoints, points, qh feasible_point);
                      new_ismalloc= True;
                      if (ismalloc)
                        qh_free(points);
                    } else {
                      hulldim= dim;
                      new_points= points;
                      new_ismalloc= ismalloc;
                 }
                 qh_init_B(new_points, numpoints, hulldim, new_ismalloc);
                 qh_qhull();
                 qh_check_output();
                 if (outfile) {
                    qh_produce_output();
                 } else {
                    qh_prepare_output();
                 }

                 if (qh VERIFYoutput && !qh STOPpoint && !qh STOPcone)
                      qh_check_points();
                  }
                  qh NOerrexit = True;
                  return exitcode;
                }


				////////////////////////////////////////////////////////////////
				/// qh_errexit
				////////////////////////////////////////////////////////////////
                void qh_errexit(int exitcode, facetT *facet, ridgeT *ridge) {
                  if (qh ERREXITcalled) {
                    qh_fprintf(qh ferr, 8126, "\nqhull error while processing previous error.  Exit program\n");
                    qh_exit(qh_ERRqhull);
                  }

                  qh ERREXITcalled= True;
                  if (!qh QHULLfinished)
                    qh hulltime= qh_CPUclock - qh hulltime;
                  qh_errprint("ERRONEOUS", facet, NULL, ridge, NULL);
                  qh_fprintf(qh ferr, 8127, "\nWhile executing: %s | %s\n", qh rbox_command, qh qhull_command);
                  qh_fprintf(qh ferr, 8128, "Options selected for Qhull %s:\n%s\n", qh_version, qh qhull_options);
                  if (qh furthest_id >= 0) {
                    qh_fprintf(qh ferr, 8129, "Last point added to hull was p%d.", qh furthest_id);
                    if (zzval_(Ztotmerge))
                      qh_fprintf(qh ferr, 8130, "  Last merge was #%d.", zzval_(Ztotmerge));
                    if (qh QHULLfinished)
                      qh_fprintf(qh ferr, 8131, "\nQhull has finished constructing the hull.");
                    else if (qh POSTmerging)
                      qh_fprintf(qh ferr, 8132, "\nQhull has started post-merging.");
                    qh_fprintf(qh ferr, 8133, "\n");
                  }

                  if (qh FORCEoutput && (qh QHULLfinished || (!facet && !ridge)))
                    qh_produce_output();
                  else if (exitcode != qh_ERRinput) {
                    if (exitcode != qh_ERRsingular && zzval_(Zsetplane) > qh hull_dim+1) {
                      qh_fprintf(qh ferr, 8134, "\nAt error exit:\n");
                      qh_printsummary(qh ferr);
                      if (qh PRINTstatistics) {
                        qh_collectstatistics();
                        qh_printstatistics(qh ferr, "at error exit");
                        qh_memstatistics(qh ferr);
                      }
                    }
                    if (qh PRINTprecision)
                      qh_printstats(qh ferr, qhstat precision, NULL);
                  }
                  if (!exitcode)
                    exitcode= qh_ERRqhull;
                  else if (exitcode == qh_ERRsingular)
                    qh_printhelp_singular(qh ferr);
                  else if (exitcode == qh_ERRprec && !qh PREmerge)
                    qh_printhelp_degenerate(qh ferr);
                  if (qh NOerrexit) {
                    qh_fprintf(qh ferr, 6187, "qhull error while ending program, or qh->NOerrexit not cleared after setjmp(). Exit program with error.\n");
                    qh_exit(qh_ERRqhull);
                  }
                  qh ERREXITcalled= False;
                  qh NOerrexit= True;
                  qh ALLOWrestart= False;  // longjmp will undo qh_build_withrestart
                  longjmp(qh errexit, exitcode);
                }

				////////////////////////////////////////////////////////////////
				/// qh_errprint
				////////////////////////////////////////////////////////////////
                void qh_errprint(const char *string, facetT *atfacet, facetT *otherfacet, ridgeT *atridge, vertexT *atvertex) {
                  int i;

                  if (atfacet) {
                    qh_fprintf(qh ferr, 8135, "%s FACET:\n", string);
                    qh_printfacet(qh ferr, atfacet);
                  }
                  if (otherfacet) {
                    qh_fprintf(qh ferr, 8136, "%s OTHER FACET:\n", string);
                    qh_printfacet(qh ferr, otherfacet);
                  }
                  if (atridge) {
                    qh_fprintf(qh ferr, 8137, "%s RIDGE:\n", string);
                    qh_printridge(qh ferr, atridge);
                    if (atridge->top && atridge->top != atfacet && atridge->top != otherfacet)
                      qh_printfacet(qh ferr, atridge->top);
                    if (atridge->bottom
                        && atridge->bottom != atfacet && atridge->bottom != otherfacet)
                      qh_printfacet(qh ferr, atridge->bottom);
                    if (!atfacet)
                      atfacet= atridge->top;
                    if (!otherfacet)
                      otherfacet= otherfacet_(atridge, atfacet);
                  }
                  if (atvertex) {
                    qh_fprintf(qh ferr, 8138, "%s VERTEX:\n", string);
                    qh_printvertex(qh ferr, atvertex);
                  }
                  if (qh fout && qh FORCEoutput && atfacet && !qh QHULLfinished && !qh IStracing) {
                    qh_fprintf(qh ferr, 8139, "ERRONEOUS and NEIGHBORING FACETS to output\n");
                    for (i=0; i < qh_PRINTEND; i++)  /* use fout for geomview output */
                      qh_printneighborhood(qh fout, qh PRINTout[i], atfacet, otherfacet,
                                            !qh_ALL);
                  }
                }

                ////////////////////////////////////////////////////////////////////
                /// qh_printhelp_degenerate
                ////////////////////////////////////////////////////////////////////
                void qh_printhelp_degenerate(FILE *fp) {

                  if (qh MERGEexact || qh PREmerge || qh JOGGLEmax < REALmax/2)
                    qh_fprintf(fp, 9368, "\n\
                A Qhull error has occurred.  Qhull should have corrected the above\n\
                precision error.  Please send the input and all of the output to\n\
                qhull_bug@qhull.org\n");
                  else if (!qh_QUICKhelp) {
                    qh_fprintf(fp, 9369, "\n\
                Precision problems were detected during construction of the convex hull.\n\
                This occurs because convex hull algorithms assume that calculations are\n\
                exact, but floating-point arithmetic has roundoff errors.\n\
                \n\
                To correct for precision problems, do not use 'Q0'.  By default, Qhull\n\
                selects 'C-0' or 'Qx' and merges non-convex facets.  With option 'QJ',\n\
                Qhull joggles the input to prevent precision problems.  See \"Imprecision\n\
                in Qhull\" (qh-impre.htm).\n\
                \n\
                If you use 'Q0', the output may include\n\
                coplanar ridges, concave ridges, and flipped facets.  In 4-d and higher,\n\
                Qhull may produce a ridge with four neighbors or two facets with the same \n\
                vertices.  Qhull reports these events when they occur.  It stops when a\n\
                concave ridge, flipped facet, or duplicate facet occurs.\n");
                #if REALfloat
                    qh_fprintf(fp, 9370, "\
                \n\
                Qhull is currently using single precision arithmetic.  The following\n\
                will probably remove the precision problems:\n\
                  - recompile qhull for realT precision(#define REALfloat 0 in user.h).\n");
                #endif
                    if (qh DELAUNAY && !qh SCALElast && qh MAXabs_coord > 1e4)
                      qh_fprintf(fp, 9371, "\
                \n\
                When computing the Delaunay triangulation of coordinates > 1.0,\n\
                  - use 'Qbb' to scale the last coordinate to [0,m] (max previous coordinate)\n");
                    if (qh DELAUNAY && !qh ATinfinity)
                      qh_fprintf(fp, 9372, "\
                When computing the Delaunay triangulation:\n\
                  - use 'Qz' to add a point at-infinity.  This reduces precision problems.\n");

                    qh_fprintf(fp, 9373, "\
                \n\
                If you need triangular output:\n\
                  - use option 'Qt' to triangulate the output\n\
                  - use option 'QJ' to joggle the input points and remove precision errors\n\
                  - use option 'Ft'.  It triangulates non-simplicial facets with added points.\n\
                \n\
                If you must use 'Q0',\n\
                try one or more of the following options.  They can not guarantee an output.\n\
                  - use 'QbB' to scale the input to a cube.\n\
                  - use 'Po' to produce output and prevent partitioning for flipped facets\n\
                  - use 'V0' to set min. distance to visible facet as 0 instead of roundoff\n\
                  - use 'En' to specify a maximum roundoff error less than %2.2g.\n\
                  - options 'Qf', 'Qbb', and 'QR0' may also help\n",
                               qh DISTround);
                    qh_fprintf(fp, 9374, "\
                \n\
                To guarantee simplicial output:\n\
                  - use option 'Qt' to triangulate the output\n\
                  - use option 'QJ' to joggle the input points and remove precision errors\n\
                  - use option 'Ft' to triangulate the output by adding points\n\
                  - use exact arithmetic (see \"Imprecision in Qhull\", qh-impre.htm)\n\
                ");
                  }
                }


                ////////////////////////////////////////////////////////////////////
                /// qh_printhelp_narrowhull
                ////////////////////////////////////////////////////////////////////
                void qh_printhelp_narrowhull(FILE *fp, realT minangle) {

                    qh_fprintf(fp, 9375, "qhull precision warning: \n\
                The initial hull is narrow (cosine of min. angle is %.16f).\n\
                Is the input lower dimensional (e.g., on a plane in 3-d)?  Qhull may\n\
                produce a wide facet.  Options 'QbB' (scale to unit box) or 'Qbb' (scale\n\
                last coordinate) may remove this warning.  Use 'Pp' to skip this warning.\n\
                See 'Limitations' in qh-impre.htm.\n",
                          -minangle); // convert from angle between normals to angle between facets
                }


                ////////////////////////////////////////////////////////////////////
                /// qh_printhelp_singular
                ////////////////////////////////////////////////////////////////////
                void qh_printhelp_singular(FILE *fp) {
                  facetT *facet;
                  vertexT *vertex, **vertexp;
                  realT min, max, *coord, dist;
                  int i,k;

                  qh_fprintf(fp, 9376, "\n\
                The input to qhull appears to be less than %d dimensional, or a\n\
                computation has overflowed.\n\n\
                Qhull could not construct a clearly convex simplex from points:\n",
                           qh hull_dim);
                  qh_printvertexlist(fp, "", qh facet_list, NULL, qh_ALL);
                  if (!qh_QUICKhelp)
                    qh_fprintf(fp, 9377, "\n\
                The center point is coplanar with a facet, or a vertex is coplanar\n\
                with a neighboring facet.  The maximum round off error for\n\
                computing distances is %2.2g.  The center point, facets and distances\n\
                to the center point are as follows:\n\n", qh DISTround);
                  qh_printpointid(fp, "center point", qh hull_dim, qh interior_point, qh_IDunknown);
                  qh_fprintf(fp, 9378, "\n");
                  FORALLfacets {
                    qh_fprintf(fp, 9379, "facet");
                    FOREACHvertex_(facet->vertices)
                      qh_fprintf(fp, 9380, " p%d", qh_pointid(vertex->point));
                    zinc_(Zdistio);
                    qh_distplane(qh interior_point, facet, &dist);
                    qh_fprintf(fp, 9381, " distance= %4.2g\n", dist);
                  }
                  if (!qh_QUICKhelp) {
                    if (qh HALFspace)
                      qh_fprintf(fp, 9382, "\n\
                These points are the dual of the given halfspaces.  They indicate that\n\
                the intersection is degenerate.\n");
                    qh_fprintf(fp, 9383,"\n\
                These points either have a maximum or minimum x-coordinate, or\n\
                they maximize the determinant for k coordinates.  Trial points\n\
                are first selected from points that maximize a coordinate.\n");
                    if (qh hull_dim >= qh_INITIALmax)
                      qh_fprintf(fp, 9384, "\n\
                Because of the high dimension, the min x-coordinate and max-coordinate\n\
                points are used if the determinant is non-zero.  Option 'Qs' will\n\
                do a better, though much slower, job.  Instead of 'Qs', you can change\n\
                the points by randomly rotating the input with 'QR0'.\n");
                  }
                  qh_fprintf(fp, 9385, "\nThe min and max coordinates for each dimension are:\n");
                  for (k=0; k < qh hull_dim; k++) {
                    min= REALmax;
                    max= -REALmin;
                    for (i=qh num_points, coord= qh first_point+k; i--; coord += qh hull_dim) {
                      maximize_(max, *coord);
                      minimize_(min, *coord);
                    }
                    qh_fprintf(fp, 9386, "  %d:  %8.4g  %8.4g  difference= %4.4g\n", k, min, max, max-min);
                  }
                  if (!qh_QUICKhelp) {
                    qh_fprintf(fp, 9387, "\n\
                If the input should be full dimensional, you have several options that\n\
                may determine an initial simplex:\n\
                  - use 'QJ'  to joggle the input and make it full dimensional\n\
                  - use 'QbB' to scale the points to the unit cube\n\
                  - use 'QR0' to randomly rotate the input for different maximum points\n\
                  - use 'Qs'  to search all points for the initial simplex\n\
                  - use 'En'  to specify a maximum roundoff error less than %2.2g.\n\
                  - trace execution with 'T3' to see the determinant for each point.\n",
                                     qh DISTround);
                #if REALfloat
                    qh_fprintf(fp, 9388, "\
                  - recompile qhull for realT precision(#define REALfloat 0 in libqhull.h).\n");
                #endif
                    qh_fprintf(fp, 9389, "\n\
                If the input is lower dimensional:\n\
                  - use 'QJ' to joggle the input and make it full dimensional\n\
                  - use 'Qbk:0Bk:0' to delete coordinate k from the input.  You should\n\
                    pick the coordinate with the least range.  The hull will have the\n\
                    correct topology.\n\
                  - determine the flat containing the points, rotate the points\n\
                    into a coordinate plane, and delete the other coordinates.\n\
                  - add one or more points to make the input full dimensional.\n\
                ");
                  }
                }

                void erase_face
                (
                	Grid& g,
                    SubsetHandler& sh,
                    size_t si,
                    Grid::VertexAttachmentAccessor<APosition>& aaPos,
                    std::vector<ug::vector3>& verts
                )
                {
                	Selector sel(g);
                	for (TriangleIterator iter = sh.begin<ug::Triangle>(si+1000); iter != sh.end<ug::Triangle>(si+1000); ++iter) {
                		UG_LOGN("triangle!");
                		bool atEnd = true;
                		for (size_t i = 0; i < (*iter)->num_vertices(); i++) {
                			bool current = std::find(verts.begin(), verts.end(), aaPos[(*iter)->vertex(i)]) != verts.end();
                			atEnd = atEnd && current;
                			UG_LOGN("current: " << current);
                		}

                		if (atEnd) {
                			UG_LOGN("Erasing face!");
                			sel.select(*iter);
                		}
                	}
                	SelectAssociatedEdges(sel, sel.begin<ug::Triangle>(), sel.end<ug::Triangle>(), true);
                	EraseSelectedObjects(sel);
                }

                void gen_face
                (
                        const std::vector<ug::vector3>& points,
                        Grid& g,
                        SubsetHandler& sh,
                        size_t si,
                        Grid::VertexAttachmentAccessor<APosition>& aaPos
                )
                {
                     int numPoints= points.size();
                     int dim = 3;
                     coordT ps[numPoints*dim];
                     for (size_t i = 0; i < numPoints; i++) {
                             for (size_t j = 0; j < dim; j++) {
                                     ps[i*dim+j] = points[i][j];
                             }
                     }

                    char flags[25];
                    sprintf (flags, "qhull s FA Qt");
                    qh_new_qhull(dim, numPoints, ps, 0, flags, NULL, NULL);
                    facetT* facet;
                    for (facet=qh facet_list;facet && facet->next;facet=facet->next) {
                        setT *vertices = facet->vertices;
                        vertexT *vertex, **vertexp;
                        int count = qh_setsize(vertices);
                        std::vector<ug::vector3> verts;
                        FOREACHvertex_(facet->vertices) {
                          ug::vector3 temp;
                          pointT *point = vertex->point;
                          int i;
                          for (i = 0; i < dim; i++) {
                        	  double coord = vertex->point[i];
                              temp[i] = coord;
                          }
                          verts.push_back(temp);
                        }

                        Vertex* v1 = *g.create<RegularVertex>();
                        Vertex* v2 = *g.create<RegularVertex>();
                        Vertex* v3 = *g.create<RegularVertex>();
                        aaPos[v1] = verts[0];
                        aaPos[v2] = verts[1];
                        aaPos[v3] = verts[2];
                        Triangle* t2 = *g.create<Triangle>(TriangleDescriptor(v1, v2, v3));
                        Selector sel(g);
                        sel.select(t2);
                        sel.select(v1);
                        sel.select(v2);
                        sel.select(v3);
                        CloseSelection(sel);
                        AssignSelectionToSubset(sel, sh, si + 1000);
                        sel.clear();
                  }
                 qh_freeqhull(!qh_ALL);
             }
		}
	}
}