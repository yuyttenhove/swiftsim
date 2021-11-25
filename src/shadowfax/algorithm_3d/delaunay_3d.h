//
// Created by yuyttenh on 02/07/2021.
//

/**
 * @file delaunay3d.h
 *
 * @brief 3D Delaunay struct and construction algorithm.
 */

#ifndef SWIFTSIM_DELAUNAY3D_H
#define SWIFTSIM_DELAUNAY3D_H

#include "geometry_3d.h"
#include "src/shadowfax/queues.h"
#include "tetrahedron.h"

struct delaunay {

  /* Activity flag, useful for debugging. */
  int active;

  /*! @brief Anchor of the simulation volume. */
  double anchor[3];

  /*! @brief Side length of the simulation volume. */
  double side;

  /*! @brief Inverse side length of the simulation volume. */
  double inverse_side;

  /*! @brief Vertex positions. This array is a copy of the array defined in
   *  main() and we probably want to get rid of it in a SWIFT implementation. */
  double* vertices;

  /*! @brief Vertex positions, rescaled to the range 1-2. Kept here in case we
   *  want to adopt hybrid geometrical checks (floating point checks when safe,
   *  integer checks when there is a risk of numerical error leading to the
   *  wrong result) to speed things up. */
  double* rescaled_vertices;

  /*! @brief Integer vertex_indices. These are the vertex coordinates that are
   *  actually used during the incremental construction. */
  unsigned long int* integer_vertices;

  /*! @brief Vertex-tetrahedron connections. For every vertex in the
   * tessellation, this array stores the index of a tetrahedron that contains
   * this vertex (usually one of the last tetrahedra that was constructed with
   * this vertex as a member vertex). This array is not required for the
   * incremental construction algorithm itself, but is indispensable for the
   * conversion from Delaunay tessellation to Voronoi grid, since it links each
   * input vertex to one of the tetrahedra it is connected to. Without this
   * array, we would have no efficient way of finding a tetrahedron that
   * contains a given vertex. */
  int* vertex_tetrahedron_links;

  /*! @brief Vertex-tetrahedron connection indices. For every vertex-tetrahedron
   * pair stored in vertex_tetrahedron_links, this array contains the index of
   * the vertex within the vertex list of the tetrahedron. This saves us from
   * having to loop through the tetrahedron vertex_indices during Voronoi grid
   * construction. */
  int* vertex_tetrahedron_index;

  /*! @brief Vertex search radii. For every vertex, this array contains twice
   *  the radius of the largest circumsphere of the tetrahedra that vertex is
   *  part of. */
  double* search_radii;

  /*! @brief Direct pointers to particles corresponding to vertices */
  struct part** part_pointers;

  /*! @brief Next available index within the vertex array. Corresponds to the
   *  actual size of the vertex array. */
  int vertex_index;

  /*! @brief Current size of the vertex array in memory. If vertex_size matches
   *  n_vertices, the memory buffer is full and needs to be expanded. */
  int vertex_size;

  /*! @brief Begin index of the normal vertex_indices. This skips the 3
   * auxiliary vertex_indices required for the incremental construction
   * algorithm. */
  int vertex_start;

  /*! @brief End index of the normal vertex_indices. This variable is set by
   * calling delaunay_consolidate() and contains the offset of the ghost
   * vertex_indices within the vertex array. */
  int vertex_end;

  /*! @brief Tetrahedra that make up the tessellation. */
  struct tetrahedron* tetrahedra;

  /*! @brief Next available index within the tetrahedron array. Corresponds to
   * the actual size of the tetrahedron array. */
  int tetrahedron_index;

  /*! @brief Current size of the tetrahedron array in memory. If
   * tetrahedron_size matches tetrahedron_index, the memory buffer is full and
   * needs to be expanded. */
  int tetrahedron_size;

  /*! @brief Index of the last tetrahedron that was created or modified. Used as
   * initial guess for the tetrahedron that contains the next vertex that will
   * be added. If vertex_indices are added in some sensible order (e.g. in
   * Peano-Hilbert curve order) then this will greatly speed up the algorithm.
   */
  int last_tetrahedron;

  /*! @brief Lifo queue of tetrahedra that need checking during the incremental
   *  construction algorithm. After a new vertex has been added, all new
   *  tetrahedra are added to this queue and are tested to see if the
   *  Delaunay criterion (empty circumcircles) still holds. New tetrahedra
   *  created when flipping invalid tetrahedra are also added to the queue. */
  struct int_lifo_queue tetrahedra_to_check;

  /*! @brief Lifo queue of free spots in the tetrahedra array. Sometimes 3
   * tetrahedra can be split into 2 new ones. This leaves a free spot in the
   * array.
   */
  struct int_lifo_queue free_tetrahedron_indices;

  /*! @brief Array of tetrahedra containing the current vertex */
  struct int_lifo_queue tetrahedra_containing_vertex;

  /*! @brief Queue to store neighbouring vertices of the current vertex when
   * looping over all the tetrahedra containing the current vertex to calculate
   * its search radius */
  struct int3_fifo_queue get_radius_neighbour_info_queue;

  /*! @brief Array to indicate which neighbours have already been added to the
   * get_radius_neighbour_info_queue during the search radius calculation. */
  int* get_radius_neighbour_flags;

  /*! @brief Geometry variables. Auxiliary variables used by the exact integer
   *  geometry3d tests that need to be stored in between tests, since allocating
   *  and deallocating them for every test is too expensive. */
  struct geometry3d geometry;

  /*! @brief Cell neighbour sids keeping track of which neighbouring cells
   *  contains a specific neighbouring vertex. */
  int* ngb_cell_sids;

  /*! @brief Pointers to neighbour cells containing corresponding particles */
  struct cell** ngb_cell_ptrs;

  /*! @brief Current used size of the neighbouring vertex bookkeeping arrays
   *  (and next valid index in this array). */
  int ngb_index;

  /*! @brief Current size in memory of the neighbouring vertex bookkeeping
   *  arrays. More memory needs to be allocated if ngb_index reaches this
   *  value. */
  int ngb_size;

  /*! @brief Offset of the neighbouring vertices within the vertex array (so
   *  that neighbouring information for vertex v is stored in
   *  ngb_cell_sids[v-ngb_offset]). */
  int ngb_offset;

  /*! @brief Array of booleans indicating whether or not neighbouring particles
   * have been tried to be added for a given sid. If this is 0 for a given sid,
   * this means that this cell should get the reflective boundary condition
   * applied for that sid. */
  int sid_is_inside_face[27];
};

#endif  // SWIFTSIM_DELAUNAY3D_H
