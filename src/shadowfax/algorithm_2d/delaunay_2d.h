/*******************************************************************************
 * This file is part of cVoronoi.
 * Copyright (c) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/**
 * @file delaunay.h
 *
 * @brief Delaunay tessellation.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#ifndef SWIFT_DELAUNAY_STRUCT_H
#define SWIFT_DELAUNAY_STRUCT_H

#include "geometry_2d.h"
#include "triangle.h"

/**
 * @brief Delaunay tessellation.
 *
 * The tessellation contains all the triangles that make up the tessellation,
 * their connectivity is stored implicitly within the triangles themselves.
 */
struct delaunay {

  /* Activity flag, useful for debugging. */
  int active;

  /* Box geometry: used to set up the initial vertices and triangles and to
   * convert input coordinates to integer coordinates. */

  /*! @brief Anchor of the simulation volume. */
  double anchor[2];

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

  /*! @brief Integer vertices. These are the vertex coordinates that are
   *  actually used during the incremental construction. */
  unsigned long int* integer_vertices;

  /*! @brief Vertex-triangle connections. For every vertex in the tessellation,
   *  this array stores the index of a triangle that contains this vertex
   *  (usually one of the last triangles that was constructed with this vertex
   *  as a member vertex). This array is not required for the incremental
   *  construction algorithm itself, but is indispensable for the conversion
   *  from Delaunay tessellation to Voronoi grid, since it links each input
   *  vertex to one of the triangles it is connected to. Without this array, we
   *  would have no efficient way of finding a triangle that contains a given
   *  vertex. */
  int* vertex_triangles;

  /*! @brief Vertex-triangle connection indices. For every vertex-triangle pair
   *  stored in vertex_triangles, this array contains the index of the vertex
   *  within the vertex list of the triangle. This saves us from having to
   *  loop through the triangle vertices during Voronoi grid construction. */
  int* vertex_triangle_index;

  /*! @brief Vertex search radii. For every vertex, this array contains twice
   *  the radius of the largest circumcircle of the triangles that vertex is
   *  part of. */
  double* search_radii;

  /*! @brief Direct pointers to particles corresponding to vertices */
  struct part** part_pointers;

  /*! @brief Next available index within the vertex array. Corresponds to the
   *  actual size of the vertex array. */
  int vertex_index;

  /*! @brief Current size of the vertex array in memory. If vertex_size matches
   *  vertex_index, the memory buffer is full and needs to be expanded. */
  int vertex_size;

  /*! @brief Begin index of the normal vertices. This skips the 3 auxiliary
   *  vertices required for the incremental construction algorithm. */
  int vertex_start;

  /*! @brief End index of the normal vertices. This variable is set by calling
   *  delaunay_consolidate() and contains the offset of the ghost vertices
   *  within the vertex array. */
  int vertex_end;

  /*! @brief Triangles that make up the tessellation. */
  struct triangle* triangles;

  /*! @brief Next available index within the triangle array. Corresponds to the
   *  actual size of the triangle array. */
  int triangle_index;

  /*! @brief Current size of the triangle array in memory. If triangle_size
   *  matches triangle_index, the memory buffer is full and needs to be
   *  expanded. */
  int triangle_size;

  /*! @brief Queue of triangles that need checking during the incremental
   *  construction algorithm. After a new vertex has been added, all new
   *  triangles are added to this queue and are tested to see if the
   *  Delaunay criterion (empty circumcircles) still holds. New triangles
   *  created when flipping invalid triangles are also added to the queue. */
  int* queue;

  /*! @brief Next available index in the queue. Determines both the actual size
   *  of the queue as the first element that will be popped from the queue. */
  int queue_index;

  /*! @brief Current size of the queue in memory. If queue_size matches
   *  queue_index, the memory buffer is full and needs to be expanded. */
  int queue_size;

  /*! @brief Index of the last triangle that was accessed. Used as initial
   *  guess for the triangle that contains the next vertex that will be added.
   *  If vertices are added in some sensible order (e.g. in Peano-Hilbert curve
   *  order) then this will greatly speed up the algorithm. */
  int last_triangle;

  /*! @brief Geometry variables. Auxiliary variables used by the exact integer
   *  geometry tests that need to be stored in between tests, since allocating
   *  and deallocating them for every test is too expensive. */
  struct geometry2d geometry;

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

#endif /* SWIFT_DELAUNAY_STRUCT_H */
