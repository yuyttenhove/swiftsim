//
// Created by yuyttenh on 02/07/2021.
//

/**
 * @file voronoi3d.h
 *
 * @brief 3D Voronoi grid and functions.
 */

#ifndef SWIFTSIM_VORONOI_3D_H
#define SWIFTSIM_VORONOI_3D_H

#include "part.h"
#include "shadowfax/queues.h"

/**
 * @brief Voronoi interface.
 *
 * An interface is a connection between two neighbouring Voronoi cells. It is
 * completely defined by the indices of the generators that generate the two
 * neighbouring cells, a surface area and a midpoint position.
 */
struct voronoi_pair {
  /*! idx of the particle on the left of this pair in its respective swift
   * cell. Since the left particle is always local this is also the index of the
   * corresponding cell in this voronoi tesselation. */
  int left_idx;

  /*! idx of the particle on the right of this pair in its respective swift cell
   * if that cell is the same as the cell holding this Voronoi tesselation (i.e.
   * the particle is local) or in the super cell of its respective swift cell if
   * that swift cell is foreign. For local particles, this is also the index of
   * the corresponding cell in this voronoi tesselation. */
  int right_idx;

  int right_del_idx;

  struct cell *right_cell;

  int right_nodeID;

  /*! Surface area of the interface. */
  double surface_area;

  /*! Midpoint of the interface. */
  double midpoint[3];

#ifdef VORONOI_STORE_FACES
  /*! Vertices of the interface. */
  double *vertices;

  /*! Number of vertices of this face. */
  int n_vertices;
#endif
};

struct voronoi {
  /*! @brief Voronoi cells. */
  struct voronoi_cell_new *cells;

  /*! @brief Number of cells. */
  int number_of_cells;

  /*! @brief The allocated size of the cells array */
  int cells_size;

  /*! @brief Voronoi cell pairs. We store these per (SWIFT) cell, i.e. pairs[0]
   *  contains all pairs that are completely contained within this cell, while
   *  pairs[1] corresponds to pairs crossing the boundary between this cell and
   *  the cell with coordinates that are lower in all coordinate directions (the
   *  cell to the left, front, bottom, sid=0), and so on. */
  struct voronoi_pair *pairs[28];

  /*! @brief Current number of pairs per cell index. */
  int pair_index[28];

  /*! @brief Allocated number of pairs per cell index. */
  int pair_size[28];

  /*! @brief cell pair connection. Queue of 2-tuples containing the index of
   * the pair and the sid of the pair */
  struct int2_lifo_queue cell_pair_connections;

  /*! @brief Flag indicating whether this voronoi struct is active (has memory
   * allocated) */
  int active;

  /*! @brief The absolute minimal surface area of faces in this voronoi
   * tessellation */
  double min_surface_area;

#ifdef VORONOI_STORE_FACES
  /*! @brief The face vertices of all the faces concatenated in one big array.
   * Every face has a pointer to the start of its vertices */
   double *face_vertices;
   int face_vertices_size;
   int face_vertices_index;
#endif

  /*! Pointer to swift cell containing this tesselation */
  struct cell *swift_cell;
};

#endif  // SWIFTSIM_VORONOI_3D_H
