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

#include <part.h>

/**
 * @brief Voronoi interface.
 *
 * An interface is a connection between two neighbouring Voronoi cells. It is
 * completely defined by the indices of the generators that generate the two
 * neighbouring cells, a surface area and a midpoint position.
 */
struct voronoi_pair {
  /*! Pointer to particle corresponding to the generator on the left of the
   * interface (always a particle within the local cell). */
  struct part *left;

  /*! idx of cell on the left in this voronoi tesselation. */
  int left_idx;

  /*! Pointer to particle corresponding to the generator on the right of the
   * interface (can be a local particle, but also a particle in a
   * neighbouring cell). */
  struct part *right;

  /*! idx of cell on the right in this voronoi tesselation If it is not a ghost
   * particle, else -1. */
  int right_idx;

  struct cell *right_cell;

  /*! Surface area of the interface. */
  double surface_area;

  /*! Midpoint of the interface. */
  double midpoint[3];

#ifdef VORONOI_STORE_CONNECTIONS
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
  struct voronoi_pair *pairs[27];

  /*! @brief Current number of pairs per cell index. */
  int pair_index[27];

  /*! @brief Allocated number of pairs per cell index. */
  int pair_size[27];

  /*! @brief Flag indicating whether this voronoi struct is active (has memory
   * allocated)
   */
  int active;

  /*! @brief The absolute minimal surface area of faces in this voronoi
   * tessellation */
  double min_surface_area;

  /*! Pointer to swift cell containing this tesselation */
  struct cell *swift_cell;
};

#endif  // SWIFTSIM_VORONOI_3D_H
