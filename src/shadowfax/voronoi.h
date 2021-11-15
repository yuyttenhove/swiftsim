//
// Created by yuyttenh on 02/07/2021.
//

#ifndef SWIFTSIM_VORONOI_H
#define SWIFTSIM_VORONOI_H

/*! @brief Store the edges of faces (so that the actual Voronoi grid can be
 *  reconstructed). */
#define VORONOI_STORE_CONNECTIONS

/*! @brief Store cell generators. */
#define VORONOI_STORE_GENERATORS

/*! @brief Activate runtime assertions. */
#define VORONOI_DO_ASSERTIONS

/*! @brief Activate extra checks */
//#define VORONOI_CHECKS

/**
 *@brief Evaluate the given condition and abort if it evaluates to true.
 *
 * This macro is similar to the standard assert() macro.
 * This macro is only defined when VORONOI_DO_ASSERTIONS is active.
 */
#if defined(VORONOI_DO_ASSERTIONS) || defined(VORONOI_CHECKS)
#define voronoi_assert(condition)                                     \
  if (!(condition)) {                                                 \
    fprintf(stderr, "%s:%s():%i: Condition failed: " #condition "\n", \
            __FILE__, __FUNCTION__, __LINE__);                        \
    abort();                                                          \
  }
#else
#define voronoi_assert(condition)
#endif

#define voronoi_error(s, ...) \
  fprintf(stderr, s, ##__VA_ARGS__); \
  abort();

/**
 * @brief minimal 1d relative size of voronoi faces
 */
extern double min_rel_voronoi_face_size;

/**
 * @brief Voronoi cell.
 *
 * A cell stores geometrical information about a Voronoi cell: its volume and
 * the location of its centroid (for compatibility reasons, these must always be
 * 3D vectors (also for the 2D voronoi algorithm).
 */
struct voronoi_cell_new {
  /*! Cell volume. */
  double volume;

  /*! Cell centroid. */
  double centroid[3];

#ifdef VORONOI_STORE_GENERATORS
  /*! Position of the cell generator. */
  struct part *generator;
#endif

  /*! Number of faces of this cell. */
  int nface;

#ifdef VORONOI_STORE_CONNECTIONS
  /*! cell_pair_connections offset */
  int pair_connections_offset;
#endif
};

#if defined(HYDRO_DIMENSION_3D)
#include "algorithm_3d/voronoi_3d.h"
#elif defined(HYDRO_DIMENSION_2D)
#include "algorithm_2d/voronoi_2d.h"
#else
#error "Only 2D or 3D moving mesh are supported!"
#endif

#endif  // SWIFTSIM_VORONOI_H
