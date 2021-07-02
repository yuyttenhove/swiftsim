//
// Created by yuyttenh on 02/07/2021.
//

#ifndef SWIFTSIM_VORONOI_H
#define SWIFTSIM_VORONOI_H

/*! @brief Store the edges of faces (so that the actual Voronoi grid can be
 *  reconstructed). */
#define VORONOI_STORE_CONNECTIONS

/*! @brief Store information about the number of faces per cell. */
#define VORONOI_STORE_CELL_STATS

/*! @brief Store cell generators. */
#define VORONOI_STORE_GENERATORS

/*! @brief Activate runtime assertions. */
#define VORONOI_DO_ASSERTIONS

#define VORONOI_CHECKS

/**
 *@brief Evaluate the given condition and abort if it evaluates to true.
 *
 * This macro is similar to the standard assert() macro.
 * This macro is only defined when VORONOI_DO_ASSERTIONS is active.
 */
#ifdef VORONOI_DO_ASSERTIONS
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

#if defined(HYDRO_DIMENSION_3D)
#include "voronoi_3d/voronoi_3d.h"
#elif defined(HYDRO_DIMENSION_2D)
#include "voronoi_2d/voronoi_2d.h"
#else
#error "Only 2D or 3D moving mesh are supported!"
#endif

#endif  // SWIFTSIM_VORONOI_H
