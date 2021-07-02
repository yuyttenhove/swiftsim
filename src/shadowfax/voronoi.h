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

#if defined(HYDRO_DIMENSION_3D)
#include "voronoi_3d/voronoi_3d.h"
#elif defined(HYDRO_DIMENSION_2D)
#include "voronoi_2d/voronoi_2d.h"
#else
#error "Only 2D or 3D moving mesh are supported!"
#endif

#endif  // SWIFTSIM_VORONOI_H
