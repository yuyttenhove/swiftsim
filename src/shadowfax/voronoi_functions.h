//
// Created by yuyttenh on 02/07/2021.
//

#ifndef SWIFTSIM_VORONOI_FUNCTIONS_H
#define SWIFTSIM_VORONOI_FUNCTIONS_H

#include "voronoi.h"

#if defined(HYDRO_DIMENSION_3D)
#include "voronoi_3d/voronoi_functions_3d.h"
#elif defined(HYDRO_DIMENSION_2D)
#include "voronoi_2d/voronoi_functions_2d.h"
#else
#error "Only 2D or 3D moving mesh are supported!"
#endif

#endif  // SWIFTSIM_VORONOI_FUNCTIONS_H
