//
// Created by yuyttenh on 02/07/2021.
//

#ifndef SWIFTSIM_DELAUNAY_FUNCTIONS_H
#define SWIFTSIM_DELAUNAY_FUNCTIONS_H

#include "delaunay.h"

#if defined(HYDRO_DIMENSION_3D)
#include "voronoi_3d/delaunay_functions_3d.h"
#elif defined(HYDRO_DIMENSION_2D)
#include "voronoi_2d/delaunay_functions_2d.h"
#else
#error "Only 2D or 3D moving mesh are supported!"
#endif

#endif  // SWIFTSIM_DELAUNAY_FUNCTIONS_H
