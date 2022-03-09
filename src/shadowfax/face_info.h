//
// Created by yuyttenh on 8/03/22.
//

#ifndef SWIFTSIM_FACE_INFO_H
#define SWIFTSIM_FACE_INFO_H

#include "voronoi.h"

struct face_info {
  int face_counts[27];
  struct voronoi_pair faces[];
};

#endif  // SWIFTSIM_FACE_INFO_H
