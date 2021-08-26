//
// Created by yuyttenh on 26/08/2021.
//

#ifndef SWIFTSIM_UTILS_H
#define SWIFTSIM_UTILS_H

static inline int sgn(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

#endif  // SWIFTSIM_UTILS_H
