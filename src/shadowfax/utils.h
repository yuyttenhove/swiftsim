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

/**
 * @brief Check whether two doubles are equal up to the given precision.
 *
 * @param double1
 * @param double2
 * @param precision
 * @return 1 for equality 0 else.
 */
inline static int double_cmp(double double1, double double2,
                             unsigned long precision) {
  long long1, long2;
  if (double1 > 0)
    long1 = (long)(double1 * precision + .5);
  else
    long1 = (long)(double1 * precision - .5);
  if (double2 > 0)
    long2 = (long)(double2 * precision + .5);
  else
    long2 = (long)(double2 * precision - .5);
  return (long1 == long2);
}

#endif  // SWIFTSIM_UTILS_H
