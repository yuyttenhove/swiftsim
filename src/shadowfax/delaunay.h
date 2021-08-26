//
// Created by yuyttenh on 02/07/2021.
//

#ifndef SWIFTSIM_DELAUNAY_H
#define SWIFTSIM_DELAUNAY_H

#include <stdio.h>
#include <stdlib.h>

/*! @brief Activate extensive log output. */
/*#define DELAUNAY_LOG_OUTPUT*/
/*! @brief Activate runtime assertions. */
//#define DELAUNAY_DO_ASSERTIONS
/*! @brief Use and output non-exact floating point geometrical tests as well as
 *  the default exact integer tests. This is especially helpful when trying to
 *  visualise the geometry, since the integer coordinates are very hard to
 *  interpret. */
#define DELAUNAY_NONEXACT
/*! @brief Check the validity of the Delaunay tessellation after every addition
 *  of a new vertex. This feature is very helpful when debugging to catch
 *  problems as they happen, but adds a very significant runtime cost. It should
 *  never be activated for production runs! */
//#define DELAUNAY_CHECKS

/**
 * @brief Print the given message to the standard output.
 *
 * This macro behaves the same as printf(), but prepends a file name, function
 * name and line number to each output line.
 *
 * This macro is only defined when DELAUNAY_LOG_OUTPUT is active.
 */
#ifdef DELAUNAY_LOG_OUTPUT
#define delaunay_log(s, ...)                                      \
  printf("%s:%s():%i: " s "\n", __FILE__, __FUNCTION__, __LINE__, \
         ##__VA_ARGS__);
#else
#define delaunay_log(s, ...)
#endif

/**
 *@brief Evaluate the given condition and abort if it evaluates to true.
 *
 * This macro is similar to the standard assert() macro.
 *
 * This macro is only defined when DELAUNAY_DO_ASSERTIONS is active.
 */
#ifdef DELAUNAY_DO_ASSERTIONS
#define delaunay_assert(condition)                                    \
  if (!(condition)) {                                                 \
    fprintf(stderr, "%s:%s():%i: Condition failed: " #condition "\n", \
            __FILE__, __FUNCTION__, __LINE__);                        \
    abort();                                                          \
  }
#else
#define delaunay_assert(condition)
#endif

/**
 * @brief Convert the given double precision floating point value to an integer,
 * by reading out its 52-bit mantissa.
 *
 * A floating point variable consists of a mantissa and an exponent, and can be
 * thought of as the base 2 equivalent of scientific notation:
 * @f[
 *    V = M \times[} 2^E
 * @f]
 * The sign of the mantissa (highest bit of the mantissa) determines the sign
 * of the floating point value.
 *
 * This code was taken from the AREPO-code with some small adaptations.
 *
 * @param d Input double precision floating point value.
 * @return Integer value of the 52-bit mantissa.
 */
static inline unsigned long int delaunay_double_to_int(double d) {
  /* the idea here is pretty simple: we set up a union consisting of a 64-bit
     double precision floating point value and a 64-bit unsigned long integer
     that occupy the same 64-bits in memory.
     We then copy the value we want to convert into the double precision
     variable and access its individual bits through the 64-bit unsigned long
     integer variable. */
  union {
    double d;
    unsigned long int ull;
  } u;
  u.d = d;
  /* the mask filters out the lowest 52 bits of the binary sequence, which
     correspond to the mantissa of the floating point variable */
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}

#if defined(HYDRO_DIMENSION_3D)
#include "algorithm_3d/delaunay_3d.h"
#elif defined(HYDRO_DIMENSION_2D)
#include "algorithm_2d/delaunay_2d.h"
#else
#error "Only 2D or 3D moving mesh are supported!"
#endif


#endif  // SWIFTSIM_DELAUNAY_H
