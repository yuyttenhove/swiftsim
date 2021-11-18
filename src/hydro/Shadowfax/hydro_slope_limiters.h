//
// Created by yuyttenh on 26/04/2021.
//

#ifndef SWIFT_HYDRO_SLOPE_LIMITERS_H
#define SWIFT_HYDRO_SLOPE_LIMITERS_H

#include "dimension.h"
#include "kernel_hydro.h"

#ifdef SHADOWFAX_SLOPE_LIMITER_PER_FACE

#define HYDRO_SLOPE_LIMITER_FACE_IMPLEMENTATION \
  "GIZMO piecewise slope limiter (Hopkins 2015)"
#include "hydro_slope_limiters_face.h"

#else

#define HYDRO_SLOPE_LIMITER_FACE_IMPLEMENTATION "No piecewise slope limiter"

/**
 * @brief Slope limit the slopes at the interface between two particles
 *
 * @param Wi Hydrodynamic variables of particle i.
 * @param Wj Hydrodynamic variables of particle j.
 * @param dWi Difference between the hydrodynamic variables of particle i at the
 * position of particle i and at the interface position.
 * @param dWj Difference between the hydrodynamic variables of particle j at the
 * position of particle j and at the interface position.
 * @param xij_i Relative position vector of the interface w.r.t. particle i.
 * @param xij_j Relative position vector of the interface w.r.t. partilce j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_face(
    double *Wi, double *Wj, double *dWi, double *dWj, const double *xij_i, const double *xij_j,
    double r) {}

#endif

#if defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE) || defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE_EXACT)

#define HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION \
  "Cell wide slope limiter (Springel 2010)"
#include "hydro_slope_limiters_cell.h"

#else

#define HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION "No cell wide slope limiter"

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell_init(
    struct part *p) {}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * neighbour loop
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect(struct part *pi, struct part *pj, float r) {}

/**
 * @brief Slope limit cell gradients
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part *p) {}

#endif

#endif  // SWIFT_HYDRO_SLOPE_LIMITERS_H
