//
// Created by yuyttenh on 26/04/2021.
//

#ifndef SWIFTSIM_HYDRO_GRADIENTS_H
#define SWIFTSIM_HYDRO_GRADIENTS_H

#include "hydro_slope_limiters.h"

#ifdef SHADOWFAX_GRADIENTS

#define HYDRO_GRADIENT_IMPLEMENTATION "Shadowfax gradients (Springel 2010)"
#include "hydro_gradients_shadowfax.h"

#else

/* No gradients. Perfectly acceptable, but we have to provide empty functions */
#define HYDRO_GRADIENT_IMPLEMENTATION "No gradients (first order scheme)"

/**
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part* p) {}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_collect(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj) {}

/**
 * @brief Gradient calculations done during the neighbour loop: non-symmetric
 * version
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_nonsym_collect(float r2, const float* dx, float hi, float hj,
                               struct part* restrict pi,
                               const struct part* restrict pj) {}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part* p) {}

/**
 * @param qL Value of the quantity on the left.
 * @param qR Value of the quantity on the right.
 * @param cLR Vector pointing from the midpoint of the particle pair to the
 * geometrical centroid of the face in between the particles.
 * @param xLR Vector pointing from the right particle to the left particle.
 * @param A Surface area of the face in between the particles.
 * @param grad Current value of the gradient for the quantity (is updated).
 */
__attribute__((always_inline)) INLINE void hydro_gradients_single_quantity(
    float qL, float qR, double *cLR, const double *xLR, double rLR, double A,
    double *grad) {}

#endif

/**
 * @brief Gradients reconstruction. Is the same for all gradient types (although
 * gradients_none does nothing, since all gradients are zero -- are they?).
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    const struct part* pi, const struct part* pj, float hi, float hj, const float* dx,
    float r, float* xij_i, float* Wi, float* Wj) {

  float dWi[5], dWj[5];
  float xij_j[3];

  /* xij_j = real_centroid - pj->x
           = xij_i + pi->x - pj->x
           = xij_i + dx */
  xij_j[0] = xij_i[0] + dx[0];
  xij_j[1] = xij_i[1] + dx[1];
  xij_j[2] = xij_i[2] + dx[2];

  dWi[0] = pi->primitives.gradients.rho[0] * xij_i[0] +
           pi->primitives.gradients.rho[1] * xij_i[1] +
           pi->primitives.gradients.rho[2] * xij_i[2];
  dWi[1] = pi->primitives.gradients.v[0][0] * xij_i[0] +
           pi->primitives.gradients.v[0][1] * xij_i[1] +
           pi->primitives.gradients.v[0][2] * xij_i[2];
  dWi[2] = pi->primitives.gradients.v[1][0] * xij_i[0] +
           pi->primitives.gradients.v[1][1] * xij_i[1] +
           pi->primitives.gradients.v[1][2] * xij_i[2];
  dWi[3] = pi->primitives.gradients.v[2][0] * xij_i[0] +
           pi->primitives.gradients.v[2][1] * xij_i[1] +
           pi->primitives.gradients.v[2][2] * xij_i[2];
  dWi[4] = pi->primitives.gradients.P[0] * xij_i[0] +
           pi->primitives.gradients.P[1] * xij_i[1] +
           pi->primitives.gradients.P[2] * xij_i[2];

  dWj[0] = pj->primitives.gradients.rho[0] * xij_j[0] +
           pj->primitives.gradients.rho[1] * xij_j[1] +
           pj->primitives.gradients.rho[2] * xij_j[2];
  dWj[1] = pj->primitives.gradients.v[0][0] * xij_j[0] +
           pj->primitives.gradients.v[0][1] * xij_j[1] +
           pj->primitives.gradients.v[0][2] * xij_j[2];
  dWj[2] = pj->primitives.gradients.v[1][0] * xij_j[0] +
           pj->primitives.gradients.v[1][1] * xij_j[1] +
           pj->primitives.gradients.v[1][2] * xij_j[2];
  dWj[3] = pj->primitives.gradients.v[2][0] * xij_j[0] +
           pj->primitives.gradients.v[2][1] * xij_j[1] +
           pj->primitives.gradients.v[2][2] * xij_j[2];
  dWj[4] = pj->primitives.gradients.P[0] * xij_j[0] +
           pj->primitives.gradients.P[1] * xij_j[1] +
           pj->primitives.gradients.P[2] * xij_j[2];

  hydro_slope_limit_face(Wi, Wj, dWi, dWj, xij_i, xij_j, r);

  Wi[0] += dWi[0];
  Wi[1] += dWi[1];
  Wi[2] += dWi[2];
  Wi[3] += dWi[3];
  Wi[4] += dWi[4];

  Wj[0] += dWj[0];
  Wj[1] += dWj[1];
  Wj[2] += dWj[2];
  Wj[3] += dWj[3];
  Wj[4] += dWj[4];

  /* Sanity check: if density or pressure becomes negative after the
     interpolation, just reset them */
  if (Wi[0] < 0.0f) {
    Wi[0] -= dWi[0];
  }
  if (Wi[4] < 0.0f) {
    Wi[4] -= dWi[4];
  }
  if (Wj[0] < 0.0f) {
    Wj[0] -= dWj[0];
  }
  if (Wj[4] < 0.0f) {
    Wj[4] -= dWj[4];
  }
}

#endif  // SWIFTSIM_HYDRO_GRADIENTS_H
