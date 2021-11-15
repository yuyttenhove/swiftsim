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
    float qL, float qR, double* cLR, const double* xLR, double rLR, double A,
    double* grad) {}

/** @brief No time extrapolation when gradients are disabled.
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_extrapolate_in_time(const struct part* p, const double* W,
                                    double dt, double* dW) {
  dW[0] = 0.;
  dW[1] = 0.;
  dW[2] = 0.;
  dW[3] = 0.;
  dW[4] = 0.;
}

#endif

/**
 * @brief Gradients reconstruction. Is the same for all gradient types (although
 * gradients_none does nothing, since all gradients are zero -- are they?).
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    const struct part* pi, const struct part* pj, float hi, float hj,
    const double* dx, double r, double* xij_i, double dt, double* Wi, double* Wj) {

  double xij_j[3], dWi[5], dWj[5], dWi_time[5], dWj_time[5];

  /* xij_j = real_centroid - pj->x
           = xij_i + pi->x - pj->x
           = xij_i + dx */
  xij_j[0] = xij_i[0] + dx[0];
  xij_j[1] = xij_i[1] + dx[1];
  xij_j[2] = xij_i[2] + dx[2];

  hydro_gradients_extrapolate_in_time(pi, Wi, dt, dWi_time);
  hydro_gradients_extrapolate_in_time(pj, Wj, dt, dWj_time);

  dWi[0] = pi->gradients.rho[0] * xij_i[0] + pi->gradients.rho[1] * xij_i[1] +
           pi->gradients.rho[2] * xij_i[2] + dWi_time[0];
  dWi[1] = pi->gradients.v[0][0] * xij_i[0] + pi->gradients.v[0][1] * xij_i[1] +
           pi->gradients.v[0][2] * xij_i[2] + dWi_time[1];
  dWi[2] = pi->gradients.v[1][0] * xij_i[0] + pi->gradients.v[1][1] * xij_i[1] +
           pi->gradients.v[1][2] * xij_i[2] + dWi_time[2];
  dWi[3] = pi->gradients.v[2][0] * xij_i[0] + pi->gradients.v[2][1] * xij_i[1] +
           pi->gradients.v[2][2] * xij_i[2] + dWi_time[3];
  dWi[4] = pi->gradients.P[0] * xij_i[0] + pi->gradients.P[1] * xij_i[1] +
           pi->gradients.P[2] * xij_i[2] + dWi_time[4];

  dWj[0] = pj->gradients.rho[0] * xij_j[0] + pj->gradients.rho[1] * xij_j[1] +
           pj->gradients.rho[2] * xij_j[2] + dWj_time[0];
  dWj[1] = pj->gradients.v[0][0] * xij_j[0] + pj->gradients.v[0][1] * xij_j[1] +
           pj->gradients.v[0][2] * xij_j[2] + dWj_time[1];
  dWj[2] = pj->gradients.v[1][0] * xij_j[0] + pj->gradients.v[1][1] * xij_j[1] +
           pj->gradients.v[1][2] * xij_j[2] + dWj_time[2];
  dWj[3] = pj->gradients.v[2][0] * xij_j[0] + pj->gradients.v[2][1] * xij_j[1] +
           pj->gradients.v[2][2] * xij_j[2] + dWj_time[3];
  dWj[4] = pj->gradients.P[0] * xij_j[0] + pj->gradients.P[1] * xij_j[1] +
           pj->gradients.P[2] * xij_j[2] + dWj_time[4];

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
