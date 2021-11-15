/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "hydro_slope_limiters.h"
#include "shadowfax/voronoi.h"
/**
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part *p) {

  p->gradients.rho[0] = 0.0f;
  p->gradients.rho[1] = 0.0f;
  p->gradients.rho[2] = 0.0f;

  p->gradients.v[0][0] = 0.0f;
  p->gradients.v[0][1] = 0.0f;
  p->gradients.v[0][2] = 0.0f;

  p->gradients.v[1][0] = 0.0f;
  p->gradients.v[1][1] = 0.0f;
  p->gradients.v[1][2] = 0.0f;

  p->gradients.v[2][0] = 0.0f;
  p->gradients.v[2][1] = 0.0f;
  p->gradients.v[2][2] = 0.0f;

  p->gradients.P[0] = 0.0f;
  p->gradients.P[1] = 0.0f;
  p->gradients.P[2] = 0.0f;

  hydro_slope_limit_cell_init(p);
}

/**
 * @brief Add the gradient estimate for a single quantity due to a particle pair
 * to the total gradient for that quantity
 *
 * This corresponds to one term of equation (21) in Springel (2010).
 *
 * @param qL Value of the quantity on the left.
 * @param qR Value of the quantity on the right.
 * @param cLR Vector pointing from the midpoint of the particle pair to the
 * geometrical centroid of the face in between the particles.
 * @param xLR Vector pointing from the right particle to the left particle.
 * @param rLR Distance between two particles.
 * @param A Surface area of the face in between the particles.
 * @param grad Current value of the gradient for the quantity (is updated).
 */
__attribute__((always_inline)) INLINE void hydro_gradients_single_quantity(
    double qL, double qR, const double *cLR, const double *xLR, double rLR,
    double A, double *grad) {

  grad[0] += A * ((qR - qL) * cLR[0] / rLR - 0.5f * (qL + qR) * xLR[0] / rLR);
  grad[1] += A * ((qR - qL) * cLR[1] / rLR - 0.5f * (qL + qR) * xLR[1] / rLR);
  grad[2] += A * ((qR - qL) * cLR[2] / rLR - 0.5f * (qL + qR) * xLR[2] / rLR);
}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * moved to cell shadowfax
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_collect(
    float r2, const float *dx, float hi, float hj, struct part *pi,
    struct part *pj) {}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * Moved to cell_shadowfax
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_nonsym_collect(float r2, const float *dx, float hi, float hj,
                               struct part *pi, const struct part *pj) {}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part *p) {

  double volume = p->voronoi.volume;

  p->gradients.rho[0] /= volume;
  p->gradients.rho[1] /= volume;
  p->gradients.rho[2] /= volume;

  p->gradients.v[0][0] /= volume;
  p->gradients.v[0][1] /= volume;
  p->gradients.v[0][2] /= volume;
  p->gradients.v[1][0] /= volume;
  p->gradients.v[1][1] /= volume;
  p->gradients.v[1][2] /= volume;
  p->gradients.v[2][0] /= volume;
  p->gradients.v[2][1] /= volume;
  p->gradients.v[2][2] /= volume;

  p->gradients.P[0] /= volume;
  p->gradients.P[1] /= volume;
  p->gradients.P[2] /= volume;

#ifdef SHADOWFAX_SLOPE_LIMITER_CELL_WIDE_EXACT
  /* Uncorrected gradients were used for the extrapolation... */
  p->limiter.extrapolations.rho[0] /= volume;
  p->limiter.extrapolations.rho[1] /= volume;

  p->limiter.extrapolations.v[0][0] /= volume;
  p->limiter.extrapolations.v[0][1] /= volume;
  p->limiter.extrapolations.v[1][0] /= volume;
  p->limiter.extrapolations.v[1][1] /= volume;
  p->limiter.extrapolations.v[2][0] /= volume;
  p->limiter.extrapolations.v[2][1] /= volume;

  p->limiter.extrapolations.P[0] /= volume;
  p->limiter.extrapolations.P[1] /= volume;
#endif

  hydro_slope_limit_cell(p);
}

/** @brief Gradients time extrapolation (makes scheme second order in time).
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_extrapolate_in_time(const struct part *p, const double *W,
                                    double dt, double *dW) {
  const double div_v =
      p->gradients.v[0][0] + p->gradients.v[1][1] + p->gradients.v[2][2];

  dW[0] = -0.5 * dt *
          (W[0] * div_v + W[1] * p->gradients.rho[0] +
           W[2] * p->gradients.rho[1] + W[3] * p->gradients.rho[2]);

  if (W[0] != 0.0f) {
    const double rho_inv = 1. / W[0];
    dW[1] = -0.5 * dt * (W[1] * div_v + rho_inv * p->gradients.P[0]);
    dW[2] = -0.5 * dt * (W[2] * div_v + rho_inv * p->gradients.P[1]);
    dW[3] = -0.5 * dt * (W[3] * div_v + rho_inv * p->gradients.P[2]);
  } else {
    dW[1] = 0.0f;
    dW[2] = 0.0f;
    dW[3] = 0.0f;
  }
  dW[4] = -0.5 * dt *
          (hydro_gamma * W[4] * div_v + W[1] * p->gradients.P[0] +
           W[2] * p->gradients.P[1] + W[3] * p->gradients.P[2]);

  /* Sanity check */
  if (W[0] + dW[0] < 0) {
    dW[0] = 0.;
  }
  if (W[4] + dW[4] < 0) {
    dW[4] = 0.;
  }
}
