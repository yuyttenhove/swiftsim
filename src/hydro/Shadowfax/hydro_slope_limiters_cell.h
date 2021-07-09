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

#include <float.h>

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell_init(
    struct part* p) {

  p->primitives.limiter.rho[0] = FLT_MAX;
  p->primitives.limiter.rho[1] = -FLT_MAX;
  p->primitives.limiter.v[0][0] = FLT_MAX;
  p->primitives.limiter.v[0][1] = -FLT_MAX;
  p->primitives.limiter.v[1][0] = FLT_MAX;
  p->primitives.limiter.v[1][1] = -FLT_MAX;
  p->primitives.limiter.v[2][0] = FLT_MAX;
  p->primitives.limiter.v[2][1] = -FLT_MAX;
  p->primitives.limiter.P[0] = FLT_MAX;
  p->primitives.limiter.P[1] = -FLT_MAX;

#if defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE)
  p->primitives.limiter.maxr = -FLT_MAX;
#elif defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE_EXACT)
  p->primitives.limiter.extrapolations.rho[0] = FLT_MAX;
  p->primitives.limiter.extrapolations.rho[1] = -FLT_MAX;
  p->primitives.limiter.extrapolations.v[0][0] = FLT_MAX;
  p->primitives.limiter.extrapolations.v[0][1] = -FLT_MAX;
  p->primitives.limiter.extrapolations.v[1][0] = FLT_MAX;
  p->primitives.limiter.extrapolations.v[1][1] = -FLT_MAX;
  p->primitives.limiter.extrapolations.v[2][0] = FLT_MAX;
  p->primitives.limiter.extrapolations.v[2][1] = -FLT_MAX;
  p->primitives.limiter.extrapolations.P[0] = FLT_MAX;
  p->primitives.limiter.extrapolations.P[1] = -FLT_MAX;
#endif
}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * neighbour loop
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect(struct part* pi, const struct part* pj,
                               float r) {

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  pi->primitives.limiter.rho[0] =
      fminf(pj->primitives.rho, pi->primitives.limiter.rho[0]);
  pi->primitives.limiter.rho[1] =
      fmaxf(pj->primitives.rho, pi->primitives.limiter.rho[1]);

  pi->primitives.limiter.v[0][0] =
      fminf(pj->primitives.v[0], pi->primitives.limiter.v[0][0]);
  pi->primitives.limiter.v[0][1] =
      fmaxf(pj->primitives.v[0], pi->primitives.limiter.v[0][1]);
  pi->primitives.limiter.v[1][0] =
      fminf(pj->primitives.v[1], pi->primitives.limiter.v[1][0]);
  pi->primitives.limiter.v[1][1] =
      fmaxf(pj->primitives.v[1], pi->primitives.limiter.v[1][1]);
  pi->primitives.limiter.v[2][0] =
      fminf(pj->primitives.v[2], pi->primitives.limiter.v[2][0]);
  pi->primitives.limiter.v[2][1] =
      fmaxf(pj->primitives.v[2], pi->primitives.limiter.v[2][1]);

  pi->primitives.limiter.P[0] =
      fminf(pj->primitives.P, pi->primitives.limiter.P[0]);
  pi->primitives.limiter.P[1] =
      fmaxf(pj->primitives.P, pi->primitives.limiter.P[1]);

#if defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE)
  pi->primitives.limiter.maxr = fmaxf(r, pi->primitives.limiter.maxr);
#endif
}

/*!
 * @brief Collect information about extrapolated primitives for the cell wide
 * limiter.
 *
 * @param p Particle to extrapolate quantities for
 * @param midpoint Position of the midpoint of a face of the voronoi cell
 * associated with p.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect_extrapolations(struct part* p,
                                              const double* midpoint) {
  float dx[3] = {(float)midpoint[0] - (float)p->x[0],
                 (float)midpoint[1] - (float)p->x[1],
                 (float)midpoint[2] - (float)p->x[2]};

  float extrapolated_rho = (p->primitives.gradients.rho[0] * dx[0] +
                            p->primitives.gradients.rho[1] * dx[1] +
                            p->primitives.gradients.rho[2] * dx[2]);
  p->primitives.limiter.extrapolations.rho[0] =
      fminf(p->primitives.limiter.extrapolations.rho[0], extrapolated_rho);
  p->primitives.limiter.extrapolations.rho[1] =
      fmaxf(p->primitives.limiter.extrapolations.rho[1], extrapolated_rho);

  float extrapolated_v[3] = {
      (p->primitives.gradients.v[0][0] * dx[0] +
       p->primitives.gradients.v[0][1] * dx[1] +
       p->primitives.gradients.v[0][2] * dx[2]),
      (p->primitives.gradients.v[1][0] * dx[0] +
       p->primitives.gradients.v[1][1] * dx[1] +
       p->primitives.gradients.v[1][2] * dx[2]),
      (p->primitives.gradients.v[2][0] * dx[0] +
       p->primitives.gradients.v[2][1] * dx[1] +
       p->primitives.gradients.v[2][2] * dx[2]),
  };
  p->primitives.limiter.extrapolations.v[0][0] =
      fminf(extrapolated_v[0], p->primitives.limiter.extrapolations.v[0][0]);
  p->primitives.limiter.extrapolations.v[0][1] =
      fmaxf(extrapolated_v[0], p->primitives.limiter.extrapolations.v[0][1]);
  p->primitives.limiter.extrapolations.v[1][0] =
      fminf(extrapolated_v[1], p->primitives.limiter.extrapolations.v[1][0]);
  p->primitives.limiter.extrapolations.v[1][1] =
      fmaxf(extrapolated_v[1], p->primitives.limiter.extrapolations.v[1][1]);
  p->primitives.limiter.extrapolations.v[2][0] =
      fminf(extrapolated_v[2], p->primitives.limiter.extrapolations.v[2][0]);
  p->primitives.limiter.extrapolations.v[2][1] =
      fmaxf(extrapolated_v[2], p->primitives.limiter.extrapolations.v[2][1]);

  float extrapolated_P = (p->primitives.gradients.P[0] * dx[0] +
                          p->primitives.gradients.P[1] * dx[1] +
                          p->primitives.gradients.P[2] * dx[2]);
  p->primitives.limiter.extrapolations.P[0] =
      fminf(p->primitives.limiter.extrapolations.P[0], extrapolated_P);
  p->primitives.limiter.extrapolations.P[1] =
      fmaxf(p->primitives.limiter.extrapolations.P[1], extrapolated_P);
}

#ifdef SHADOWFAX_SLOPE_LIMITER_CELL_WIDE
/**
 * @brief Apply the cell wide slope limiter to the gradient of a single quantity
 *
 * This corresponds to equation (B2) in Hopkins (2015).
 *
 * @param grad Gradient to slope limit
 * @param qval Value of the quantity at the cell generator
 * @param qmin Minimal value of the quantity among all cell neighbours
 * @param qmax Maximal value of the quantity among all cell neighbours
 * @param maxr Maximal distance between the generator and all of its neighbours
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_quantity(float* grad, float qval, float qmin, float qmax,
                                float maxr) {

  float gradtrue, gradmax, gradmin, alpha;

  gradtrue = sqrtf(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
  if (gradtrue) {
    gradtrue *= maxr;
    gradmax = qmax - qval;
    gradmin = qval - qmin;
    alpha = fminf(1.0f, fminf(gradmax / gradtrue, gradmin / gradtrue));
    grad[0] *= alpha;
    grad[1] *= alpha;
    grad[2] *= alpha;
  }
}

/**
 * @brief Slope limit cell gradients
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part* p) {

  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.rho, p->primitives.rho,
      p->primitives.limiter.rho[0], p->primitives.limiter.rho[1],
      p->primitives.limiter.maxr);

  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[0], p->primitives.v[0],
      p->primitives.limiter.v[0][0], p->primitives.limiter.v[0][1],
      p->primitives.limiter.maxr);
  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[1], p->primitives.v[1],
      p->primitives.limiter.v[1][0], p->primitives.limiter.v[1][1],
      p->primitives.limiter.maxr);
  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[2], p->primitives.v[2],
      p->primitives.limiter.v[2][0], p->primitives.limiter.v[2][1],
      p->primitives.limiter.maxr);

  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.P, p->primitives.P, p->primitives.limiter.P[0],
      p->primitives.limiter.P[1], p->primitives.limiter.maxr);
}

#elif defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE_EXACT)
/**
 * @brief Apply the cell wide slope limiter to the gradient of a single quantity
 *
 * This corresponds to equation (B2) in Hopkins (2015).
 *
 * @param grad Gradient to slope limit
 * @param qval Value of the quantity at the cell generator
 * @param qmin Minimal value of the quantity among all cell neighbours
 * @param qmax Maximal value of the quantity among all cell neighbours
 * @param emax Maximal extrapolated value of the quantity among all cell
 *             neighbours
 * @param emin Minimal extrapolated value of the quantity among all cell
 *             neighbours
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_quantity(float* grad, float qval, float qmin, float qmax,
                                float emax, float emin) {
  float gradmax = qmax - qval;
  float gradmin = qval - qmin;
  float alpha = fminf(1.0f, fminf(gradmax / emax, gradmin / emin));
  grad[0] *= alpha;
  grad[1] *= alpha;
  grad[2] *= alpha;
}

/**
 * @brief Slope limit cell gradients
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part* p) {
  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.rho, p->primitives.rho,
      p->primitives.limiter.rho[0], p->primitives.limiter.rho[1],
      p->primitives.limiter.extrapolations.rho[0],
      p->primitives.limiter.extrapolations.rho[1]);

  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[0], p->primitives.v[0],
      p->primitives.limiter.v[0][0], p->primitives.limiter.v[0][1],
      p->primitives.limiter.extrapolations.v[0][0],
      p->primitives.limiter.extrapolations.v[0][1]);
  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[1], p->primitives.v[1],
      p->primitives.limiter.v[1][0], p->primitives.limiter.v[1][1],
      p->primitives.limiter.extrapolations.v[1][0],
      p->primitives.limiter.extrapolations.v[1][1]);
  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[2], p->primitives.v[2],
      p->primitives.limiter.v[2][0], p->primitives.limiter.v[2][1],
      p->primitives.limiter.extrapolations.v[2][0],
      p->primitives.limiter.extrapolations.v[2][1]);

  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.P, p->primitives.P, p->primitives.limiter.P[0],
      p->primitives.limiter.P[1], p->primitives.limiter.extrapolations.P[0],
      p->primitives.limiter.extrapolations.P[1]);
}
#endif
