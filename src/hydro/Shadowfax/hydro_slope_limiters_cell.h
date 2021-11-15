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

  p->limiter.rho[0] = DBL_MAX;
  p->limiter.rho[1] = -DBL_MAX;
  p->limiter.v[0][0] = DBL_MAX;
  p->limiter.v[0][1] = -DBL_MAX;
  p->limiter.v[1][0] = DBL_MAX;
  p->limiter.v[1][1] = -DBL_MAX;
  p->limiter.v[2][0] = DBL_MAX;
  p->limiter.v[2][1] = -DBL_MAX;
  p->limiter.P[0] = DBL_MAX;
  p->limiter.P[1] = -DBL_MAX;

#if defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE)
  p->limiter.maxr = -DBL_MAX;
#elif defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE_EXACT)
  p->limiter.extrapolations.rho[0] = DBL_MAX;
  p->limiter.extrapolations.rho[1] = -DBL_MAX;
  p->limiter.extrapolations.v[0][0] = DBL_MAX;
  p->limiter.extrapolations.v[0][1] = -DBL_MAX;
  p->limiter.extrapolations.v[1][0] = DBL_MAX;
  p->limiter.extrapolations.v[1][1] = -DBL_MAX;
  p->limiter.extrapolations.v[2][0] = DBL_MAX;
  p->limiter.extrapolations.v[2][1] = -DBL_MAX;
  p->limiter.extrapolations.P[0] = DBL_MAX;
  p->limiter.extrapolations.P[1] = -DBL_MAX;
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
                               double r) {

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  pi->limiter.rho[0] = fmin(pj->rho, pi->limiter.rho[0]);
  pi->limiter.rho[1] = fmax(pj->rho, pi->limiter.rho[1]);

  pi->limiter.v[0][0] = fmin(pj->fluid_v[0], pi->limiter.v[0][0]);
  pi->limiter.v[0][1] = fmax(pj->fluid_v[0], pi->limiter.v[0][1]);
  pi->limiter.v[1][0] = fmin(pj->fluid_v[1], pi->limiter.v[1][0]);
  pi->limiter.v[1][1] = fmax(pj->fluid_v[1], pi->limiter.v[1][1]);
  pi->limiter.v[2][0] = fmin(pj->fluid_v[2], pi->limiter.v[2][0]);
  pi->limiter.v[2][1] = fmax(pj->fluid_v[2], pi->limiter.v[2][1]);

  pi->limiter.P[0] = fmin(pj->P, pi->limiter.P[0]);
  pi->limiter.P[1] = fmax(pj->P, pi->limiter.P[1]);

#if defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE)
  pi->limiter.maxr = fmax(r, pi->limiter.maxr);
#endif
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
hydro_slope_limit_cell_quantity(double* grad, double qval, double qmin,
                                double qmax, double maxr) {

  double gradmax, gradmin, alpha, gradtrue;

  gradtrue = sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
  if (gradtrue) {
    gradtrue *= maxr;
    gradmax = qmax - qval;
    gradmin = qval - qmin;
    alpha = fmin(1.0, fmin(gradmax / gradtrue, gradmin / gradtrue));
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

  hydro_slope_limit_cell_quantity(p->gradients.rho, p->rho, p->limiter.rho[0],
                                  p->limiter.rho[1], p->limiter.maxr);

  hydro_slope_limit_cell_quantity(p->gradients.v[0], p->fluid_v[0],
                                  p->limiter.v[0][0], p->limiter.v[0][1],
                                  p->limiter.maxr);
  hydro_slope_limit_cell_quantity(p->gradients.v[1], p->fluid_v[1],
                                  p->limiter.v[1][0], p->limiter.v[1][1],
                                  p->limiter.maxr);
  hydro_slope_limit_cell_quantity(p->gradients.v[2], p->fluid_v[2],
                                  p->limiter.v[2][0], p->limiter.v[2][1],
                                  p->limiter.maxr);

  hydro_slope_limit_cell_quantity(p->gradients.P, p->P, p->limiter.P[0],
                                  p->limiter.P[1], p->limiter.maxr);
}

#elif defined(SHADOWFAX_SLOPE_LIMITER_CELL_WIDE_EXACT)
/*!
 * @brief Collect information about extrapolated for the cell wide
 * limiter.
 *
 * @param p Particle to extrapolate quantities for
 * @param midpoint Position of the midpoint of a face of the voronoi cell
 * associated with p.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect_extrapolations(struct part* p,
                                              const double* midpoint) {
  double dx[3] = {midpoint[0] - p->x[0], midpoint[1] - p->x[1],
                  midpoint[2] - p->x[2]};

  double extrapolated_rho =
      (p->gradients.rho[0] * dx[0] + p->gradients.rho[1] * dx[1] +
       p->gradients.rho[2] * dx[2]);
  p->limiter.extrapolations.rho[0] =
      fmin(p->limiter.extrapolations.rho[0], extrapolated_rho);
  p->limiter.extrapolations.rho[1] =
      fmax(p->limiter.extrapolations.rho[1], extrapolated_rho);

  double extrapolated_v[3] = {
      (p->gradients.v[0][0] * dx[0] + p->gradients.v[0][1] * dx[1] +
       p->gradients.v[0][2] * dx[2]),
      (p->gradients.v[1][0] * dx[0] + p->gradients.v[1][1] * dx[1] +
       p->gradients.v[1][2] * dx[2]),
      (p->gradients.v[2][0] * dx[0] + p->gradients.v[2][1] * dx[1] +
       p->gradients.v[2][2] * dx[2]),
  };
  p->limiter.extrapolations.v[0][0] =
      fmin(extrapolated_v[0], p->limiter.extrapolations.v[0][0]);
  p->limiter.extrapolations.v[0][1] =
      fmax(extrapolated_v[0], p->limiter.extrapolations.v[0][1]);
  p->limiter.extrapolations.v[1][0] =
      fmin(extrapolated_v[1], p->limiter.extrapolations.v[1][0]);
  p->limiter.extrapolations.v[1][1] =
      fmax(extrapolated_v[1], p->limiter.extrapolations.v[1][1]);
  p->limiter.extrapolations.v[2][0] =
      fmin(extrapolated_v[2], p->limiter.extrapolations.v[2][0]);
  p->limiter.extrapolations.v[2][1] =
      fmax(extrapolated_v[2], p->limiter.extrapolations.v[2][1]);

  double extrapolated_P =
      (p->gradients.P[0] * dx[0] + p->gradients.P[1] * dx[1] +
       p->gradients.P[2] * dx[2]);
  p->limiter.extrapolations.P[0] =
      fmin(p->limiter.extrapolations.P[0], extrapolated_P);
  p->limiter.extrapolations.P[1] =
      fmax(p->limiter.extrapolations.P[1], extrapolated_P);
}

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
hydro_slope_limit_cell_quantity(double* grad, double qval, double qmin,
                                double qmax, double emin, double emax) {
  double delta_max = qmax - qval;
  double delta_min = qmin - qval;
  double alpha = 1.;
  if (emin != 0 && emax != 0) {
    alpha = fmin(1.0, fmin(delta_max / emax, delta_min / emin));
  } else if (emin != 0) {
    alpha = fmin(1.0, delta_min / emin);
  } else if (emax != 0) {
    alpha = fmin(1.0, delta_max / emax);
  }
  if (alpha != 1.f) {
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
      p->gradients.rho, p->rho, p->limiter.rho[0], p->limiter.rho[1],
      p->limiter.extrapolations.rho[0], p->limiter.extrapolations.rho[1]);

  hydro_slope_limit_cell_quantity(
      p->gradients.v[0], p->fluid_v[0], p->limiter.v[0][0], p->limiter.v[0][1],
      p->limiter.extrapolations.v[0][0], p->limiter.extrapolations.v[0][1]);
  hydro_slope_limit_cell_quantity(
      p->gradients.v[1], p->fluid_v[1], p->limiter.v[1][0], p->limiter.v[1][1],
      p->limiter.extrapolations.v[1][0], p->limiter.extrapolations.v[1][1]);
  hydro_slope_limit_cell_quantity(
      p->gradients.v[2], p->fluid_v[2], p->limiter.v[2][0], p->limiter.v[2][1],
      p->limiter.extrapolations.v[2][0], p->limiter.extrapolations.v[2][1]);

  hydro_slope_limit_cell_quantity(
      p->gradients.P, p->P, p->limiter.P[0], p->limiter.P[1],
      p->limiter.extrapolations.P[0], p->limiter.extrapolations.P[1]);
}
#endif
