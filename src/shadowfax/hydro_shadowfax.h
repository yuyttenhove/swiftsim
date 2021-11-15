#ifndef SWIFT_HYDRO_SHADOWFAX_H
#define SWIFT_HYDRO_SHADOWFAX_H

/* Local headers. */
#include "equation_of_state.h"
#include "hydro/Shadowfax/hydro_flux.h"

/**
 * @brief Convert conserved variables into primitive variables.
 *
 * This method also initializes the gradient variables (if gradients are used).
 *
 * @param p The particle to act upon.
 * @param volume The volume of the particle's associated voronoi cell
 */
__attribute__((always_inline)) INLINE static void
hydro_shadowfax_convert_conserved_to_primitive(struct part *restrict p) {
  float m = p->conserved.mass;
  const float inv_m = 1.f / m;
  float energy;
  if (m > 0.) {
    p->rho = (float)(m / p->voronoi.volume);
    p->fluid_v[0] = p->conserved.momentum[0] * inv_m;
    p->fluid_v[1] = p->conserved.momentum[1] * inv_m;
    p->fluid_v[2] = p->conserved.momentum[2] * inv_m;

    energy = p->conserved.energy;

#ifdef SHADOWFAX_TOTAL_ENERGY
    energy -=
        0.5f * (momentum[0] * p->fluid_v[0] + momentum[1] * p->fluid_v[1] +
                momentum[2] * p->fluid_v[2]);
#endif

    energy /= m;

    p->P = gas_pressure_from_internal_energy(p->rho, energy);
  } else {
    p->rho = 0.f;
    p->fluid_v[0] = 0.f;
    p->fluid_v[1] = 0.f;
    p->fluid_v[2] = 0.f;
    p->P = 0.f;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (p->rho < 0.) {
    error("Negative density!");
  }

  if (p->P < 0.) {
    error("Negative pressure!");
  }
#endif
}

__attribute__((always_inline)) INLINE static void hydro_shadowfax_flux_exchange(
    struct part *pi, struct part *pj, double const *centroid,
    double surface_area, const double *shift, const int symmetric) {

  if (pi->x[0] > 0.8) {
    pi->x[0] = pi->x[0];
  }

  /* Initialize local variables */
  /* Vector from pj to pi */
  float dx[3];
  for (int k = 0; k < 3; k++) {
    dx[k] = (float)pi->x[k] - (float)pj->x[k] - (float)shift[k];
  }
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);

  /* Midpoint between pj and pi */
  float midpoint[3];
  for (int k = 0; k < 3; k++) {
    midpoint[k] = 0.5f * (float)(pi->x[k] + pj->x[k] + shift[k]);
  }

  /* Primitive quantities */
  double Wi[5], Wj[5];
  hydro_get_primitives(pi, Wi);
  hydro_get_primitives(pj, Wj);

  /* particle velocities */
  float vi[3], vj[3];
  for (int k = 0; k < 3; k++) {
    vi[k] = pi->v[k];
    vj[k] = pj->v[k];
#if defined(SWIFT_DEBUG_CHECKS) && !defined(SHADOWFAX_STEER_CELL_MOTION) && \
    !defined(SHADOWFAX_FIX_CELLS)
    assert(pi->fluid_v[k] == pi->v[k]);
    assert(pj->fluid_v[k] == pj->v[k]);
#endif
  }

  /* calculate the maximal signal velocity */
  double vmax = 0.0f;
  if (Wi[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(pi->rho, pi->P);
  }
  if (Wj[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(pj->rho, pj->P);
  }
  double dvdotdx = (Wi[1] - Wj[1]) * dx[0] + (Wi[2] - Wj[2]) * dx[1] +
                   (Wi[3] - Wj[3]) * dx[2];
  if (dvdotdx > 0.) {
    vmax -= dvdotdx / r;
  }
  pi->timestepvars.vmax = (float)fmax(pi->timestepvars.vmax, vmax);
  pj->timestepvars.vmax = (float)fmax(pj->timestepvars.vmax, vmax);

  /* Compute interface velocity, see Springel 2010 (33) */
  double vij[3];
  double fac = ((vj[0] - vi[0]) * (centroid[0] - midpoint[0]) +
                (vj[1] - vi[1]) * (centroid[1] - midpoint[1]) +
                (vj[2] - vi[2]) * (centroid[2] - midpoint[2])) /
               r2;
  vij[0] = 0.5f * (vi[0] + vj[0]) + fac * dx[0];
  vij[1] = 0.5f * (vi[1] + vj[1]) + fac * dx[1];
  vij[2] = 0.5f * (vi[2] + vj[2]) + fac * dx[2];
#if defined(SWIFT_DEBUG_CHECKS) && defined(SHADOWFAX_FIX_CELLS)
  assert(vij[0] == 0.f && vij[1] == 0.f && vij[2] == 0.);
#endif

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

  /* get the time step for the flux exchange. This is always the smallest time
     step among the two particles */
  const float min_dt = (pj->conserved.flux.dt > 0.0f)
                           ? fminf(pi->conserved.flux.dt, pj->conserved.flux.dt)
                           : pi->conserved.flux.dt;

  float xij_i[3];
  for (int k = 0; k < 3; k++) {
    xij_i[k] = (float)(centroid[k] - pi->x[k]);
  }
  hydro_gradients_predict(pi, pj, pi->h, pj->h, dx, r, xij_i, min_dt, Wi, Wj);

  /* compute the normal vector of the interface */
  float n_unit[3];
  for (int k = 0; k < 3; ++k) {
    n_unit[k] = -dx[k] / r;
  }

#ifdef SWIFT_DEBUG_CHECKS
  assert(pi->conserved.flux.dt >= 0);
  assert(min_dt >= 0);
#endif

  float totflux[5];
  hydro_compute_flux(Wi, Wj, n_unit, vij, (float)surface_area, min_dt, totflux);

  /* Update conserved variables */
  /* eqn. (16) */
  pi->conserved.flux.mass -= totflux[0];
  pi->conserved.flux.momentum[0] -= totflux[1];
  pi->conserved.flux.momentum[1] -= totflux[2];
  pi->conserved.flux.momentum[2] -= totflux[3];
  pi->conserved.flux.energy -= totflux[4];

#ifndef SHADOWFAX_TOTAL_ENERGY
  float ekin = 0.5f * (pi->fluid_v[0] * pi->fluid_v[0] +
                       pi->fluid_v[1] * pi->fluid_v[1] +
                       pi->fluid_v[2] * pi->fluid_v[2]);
  pi->conserved.flux.energy += totflux[1] * pi->fluid_v[0];
  pi->conserved.flux.energy += totflux[2] * pi->fluid_v[1];
  pi->conserved.flux.energy += totflux[3] * pi->fluid_v[2];
  pi->conserved.flux.energy -= totflux[0] * ekin;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  ++pi->voronoi.nfluxes;
#endif

  if (symmetric || (pj->conserved.flux.dt < 0.0f)) {
    pj->conserved.flux.mass += totflux[0];
    pj->conserved.flux.momentum[0] += totflux[1];
    pj->conserved.flux.momentum[1] += totflux[2];
    pj->conserved.flux.momentum[2] += totflux[3];
    pj->conserved.flux.energy += totflux[4];

#ifndef SHADOWFAX_TOTAL_ENERGY
    ekin = 0.5f *
           (pj->fluid_v[0] * pj->fluid_v[0] + pj->fluid_v[1] * pj->fluid_v[1] +
            pj->fluid_v[2] * pj->fluid_v[2]);
    pj->conserved.flux.energy -= totflux[1] * pj->fluid_v[0];
    pj->conserved.flux.energy -= totflux[2] * pj->fluid_v[1];
    pj->conserved.flux.energy -= totflux[3] * pj->fluid_v[2];
    pj->conserved.flux.energy += totflux[0] * ekin;
#endif

#ifdef SWIFT_DEBUG_CHECKS
    ++pj->voronoi.nfluxes;
#endif
  }
}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * @param pi Particle i (left particle).
 * @param pj Particle j (right particle).
 * @param midpoint Midpoint of the interface between pi's and pj's cell.
 * @param surface_area Surface area of the interface between pi's and pj's cell.
 * @param shift The shift vector to apply to pj.
 * @param symmetric Whether or not to update pj also
 */
__attribute__((always_inline)) INLINE static void
hydro_shadowfax_gradients_collect(struct part *pi, struct part *pj,
                                  double const *midpoint, double surface_area,
                                  const double *shift, int symmetric) {
  if (!surface_area) {
    /* particle is not a cell neighbour: do nothing */
    return;
  }

  /* Initialize local variables */
  const double dx[3] = {pi->x[0] - pj->x[0] - shift[0],
                        pi->x[1] - pj->x[1] - shift[1],
                        pi->x[2] - pj->x[2] - shift[2]};
  const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* c is supposed to be the vector pointing from the midpoint of pi and pj to
     the centroid of the face between pi and pj.
     The coordinates of the centroid of the face of the voronoi cell of particle
     pi are given in the case of periodic boundary conditions. */
  double c[3] = {midpoint[0] - 0.5 * (pi->x[0] + pj->x[0] + shift[0]),
                 midpoint[1] - 0.5 * (pi->x[1] + pj->x[1] + shift[1]),
                 midpoint[2] - 0.5 * (pi->x[2] + pj->x[2] + shift[2])};

  double r = sqrt(r2);
  hydro_gradients_single_quantity(pi->rho, pj->rho, c, dx, r, surface_area,
                                  pi->gradients.rho);
  hydro_gradients_single_quantity(pi->fluid_v[0], pj->fluid_v[0], c, dx, r,
                                  surface_area, pi->gradients.v[0]);
  hydro_gradients_single_quantity(pi->fluid_v[1], pj->fluid_v[1], c, dx, r,
                                  surface_area, pi->gradients.v[1]);
  hydro_gradients_single_quantity(pi->fluid_v[2], pj->fluid_v[2], c, dx, r,
                                  surface_area, pi->gradients.v[2]);
  hydro_gradients_single_quantity(pi->P, pj->P, c, dx, r, surface_area,
                                  pi->gradients.P);

  double f_ij[3] = {midpoint[0] - pi->x[0], midpoint[1] - pi->x[1],
                    midpoint[2] - pi->x[2]};
  double r_ij = sqrt(f_ij[0] * f_ij[0] + f_ij[1] * f_ij[1] + f_ij[2] * f_ij[2]);
  hydro_slope_limit_cell_collect(pi, pj, r_ij);

  if (symmetric || (pj->conserved.flux.dt < 0.0f)) {
    double mindx[3];
    mindx[0] = -dx[0];
    mindx[1] = -dx[1];
    mindx[2] = -dx[2];
    hydro_gradients_single_quantity(pj->rho, pi->rho, c, mindx, r, surface_area,
                                    pj->gradients.rho);
    hydro_gradients_single_quantity(pj->fluid_v[0], pi->fluid_v[0], c, mindx, r,
                                    surface_area, pj->gradients.v[0]);
    hydro_gradients_single_quantity(pj->fluid_v[1], pi->fluid_v[1], c, mindx, r,
                                    surface_area, pj->gradients.v[1]);
    hydro_gradients_single_quantity(pj->fluid_v[2], pi->fluid_v[2], c, mindx, r,
                                    surface_area, pj->gradients.v[2]);
    hydro_gradients_single_quantity(pj->P, pi->P, c, mindx, r, surface_area,
                                    pj->gradients.P);

    double f_ji[3] = {midpoint[0] - pj->x[0] - shift[0],
                      midpoint[1] - pj->x[1] - shift[1],
                      midpoint[2] - pj->x[2] - shift[2]};
    double r_ji =
        sqrt(f_ji[0] * f_ji[0] + f_ji[1] * f_ji[1] + f_ji[2] * f_ji[2]);
    hydro_slope_limit_cell_collect(pj, pi, r_ji);
  }
}

#endif /* SWIFT_HYDRO_SHADOWFAX_H */
