#ifndef SWIFT_HYDRO_SHADOWFAX_H
#define SWIFT_HYDRO_SHADOWFAX_H

/* Local headers. */
#include "riemann.h"
#include "equation_of_state.h"


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
  float energy, momentum[3];
  if (m > 0.) {
    momentum[0] = p->conserved.momentum[0];
    momentum[1] = p->conserved.momentum[1];
    momentum[2] = p->conserved.momentum[2];
    p->primitives.rho = m / p->voronoi.volume;
    p->primitives.v[0] = momentum[0] / m;
    p->primitives.v[1] = momentum[1] / m;
    p->primitives.v[2] = momentum[2] / m;

    energy = p->conserved.energy;

#ifdef SHADOWFAX_TOTAL_ENERGY
    energy -= 0.5f * (momentum[0] * p->primitives.v[0] +
                      momentum[1] * p->primitives.v[1] +
                      momentum[2] * p->primitives.v[2]);
#endif

    energy /= m;

    p->primitives.P =
        gas_pressure_from_internal_energy(p->primitives.rho, energy);
  } else {
    p->primitives.rho = 0.;
    p->primitives.v[0] = 0.;
    p->primitives.v[1] = 0.;
    p->primitives.v[2] = 0.;
    p->primitives.P = 0.;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (p->primitives.rho < 0.) {
    error("Negative density!");
  }

  if (p->primitives.P < 0.) {
    error("Negative pressure!");
  }
#endif
}

__attribute__((always_inline)) INLINE static void hydro_shadowfax_flux_exchange(
    struct part *restrict pi, struct part *restrict pj, double const *midpoint,
    double surface_area, const double *shift) {

  /* Initialize local variables */
  float dx[3];
  for (int k = 0; k < 3; k++) {
    dx[k] = (float)pi->x[k] - (float)pj->x[k] - (float)shift[k];
  }
  const float r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
  const float r = sqrtf(r2);

  /* Primitive quantities */
  float Wi[5] = {
      pi->primitives.rho,
      pi->primitives.v[0],
      pi->primitives.v[1],
      pi->primitives.v[2],
      pi->primitives.P
  };
  float Wj[5] = {
      pj->primitives.rho,
      pj->primitives.v[0],
      pj->primitives.v[1],
      pj->primitives.v[2],
      pj->primitives.P
  };

  /* particle velocities */
  float vi[3], vj[3];
  for (int k = 0; k < 3; k++) {
    vi[k] = pi->force.v_full[k];
    vj[k] = pj->force.v_full[k];
  }

  /* calculate the maximal signal velocity */
  float vmax = 0.0f;
  if (Wi[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(Wi[0], Wi[4]);
  }
  if (Wj[0] > 0.) {
    vmax += gas_soundspeed_from_pressure(Wj[0], Wj[4]);
  }
  float dvdotdx = (Wi[1] - Wj[1]) * dx[0] + (Wi[2] - Wj[2]) * dx[1] +
                  (Wi[3] - Wj[3]) * dx[2];
  if (dvdotdx > 0.) {
    vmax -= dvdotdx / r;
  }
  pi->timestepvars.vmax = fmaxf(pi->timestepvars.vmax, vmax);
  pj->timestepvars.vmax = fmaxf(pj->timestepvars.vmax, vmax);

  /* Compute interface velocity */
  float vij[3];
  float fac = (vi[0] - vj[0]) * (midpoint[0] + 0.5f * dx[0]) +
              (vi[1] - vj[1]) * (midpoint[1] + 0.5f * dx[1]) +
              (vi[2] - vj[2]) * (midpoint[2] + 0.5f * dx[2]);
  fac /= r;
  vij[0] = 0.5f * (vi[0] + vj[0]) - fac * dx[0];
  vij[1] = 0.5f * (vi[1] + vj[1]) - fac * dx[1];
  vij[2] = 0.5f * (vi[2] + vj[2]) - fac * dx[2];

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

  /* TODO add gradients */
//  hydro_gradients_predict(pi, pj, hi, hj, dx, r, xij_i, Wi, Wj);

  /* compute the normal vector of the interface */
  float n_unit[3];
  for (int k = 0; k < 3; ++k) {
    n_unit[k] = -dx[k] / r;
  }
  /* we don't need to rotate, we can use the unit vector in the Riemann problem
   * itself (see GIZMO) */
  float totflux[5];
  riemann_solve_for_flux(Wi, Wj, n_unit, vij, totflux);

  /* Update conserved variables */
  /* eqn. (16) */
  pi->conserved.flux.mass -= surface_area * totflux[0];
  pi->conserved.flux.momentum[0] -= surface_area * totflux[1];
  pi->conserved.flux.momentum[1] -= surface_area * totflux[2];
  pi->conserved.flux.momentum[2] -= surface_area * totflux[3];
  pi->conserved.flux.energy -= surface_area * totflux[4];

  pj->conserved.flux.mass += surface_area * totflux[0];
  pj->conserved.flux.momentum[0] += surface_area * totflux[1];
  pj->conserved.flux.momentum[1] += surface_area * totflux[2];
  pj->conserved.flux.momentum[2] += surface_area * totflux[3];
  pj->conserved.flux.energy += surface_area * totflux[4];

#ifndef SHADOWFAX_TOTAL_ENERGY
  float ekin = 0.5f * (pi->primitives.v[0] * pi->primitives.v[0] +
                       pi->primitives.v[1] * pi->primitives.v[1] +
                       pi->primitives.v[2] * pi->primitives.v[2]);
  pi->conserved.flux.energy += surface_area * totflux[1] * pi->primitives.v[0];
  pi->conserved.flux.energy += surface_area * totflux[2] * pi->primitives.v[1];
  pi->conserved.flux.energy += surface_area * totflux[3] * pi->primitives.v[2];
  pi->conserved.flux.energy -= surface_area * totflux[0] * ekin;

  ekin = 0.5f * (pj->primitives.v[0] * pj->primitives.v[0] +
                 pj->primitives.v[1] * pj->primitives.v[1] +
                 pj->primitives.v[2] * pj->primitives.v[2]);
  pj->conserved.flux.energy -= surface_area * totflux[1] * pj->primitives.v[0];
  pj->conserved.flux.energy -= surface_area * totflux[2] * pj->primitives.v[1];
  pj->conserved.flux.energy -= surface_area * totflux[3] * pj->primitives.v[2];
  pj->conserved.flux.energy += surface_area * totflux[0] * ekin;
#endif

  ++pi->voronoi.nface;
  ++pj->voronoi.nface;
}

#endif /* SWIFT_HYDRO_SHADOWFAX_H */
