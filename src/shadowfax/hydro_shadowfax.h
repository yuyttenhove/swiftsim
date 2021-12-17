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
hydro_shadowfax_convert_conserved_to_primitive(struct part *restrict p,
                                               struct xpart *restrict xp) {
  double m = p->conserved.mass;
  double energy;
  if (m > 0.) {
    p->rho = m / p->voronoi.volume;

    hydro_velocities_from_momentum(p, p->fluid_v);

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

  hydro_gravity_velocity_drift(p->fluid_v, p->v, xp->v_full);

  if (m == 0. &&
      (p->fluid_v[0] != 0. || p->fluid_v[1] != 0. || p->fluid_v[2] != 0.)) {
    error("Nonzero fluid_v for particle with zero mass!");
  }

  if (p->rho < 0.) {
    error("Negative density!");
  }

  if (p->P < 0.) {
    error("Negative pressure!");
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
