//
// Created by yuyttenh on 14/10/2021.
//

#ifndef SWIFTSIM_HYDRO_FLUX_H
#define SWIFTSIM_HYDRO_FLUX_H

#include "riemann.h"

/**
 * @brief Reset the hydrodynamical fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_reset_hydro_fluxes(
    struct part* restrict p) {

  p->conserved.flux.mass = 0.0f;
  p->conserved.flux.momentum[0] = 0.0f;
  p->conserved.flux.momentum[1] = 0.0f;
  p->conserved.flux.momentum[2] = 0.0f;
  p->conserved.flux.energy = 0.0f;
}

/**
 * @brief Reset the gravity fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_reset_gravity_fluxes(struct part* restrict p) {

  p->gravity.mflux[0] = 0.0f;
  p->gravity.mflux[1] = 0.0f;
  p->gravity.mflux[2] = 0.0f;
}

/**
 * @brief Compute the time and surface integrated flux for the Riemann problem
 * with the given left and right state, and interface normal, surface area,
 * time step and velocity.
 *
 * @param WL Left state variables.
 * @param WR Right state variables.
 * @param n_unit Unit vector of the interface.
 * @param vLR Velocity of the interface.
 * @param surface_area Surface area of the interface.
 * @param dt Time step
 * @param fluxes Array to store the result in (of size 5 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_compute_flux(
    const double* WL, const double* WR, const float* n_unit, const double* vLR,
    const double surface_area, const float dt, float* fluxes) {

  /* we don't need to rotate, we can use the unit vector in the Riemann problem
   * itself (see GIZMO) */
  float WLf[5] = {(float)WL[0], (float)WL[1], (float)WL[2], (float)WL[3],
                  (float)WL[4]};
  float WRf[5] = {(float)WR[0], (float)WR[1], (float)WR[2], (float)WR[3],
                  (float)WR[4]};
  float vLRf[3] = {(float)vLR[0], (float)vLR[1], (float)vLR[2]};
  riemann_solve_for_flux(WLf, WRf, n_unit, vLRf, fluxes);

  const float Adt = (float)(surface_area * dt);
  fluxes[0] *= Adt;
  fluxes[1] *= Adt;
  fluxes[2] *= Adt;
  fluxes[3] *= Adt;
  fluxes[4] *= Adt;
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the left of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Fluxes accross the interface.
 * @param dx Distance between the particles that share the interface.
 * @param dt Time step for the flux exchange.
 */
__attribute__((always_inline)) INLINE static void hydro_flux_update_fluxes_left(
    struct part* restrict p, const float* fluxes, const double* dx) {

  /* Update mass flux vector (used in gravity-hydro coupling) */
  p->gravity.mflux[0] += fluxes[0] * dx[0];
  p->gravity.mflux[1] += fluxes[0] * dx[1];
  p->gravity.mflux[2] += fluxes[0] * dx[2];

  /* Update conserved variables */
  /* eqn. (16) */
  p->conserved.flux.mass -= fluxes[0];
  p->conserved.flux.momentum[0] -= fluxes[1];
  p->conserved.flux.momentum[1] -= fluxes[2];
  p->conserved.flux.momentum[2] -= fluxes[3];
  p->conserved.flux.energy -= fluxes[4];

#ifndef SHADOWFAX_TOTAL_ENERGY
  double ekin =
      0.5 * (p->fluid_v[0] * p->fluid_v[0] + p->fluid_v[1] * p->fluid_v[1] +
             p->fluid_v[2] * p->fluid_v[2]);
  p->conserved.flux.energy += fluxes[1] * p->fluid_v[0];
  p->conserved.flux.energy += fluxes[2] * p->fluid_v[1];
  p->conserved.flux.energy += fluxes[3] * p->fluid_v[2];
  p->conserved.flux.energy -= fluxes[0] * ekin;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  ++p->voronoi.nfluxes;
#endif
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the right of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Fluxes accross the interface.
 * @param dx Distance between the particles that share the interface.
 * @param dt Time step for the flux exchange.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_update_fluxes_right(struct part* restrict p, const float* fluxes,
                               const double* dx) {

  /* Update mass flux vector (used in gravity-hydro coupling) */
  p->gravity.mflux[0] += fluxes[0] * dx[0];
  p->gravity.mflux[1] += fluxes[0] * dx[1];
  p->gravity.mflux[2] += fluxes[0] * dx[2];

  /* Update conserved variables */
  /* eqn. (16) */
  p->conserved.flux.mass += fluxes[0];
  p->conserved.flux.momentum[0] += fluxes[1];
  p->conserved.flux.momentum[1] += fluxes[2];
  p->conserved.flux.momentum[2] += fluxes[3];
  p->conserved.flux.energy += fluxes[4];

#ifndef SHADOWFAX_TOTAL_ENERGY
  double ekin =
      0.5 * (p->fluid_v[0] * p->fluid_v[0] + p->fluid_v[1] * p->fluid_v[1] +
             p->fluid_v[2] * p->fluid_v[2]);
  p->conserved.flux.energy -= fluxes[1] * p->fluid_v[0];
  p->conserved.flux.energy -= fluxes[2] * p->fluid_v[1];
  p->conserved.flux.energy -= fluxes[3] * p->fluid_v[2];
  p->conserved.flux.energy += fluxes[0] * ekin;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  ++p->voronoi.nfluxes;
#endif
}

#endif  // SWIFTSIM_HYDRO_FLUX_H
