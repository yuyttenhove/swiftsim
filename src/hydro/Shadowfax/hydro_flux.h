//
// Created by yuyttenh on 14/10/2021.
//

#ifndef SWIFTSIM_HYDRO_FLUX_H
#define SWIFTSIM_HYDRO_FLUX_H

#include "riemann.h"

/* Forward declaration */
__attribute__((always_inline)) static void hydro_convert_conserved_to_primitive(
    struct part* restrict p);

/**
 * @brief Reset the hydrodynamical fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_flux_reset(
    struct part* restrict p) {

  p->conserved.flux.mass = 0.0;
  p->conserved.flux.momentum[0] = 0.0;
  p->conserved.flux.momentum[1] = 0.0;
  p->conserved.flux.momentum[2] = 0.0;
  p->conserved.flux.energy = 0.0;
  p->conserved.flux.dt = -1.0f;
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
__attribute__((always_inline)) INLINE static void hydro_flux_compute(
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

__attribute__((always_inline)) INLINE static void hydro_flux_apply(
    struct part* p) {
  /* Update conserved quantities */
  p->conserved.mass += p->conserved.flux.mass;
  p->conserved.momentum[0] += p->conserved.flux.momentum[0];
  p->conserved.momentum[1] += p->conserved.flux.momentum[1];
  p->conserved.momentum[2] += p->conserved.flux.momentum[2];
#ifdef EOS_ISOTHERMAL_GAS
  /* reset the thermal energy */
  p->conserved.energy =
      p->conserved.mass * gas_internal_energy_from_entropy(0.f, 0.f);
#else
  p->conserved.energy += p->conserved.flux.energy;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  if (p->conserved.mass < 0.0) {
    error("Negative mass!");
  }
  if (p->conserved.energy < 0.0) {
    error("Negative energy!");
  }
#endif

  /* reset fluxes */
  hydro_flux_reset(p);
}

#endif  // SWIFTSIM_HYDRO_FLUX_H
