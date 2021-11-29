//
// Created by yuyttenh on 29/11/2021.
//

#ifndef SWIFTSIM_SHADOWFAX_HYDRO_GRAVITY_H
#define SWIFTSIM_SHADOWFAX_HYDRO_GRAVITY_H

/**
 * @brief Add the gravitational contribution to the fluid velocity drift.
 *
 * @param fluid_v Fluid velocity.
 * @param v (Drifted) particle velocity.
 * @param v_full (Undrifted) particle velocity.
 */
__attribute__((always_inline)) INLINE static void hydro_gravity_velocity_drift(
    double* fluid_v, const double* v, const double* v_full) {

  fluid_v[0] += v[0] - v_full[0];
  fluid_v[1] += v[1] - v_full[1];
  fluid_v[2] += v[2] - v_full[2];
}

/**
 * @brief Update the mass of the gpart associated with the given particle after
 * the mass has been updated with the hydrodynamical mass flux.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
hydro_gravity_update_gpart_mass(struct part* restrict p) {

  if (p->gpart) {
    /* Make sure the gpart knows the mass has changed. */
    p->gpart->mass = p->conserved.mass;
  }
}

/**
 * @brief Get the term required to update the MFV energy due to the change in
 * gravitational energy.
 *
 * @param dt_kick_corr Time step for the potential energy correction.
 * @param dt_grav Time step for the (optional) kinetic energy correction.
 * @param p Particle.
 * @param momentum Momentum of the particle, explicitly requested so that it is
 * clear from the code that the momentum needs to be updated after the call to
 * this function.
 * @param a_grav Gravitational acceleration.
 * @return Term used to update the energy variable.
 */
__attribute__((always_inline)) INLINE static double
hydro_gravity_energy_update_term(const float dt_kick_corr, const float dt_grav,
                                 const struct part* restrict p,
                                 const double* momentum, const float* a_grav) {

  double dE =
      -0.5f * dt_kick_corr *
      (p->gravity.mflux[0] * a_grav[0] + p->gravity.mflux[1] * a_grav[1] +
       p->gravity.mflux[2] * a_grav[2]);
  // TODO: subtract this always?
#if defined(SHADOWFAX_TOTAL_ENERGY)
  dE += dt_grav * (momentum[0] * a_grav[0] + momentum[1] * a_grav[1] +
                   momentum[2] * a_grav[2]);
#endif
  return dE;
}

#endif  // SWIFTSIM_SHADOWFAX_HYDRO_GRAVITY_H
