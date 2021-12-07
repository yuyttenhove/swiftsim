//
// Created by yuyttenh on 29/11/2021.
//

#ifndef SWIFTSIM_SHADOWFAX_HYDRO_VELOCITIES_H
#define SWIFTSIM_SHADOWFAX_HYDRO_VELOCITIES_H

/**
 * @brief Initialize the GIZMO particle velocities before the start of the
 * actual run based on the initial value of the primitive velocity.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_init(
    struct part* restrict p, struct xpart* restrict xp) {
#if defined(SHADOWFAX_FIX_CELLS)
  p->v[0] = 0.;
  p->v[1] = 0.;
  p->v[2] = 0.;
#endif

  /* set the initial velocity of the cells */
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
}

__attribute__((always_inline)) INLINE static void hydro_velocities_from_momentum(struct part* p, double* ret) {
    if (p->conserved.mass > 0.) {

      const double inverse_mass = 1.f / p->conserved.mass;

      const double rho = p->conserved.mass / p->voronoi.volume;
      if (rho < 1e-10) {
        /* Suppress velocity linearly near vacuum */
        const double fac = rho * 1e10;
        ret[0] = fac * p->conserved.momentum[0] * inverse_mass;
        ret[1] = fac * p->conserved.momentum[1] * inverse_mass;
        ret[2] = fac * p->conserved.momentum[2] * inverse_mass;
      } else {
        /* Normal case: update fluid velocity and set particle velocity accordingly. */
        ret[0] = p->conserved.momentum[0] * inverse_mass;
        ret[1] = p->conserved.momentum[1] * inverse_mass;
        ret[2] = p->conserved.momentum[2] * inverse_mass;
      }
    } else {
      ret[0] = 0.;
      ret[1] = 0.;
      ret[2] = 0.;
    }
}

/**
 * @brief Set the velocity of a GIZMO particle, based on the values of its
 * primitive variables and the geometry of its mesh-free "cell".
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_set(
    struct part* restrict p, struct xpart* restrict xp) {

#if defined(SHADOWFAX_FIX_CELLS)
  p->v[0] = 0.0f;
  p->v[1] = 0.0f;
  p->v[2] = 0.0f;
#else

  hydro_velocities_from_momentum(p, p->fluid_v);

  p->v[0] = p->fluid_v[0];
  p->v[1] = p->fluid_v[1];
  p->v[2] = p->fluid_v[2];

#ifdef SHADOWFAX_STEER_CELL_MOTION
  if (p->conserved.mass > 0.) {
    double d[3], vfac;
    d[0] = p->voronoi.centroid[0] - p->x[0];
    d[1] = p->voronoi.centroid[1] - p->x[1];
    d[2] = p->voronoi.centroid[2] - p->x[2];

    double d_norm = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    double R = get_radius_dimension_sphere((float)p->voronoi.volume);
    double fac = 4.0f * d_norm / R;
    if (fac > 0.9f) {
      float sound_speed = gas_soundspeed_from_pressure((float)p->rho, (float)p->P);
      if (fac < 1.1f) {
        vfac = sound_speed * (d_norm - 0.225f * R) / d_norm / (0.05f * R);
      } else {
        vfac = sound_speed / d_norm;
      }
      p->v[0] += vfac * d[0];
      p->v[1] += vfac * d[1];
      p->v[2] += vfac * d[2];
    }
  }
#endif
#endif

  /* Now make sure all velocity variables are up to date. */
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];

  if (p->gpart) {
    p->gpart->v_full[0] = (float)p->v[0];
    p->gpart->v_full[1] = (float)p->v[1];
    p->gpart->v_full[2] = (float)p->v[2];
  }
}

#endif // SWIFTSIM_SHADOWFAX_HYDRO_VELOCITIES_H
