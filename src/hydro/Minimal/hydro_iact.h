/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_MINIMAL_HYDRO_IACT_H
#define SWIFT_MINIMAL_HYDRO_IACT_H

/**
 * @file Minimal/hydro_iact.h
 * @brief Minimal conservative implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term without switches is implemented. No thermal conduction
 * term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

#include "adiabatic_index.h"
#include "minmax.h"
#include "hydro_cache.h"

/**
 * @brief Density loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wj, wi_dx, wj_dx;

  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);
}

/**
 * @brief Density loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);
}

#ifdef WITH_VECTORIZATION

/**
 * @brief Density interaction computed using 1 vector
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_1_vec_density(vector *r2, vector *dx, vector *dy, vector *dz,
                                 const struct input_params_density *params, const struct cache *cell_cache, const int cache_idx, struct update_cache_density *sum_cache,
                                 mask_t mask) {

  vector r, ri, ui, wi, wi_dx;

  /* Fill the vectors. */
  const vector mj = vector_load(&cell_cache->m[cache_idx]);

  /* Get the radius and inverse radius. */
  ri = vec_reciprocal_sqrt(*r2);
  r.v = vec_mul(r2->v, ri.v);

  ui.v = vec_mul(r.v, params->v_hi_inv.v);

  /* Calculate the kernel for two particles. */
  kernel_deval_1_vec(&ui, &wi, &wi_dx);

  vector wcount_dh_update;
  wcount_dh_update.v =
      vec_fma(vec_set1(hydro_dimension), wi.v, vec_mul(ui.v, wi_dx.v));

  /* Mask updates to intermediate vector sums for particle pi. */
  sum_cache->v_rhoSum.v = vec_mask_add(sum_cache->v_rhoSum.v, vec_mul(mj.v, wi.v), mask);
  sum_cache->v_rho_dhSum.v =
      vec_mask_sub(sum_cache->v_rho_dhSum.v, vec_mul(mj.v, wcount_dh_update.v), mask);
  sum_cache->v_wcountSum.v = vec_mask_add(sum_cache->v_wcountSum.v, wi.v, mask);
  sum_cache->v_wcount_dhSum.v = vec_mask_sub(sum_cache->v_wcount_dhSum.v, wcount_dh_update.v, mask);
}

/**
 * @brief Density interaction computed using 2 interleaved vectors
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_2_vec_density(struct c2_cache *int_cache, const int int_cache_idx,
                                 const struct input_params_density *params, 
                                 struct update_cache_density *sum_cache,
                                 mask_t mask, mask_t mask2, short mask_cond) {

  vector r, ri, ui, wi, wi_dx;
  vector r_2, ri2, ui2, wi2, wi_dx2;

  /* Fill the vectors. */
  const vector mj = vector_load(&int_cache->mq[int_cache_idx]);
  const vector mj2 = vector_load(&int_cache->mq[int_cache_idx + VEC_SIZE]);

  /* Get the radius and inverse radius. */
  const vector r2 = vector_load(&int_cache->r2q[int_cache_idx]);
  const vector r2_2 = vector_load(&int_cache->r2q[int_cache_idx + VEC_SIZE]);
  ri = vec_reciprocal_sqrt(r2);
  ri2 = vec_reciprocal_sqrt(r2_2);
  r.v = vec_mul(r2.v, ri.v);
  r_2.v = vec_mul(r2_2.v, ri2.v);

  ui.v = vec_mul(r.v, params->v_hi_inv.v);
  ui2.v = vec_mul(r_2.v, params->v_hi_inv.v);

  /* Calculate the kernel for two particles. */
  kernel_deval_2_vec(&ui, &wi, &wi_dx, &ui2, &wi2, &wi_dx2);

  vector wcount_dh_update, wcount_dh_update2;
  wcount_dh_update.v =
      vec_fma(vec_set1(hydro_dimension), wi.v, vec_mul(ui.v, wi_dx.v));
  wcount_dh_update2.v =
      vec_fma(vec_set1(hydro_dimension), wi2.v, vec_mul(ui2.v, wi_dx2.v));

  /* Mask updates to intermediate vector sums for particle pi. */
  /* Mask only when needed. */
  if (mask_cond) {
    sum_cache->v_rhoSum.v = vec_mask_add(sum_cache->v_rhoSum.v, vec_mul(mj.v, wi.v), mask);
    sum_cache->v_rhoSum.v = vec_mask_add(sum_cache->v_rhoSum.v, vec_mul(mj2.v, wi2.v), mask2);
    sum_cache->v_rho_dhSum.v =
        vec_mask_sub(sum_cache->v_rho_dhSum.v, vec_mul(mj.v, wcount_dh_update.v), mask);
    sum_cache->v_rho_dhSum.v =
        vec_mask_sub(sum_cache->v_rho_dhSum.v, vec_mul(mj2.v, wcount_dh_update2.v), mask2);
    sum_cache->v_wcountSum.v = vec_mask_add(sum_cache->v_wcountSum.v, wi.v, mask);
    sum_cache->v_wcountSum.v = vec_mask_add(sum_cache->v_wcountSum.v, wi2.v, mask2);
    sum_cache->v_wcount_dhSum.v = vec_mask_sub(sum_cache->v_wcount_dhSum.v, wcount_dh_update.v, mask);
    sum_cache->v_wcount_dhSum.v = vec_mask_sub(sum_cache->v_wcount_dhSum.v, wcount_dh_update2.v, mask2);
  } else {
    sum_cache->v_rhoSum.v = vec_add(sum_cache->v_rhoSum.v, vec_mul(mj.v, wi.v));
    sum_cache->v_rhoSum.v = vec_add(sum_cache->v_rhoSum.v, vec_mul(mj2.v, wi2.v));
    sum_cache->v_rho_dhSum.v = vec_sub(sum_cache->v_rho_dhSum.v, vec_mul(mj.v, wcount_dh_update.v));
    sum_cache->v_rho_dhSum.v = vec_sub(sum_cache->v_rho_dhSum.v, vec_mul(mj2.v, wcount_dh_update2.v));
    sum_cache->v_wcountSum.v = vec_add(sum_cache->v_wcountSum.v, wi.v);
    sum_cache->v_wcountSum.v = vec_add(sum_cache->v_wcountSum.v, wi2.v);
    sum_cache->v_wcount_dhSum.v = vec_sub(sum_cache->v_wcount_dhSum.v, wcount_dh_update.v);
    sum_cache->v_wcount_dhSum.v = vec_sub(sum_cache->v_wcount_dhSum.v, wcount_dh_update2.v);
  }
}
#endif

/**
 * @brief Force loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float P_over_rho2_i = pi->force.P_over_rho2;
  const float P_over_rho2_j = pj->force.P_over_rho2;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.5f * const_viscosity_alpha * v_sig * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;
  const float sph_du_term_j = P_over_rho2_j * dvdr * r_inv * wj_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);
}

/**
 * @brief Force loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Recover some data */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float P_over_rho2_i = pi->force.P_over_rho2;
  const float P_over_rho2_j = pj->force.P_over_rho2;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.5f * const_viscosity_alpha * v_sig * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
