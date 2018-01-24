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

/**
 * @brief Density loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density_scalar(
    const float r2, const float hi_inv, float *restrict rhoSum, float *restrict rho_dhSum, float *restrict wcountSum, float *restrict wcount_dhSum, const float mj) {
 
 float wi, wi_dx;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  *rhoSum += mj * wi;
  *rho_dhSum -= mj * (hydro_dimension * wi + ui * wi_dx);
  *wcountSum += wi;
  *wcount_dhSum -= (hydro_dimension * wi + ui * wi_dx);
}

#ifdef WITH_VECTORIZATION

/**
 * @brief Density interaction computed using 1 vector
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_1_vec_density(vector *r2, vector *dx, vector *dy, vector *dz,
                                 const struct input_params_density *params, const struct cache *cell_cache, const int cache_idx, struct update_cache_density *sum_cache, mask_t mask) {

  const int int_mask = vec_is_mask_true(mask);

  for(int i=0; i<VEC_SIZE; i++) {
  float rho = 0.f, rho_dh = 0.f, wcount = 0.f, wcount_dh = 0.f;
    if (int_mask & (1 << i)) {
      runner_iact_nonsym_density_scalar(r2->f[i], params->input[input_params_density_hi_inv].f[0], &rho, &rho_dh, &wcount, &wcount_dh, cell_cache->m[cache_idx + i]);
    }
    sum_cache->updates[update_cache_density_rho].f[i] += rho;
    sum_cache->updates[update_cache_density_rho_dh].f[i] += rho_dh;
    sum_cache->updates[update_cache_density_wcount].f[i] += wcount;
    sum_cache->updates[update_cache_density_wcount_dh].f[i] += wcount_dh;

  }

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

  ui.v = vec_mul(r.v, params->input[input_params_density_hi_inv].v);
  ui2.v = vec_mul(r_2.v, params->input[input_params_density_hi_inv].v);

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
    sum_cache->updates[update_cache_density_rho].v = vec_mask_add(sum_cache->updates[update_cache_density_rho].v, vec_mul(mj.v, wi.v), mask);
    sum_cache->updates[update_cache_density_rho].v = vec_mask_add(sum_cache->updates[update_cache_density_rho].v, vec_mul(mj2.v, wi2.v), mask2);
    sum_cache->updates[update_cache_density_rho_dh].v =
        vec_mask_sub(sum_cache->updates[update_cache_density_rho_dh].v, vec_mul(mj.v, wcount_dh_update.v), mask);
    sum_cache->updates[update_cache_density_rho_dh].v =
        vec_mask_sub(sum_cache->updates[update_cache_density_rho_dh].v, vec_mul(mj2.v, wcount_dh_update2.v), mask2);
    sum_cache->updates[update_cache_density_wcount].v = vec_mask_add(sum_cache->updates[update_cache_density_wcount].v, wi.v, mask);
    sum_cache->updates[update_cache_density_wcount].v = vec_mask_add(sum_cache->updates[update_cache_density_wcount].v, wi2.v, mask2);
    sum_cache->updates[update_cache_density_wcount_dh].v = vec_mask_sub(sum_cache->updates[update_cache_density_wcount_dh].v, wcount_dh_update.v, mask);
    sum_cache->updates[update_cache_density_wcount_dh].v = vec_mask_sub(sum_cache->updates[update_cache_density_wcount_dh].v, wcount_dh_update2.v, mask2);
  } else {
    sum_cache->updates[update_cache_density_rho].v = vec_add(sum_cache->updates[update_cache_density_rho].v, vec_mul(mj.v, wi.v));
    sum_cache->updates[update_cache_density_rho].v = vec_add(sum_cache->updates[update_cache_density_rho].v, vec_mul(mj2.v, wi2.v));
    sum_cache->updates[update_cache_density_rho_dh].v = vec_sub(sum_cache->updates[update_cache_density_rho_dh].v, vec_mul(mj.v, wcount_dh_update.v));
    sum_cache->updates[update_cache_density_rho_dh].v = vec_sub(sum_cache->updates[update_cache_density_rho_dh].v, vec_mul(mj2.v, wcount_dh_update2.v));
    sum_cache->updates[update_cache_density_wcount].v = vec_add(sum_cache->updates[update_cache_density_wcount].v, wi.v);
    sum_cache->updates[update_cache_density_wcount].v = vec_add(sum_cache->updates[update_cache_density_wcount].v, wi2.v);
    sum_cache->updates[update_cache_density_wcount_dh].v = vec_sub(sum_cache->updates[update_cache_density_wcount_dh].v, wcount_dh_update.v);
    sum_cache->updates[update_cache_density_wcount_dh].v = vec_sub(sum_cache->updates[update_cache_density_wcount_dh].v, wcount_dh_update2.v);
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

/**
 * @brief Force loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force_scalar(
     float r2, float dx, float dy, float dz, float hj_inv, const struct input_params_force *params, const struct cache *cell_cache, const int cache_idx, float *restrict a_hydro_xSum, float *restrict a_hydro_ySum, float *restrict a_hydro_zSum, float *restrict u_dtSum, float *restrict h_dtSum, float *restrict sigSum) {

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Recover some data */
  const float mj = cell_cache->m[cache_idx];
  const float rhoi = params->input[input_params_force_rhoi].f[0];
  const float rhoj = cell_cache->rho[cache_idx];
  const float P_over_rho2_i = params->input[input_params_force_pOrhoi2].f[0];
  const float P_over_rho2_j = cell_cache->pOrho2[cache_idx];

  /* Get the kernel for hi. */
  const float hi_inv = params->input[input_params_force_hi_inv].f[0];
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  //const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (params->input[input_params_force_vix].f[0] - cell_cache->vx[cache_idx]) * dx +
                     (params->input[input_params_force_viy].f[0] - cell_cache->vy[cache_idx]) * dy +
                     (params->input[input_params_force_viz].f[0] - cell_cache->vz[cache_idx]) * dz;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = params->input[input_params_force_ci].f[0];
  const float cj = cell_cache->soundspeed[cache_idx];
  *sigSum = ci + cj - 3.f * mu_ij;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.5f * const_viscosity_alpha * (*sigSum) * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  *a_hydro_xSum -= mj * acc * dx;
  *a_hydro_ySum -= mj * acc * dy;
  *a_hydro_zSum -= mj * acc * dz;

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;

  /* Internal energy time derivatibe */
  *u_dtSum += du_dt_i * mj;

  /* Get the time derivative for h. */
  *h_dtSum -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  //*sigSum = max(pi->force.v_sig, v_sig);
  //*sigSum = v_sig;
}

#ifdef WITH_VECTORIZATION
static const vector const_viscosity_alpha_fac =
    FILL_VEC(-0.5f * const_viscosity_alpha);

/**
 * @brief Force interaction computed using 1 vector
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_1_vec_force(
    vector *r2, vector *dx, vector *dy, vector *dz, const struct input_params_force *params, const struct cache *cell_cache, const int cache_idx, vector hj_inv, struct update_cache_force *sum_cache, mask_t mask) {

  const int int_mask = vec_is_mask_true(mask);

  //for(int i=0; i<VEC_SIZE; i++) {
  //    message("r2: %f", r2->f[i]);
  //    message("rhoi: %f, rhoj: %f", params->input[input_params_force_rhoi].f[0], cell_cache->rho[cache_idx + i]);
  //    message("hj_inv: %f", hj_inv.f[i]);
  //}
  for(int i=0; i<VEC_SIZE; i++) {
  float a_hydro_x = 0.f, a_hydro_y = 0.f, a_hydro_z = 0.f, u_dt = 0.f, h_dt = 0.f, sig = 0.f;
    if (int_mask & (1 << i)) {
      runner_iact_nonsym_force_scalar(r2->f[i], dx->f[i], dy->f[i], dz->f[i], hj_inv.f[i], params, cell_cache, cache_idx + i, &a_hydro_x, &a_hydro_y, &a_hydro_z, &u_dt, &h_dt, &sig);
    }
    sum_cache->updates[update_cache_force_a_hydro_x].f[i] += a_hydro_x;
    sum_cache->updates[update_cache_force_a_hydro_y].f[i] += a_hydro_y;
    sum_cache->updates[update_cache_force_a_hydro_z].f[i] += a_hydro_z;
    sum_cache->updates[update_cache_force_u_dt].f[i] += u_dt;
    sum_cache->updates[update_cache_force_h_dt].f[i] += h_dt;
    sum_cache->updates[update_cache_force_sig].f[i] = max(sum_cache->updates[update_cache_force_sig].f[i], sig);

  }

//#ifdef WITH_VECTORIZATION
//
//  vector r, ri;
//  vector dvx, dvy, dvz;
//  vector xi, xj;
//  vector hid_inv, hjd_inv;
//  vector wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
//  vector piax, piay, piaz;
//  vector pih_dt;
//  vector v_sig;
//  vector omega_ij, mu_ij;
//  vector rho_ij, visc, visc_term, sph_term, acc;
//  vector sph_du_term_i, visc_du_term, du_dt_i;
//
//  /* Fill vectors. */
//  const vector vjx = vector_load(&cell_cache->vx[cache_idx]);
//  const vector vjy = vector_load(&cell_cache->vy[cache_idx]);
//  const vector vjz = vector_load(&cell_cache->vz[cache_idx]);
//  const vector mj = vector_load(&cell_cache->m[cache_idx]);
//  const vector pjrho = vector_load(&cell_cache->rho[cache_idx]);
//  const vector pjPOrho2 = vector_load(&cell_cache->pOrho2[cache_idx]);
//  const vector cj = vector_load(&cell_cache->soundspeed[cache_idx]);
//
//  const vector fac_mu =
//      vector_set1(1.f); /* Will change with cosmological integration */
//
//  /* Get the radius and inverse radius. */
//  ri = vec_reciprocal_sqrt(*r2);
//  r.v = vec_mul(r2->v, ri.v);
//
//  /* Get the kernel for hi. */
//  hid_inv = pow_dimension_plus_one_vec(params->input[input_params_force_hi_inv]);
//  xi.v = vec_mul(r.v, params->input[input_params_force_hi_inv].v);
//  kernel_eval_dWdx_force_vec(&xi, &wi_dx);
//  wi_dr.v = vec_mul(hid_inv.v, wi_dx.v);
//
//  /* Get the kernel for hj. */
//  hjd_inv = pow_dimension_plus_one_vec(hj_inv);
//  xj.v = vec_mul(r.v, hj_inv.v);
//
//  /* Calculate the kernel. */
//  kernel_eval_dWdx_force_vec(&xj, &wj_dx);
//
//  wj_dr.v = vec_mul(hjd_inv.v, wj_dx.v);
//
//  /* Compute dv. */
//  dvx.v = vec_sub(params->input[input_params_force_vix].v, vjx.v);
//  dvy.v = vec_sub(params->input[input_params_force_viy].v, vjy.v);
//  dvz.v = vec_sub(params->input[input_params_force_viz].v, vjz.v);
//
//  /* Compute dv dot r. */
//  dvdr.v = vec_fma(dvx.v, dx->v, vec_fma(dvy.v, dy->v, vec_mul(dvz.v, dz->v)));
//
//  /* Compute the relative velocity. (This is 0 if the particles move away from
//   * each other and negative otherwise) */
//  omega_ij.v = vec_fmin(dvdr.v, vec_setzero());
//  mu_ij.v =
//      vec_mul(fac_mu.v, vec_mul(ri.v, omega_ij.v)); /* This is 0 or negative */
//
//  /* Compute signal velocity */
//  v_sig.v = vec_fnma(vec_set1(3.f), mu_ij.v, vec_add(params->input[input_params_force_ci].v, cj.v));
//
//  /* Now construct the full viscosity term */
//  rho_ij.v = vec_mul(vec_set1(0.5f), vec_add(params->input[input_params_force_rhoi].v, pjrho.v));
//  visc.v = vec_div(vec_mul(const_viscosity_alpha_fac.v,
//                           vec_mul(v_sig.v, mu_ij.v)),
//                   rho_ij.v);
//
//  /* Now, convolve with the kernel */
//  visc_term.v =
//      vec_mul(vec_set1(0.5f),
//              vec_mul(visc.v, vec_mul(vec_add(wi_dr.v, wj_dr.v), ri.v)));
//
//  sph_term.v =
//      vec_mul(vec_fma(params->input[input_params_force_pOrhoi2].v, wi_dr.v,
//                      vec_mul(pjPOrho2.v, wj_dr.v)),
//              ri.v);
//
//  /* Eventually get the acceleration */
//  acc.v = vec_add(visc_term.v, sph_term.v);
//
//  /* Use the force, Luke! */
//  piax.v = vec_mul(mj.v, vec_mul(dx->v, acc.v));
//  piay.v = vec_mul(mj.v, vec_mul(dy->v, acc.v));
//  piaz.v = vec_mul(mj.v, vec_mul(dz->v, acc.v));
//
//  /* Get the time derivative for u. */
//  sph_du_term_i.v = vec_mul(vec_mul(params->input[input_params_force_pOrhoi2].v, dvdr.v), vec_mul(ri.v, wi_dr.v));
//
//  /* Viscosity term */
//  visc_du_term.v = vec_mul(vec_set1(0.5f), vec_mul(visc_term.v, dvdr.v));
//
//  /* Assemble the energy equation term */
//  du_dt_i.v = vec_mul(mj.v, vec_add(sph_du_term_i.v, visc_du_term.v));
//
//  /* Get the time derivative for h. */
//  pih_dt.v =
//      vec_div(vec_mul(mj.v, vec_mul(dvdr.v, vec_mul(ri.v, wi_dr.v))), pjrho.v);
//
//  /* Store the forces back on the particles. */
//  sum_cache->updates[update_cache_force_a_hydro_x].v = vec_mask_sub(sum_cache->updates[update_cache_force_a_hydro_x].v, piax.v, mask);
//  sum_cache->updates[update_cache_force_a_hydro_y].v = vec_mask_sub(sum_cache->updates[update_cache_force_a_hydro_y].v, piay.v, mask);
//  sum_cache->updates[update_cache_force_a_hydro_z].v = vec_mask_sub(sum_cache->updates[update_cache_force_a_hydro_z].v, piaz.v, mask);
//  sum_cache->updates[update_cache_force_u_dt].v = vec_mask_add(sum_cache->updates[update_cache_force_u_dt].v, du_dt_i.v, mask);
//  sum_cache->updates[update_cache_force_h_dt].v = vec_mask_sub(sum_cache->updates[update_cache_force_h_dt].v, pih_dt.v, mask);
//  sum_cache->updates[update_cache_force_sig].v = vec_fmax(sum_cache->updates[update_cache_force_sig].v, vec_and_mask(v_sig.v, mask));
//
//#else
//
//  error(
//      "The Gadget2 serial version of runner_iact_nonsym_force was called when "
//      "the vectorised version should have been used.");
//
//#endif
}

/**
 * @brief Force interaction computed using 2 interleaved vectors
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_2_vec_force(
    float *R2, float *Dx, float *Dy, float *Dz, vector vix, vector viy,
    vector viz, vector pirho, vector grad_hi, vector piPOrho2, vector balsara_i,
    vector ci, float *Vjx, float *Vjy, float *Vjz, float *Pjrho, float *Grad_hj,
    float *PjPOrho2, float *Balsara_j, float *Cj, float *Mj, vector hi_inv,
    float *Hj_inv, vector *a_hydro_xSum, vector *a_hydro_ySum,
    vector *a_hydro_zSum, vector *h_dtSum, vector *v_sigSum,
    vector *entropy_dtSum, mask_t mask, mask_t mask_2, short mask_cond) {

#ifdef WITH_VECTORIZATION

  vector r, ri;
  vector dvx, dvy, dvz;
  vector ui, uj;
  vector hid_inv, hjd_inv;
  vector wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector piax, piay, piaz;
  vector pih_dt;
  vector v_sig;
  vector omega_ij, mu_ij, balsara;
  vector rho_ij, visc, visc_term, sph_term, acc, entropy_dt;

  vector r_2, ri_2;
  vector dvx_2, dvy_2, dvz_2;
  vector ui_2, uj_2;
  vector hjd_inv_2;
  vector wi_dx_2, wj_dx_2, wi_dr_2, wj_dr_2, dvdr_2;
  vector piax_2, piay_2, piaz_2;
  vector pih_dt_2;
  vector v_sig_2;
  vector omega_ij_2, mu_ij_2, balsara_2;
  vector rho_ij_2, visc_2, visc_term_2, sph_term_2, acc_2, entropy_dt_2;

  /* Fill vectors. */
  const vector mj = vector_load(Mj);
  const vector mj_2 = vector_load(&Mj[VEC_SIZE]);
  const vector vjx = vector_load(Vjx);
  const vector vjx_2 = vector_load(&Vjx[VEC_SIZE]);
  const vector vjy = vector_load(Vjy);
  const vector vjy_2 = vector_load(&Vjy[VEC_SIZE]);
  const vector vjz = vector_load(Vjz);
  const vector vjz_2 = vector_load(&Vjz[VEC_SIZE]);
  const vector dx = vector_load(Dx);
  const vector dx_2 = vector_load(&Dx[VEC_SIZE]);
  const vector dy = vector_load(Dy);
  const vector dy_2 = vector_load(&Dy[VEC_SIZE]);
  const vector dz = vector_load(Dz);
  const vector dz_2 = vector_load(&Dz[VEC_SIZE]);

  /* Get the radius and inverse radius. */
  const vector r2 = vector_load(R2);
  const vector r2_2 = vector_load(&R2[VEC_SIZE]);
  ri = vec_reciprocal_sqrt(r2);
  ri_2 = vec_reciprocal_sqrt(r2_2);
  r.v = vec_mul(r2.v, ri.v);
  r_2.v = vec_mul(r2_2.v, ri_2.v);

  /* Get remaining properties. */
  const vector pjrho = vector_load(Pjrho);
  const vector pjrho_2 = vector_load(&Pjrho[VEC_SIZE]);
  const vector grad_hj = vector_load(Grad_hj);
  const vector grad_hj_2 = vector_load(&Grad_hj[VEC_SIZE]);
  const vector pjPOrho2 = vector_load(PjPOrho2);
  const vector pjPOrho2_2 = vector_load(&PjPOrho2[VEC_SIZE]);
  const vector balsara_j = vector_load(Balsara_j);
  const vector balsara_j_2 = vector_load(&Balsara_j[VEC_SIZE]);
  const vector cj = vector_load(Cj);
  const vector cj_2 = vector_load(&Cj[VEC_SIZE]);
  const vector hj_inv = vector_load(Hj_inv);
  const vector hj_inv_2 = vector_load(&Hj_inv[VEC_SIZE]);

  const vector fac_mu =
      vector_set1(1.f); /* Will change with cosmological integration */

  /* Find the balsara switch. */
  balsara.v = vec_add(balsara_i.v, balsara_j.v);
  balsara_2.v = vec_add(balsara_i.v, balsara_j_2.v);

  /* Get the kernel for hi. */
  hid_inv = pow_dimension_plus_one_vec(hi_inv);
  ui.v = vec_mul(r.v, hi_inv.v);
  ui_2.v = vec_mul(r_2.v, hi_inv.v);
  kernel_eval_dWdx_force_vec(&ui, &wi_dx);
  kernel_eval_dWdx_force_vec(&ui_2, &wi_dx_2);
  wi_dr.v = vec_mul(hid_inv.v, wi_dx.v);
  wi_dr_2.v = vec_mul(hid_inv.v, wi_dx_2.v);

  /* Get the kernel for hj. */
  hjd_inv = pow_dimension_plus_one_vec(hj_inv);
  hjd_inv_2 = pow_dimension_plus_one_vec(hj_inv_2);
  uj.v = vec_mul(r.v, hj_inv.v);
  uj_2.v = vec_mul(r_2.v, hj_inv_2.v);

  /* Calculate the kernel for two particles. */
  kernel_eval_dWdx_force_vec(&uj, &wj_dx);
  kernel_eval_dWdx_force_vec(&uj_2, &wj_dx_2);

  wj_dr.v = vec_mul(hjd_inv.v, wj_dx.v);
  wj_dr_2.v = vec_mul(hjd_inv_2.v, wj_dx_2.v);

  /* Compute dv. */
  dvx.v = vec_sub(vix.v, vjx.v);
  dvx_2.v = vec_sub(vix.v, vjx_2.v);
  dvy.v = vec_sub(viy.v, vjy.v);
  dvy_2.v = vec_sub(viy.v, vjy_2.v);
  dvz.v = vec_sub(viz.v, vjz.v);
  dvz_2.v = vec_sub(viz.v, vjz_2.v);

  /* Compute dv dot r. */
  dvdr.v = vec_fma(dvx.v, dx.v, vec_fma(dvy.v, dy.v, vec_mul(dvz.v, dz.v)));
  dvdr_2.v = vec_fma(dvx_2.v, dx_2.v,
                     vec_fma(dvy_2.v, dy_2.v, vec_mul(dvz_2.v, dz_2.v)));

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_setzero());
  omega_ij_2.v = vec_fmin(dvdr_2.v, vec_setzero());
  mu_ij.v =
      vec_mul(fac_mu.v, vec_mul(ri.v, omega_ij.v)); /* This is 0 or negative */
  mu_ij_2.v = vec_mul(
      fac_mu.v, vec_mul(ri_2.v, omega_ij_2.v)); /* This is 0 or negative */

  /* Compute signal velocity */
  v_sig.v = vec_fnma(vec_set1(3.f), mu_ij.v, vec_add(ci.v, cj.v));
  v_sig_2.v = vec_fnma(vec_set1(3.f), mu_ij_2.v, vec_add(ci.v, cj_2.v));

  /* Now construct the full viscosity term */
  rho_ij.v = vec_mul(vec_set1(0.5f), vec_add(pirho.v, pjrho.v));
  rho_ij_2.v = vec_mul(vec_set1(0.5f), vec_add(pirho.v, pjrho_2.v));

  visc.v = vec_div(vec_mul(const_viscosity_alpha_fac.v,
                           vec_mul(v_sig.v, vec_mul(mu_ij.v, balsara.v))),
                   rho_ij.v);
  visc_2.v =
      vec_div(vec_mul(const_viscosity_alpha_fac.v,
                      vec_mul(v_sig_2.v, vec_mul(mu_ij_2.v, balsara_2.v))),
              rho_ij_2.v);

  /* Now, convolve with the kernel */
  visc_term.v =
      vec_mul(vec_set1(0.5f),
              vec_mul(visc.v, vec_mul(vec_add(wi_dr.v, wj_dr.v), ri.v)));
  visc_term_2.v = vec_mul(
      vec_set1(0.5f),
      vec_mul(visc_2.v, vec_mul(vec_add(wi_dr_2.v, wj_dr_2.v), ri_2.v)));

  vector grad_hi_mul_piPOrho2;
  grad_hi_mul_piPOrho2.v = vec_mul(grad_hi.v, piPOrho2.v);

  sph_term.v =
      vec_mul(vec_fma(grad_hi_mul_piPOrho2.v, wi_dr.v,
                      vec_mul(grad_hj.v, vec_mul(pjPOrho2.v, wj_dr.v))),
              ri.v);
  sph_term_2.v =
      vec_mul(vec_fma(grad_hi_mul_piPOrho2.v, wi_dr_2.v,
                      vec_mul(grad_hj_2.v, vec_mul(pjPOrho2_2.v, wj_dr_2.v))),
              ri_2.v);

  /* Eventually get the acceleration */
  acc.v = vec_add(visc_term.v, sph_term.v);
  acc_2.v = vec_add(visc_term_2.v, sph_term_2.v);

  /* Use the force, Luke! */
  piax.v = vec_mul(mj.v, vec_mul(dx.v, acc.v));
  piax_2.v = vec_mul(mj_2.v, vec_mul(dx_2.v, acc_2.v));
  piay.v = vec_mul(mj.v, vec_mul(dy.v, acc.v));
  piay_2.v = vec_mul(mj_2.v, vec_mul(dy_2.v, acc_2.v));
  piaz.v = vec_mul(mj.v, vec_mul(dz.v, acc.v));
  piaz_2.v = vec_mul(mj_2.v, vec_mul(dz_2.v, acc_2.v));

  /* Get the time derivative for h. */
  pih_dt.v =
      vec_div(vec_mul(mj.v, vec_mul(dvdr.v, vec_mul(ri.v, wi_dr.v))), pjrho.v);
  pih_dt_2.v =
      vec_div(vec_mul(mj_2.v, vec_mul(dvdr_2.v, vec_mul(ri_2.v, wi_dr_2.v))),
              pjrho_2.v);

  /* Change in entropy */
  entropy_dt.v = vec_mul(mj.v, vec_mul(visc_term.v, dvdr.v));
  entropy_dt_2.v = vec_mul(mj_2.v, vec_mul(visc_term_2.v, dvdr_2.v));

  /* Store the forces back on the particles. */
  if (mask_cond) {
    a_hydro_xSum->v = vec_mask_sub(a_hydro_xSum->v, piax.v, mask);
    a_hydro_xSum->v = vec_mask_sub(a_hydro_xSum->v, piax_2.v, mask_2);
    a_hydro_ySum->v = vec_mask_sub(a_hydro_ySum->v, piay.v, mask);
    a_hydro_ySum->v = vec_mask_sub(a_hydro_ySum->v, piay_2.v, mask_2);
    a_hydro_zSum->v = vec_mask_sub(a_hydro_zSum->v, piaz.v, mask);
    a_hydro_zSum->v = vec_mask_sub(a_hydro_zSum->v, piaz_2.v, mask_2);
    h_dtSum->v = vec_mask_sub(h_dtSum->v, pih_dt.v, mask);
    h_dtSum->v = vec_mask_sub(h_dtSum->v, pih_dt_2.v, mask_2);
    v_sigSum->v = vec_fmax(v_sigSum->v, vec_and_mask(v_sig.v, mask));
    v_sigSum->v = vec_fmax(v_sigSum->v, vec_and_mask(v_sig_2.v, mask_2));
    entropy_dtSum->v = vec_mask_add(entropy_dtSum->v, entropy_dt.v, mask);
    entropy_dtSum->v = vec_mask_add(entropy_dtSum->v, entropy_dt_2.v, mask_2);
  } else {
    a_hydro_xSum->v = vec_sub(a_hydro_xSum->v, piax.v);
    a_hydro_xSum->v = vec_sub(a_hydro_xSum->v, piax_2.v);
    a_hydro_ySum->v = vec_sub(a_hydro_ySum->v, piay.v);
    a_hydro_ySum->v = vec_sub(a_hydro_ySum->v, piay_2.v);
    a_hydro_zSum->v = vec_sub(a_hydro_zSum->v, piaz.v);
    a_hydro_zSum->v = vec_sub(a_hydro_zSum->v, piaz_2.v);
    h_dtSum->v = vec_sub(h_dtSum->v, pih_dt.v);
    h_dtSum->v = vec_sub(h_dtSum->v, pih_dt_2.v);
    v_sigSum->v = vec_fmax(v_sigSum->v, v_sig.v);
    v_sigSum->v = vec_fmax(v_sigSum->v, v_sig_2.v);
    entropy_dtSum->v = vec_add(entropy_dtSum->v, entropy_dt.v);
    entropy_dtSum->v = vec_add(entropy_dtSum->v, entropy_dt_2.v);
  }
#else

  error(
      "The Gadget2 serial version of runner_iact_nonsym_force was called when "
      "the vectorised version should have been used.");

#endif
}

#endif

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
