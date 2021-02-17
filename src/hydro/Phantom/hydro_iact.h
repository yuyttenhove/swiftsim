/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
 *                    Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PHANTOM_HYDRO_IACT_H
#define SWIFT_PHANTOM_HYDRO_IACT_H

/**
 * @file Phantom/hydro_iact.h
 * @brief Density-Energy conservative implementation of SPH based
 *        on Price 2017 (PHANTOM) (interaction routines)
 */

#include "adiabatic_index.h"
#include "const.h"
#include "dimension.h"
#include "minmax.h"

#include "./hydro_parameters.h"
#include "hydro_part.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, const float* dx, float hi, float hj, struct part* pi,
    struct part* pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;
  float dv[3], curlvr[3];

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

  /* Now we need to compute the div terms */
  const float r_inv = 1.f / r;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->viscosity.div_v -= faci * dvdr;
  pj->viscosity.div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

  /* Negative because of the change in sign of dx & dv. */
  pj->density.rot_v[0] += facj * curlvr[0];
  pj->density.rot_v[1] += facj * curlvr[1];
  pj->density.rot_v[2] += facj * curlvr[2];

  /* Compute the magnetic field divergence */
  const float dB[3] = {
                       pi->mhd.B_pred[0] - pj->mhd.B_pred[0],
                       pi->mhd.B_pred[1] - pj->mhd.B_pred[1],
                       pi->mhd.B_pred[2] - pj->mhd.B_pred[2],
  };
  const float dBdx = dB[0] * dx[0]
    + dB[1] * dx[1]
    + dB[2] * dx[2];

  pi->mhd.divB -= mj * dBdx * wi_dx * r_inv;
  pj->mhd.divB -= mi * dBdx * wj_dx * r_inv;
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, const float* dx, float hi, float hj, struct part* pi,
    const struct part* pj, float a, float H) {

  float wi, wi_dx;
  float dv[3], curlvr[3];

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

  const float r_inv = 1.f / r;
  const float faci = mj * wi_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->viscosity.div_v -= faci * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

  /* Compute the magnetic field divergence */
  const float dB[3] = {
                       pi->mhd.B_pred[0] - pj->mhd.B_pred[0],
                       pi->mhd.B_pred[1] - pj->mhd.B_pred[1],
                       pi->mhd.B_pred[2] - pj->mhd.B_pred[2],
  };
  const float dBdx = dB[0] * dx[0]
    + dB[1] * dx[1]
    + dB[2] * dx[2];

  pi->mhd.divB -= mj * dBdx * wi_dx * r_inv;
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj, float a, float H) {

  /* /\* We need to construct the maximal signal velocity between our particle */
  /*  * and all of it's neighbours *\/ */

  /* const float r = sqrtf(r2); */
  /* const float r_inv = 1.f / r; */
  /* const float ci = pi->mhd.force.soundspeed; */
  /* const float cj = pj->mhd.force.soundspeed; */

  /* /\* Cosmology terms for the signal velocity *\/ */
  /* const float fac_mu = pow_three_gamma_minus_five_over_two(a); */
  /* const float a2_Hubble = a * a * H; */

  /* const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] + */
  /*                    (pi->v[1] - pj->v[1]) * dx[1] + */
  /*                    (pi->v[2] - pj->v[2]) * dx[2]; */

  /* /\* Add Hubble flow *\/ */

  /* const float dvdr_Hubble = dvdr + a2_Hubble * r2; */
  /* /\* Are the particles moving towards each others ? *\/ */
  /* const float omega_ij = min(dvdr_Hubble, 0.f); */
  /* const float mu_ij = fac_mu * r_inv * omega_ij; /\* This is 0 or negative *\/ */

  /* /\* Signal velocity *\/ */
  /* const float new_v_sig = ci + cj - const_viscosity_beta * mu_ij; */

  /* /\* Update if we need to *\/ */
  /* pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig); */
  /* pj->viscosity.v_sig = max(pj->viscosity.v_sig, new_v_sig); */
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj, float a, float H) {

  /* /\* We need to construct the maximal signal velocity between our particle */
  /*  * and all of it's neighbours *\/ */

  /* const float r = sqrtf(r2); */
  /* const float r_inv = 1.f / r; */
  /* const float ci = pi->mhd.force.soundspeed; */
  /* const float cj = pj->mhd.force.soundspeed; */

  /* /\* Cosmology terms for the signal velocity *\/ */
  /* const float fac_mu = pow_three_gamma_minus_five_over_two(a); */
  /* const float a2_Hubble = a * a * H; */

  /* const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] + */
  /*                    (pi->v[1] - pj->v[1]) * dx[1] + */
  /*                    (pi->v[2] - pj->v[2]) * dx[2]; */

  /* /\* Add Hubble flow *\/ */

  /* const float dvdr_Hubble = dvdr + a2_Hubble * r2; */
  /* /\* Are the particles moving towards each others ? *\/ */
  /* const float omega_ij = min(dvdr_Hubble, 0.f); */
  /* const float mu_ij = fac_mu * r_inv * omega_ij; /\* This is 0 or negative *\/ */

  /* /\* Signal velocity *\/ */
  /* const float new_v_sig = ci + cj - const_viscosity_beta * mu_ij; */

  /* /\* Update if we need to *\/ */
  /* pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig); */
}

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
                                                                           float r2, const float* dx, float hi, float hj, struct part* pi,
                                                                           const struct part* pj, float a, float H);

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, const float* dx, float hi, float hj, struct part* pi,
    struct part* pj, float a, float H) {

  runner_iact_nonsym_force(r2, dx, hi, hj, pi, pj, a, H);

  float inv_dx[3] = {-dx[0], -dx[1], -dx[2]};

  runner_iact_nonsym_force(r2, inv_dx, hj, hi, pj, pi, a, H);

    /* /\* Cosmological factors entering the EoMs *\/ */
  /* const float fac_mu = pow_three_gamma_minus_five_over_two(a); */
  /* const float a2_Hubble = a * a * H; */

  /* const float r = sqrtf(r2); */
  /* const float r_inv = 1.0f / r; */

  /* /\* Recover some data *\/ */
  /* const float mj = pj->mass; */
  /* const float mi = pi->mass; */

  /* const float rhoi = pi->rho; */
  /* const float rhoj = pj->rho; */

  /* const float pressurei = pi->force.pressure; */
  /* const float pressurej = pj->force.pressure; */

  /* /\* Get the kernel for hi. *\/ */
  /* const float hi_inv = 1.0f / hi; */
  /* const float hid_inv = pow_dimension_plus_one(hi_inv); /\* 1/h^(d+1) *\/ */
  /* const float xi = r * hi_inv; */
  /* float wi, wi_dx; */
  /* kernel_deval(xi, &wi, &wi_dx); */
  /* const float wi_dr = hid_inv * wi_dx; */

  /* /\* Get the kernel for hj. *\/ */
  /* const float hj_inv = 1.0f / hj; */
  /* const float hjd_inv = pow_dimension_plus_one(hj_inv); /\* 1/h^(d+1) *\/ */
  /* const float xj = r * hj_inv; */
  /* float wj, wj_dx; */
  /* kernel_deval(xj, &wj, &wj_dx); */
  /* const float wj_dr = hjd_inv * wj_dx; */

  /* /\* Compute B * gradW *\/ */
  /* const float B_gradW_i = (dx[0] * pi->mhd.B_pred[0] + */
  /*                          dx[1] * pi->mhd.B_pred[1] + */
  /*                          dx[2] * pi->mhd.B_pred[2]) * wi_dr * r_inv; */
  /* const float B_gradW_j = (dx[0] * pj->mhd.B_pred[0] + */
  /*                          dx[1] * pj->mhd.B_pred[1] + */
  /*                          dx[2] * pj->mhd.B_pred[2]) * wj_dr * r_inv; */

  /* /\* Compute dv dot r. *\/ */
  /* const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] + */
  /*                    (pi->v[1] - pj->v[1]) * dx[1] + */
  /*                    (pi->v[2] - pj->v[2]) * dx[2]; */

  /* /\* Includes the hubble flow term; not used for du/dt *\/ */
  /* const float dvdr_Hubble = dvdr + a2_Hubble * r2; */

  /* /\* Are the particles moving towards each others ? *\/ */
  /* const float omega_ij = min(dvdr_Hubble, 0.f); */
  /* const float mu_ij = fac_mu * r_inv * omega_ij; /\* This is 0 or negative *\/ */

  /* /\* Compute sound speeds and signal velocity *\/ */
  /* const float v_sig = pi->force.soundspeed + pj->force.soundspeed - */
  /*                     const_viscosity_beta * mu_ij; */

  /* /\* Balsara term *\/ */
  /* const float balsara_i = pi->force.balsara; */
  /* const float balsara_j = pj->force.balsara; */

  /* /\* Construct the full viscosity term *\/ */
  /* /\* Balsara includes the alphas *\/ */
  /* const float visc = -0.125f * v_sig * mu_ij * (balsara_i + balsara_j); */

  /* /\* Convolve with the kernel *\/ */
  /* const float visc_acc_term = */
  /*     0.5f * visc * (wi_dr * pi->force.f / rhoi + wj_dr * pj->force.f / rhoj) * */
  /*     r_inv; */

  /* /\* Compute the scalar terms in the acceleration *\/ */
  /* const float one_over_rho2_i = 1.f / (rhoi * rhoi) * pi->force.f; */
  /* const float one_over_rho2_j = 1.f / (rhoj * rhoj) * pj->force.f; */

  /* /\* SPH acceleration term *\/ */
  /* const float sph_acc_term_i = one_over_rho2_i * wi_dr * r_inv; */
  /* const float sph_acc_term_j = one_over_rho2_j * wj_dr * r_inv; */

  /* for(int i = 0; i < 3; i++) { */
  /*   /\* Compute the matrix product between the maxwell tensor and dx  *\/ */
  /*   float Mdx_i = 0.f; */
  /*   float Mdx_j = 0.f; */

  /*   for(int j = 0; j < 3; j++) { */
  /*     Mdx_i += pi->mhd.maxwell_stress[i][j] * dx[j]; */
  /*     Mdx_j += pj->mhd.maxwell_stress[i][j] * dx[j]; */
  /*   } */

  /*   /\* Compute the stabilizing term of the magnetic divergence *\/ */
  /*   const float beta_i = pi->mhd.force.beta; */
  /*   const float beta_j = pj->mhd.force.beta; */
  /*   const float B_hat_i = beta_i > 10.f ? 0 : */
  /*     (beta_i > 2.f ? (10.f - beta_i) * pi->mhd.B_pred[i]: pi->mhd.B_pred[i]); */
  /*   const float B_hat_j = beta_j > 10.f ? 0 : */
  /*     (beta_j > 2.f ? (10.f - beta_j) * pj->mhd.B_pred[i]: pj->mhd.B_pred[i]); */

  /*   /\* Compute the acceleration due to the magnetic divergence *\/ */
  /*   const float div_b_acc = B_gradW_i * pi->force.f / (rhoi * rhoi) + */
  /*     B_gradW_j * pj->force.f / (rhoj * rhoj); */

  /*   /\* Assemble the acceleration *\/ */
  /*   const float acc = sph_acc_term_i * Mdx_i + */
  /*     sph_acc_term_j * Mdx_j + visc_acc_term * dx[i]; */


  /*   /\* Use the force Luke ! *\/ */
  /*   pi->a_hydro[i] -= mj * (acc + B_hat_i * div_b_acc); */

  /*   pj->a_hydro[i] += mi * (acc + B_hat_j * div_b_acc); */
  /* } */

  /* /\* Get the time derivative for u. *\/ */
  /* const float sph_du_term_i = pressurei * one_over_rho2_i * dvdr * r_inv * wi_dr; */
  /* const float sph_du_term_j = pressurej * one_over_rho2_j * dvdr * r_inv * wj_dr; */

  /* /\* Viscosity term *\/ */
  /* const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble; */

  /* /\* Diffusion term *\/ */
  /* const float v_diff = */
  /*     sqrtf(2.0 * fabsf(pressurei - pressurej) / (rhoi + rhoj)) + */
  /*     dvdr_Hubble * r_inv; */

  /* const float alpha_diff = 0.5f * (pi->diffusion.alpha + pj->diffusion.alpha); */

  /* /\* wi_dx + wj_dx / 2 is F_ij *\/ */
  /* const float diff_du_term = */
  /*     alpha_diff * v_diff * (pi->u - pj->u) * 0.5f * */
  /*     (wi_dr * pi->force.f / pi->rho + wj_dr * pi->force.f / pi->rho); */

  /* /\* Assemble the energy equation term *\/ */
  /* const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term; */
  /* const float du_dt_j = sph_du_term_j + visc_du_term - diff_du_term; */

  /* /\* Internal energy time derivative *\/ */
  /* pi->u_dt += du_dt_i * mj; */
  /* pj->u_dt += du_dt_j * mi; */

  /* /\* Get the time derivative for h. *\/ */
  /* pi->force.h_dt -= mj * dvdr * pi->force.f * r_inv / rhoj * wi_dr; */
  /* pj->force.h_dt -= mi * dvdr * pj->force.f * r_inv / rhoi * wj_dr; */

  /* /\* Compute the evolution of the magnetic field *\/ */
  /* const float c_mag_i = over_clean_fac * pi->mhd.force.soundspeed; */
  /* const float c_mag_j = over_clean_fac * pj->mhd.force.soundspeed; */

  /* const float psi_i = pi->mhd.psi_pred; */
  /* const float psi_j = pj->mhd.psi_pred; */

  /* /\* Compute part of the magnetic field evolution *\/ */
  /* const float f_over_rho2_i = pi->force.f / (rhoi * rhoi); */
  /* const float f_over_rho2_j = pj->force.f / (rhoj * rhoj); */
  /* const float B_term_i = mj * f_over_rho2_i * B_gradW_i; */
  /* const float B_term_j = mi * f_over_rho2_j * B_gradW_j; */

  /* const float dv[3] = { */
  /*   pi->v[0] - pj->v[0], */
  /*   pi->v[1] - pj->v[1], */
  /*   pi->v[2] - pj->v[2] */
  /* }; */

  /* /\* Magnetic field evolution *\/ */
  /* pi->mhd.force.B_rho_dt[0] -= B_term_i * dv[0]; */
  /* pi->mhd.force.B_rho_dt[1] -= B_term_i * dv[1]; */
  /* pi->mhd.force.B_rho_dt[2] -= B_term_i * dv[2]; */

  /* pj->mhd.force.B_rho_dt[0] += B_term_j * dv[0]; */
  /* pj->mhd.force.B_rho_dt[1] += B_term_j * dv[1]; */
  /* pj->mhd.force.B_rho_dt[2] += B_term_j * dv[2]; */

  /* /\* Magnetic field evolution due to the divergence cleaning *\/ */
  /* const float divB_term = (psi_i * f_over_rho2_i * wi_dr + psi_j * f_over_rho2_j * wj_dr) * r_inv; */

  /* pi->mhd.force.B_rho_dt[0] -= mj * divB_term * dx[0]; */
  /* pi->mhd.force.B_rho_dt[1] -= mj * divB_term * dx[1]; */
  /* pi->mhd.force.B_rho_dt[2] -= mj * divB_term * dx[2]; */

  /* pj->mhd.force.B_rho_dt[0] += mi * divB_term * dx[0]; */
  /* pj->mhd.force.B_rho_dt[1] += mi * divB_term * dx[1]; */
  /* pj->mhd.force.B_rho_dt[2] += mi * divB_term * dx[2]; */

  /* /\* Compute the evolution of the divergence cleaning field *\/ */
  /* const float dB[3] = {pi->mhd.B_pred[0] - pj->mhd.B_pred[0], */
  /*                      pi->mhd.B_pred[1] - pj->mhd.B_pred[1], */
  /*                      pi->mhd.B_pred[2] - pj->mhd.B_pred[2]}; */

  /* const float dBdx = (dB[0] * dx[0] + */
  /*                     dB[1] * dx[1] + */
  /*                     dB[2] * dx[2]); */

  /* /\* Magnetic divergence contribution *\/ */
  /* pi->mhd.force.psi_c_dt += c_mag_i * f_over_rho2_i * rhoi * mj * dBdx * wi_dr * r_inv; */
  /* pj->mhd.force.psi_c_dt -= c_mag_j * f_over_rho2_j * rhoj * mi * dBdx * wj_dr * r_inv; */

  /* /\* Velocity divergence contribution *\/ */
  /* pi->mhd.force.psi_c_dt += 0.5f * psi_i * f_over_rho2_i * rhoi * mj * dvdr * wi_dr * r_inv / c_mag_i; */
  /* pj->mhd.force.psi_c_dt -= 0.5f * psi_j * f_over_rho2_j * rhoj * mi * dvdr * wj_dr * r_inv / c_mag_j; */

  /* /\* divergence damping *\/ */
  /* const float sigma_c = 1.f; */
  /* const float tau_i = hi / (c_mag_i * sigma_c); */
  /* const float tau_j = hj / (c_mag_j * sigma_c); */
  /* pi->mhd.force.psi_c_dt -= psi_i / (c_mag_i * tau_i); */
  /* pj->mhd.force.psi_c_dt -= psi_j / (c_mag_j * tau_j); */

  /* /\* Compute the divergence of B *\/ */
  /* pi->mhd.divB -= f_over_rho2_i * rhoi * mj * dBdx * wi_dr * r_inv; */
  /* pj->mhd.divB += f_over_rho2_j * rhoj * mi * dBdx * wj_dr * r_inv; */
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, const float* dx, float hi, float hj, struct part* pi,
    const struct part* pj, float a, float H) {

  /* Cosmological factors entering the EoMs */
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Recover some data */
  // const float mi = pi->mass;
  const float mj = pj->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

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

  /* Compute B * gradW */
  const float B_gradW_i = (dx[0] * pi->mhd.B_pred[0] +
                           dx[1] * pi->mhd.B_pred[1] +
                           dx[2] * pi->mhd.B_pred[2]) * wi_dr * r_inv;
  const float B_gradW_j = (dx[0] * pj->mhd.B_pred[0] +
                           dx[1] * pj->mhd.B_pred[1] +
                           dx[2] * pj->mhd.B_pred[2]) * wj_dr * r_inv;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Compute dB */
  const float dB[3] = {pi->mhd.B_pred[0] - pj->mhd.B_pred[0],
                       pi->mhd.B_pred[1] - pj->mhd.B_pred[1],
                       pi->mhd.B_pred[2] - pj->mhd.B_pred[2]};

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Compute sound speeds and signal velocity */
  const float v_sig_i = pi->viscosity.alpha * pi->force.soundspeed +
    const_viscosity_beta * fabs(dvdr_Hubble) * r_inv;

  /* Update the velocity signal */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, v_sig_i);

  /* Construct the full viscosity term */
  /* Balsara includes the alphas */
  const float q_i = -0.5f * rhoi * v_sig_i * max(dvdr_Hubble, 0.f) * r_inv;

  /* Compute gradient terms */
  const float f_over_rho2_i = pi->force.f / (rhoi * rhoi);
  const float f_over_rho2_j = pj->force.f / (rhoj * rhoj);

  /* SPH acceleration term */
  const float sph_acc_term_i = f_over_rho2_i * wi_dr * r_inv;
  const float sph_acc_term_j = f_over_rho2_j * wj_dr * r_inv;

  /* Loop over the dimensions */
  for(int i = 0; i < 3; i++) {
    /* Compute the matrix product between the maxwell tensor and dx  */
    float Mdx_i = 0;
    float Mdx_j = 0;

    for(int j = 0; j < 3; j++) {
      Mdx_i += pi->mhd.maxwell_stress[i][j] * dx[j];
      Mdx_j += pj->mhd.maxwell_stress[i][j] * dx[j];
    }
    // TODO remove
    Mdx_i = pi->force.pressure * dx[i];
    Mdx_j = pj->force.pressure * dx[i];

    /* Compute the stabilizing term of the magnetic divergence */
    const float beta_i = pi->mhd.force.beta;
    const float B_hat_i = beta_i > 10.f ? 0 :
      (beta_i > 2.f ? 0.125 * (10.f - beta_i) * pi->mhd.B_pred[i] :
       pi->mhd.B_pred[i]);

    /* Compute the acceleration due to the magnetic divergence */
    const float div_b_acc = B_gradW_i * f_over_rho2_i +
      B_gradW_j * f_over_rho2_j;

    /* Compute the viscosity term */
    const float visc_term = q_i * dx[i];

    /* Assemble the acceleration */
    const float acc = sph_acc_term_i * (Mdx_i + visc_term) +
      sph_acc_term_j * (Mdx_j + visc_term);


    /* Use the force Luke ! */
    pi->a_hydro[i] -= mj * (acc + 0 * B_hat_i * div_b_acc);
  }

  /* Get the time derivative for u. */
  const float sph_du_term_i = pressurei * f_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Diffusion term */
  const float v_diff =
      sqrtf(2.0 * fabsf(pressurei - pressurej) / (rhoi + rhoj)) +
      dvdr_Hubble * r_inv;

  /* wi_dx + wj_dx / 2 is F_ij */
  const float alpha = 1.0f;
  const float diff_du_term =
    alpha * v_diff * (pi->u - pj->u) * 0.5f *
    (wi_dr * pi->force.f / pi->rho + wj_dr * pi->force.f / pi->rho);

  /* Viscosity term */
  const float visc_du_term = - v_sig_i * 0.5 * dvdr_Hubble * dvdr_Hubble * r_inv * r_inv * wi_dr;

  /* Get the B artificial resistivity contribution to u */
  const float dv[3] = {
                       pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]
  };
  const float v_curl_r[3] = {
                             dv[1] * dx[2] - dv[2] * dx[1],
                             dv[2] * dx[0] - dv[0] * dx[2],
                             dv[0] * dx[1] - dv[1] * dx[0],
  };
  const float v_sig_B = sqrtf(v_curl_r[0] * v_curl_r[0]
                              + v_curl_r[1] * v_curl_r[1]
                              + v_curl_r[2] * v_curl_r[2]) * r_inv;

  const float res_i = v_sig_B * f_over_rho2_i * wi_dr;
  const float res_j = v_sig_B * f_over_rho2_j * wj_dr;
  const float dB2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  const float res_du_term = 0.25 * (res_i + res_j) * dB2;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + diff_du_term + visc_du_term
    - 0 * res_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * pi->force.f * r_inv / rhoj * wi_dr;

  /* /\* Compute the evolution of the magnetic field *\/ */
  /* const float c_mag_i = over_clean_fac * pi->mhd.force.soundspeed; */

  /* const float psi_i = pi->mhd.psi_pred; */
  /* const float psi_j = pj->mhd.psi_pred; */

  /* /\* Compute part of the magnetic field evolution *\/ */
  /* const float B_term_i = mj * f_over_rho2_i * B_gradW_i; */

  /* /\* Magnetic field evolution *\/ */
  /* pi->mhd.force.B_rho_dt[0] -= B_term_i * dv[0]; */
  /* pi->mhd.force.B_rho_dt[1] -= B_term_i * dv[1]; */
  /* pi->mhd.force.B_rho_dt[2] -= B_term_i * dv[2]; */

  /* /\* Magnetic field evolution due to the divergence cleaning *\/ */
  /* // TODO change name */
  /* const float divB_term = 0 * (psi_i * f_over_rho2_i * wi_dr + psi_j * f_over_rho2_j * wj_dr) * r_inv; */

  /* pi->mhd.force.B_rho_dt[0] -= mj * divB_term * dx[0]; */
  /* pi->mhd.force.B_rho_dt[1] -= mj * divB_term * dx[1]; */
  /* pi->mhd.force.B_rho_dt[2] -= mj * divB_term * dx[2]; */

  /* /\* Magnetic field evolution due to the dissipation (shock capture) *\/ */
  /* const float res = 0 * 0.5 * mj * (res_i + res_j); */
  /* pi->mhd.force.B_rho_dt[0] += res * dB[0]; */
  /* pi->mhd.force.B_rho_dt[1] += res * dB[1]; */
  /* pi->mhd.force.B_rho_dt[2] += res * dB[2]; */

  /* /\* Compute the evolution of the divergence cleaning field *\/ */

  /* /\* Velocity divergence contribution *\/ */
  /* pi->mhd.force.psi_c_dt += 0.5f * psi_i * f_over_rho2_i * rhoi * mj * dvdr * wi_dr * r_inv / c_mag_i; */

  /* /\* div B contribution done in density and hydro.h *\/ */
}

/**
 * @brief Timestep limiter loop
 *
 * @param r2 Comoving square distance between the two particles.
   * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_limiter(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj, float a, float H) {

  /* Nothing to do here if both particles are active */
}

/**
 * @brief Timestep limiter loop (non-symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_limiter(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj, float a, float H) {

  /* Wake up the neighbour? */
  if (pi->viscosity.v_sig >
      const_limiter_max_v_sig_ratio * pj->viscosity.v_sig) {

    pj->wakeup = time_bin_awake;
  }
}

#endif /* SWIFT_PHANTOM_HYDRO_IACT_H */
