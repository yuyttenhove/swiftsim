/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (jame.s.willis@durham.ac.uk)
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
#ifndef SWIFT_CACHE_H
#define SWIFT_CACHE_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "align.h"
#include "cell.h"
#include "error.h"
#include "part.h"
#include "sort_part.h"
#include "vector.h"

#define NUM_VEC_PROC 2
#define C2_CACHE_SIZE (NUM_VEC_PROC * VEC_SIZE * 6) + (NUM_VEC_PROC * VEC_SIZE)

#ifdef WITH_VECTORIZATION
/* Cache struct to hold a local copy of a cells' particle
 * properties required for density/force calculations.*/
struct cache {

  /* Particle x position. */
  float *restrict x SWIFT_CACHE_ALIGN;

  /* Particle y position. */
  float *restrict y SWIFT_CACHE_ALIGN;

  /* Particle z position. */
  float *restrict z SWIFT_CACHE_ALIGN;

  /* Particle smoothing length. */
  float *restrict h SWIFT_CACHE_ALIGN;

  /* Particle mass. */
  float *restrict m SWIFT_CACHE_ALIGN;

  /* Particle x velocity. */
  float *restrict vx SWIFT_CACHE_ALIGN;

  /* Particle y velocity. */
  float *restrict vy SWIFT_CACHE_ALIGN;

  /* Particle z velocity. */
  float *restrict vz SWIFT_CACHE_ALIGN;

  /* Maximum index into neighbouring cell for particles that are in range. */
  int *restrict max_index SWIFT_CACHE_ALIGN;

  /* Particle density. */
  float *restrict rho SWIFT_CACHE_ALIGN;

  /* Pressure. */
  float *restrict pressure SWIFT_CACHE_ALIGN;
  
  /* Particle smoothing length gradient. */
  float *restrict grad_h SWIFT_CACHE_ALIGN;

  /* Particle sound speed. */
  float *restrict soundspeed SWIFT_CACHE_ALIGN;

  /* Cache size. */
  int count;
};

/* Secondary cache struct to hold a list of interactions between two
 * particles.*/
struct c2_cache {

  /* Separation between two particles squared. */
  float r2q[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* x separation between two particles. */
  float dxq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* y separation between two particles. */
  float dyq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* z separation between two particles. */
  float dzq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

  /* Mass of particle pj. */
  float mq[C2_CACHE_SIZE] SWIFT_CACHE_ALIGN;

};

/* List which particle density parameters are required as input. */
enum input_params_density_types {
  input_params_density_hi_inv = 0,
  input_params_density_length
};

/* List which particle density parameters need updating. */
enum update_cache_density_types {
  update_cache_density_rho = 0,
  update_cache_density_rho_dh,
  update_cache_density_wcount,
  update_cache_density_wcount_dh,
  update_cache_density_length
};

/* List which particle force parameters are required as input. */
enum input_params_force_types {
  input_params_force_vix = 0,
  input_params_force_viy,
  input_params_force_viz,
  input_params_force_hi_inv,
  input_params_force_rhoi,
  input_params_force_pressure,
  input_params_force_grad_h,
  input_params_force_ci,
  input_params_force_length
};

/* List which particle force parameters need updating. */
enum update_cache_force_types {
  update_cache_force_a_hydro_x = 0,
  update_cache_force_a_hydro_y,
  update_cache_force_a_hydro_z,
  update_cache_force_u_dt,
  update_cache_force_h_dt,
  update_cache_force_sig,
  update_cache_force_length
};

/* Cache to hold a list of vectors used to update particle properties after a density interaction. */
struct update_cache_density {
  vector updates[update_cache_density_length];  
};

/* Cache to hold a list of vectors used to update particle properties after a force interaction. */
struct update_cache_force {
  vector updates[update_cache_force_length];  
};

/* Input parameters needed for computing the density interaction. */
struct input_params_density {
  vector input[input_params_density_length];
};

/* Input parameters needed for computing the force interaction. */
struct input_params_force {
  vector input[input_params_force_length];
};

/**
 * @brief Reset the density update cache to zero.
 *
 * @param c The update cache.
 */
__attribute__((always_inline)) INLINE void update_cache_density_init(struct update_cache_density *c) {

  for(int i=0; i<update_cache_density_length; i++) c->updates[i].v = vec_setzero();

}

/**
 * @brief Reset the force update cache to zero.
 *
 * @param c The update cache.
 */
__attribute__((always_inline)) INLINE void update_cache_force_init(struct update_cache_force *c) {

  for(int i=0; i<update_cache_force_length; i++) c->updates[i].v = vec_setzero();
}

/**
 * @brief Perform horizontal adds on vector sums and store result in particle pi. Density interaction.
 *
 * @param pi Particle to update.
 * @param c Update cache.
 */
__attribute__((always_inline)) INLINE void update_density_particle(struct part *restrict pi, struct update_cache_density *c) {

  VEC_HADD(c->updates[update_cache_density_rho], pi->rho);
  VEC_HADD(c->updates[update_cache_density_rho_dh], pi->density.rho_dh);
  VEC_HADD(c->updates[update_cache_density_wcount], pi->density.wcount);
  VEC_HADD(c->updates[update_cache_density_wcount_dh], pi->density.wcount_dh);
}

/**
 * @brief Perform horizontal adds on vector sums and store result in particle pi. Force interaction.
 *
 * @param pi Particle to update.
 * @param c Update cache.
 */
__attribute__((always_inline)) INLINE void update_force_particle(struct part *restrict pi, struct update_cache_force *c) {

  VEC_HADD(c->updates[update_cache_force_a_hydro_x], pi->a_hydro[0]);
  VEC_HADD(c->updates[update_cache_force_a_hydro_y], pi->a_hydro[1]);
  VEC_HADD(c->updates[update_cache_force_a_hydro_z], pi->a_hydro[2]);
  VEC_HADD(c->updates[update_cache_force_u_dt], pi->u_dt);
  VEC_HADD(c->updates[update_cache_force_h_dt], pi->force.h_dt);
  float max_v_sig = 0.f;
  VEC_HMAX(c->updates[update_cache_force_sig], max_v_sig);
  pi->force.v_sig = max(pi->force.v_sig, max_v_sig);

}

/**
 * @brief Populate the parameters used in the interaction function using a cache. Density interaction.
 *
 * @param c Particle cache.
 * @param cache_index Cache index.
 * @param params Input parameters.
 */
__attribute__((always_inline)) INLINE void populate_input_params_density_cache(const struct cache *restrict c, const int cache_index, struct input_params_density *params) {

  const float hi = c->h[cache_index];
  const float hi_inv = 1.f / hi;
  
  params->input[input_params_density_hi_inv] = vector_set1(hi_inv);
}

/**
 * @brief Populate the parameters used in the interaction function using a cache. Force interaction.
 *
 * @param c Particle cache.
 * @param cache_index Cache index.
 * @param params Input parameters.
 */
__attribute__((always_inline)) INLINE void populate_input_params_force_cache(const struct cache *restrict c, const int cache_index, struct input_params_force *params) {

  const float hi = c->h[cache_index];
  const float hi_inv = 1.f / hi;
  
  params->input[input_params_force_vix] = vector_set1(c->vx[cache_index]);
  params->input[input_params_force_viy] = vector_set1(c->vy[cache_index]);
  params->input[input_params_force_viz] = vector_set1(c->vz[cache_index]);
  params->input[input_params_force_hi_inv] = vector_set1(hi_inv); 
  params->input[input_params_force_rhoi] = vector_set1(c->rho[cache_index]);
  params->input[input_params_force_pressure] = vector_set1(c->pressure[cache_index]);
  params->input[input_params_force_grad_h] = vector_set1(c->grad_h[cache_index]);
  params->input[input_params_force_ci] = vector_set1(c->soundspeed[cache_index]);

}

/**
 * @brief Populate the parameters used in the interaction function. Density interaction.
 *
 * @param pi Particle to update.
 * @param params Input parameters.
 */
__attribute__((always_inline)) INLINE void populate_input_params_density(struct part *restrict pi, struct input_params_density *params) {

  const float hi = pi->h;
  const float hi_inv = 1.f / hi;
  
  params->input[input_params_density_hi_inv] = vector_set1(hi_inv);
}

/**
 * @brief Allocate memory and initialise cache.
 *
 * @param c The cache.
 * @param count Number of particles to allocate space for.
 */
__attribute__((always_inline)) INLINE void cache_init(struct cache *c,
                                                      size_t count) {

  /* Align cache on correct byte boundary and pad cache size to be a multiple of
   * the vector size
   * and include 2 vector lengths for remainder operations. */
  size_t pad = 2 * VEC_SIZE, rem = count % VEC_SIZE;
  if (rem > 0) pad += VEC_SIZE - rem;
  size_t sizeBytes = (count + pad) * sizeof(float);
  size_t sizeIntBytes = (count + pad) * sizeof(int);
  int error = 0;

  /* Free memory if cache has already been allocated. */
  if (c->count > 0) {
    free(c->x);
    free(c->y);
    free(c->z);
    free(c->m);
    free(c->vx);
    free(c->vy);
    free(c->vz);
    free(c->h);
    free(c->max_index);
    free(c->rho);
    free(c->pressure);
    free(c->grad_h);
    free(c->soundspeed);
  }

  error += posix_memalign((void **)&c->x, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->y, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->z, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->m, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->vx, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->vy, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->vz, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->h, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error += posix_memalign((void **)&c->max_index, SWIFT_CACHE_ALIGNMENT,
                          sizeIntBytes);
  error += posix_memalign((void **)&c->rho, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error +=
      posix_memalign((void **)&c->pressure, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error +=
      posix_memalign((void **)&c->grad_h, SWIFT_CACHE_ALIGNMENT, sizeBytes);
  error +=
      posix_memalign((void **)&c->soundspeed, SWIFT_CACHE_ALIGNMENT, sizeBytes);

  if (error != 0)
    error("Couldn't allocate cache, no. of particles: %d", (int)count);
  c->count = count;
}

/**
 * @brief Populate cache by reading in the particles in unsorted order.
 *
 * @param ci The #cell.
 * @param ci_cache The cache.
 */
__attribute__((always_inline)) INLINE void cache_read_particles(
    const struct cell *restrict const ci,
    struct cache *restrict const ci_cache, const int ci_count) {

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  
  const struct part *restrict parts = ci->parts;
  const double loc[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};

  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
  for (int i = 0; i < ci_count; i++) {
    x[i] = (float)(parts[i].x[0] - loc[0]);
    y[i] = (float)(parts[i].x[1] - loc[1]);
    z[i] = (float)(parts[i].x[2] - loc[2]);
    h[i] = parts[i].h;
    m[i] = parts[i].mass;
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  const double max_dx = ci->dx_max_part;
  const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                               -(2. * ci->width[1] + max_dx),
                               -(2. * ci->width[2] + max_dx)};
  const float h_padded = ci->parts[0].h;
  const int ci_count_padded = ci_count - (ci_count % (NUM_VEC_PROC * VEC_SIZE)) + NUM_VEC_PROC * VEC_SIZE;

  for (int i = ci_count; i < ci_count_padded; i++) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    h[i] = h_padded;
    m[i] = 1.f;
  }
}

/**
 * @brief Populate cache by only reading particles that are within range of
 * each other within the adjoining cell. Also read the particles into the cache
 * in sorted order.
 *
 * @param ci The i #cell.
 * @param ci_cache The #cache for cell ci.
 * @param sort_i The array of sorted particle indices for cell ci.
 * @param first_pi The first particle in cell ci that is in range.
 * @param last_pi The last particle in cell ci that is in range.
 * @param last_pi The last particle in cell ci that is in range.
 * @param loc The cell location to remove from the particle positions.
 * @param flipped Flag to check whether the cells have been flipped or not.
 */
__attribute__((always_inline)) INLINE void cache_read_particles_subset(
    const struct cell *restrict const ci, struct cache *restrict const ci_cache,
    const struct entry *restrict sort_i, int *first_pi, int *last_pi,
    const double *loc, const int flipped) {

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  
  const struct part *restrict parts = ci->parts;

  /* The cell is on the right so read the particles
   * into the cache from the start of the cell. */
  if (!flipped) {
    const int rem = (*last_pi + 1) % VEC_SIZE;
    if (rem != 0) {
      const int pad = VEC_SIZE - rem;

      /* Increase last_pi if there are particles in the cell left to read. */
      if (*last_pi + pad < ci->count) *last_pi += pad;
    }

    /* Shift the particles positions to a local frame so single precision can be
     * used instead of double precision. */
    for (int i = 0; i < *last_pi; i++) {
      const int idx = sort_i[i].i;
      x[i] = (float)(parts[idx].x[0] - loc[0]);
      y[i] = (float)(parts[idx].x[1] - loc[1]);
      z[i] = (float)(parts[idx].x[2] - loc[2]);
      h[i] = parts[idx].h;
      m[i] = parts[idx].mass;
     }

    /* Pad cache with fake particles that exist outside the cell so will not
     * interact. We use values of the same magnitude (but negative!) as the real
     * particles to avoid overflow problems. */
    const double max_dx = ci->dx_max_part;
    const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                                 -(2. * ci->width[1] + max_dx),
                                 -(2. * ci->width[2] + max_dx)};
    const float h_padded = ci->parts[0].h;

    for (int i = *last_pi; i < *last_pi + VEC_SIZE; i++) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];
      h[i] = h_padded;
      m[i] = 1.f;
    }
  }
  /* The cell is on the left so read the particles
   * into the cache from the end of the cell. */
  else {
    const int rem = (ci->count - *first_pi) % VEC_SIZE;
    if (rem != 0) {
      const int pad = VEC_SIZE - rem;

      /* Decrease first_pi if there are particles in the cell left to read. */
      if (*first_pi - pad >= 0) *first_pi -= pad;
    }

    const int ci_cache_count = ci->count - *first_pi;

    /* Shift the particles positions to a local frame so single precision can be
     * used instead of double precision. */
    for (int i = 0; i < ci_cache_count; i++) {
      const int idx = sort_i[i + *first_pi].i;
      x[i] = (float)(parts[idx].x[0] - loc[0]);
      y[i] = (float)(parts[idx].x[1] - loc[1]);
      z[i] = (float)(parts[idx].x[2] - loc[2]);
      h[i] = parts[idx].h;
      m[i] = parts[idx].mass;
    }

    /* Pad cache with fake particles that exist outside the cell so will not
     * interact. We use values of the same magnitude (but negative!) as the real
     * particles to avoid overflow problems. */
    const double max_dx = ci->dx_max_part;
    const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                                 -(2. * ci->width[1] + max_dx),
                                 -(2. * ci->width[2] + max_dx)};
    const float h_padded = ci->parts[0].h;

    for (int i = ci->count - *first_pi; i < ci->count - *first_pi + VEC_SIZE;
         i++) {
      x[i] = pos_padded[0];
      y[i] = pos_padded[1];
      z[i] = pos_padded[2];
      h[i] = h_padded;
      m[i] = 1.f;
    }
  }

}

/**
 * @brief Populate cache for force interactions by reading in the particles in
 * unsorted order.
 *
 * @param ci The #cell.
 * @param ci_cache The cache.
 */
__attribute__((always_inline)) INLINE void cache_read_force_particles(
    const struct cell *restrict const ci,
    struct cache *restrict const ci_cache, const int ci_count) {


  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vx, ci_cache->vx, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vy, ci_cache->vy, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vz, ci_cache->vz, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, rho, ci_cache->rho, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pressure, ci_cache->pressure,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, grad_h, ci_cache->grad_h,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, soundspeed, ci_cache->soundspeed,
                            SWIFT_CACHE_ALIGNMENT);

  const struct part *restrict parts = ci->parts;
  const double loc[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};

  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
  for (int i = 0; i < ci_count; i++) {
    x[i] = (float)(parts[i].x[0] - loc[0]);
    y[i] = (float)(parts[i].x[1] - loc[1]);
    z[i] = (float)(parts[i].x[2] - loc[2]);
    h[i] = parts[i].h;
    m[i] = parts[i].mass;
    vx[i] = parts[i].v[0];
    vy[i] = parts[i].v[1];
    vz[i] = parts[i].v[2];
    rho[i] = parts[i].rho;
    pressure[i] = parts[i].force.pressure;
    grad_h[i] = parts[i].force.f;
    soundspeed[i] = parts[i].force.soundspeed;
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  const double max_dx = ci->dx_max_part;
  const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                               -(2. * ci->width[1] + max_dx),
                               -(2. * ci->width[2] + max_dx)};
  const float h_padded = ci->parts[0].h;
  const int ci_count_padded = ci_count - (ci_count % VEC_SIZE) + VEC_SIZE;

  for (int i = ci_count; i < ci_count_padded; i++) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    h[i] = h_padded;
    m[i] = 1.f;
    vx[i] = 1.f;
    vy[i] = 1.f;
    vz[i] = 1.f;
    rho[i] = 1.f;
    pressure[i] = 1.f;
    grad_h[i] = 1.f;
    soundspeed[i] = 1.f;
  }
}

/**
 * @brief Populate caches by only reading particles that are within range of
 * each other within the adjoining cell.Also read the particles into the cache
 * in sorted order.
 *
 * @param ci The i #cell.
 * @param cj The j #cell.
 * @param ci_cache The #cache for cell ci.
 * @param cj_cache The #cache for cell cj.
 * @param sort_i The array of sorted particle indices for cell ci.
 * @param sort_j The array of sorted particle indices for cell ci.
 * @param shift The amount to shift the particle positions to account for BCs
 * @param first_pi The first particle in cell ci that is in range.
 * @param last_pj The last particle in cell cj that is in range.
 */
__attribute__((always_inline)) INLINE void cache_read_two_partial_cells_sorted(
    const struct cell *restrict const ci, const struct cell *restrict const cj,
    struct cache *restrict const ci_cache,
    struct cache *restrict const cj_cache, const struct entry *restrict sort_i,
    const struct entry *restrict sort_j, const double *restrict const shift,
    int *first_pi, int *last_pj) {

  /* Make the number of particles to be read a multiple of the vector size.
   * This eliminates serial remainder loops where possible when populating the
   * cache. */

  /* Is the number of particles to read a multiple of the vector size? */
  int rem = (ci->count - *first_pi) % VEC_SIZE;
  if (rem != 0) {
    int pad = VEC_SIZE - rem;

    /* Decrease first_pi if there are particles in the cell left to read. */
    if (*first_pi - pad >= 0) *first_pi -= pad;
  }

  rem = (*last_pj + 1) % VEC_SIZE;
  if (rem != 0) {
    int pad = VEC_SIZE - rem;

    /* Increase last_pj if there are particles in the cell left to read. */
    if (*last_pj + pad < cj->count) *last_pj += pad;
  }

  /* Get some local pointers */
  const int first_pi_align = *first_pi;
  const int last_pj_align = *last_pj;
  const struct part *restrict parts_i = ci->parts;
  const struct part *restrict parts_j = cj->parts;

  /* Shift particles to the local frame and account for boundary conditions.*/
  const double total_ci_shift[3] = {
      cj->loc[0] + shift[0], cj->loc[1] + shift[1], cj->loc[2] + shift[2]};
  const double total_cj_shift[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  
  int ci_cache_count = ci->count - first_pi_align;

  /* Shift the particles positions to a local frame (ci frame) so single
   * precision can be used instead of double precision.  */
  for (int i = 0; i < ci_cache_count; i++) {
    const int idx = sort_i[i + first_pi_align].i;
    x[i] = (float)(parts_i[idx].x[0] - total_ci_shift[0]);
    y[i] = (float)(parts_i[idx].x[1] - total_ci_shift[1]);
    z[i] = (float)(parts_i[idx].x[2] - total_ci_shift[2]);
    h[i] = parts_i[idx].h;
    m[i] = parts_i[idx].mass;
  }

#ifdef SWIFT_DEBUG_CHECKS
  const float shift_threshold_x =
      2. * ci->width[0] + 2. * max(ci->dx_max_part, cj->dx_max_part);
  const float shift_threshold_y =
      2. * ci->width[1] + 2. * max(ci->dx_max_part, cj->dx_max_part);
  const float shift_threshold_z =
      2. * ci->width[2] + 2. * max(ci->dx_max_part, cj->dx_max_part);

  /* Make sure that particle positions have been shifted correctly. */
  for (int i = 0; i < ci_cache_count; i++) {
    if (x[i] > shift_threshold_x || x[i] < -shift_threshold_x)
      error(
          "Error: ci->loc[%lf,%lf,%lf],cj->loc[%lf,%lf,%lf] Particle %d x pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. x=%f, ci->width[0]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, x[i], ci->width[0]);
    if (y[i] > shift_threshold_y || y[i] < -shift_threshold_y)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d y pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. y=%f, ci->width[1]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, y[i], ci->width[1]);
    if (z[i] > shift_threshold_z || z[i] < -shift_threshold_z)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d z pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. z=%f, ci->width[2]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, z[i], ci->width[2]);
  }
#endif

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  const double max_dx = max(ci->dx_max_part, cj->dx_max_part);
  const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                               -(2. * ci->width[1] + max_dx),
                               -(2. * ci->width[2] + max_dx)};
  const float h_padded = ci->parts[0].h;

  for (int i = ci->count - first_pi_align;
       i < ci->count - first_pi_align + VEC_SIZE; i++) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    h[i] = h_padded;
    m[i] = 1.f;
  }

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, xj, cj_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, yj, cj_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, zj, cj_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, hj, cj_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, mj, cj_cache->m, SWIFT_CACHE_ALIGNMENT);
  
  for (int i = 0; i <= last_pj_align; i++) {
    const int idx = sort_j[i].i;
    xj[i] = (float)(parts_j[idx].x[0] - total_cj_shift[0]);
    yj[i] = (float)(parts_j[idx].x[1] - total_cj_shift[1]);
    zj[i] = (float)(parts_j[idx].x[2] - total_cj_shift[2]);
    hj[i] = parts_j[idx].h;
    mj[i] = parts_j[idx].mass;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure that particle positions have been shifted correctly. */
  for (int i = 0; i <= last_pj_align; i++) {
    if (xj[i] > shift_threshold_x || xj[i] < -shift_threshold_x)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d xj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. xj=%f, ci->width[0]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, xj[i], ci->width[0]);
    if (yj[i] > shift_threshold_y || yj[i] < -shift_threshold_y)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d yj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. yj=%f, ci->width[1]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, yj[i], ci->width[1]);
    if (zj[i] > shift_threshold_z || zj[i] < -shift_threshold_z)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d zj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. zj=%f, ci->width[2]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, zj[i], ci->width[2]);
  }
#endif

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  const float pos_padded_j[3] = {-(2. * cj->width[0] + max_dx),
                                 -(2. * cj->width[1] + max_dx),
                                 -(2. * cj->width[2] + max_dx)};
  const float h_padded_j = cj->parts[0].h;

  for (int i = last_pj_align + 1; i < last_pj_align + 1 + VEC_SIZE; i++) {
    xj[i] = pos_padded_j[0];
    yj[i] = pos_padded_j[1];
    zj[i] = pos_padded_j[2];
    hj[i] = h_padded_j;
    mj[i] = 1.f;
  }
}

/**
 * @brief Populate caches by only reading particles that are within range of
 * each other within the adjoining cell.Also read the particles into the cache
 * in sorted order.
 *
 * @param ci The i #cell.
 * @param cj The j #cell.
 * @param ci_cache The #cache for cell ci.
 * @param cj_cache The #cache for cell cj.
 * @param sort_i The array of sorted particle indices for cell ci.
 * @param sort_j The array of sorted particle indices for cell ci.
 * @param shift The amount to shift the particle positions to account for BCs
 * @param first_pi The first particle in cell ci that is in range.
 * @param last_pj The last particle in cell cj that is in range.
 */
__attribute__((always_inline)) INLINE void
cache_read_two_partial_cells_sorted_force(
    const struct cell *const ci, const struct cell *const cj,
    struct cache *const ci_cache, struct cache *const cj_cache,
    const struct entry *restrict sort_i, const struct entry *restrict sort_j,
    const double *const shift, int *first_pi, int *last_pj) {

  /* Make the number of particles to be read a multiple of the vector size.
   * This eliminates serial remainder loops where possible when populating the
   * cache. */

  /* Is the number of particles to read a multiple of the vector size? */
  int rem = (ci->count - *first_pi) % VEC_SIZE;
  if (rem != 0) {
    int pad = VEC_SIZE - rem;

    /* Decrease first_pi if there are particles in the cell left to read. */
    if (*first_pi - pad >= 0) *first_pi -= pad;
  }

  rem = (*last_pj + 1) % VEC_SIZE;
  if (rem != 0) {
    int pad = VEC_SIZE - rem;

    /* Increase last_pj if there are particles in the cell left to read. */
    if (*last_pj + pad < cj->count) *last_pj += pad;
  }

  /* Get some local pointers */
  const int first_pi_align = *first_pi;
  const int last_pj_align = *last_pj;
  const struct part *restrict parts_i = ci->parts;
  const struct part *restrict parts_j = cj->parts;

  /* Shift particles to the local frame and account for boundary conditions.*/
  const double total_ci_shift[3] = {
      cj->loc[0] + shift[0], cj->loc[1] + shift[1], cj->loc[2] + shift[2]};
  const double total_cj_shift[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, m, ci_cache->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vx, ci_cache->vx, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vy, ci_cache->vy, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vz, ci_cache->vz, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, rho, ci_cache->rho, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pressure, ci_cache->pressure,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, grad_h, ci_cache->grad_h,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, soundspeed, ci_cache->soundspeed,
                            SWIFT_CACHE_ALIGNMENT);

  int ci_cache_count = ci->count - first_pi_align;
  /* Shift the particles positions to a local frame (ci frame) so single
   * precision can be  used instead of double precision.  */
  for (int i = 0; i < ci_cache_count; i++) {

    const int idx = sort_i[i + first_pi_align].i;
    x[i] = (float)(parts_i[idx].x[0] - total_ci_shift[0]);
    y[i] = (float)(parts_i[idx].x[1] - total_ci_shift[1]);
    z[i] = (float)(parts_i[idx].x[2] - total_ci_shift[2]);
    h[i] = parts_i[idx].h;
    vx[i] = parts_i[idx].v[0];
    vy[i] = parts_i[idx].v[1];
    vz[i] = parts_i[idx].v[2];
    m[i] = parts_i[idx].mass;
    rho[i] = parts_i[idx].rho;
    pressure[i] = parts_i[idx].force.pressure;
    grad_h[i] = parts_i[idx].force.f;
    soundspeed[i] = parts_i[idx].force.soundspeed;
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  const double max_dx = max(ci->dx_max_part, cj->dx_max_part);
  const float pos_padded[3] = {-(2. * ci->width[0] + max_dx),
                               -(2. * ci->width[1] + max_dx),
                               -(2. * ci->width[2] + max_dx)};
  const float h_padded = ci->parts[0].h;

  for (int i = ci->count - first_pi_align;
       i < ci->count - first_pi_align + VEC_SIZE; i++) {
    x[i] = pos_padded[0];
    y[i] = pos_padded[1];
    z[i] = pos_padded[2];
    h[i] = h_padded;
    m[i] = 1.f;
    vx[i] = 1.f;
    vy[i] = 1.f;
    vz[i] = 1.f;
    rho[i] = 1.f;
    pressure[i] = 1.f;
    grad_h[i] = 1.f;
    soundspeed[i] = 1.f;
  }

  /* Let the compiler know that the data is aligned and create pointers to the
   * arrays inside the cache. */
  swift_declare_aligned_ptr(float, xj, cj_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, yj, cj_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, zj, cj_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, hj, cj_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, mj, cj_cache->m, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vxj, cj_cache->vx, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vyj, cj_cache->vy, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, vzj, cj_cache->vz, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, rhoj, cj_cache->rho, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, pressurej, cj_cache->pressure,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, grad_hj, cj_cache->grad_h,
                            SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, soundspeedj, cj_cache->soundspeed,
                            SWIFT_CACHE_ALIGNMENT);

  for (int i = 0; i <= last_pj_align; i++) {
    const int idx = sort_j[i].i;
    xj[i] = (float)(parts_j[idx].x[0] - total_cj_shift[0]);
    yj[i] = (float)(parts_j[idx].x[1] - total_cj_shift[1]);
    zj[i] = (float)(parts_j[idx].x[2] - total_cj_shift[2]);
    hj[i] = parts_j[idx].h;
    vxj[i] = parts_j[idx].v[0];
    vyj[i] = parts_j[idx].v[1];
    vzj[i] = parts_j[idx].v[2];
    mj[i] = parts_j[idx].mass;
    rhoj[i] = parts_j[idx].rho;
    pressurej[i] = parts_j[idx].force.pressure;
    grad_hj[i] = parts_j[idx].force.f;
    soundspeedj[i] = parts_j[idx].force.soundspeed;
  }

  /* Pad cache with fake particles that exist outside the cell so will not
   * interact. We use values of the same magnitude (but negative!) as the real
   * particles to avoid overflow problems. */
  const float pos_padded_j[3] = {-(2. * cj->width[0] + max_dx),
                                 -(2. * cj->width[1] + max_dx),
                                 -(2. * cj->width[2] + max_dx)};
  const float h_padded_j = cj->parts[0].h;

  for (int i = last_pj_align + 1; i < last_pj_align + 1 + VEC_SIZE; i++) {
    xj[i] = pos_padded_j[0];
    yj[i] = pos_padded_j[1];
    zj[i] = pos_padded_j[2];
    hj[i] = h_padded_j;
    mj[i] = 1.f;
    vxj[i] = 1.f;
    vyj[i] = 1.f;
    vzj[i] = 1.f;
    rhoj[i] = 1.f;
    pressurej[i] = 1.f;
    grad_hj[i] = 1.f;
    soundspeedj[i] = 1.f;
  }
}

/* @brief Pads the secondary cache so that there are no contributions in the interaction
 * function.  
 *
 * @param int_cache (return) secondary cache of interactions between two
 * particles.
 * @param icount Interaction count.
 * @param icount_padded Interaction count padded to a multiple of two vector lengths.
*/
__attribute__((always_inline)) INLINE static void pad_c2_cache(
    struct c2_cache *const int_cache, const int icount, const int icount_padded) {

  for (int i = icount; i < icount_padded; i++) {
      int_cache->mq[i] = 0.f;
      int_cache->r2q[i] = 1.f;
      int_cache->dxq[i] = 0.f;
      int_cache->dyq[i] = 0.f;
      int_cache->dzq[i] = 0.f;
  }
}

/**
 * @brief Left-packs the values needed by an interaction into the secondary
 * cache (Supports AVX, AVX2 and AVX512 instruction sets).
 *
 * @param mask Contains which particles need to interact.
 * @param pjd Index of the particle to store into.
 * @param v_r2 #vector of the separation between two particles squared.
 * @param v_dx #vector of the x separation between two particles.
 * @param v_dy #vector of the y separation between two particles.
 * @param v_dz #vector of the z separation between two particles.
 * @param cell_cache cache of all particles in the cell.
 * @param int_cache (return) secondary cache of interactions between two
 * particles.
 * @param icount Interaction count.
 */
__attribute__((always_inline)) INLINE static void left_pack_c2_cache(
    const int mask, const int pjd, vector *v_r2, vector *v_dx, vector *v_dy,
    vector *v_dz, const struct cache *const cell_cache,
    struct c2_cache *const int_cache, int *icount) {

  /* Left-pack values needed into the secondary cache using the interaction mask.*/
#if defined(HAVE_AVX2) || defined(HAVE_AVX512_F)
  mask_t packed_mask;
  VEC_FORM_PACKED_MASK(mask, packed_mask);

  VEC_LEFT_PACK(v_r2->v, packed_mask, &int_cache->r2q[*icount]);
  VEC_LEFT_PACK(v_dx->v, packed_mask, &int_cache->dxq[*icount]);
  VEC_LEFT_PACK(v_dy->v, packed_mask, &int_cache->dyq[*icount]);
  VEC_LEFT_PACK(v_dz->v, packed_mask, &int_cache->dzq[*icount]);
  VEC_LEFT_PACK(vec_load(&cell_cache->m[pjd]), packed_mask,
                &int_cache->mq[*icount]);
  
  /* Increment interaction count by number of bits set in mask. */
  (*icount) += __builtin_popcount(mask);
#else
  /* Quicker to do it serially in AVX rather than use intrinsics. */
  for (int bit_index = 0; bit_index < VEC_SIZE; bit_index++) {
    if (mask & (1 << bit_index)) {
      /* Add this interaction to the queue. */
      int_cache->r2q[*icount] = v_r2->f[bit_index];
      int_cache->dxq[*icount] = v_dx->f[bit_index];
      int_cache->dyq[*icount] = v_dy->f[bit_index];
      int_cache->dzq[*icount] = v_dz->f[bit_index];
      int_cache->mq[*icount] = cell_cache->m[pjd + bit_index];
      
      (*icount)++;
    }
  }

#endif /* defined(HAVE_AVX2) || defined(HAVE_AVX512_F) */

}

/* @brief Clean the memory allocated by a #cache object.
 *
 * @param c The #cache to clean.
 */
static INLINE void cache_clean(struct cache *c) {
  if (c->count > 0) {
    free(c->x);
    free(c->y);
    free(c->z);
    free(c->m);
    free(c->vx);
    free(c->vy);
    free(c->vz);
    free(c->h);
    free(c->max_index);
    free(c->rho);
    free(c->pressure);
    free(c->grad_h);
    free(c->soundspeed);
  }
}

#endif /* WITH_VECTORIZATION */

#endif /* SWIFT_CACHE_H */
