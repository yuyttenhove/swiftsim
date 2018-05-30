/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James Willis (james.s.willis@durham.ac.uk)
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
#include "generic_cache.h"
#include "hydro.h"

#ifdef WITH_VECTORIZATION

/**
 * @brief Reset the update cache to zero.
 *
 * @param num_fields No. of fields to reset.
 * @param props #cache_props which points to update cache
 */
__attribute__((always_inline)) INLINE void update_cache_init(const int num_fields, struct cache_props *props) {

  for(int i=0; i<num_fields; i++) ((vector *)props[i].cache_addr)->v = vec_setzero();

}

/**
 * @brief Perform reduction operations on sum vectors and store result in particle pi.
 *
 * @param props #cache_props which points to update cache
 * @param pid Particle offset into particle array.
 * @param num_fields No. of fields to reset.
 */
__attribute__((always_inline)) INLINE void update_particle(struct cache_props *props, const int pid, const int num_fields) {

  const size_t sizePart = sizeof(struct part);

  for(int i=0; i<num_fields; i++) {
    props[i].reduction_f(*(vector*)props[i].cache_addr, (float *)&props[i].field[pid*sizePart]);
  }

}

/**
 * @brief Populate cache by reading in the particles in unsorted order.
 *
 * @param ci The #cell.
 * @param ci_cache The cache.
 * @param ci_count ci particle count.
 */
__attribute__((always_inline)) INLINE void cache_read_particles(
    const struct cell *restrict const ci,
    struct cache *restrict const ci_cache, const int ci_count) {

  const double loc[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};
  const struct part *restrict parts = ci->parts;
  const size_t sizePart = sizeof(struct part);

  /* Construct a list of particle fields to read into the cache. */
  struct cache_props props[NUM_OF_DENSITY_CACHE_FIELDS];
  cache_read_particle_fields_density(parts, props, ci_cache);

  /* Construct a list of pointers to the arrays inside the cache. */
  float *restrict fields[MAX_NUM_OF_CACHE_FIELDS];  

  for(int i=0; i<NUM_OF_DENSITY_CACHE_FIELDS; i++) {
    fields[i] = props[i].cache_addr;
    swift_align_information(float, fields[i], SWIFT_CACHE_ALIGNMENT);
  }
 
  /* Let the compiler know that the data is aligned. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_align_information_loop(fields, MAX_NUM_OF_CACHE_FIELDS);

  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
PRAGMA_IVDEP
PRAGMA_OMP_SIMD
  for (int i = 0; i < ci_count; i++) {
    x[i] = (float)(parts[i].x[0] - loc[0]);
    y[i] = (float)(parts[i].x[1] - loc[1]);
    z[i] = (float)(parts[i].x[2] - loc[2]);
    h[i] = parts[i].h;
PRAGMA_UNROLL
    for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) {
      fields[j][i] = *(float *)&(props[j].field[i*sizePart]);
    }
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
    for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) fields[j][i] = 1.f;
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

  const size_t sizePart = sizeof(struct part);
  const struct part *restrict parts = ci->parts;
  
  /* Construct a list of particle fields to read into the cache. */
  struct cache_props props[NUM_OF_DENSITY_CACHE_FIELDS];
  cache_read_particle_fields_density(parts, props, ci_cache);
  
  /* Construct a list of pointers to the arrays inside the cache. */
  float *restrict fields[NUM_OF_DENSITY_CACHE_FIELDS];  

  for(int i=0; i<NUM_OF_DENSITY_CACHE_FIELDS; i++) {
    fields[i] = props[i].cache_addr;
    swift_align_information(float, fields[i], SWIFT_CACHE_ALIGNMENT);
  }
 
  /* The cell is on the right so read the particles
   * into the cache from the start of the cell. */
  if (!flipped) {
    const int rem = (*last_pi + 1) % VEC_SIZE;
    if (rem != 0) {
      const int pad = VEC_SIZE - rem;

      /* Increase last_pi if there are particles in the cell left to read. */
      if (*last_pi + pad < ci->count) *last_pi += pad;
    }

    /* Let the compiler know that the data is aligned. */
    swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
    swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
    swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
    swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
    swift_align_information_loop(fields, MAX_NUM_OF_CACHE_FIELDS);

    /* Shift the particles positions to a local frame so single precision can be
     * used instead of double precision. */
PRAGMA_IVDEP
    for (int i = 0; i < *last_pi; i++) {
      const int idx = sort_i[i].i;
      x[i] = (float)(parts[idx].x[0] - loc[0]);
      y[i] = (float)(parts[idx].x[1] - loc[1]);
      z[i] = (float)(parts[idx].x[2] - loc[2]);
      h[i] = parts[idx].h;
      for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) {
        fields[j][i] = *(float *)&(props[j].field[idx*sizePart]);
      }
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
      for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) fields[j][i] = 1.f;
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

    /* Let the compiler know that the data is aligned. */
    swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
    swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
    swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
    swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);

    /* Shift the particles positions to a local frame so single precision can be
     * used instead of double precision. */
PRAGMA_IVDEP
    for (int i = 0; i < ci_cache_count; i++) {
      const int idx = sort_i[i + *first_pi].i;
      x[i] = (float)(parts[idx].x[0] - loc[0]);
      y[i] = (float)(parts[idx].x[1] - loc[1]);
      z[i] = (float)(parts[idx].x[2] - loc[2]);
      h[i] = parts[idx].h;
      for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) {
        fields[j][i] = *(float *)&(props[j].field[idx*sizePart]);
      }
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
      for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) fields[j][i] = 1.f;
    }
  }

}

/**
 * @brief Populate cache for force interactions by reading in the particles in
 * unsorted order.
 *
 * @param ci The #cell.
 * @param ci_cache The cache.
 * @param ci_count ci particle count.
 */
__attribute__((always_inline)) INLINE void cache_read_force_particles(
    const struct cell *restrict const ci,
    struct cache *restrict const ci_cache, const int ci_count) {

  const double loc[3] = {ci->loc[0], ci->loc[1], ci->loc[2]};
  const struct part *restrict parts = ci->parts;
  const size_t sizePart = sizeof(struct part);
  
  /* Construct a list of particle fields to read into the cache. */
  struct cache_props props[NUM_OF_FORCE_CACHE_FIELDS];
  cache_read_particle_fields_force(parts, props, ci_cache);

  /* Construct a list of particle fields to read into the cache. */
  float *restrict fields[MAX_NUM_OF_CACHE_FIELDS];  

  for(int i=0; i<NUM_OF_FORCE_CACHE_FIELDS; i++) {
    fields[i] = props[i].cache_addr;
    swift_align_information(float, fields[i], SWIFT_CACHE_ALIGNMENT);
  }

  /* Let the compiler know that the data is aligned. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_align_information_loop(fields, MAX_NUM_OF_CACHE_FIELDS);
  
  /* Shift the particles positions to a local frame so single precision can be
   * used instead of double precision. */
PRAGMA_IVDEP
  for (int i = 0; i < ci_count; i++) {
    x[i] = (float)(parts[i].x[0] - loc[0]);
    y[i] = (float)(parts[i].x[1] - loc[1]);
    z[i] = (float)(parts[i].x[2] - loc[2]);
    h[i] = parts[i].h;
PRAGMA_UNROLL
    for(int j = 0; j < NUM_OF_FORCE_CACHE_FIELDS; j++) {
      fields[j][i] = *(float *)&(props[j].field[i*sizePart]);
    }
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
    for(int j = 0; j < NUM_OF_FORCE_CACHE_FIELDS; j++) fields[j][i] = 1.f;
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
  const size_t sizePart = sizeof(struct part);

  /* Shift particles to the local frame and account for boundary conditions.*/
  const double total_ci_shift[3] = {
      cj->loc[0] + shift[0], cj->loc[1] + shift[1], cj->loc[2] + shift[2]};
  const double total_cj_shift[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  /* Construct a list of particle fields to read into the cache. */
  struct cache_props props[NUM_OF_DENSITY_CACHE_FIELDS];
  cache_read_particle_fields_density(parts_i, props, ci_cache);

  /* Construct a list of pointers to the arrays inside the cache. */
  float *restrict fields[NUM_OF_DENSITY_CACHE_FIELDS];  

  for(int i=0; i<NUM_OF_DENSITY_CACHE_FIELDS; i++) {
    fields[i] = props[i].cache_addr;
    swift_align_information(float, fields[i], SWIFT_CACHE_ALIGNMENT);
  }
  
  int ci_cache_count = ci->count - first_pi_align;

  /* Let the compiler know that the data is aligned. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_align_information_loop(fields, MAX_NUM_OF_CACHE_FIELDS);
  
  /* Shift the particles positions to a local frame (ci frame) so single
   * precision can be used instead of double precision.  */
PRAGMA_IVDEP
  for (int i = 0; i < ci_cache_count; i++) {
    const int idx = sort_i[i + first_pi_align].i;
    x[i] = (float)(parts_i[idx].x[0] - total_ci_shift[0]);
    y[i] = (float)(parts_i[idx].x[1] - total_ci_shift[1]);
    z[i] = (float)(parts_i[idx].x[2] - total_ci_shift[2]);
    h[i] = parts_i[idx].h;
    for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) {
      fields[j][i] = *(float *)&(props[j].field[idx*sizePart]);
    }
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
    if (fields[0][i] > shift_threshold_x || fields[0][i] < -shift_threshold_x)
      error(
          "Error: ci->loc[%lf,%lf,%lf],cj->loc[%lf,%lf,%lf] Particle %d x pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. x=%f, ci->width[0]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, fields[0][i], ci->width[0]);
    if (fields[1][i] > shift_threshold_y || fields[1][i] < -shift_threshold_y)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d y pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. y=%f, ci->width[1]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, fields[1][i], ci->width[1]);
    if (fields[2][i] > shift_threshold_z || fields[2][i] < -shift_threshold_z)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d z pos "
          "is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. z=%f, ci->width[2]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, fields[2][i], ci->width[2]);
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
    for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) fields[j][i] = 1.f;
  }

  /* Construct a list of particle fields to read into the cache. */
  cache_read_particle_fields_density(parts_j, props, cj_cache);
 
  /* Construct a list of pointers to the arrays inside the cache. */
  for(int i=0; i<NUM_OF_DENSITY_CACHE_FIELDS; i++) {
    fields[i] = props[i].cache_addr;
    swift_align_information(float, fields[i], SWIFT_CACHE_ALIGNMENT);
  }
  
  /* Let the compiler know that the data is aligned. */
  swift_declare_aligned_ptr(float, xj, cj_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, yj, cj_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, zj, cj_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, hj, cj_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_align_information_loop(fields, MAX_NUM_OF_CACHE_FIELDS);

PRAGMA_IVDEP
  for (int i = 0; i <= last_pj_align; i++) {
    const int idx = sort_j[i].i;
    xj[i] = (float)(parts_j[idx].x[0] - total_cj_shift[0]);
    yj[i] = (float)(parts_j[idx].x[1] - total_cj_shift[1]);
    zj[i] = (float)(parts_j[idx].x[2] - total_cj_shift[2]);
    hj[i] = parts_j[idx].h;
    for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) {
      fields[j][i] = *(float *)&(props[j].field[idx*sizePart]);
    }
 }

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure that particle positions have been shifted correctly. */
  for (int i = 0; i <= last_pj_align; i++) {
    if (fields[0][i] > shift_threshold_x || fields[0][i] < -shift_threshold_x)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d xj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. xj=%f, ci->width[0]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, fields[0][i], ci->width[0]);
    if (fields[1][i] > shift_threshold_y || fields[1][i] < -shift_threshold_y)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d yj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. yj=%f, ci->width[1]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, fields[1][i], ci->width[1]);
    if (fields[2][i] > shift_threshold_z || fields[2][i] < -shift_threshold_z)
      error(
          "Error: ci->loc[%lf,%lf,%lf], cj->loc[%lf,%lf,%lf] Particle %d zj "
          "pos is not within "
          "[-4*ci->width*(1 + 2*space_maxreldx), 4*ci->width*(1 + "
          "2*space_maxreldx)]. zj=%f, ci->width[2]=%f",
          ci->loc[0], ci->loc[1], ci->loc[2], cj->loc[0], cj->loc[1],
          cj->loc[2], i, fields[2][i], ci->width[2]);
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
    for(int j = 0; j < NUM_OF_DENSITY_CACHE_FIELDS; j++) fields[j][i] = 1.f;
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
  const size_t sizePart = sizeof(struct part);

  /* Shift particles to the local frame and account for boundary conditions.*/
  const double total_ci_shift[3] = {
      cj->loc[0] + shift[0], cj->loc[1] + shift[1], cj->loc[2] + shift[2]};
  const double total_cj_shift[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  /* Construct a list of particle fields to read into the cache. */
  struct cache_props props[NUM_OF_FORCE_CACHE_FIELDS];
  cache_read_particle_fields_force(parts_i, props, ci_cache);

  /* Construct a list of pointers to the arrays inside the cache. */
  float *restrict fields[NUM_OF_FORCE_CACHE_FIELDS];  

  for(int i=0; i<NUM_OF_FORCE_CACHE_FIELDS; i++) {
    fields[i] = props[i].cache_addr;
    swift_align_information(float, fields[i], SWIFT_CACHE_ALIGNMENT);
  }

  /* Let the compiler know that the data is aligned. */
  swift_declare_aligned_ptr(float, x, ci_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, y, ci_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, z, ci_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, h, ci_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_align_information_loop(fields, MAX_NUM_OF_CACHE_FIELDS);

  int ci_cache_count = ci->count - first_pi_align;
  
  /* Shift the particles positions to a local frame (ci frame) so single
   * precision can be  used instead of double precision.  */
PRAGMA_IVDEP
  for (int i = 0; i < ci_cache_count; i++) {

    const int idx = sort_i[i + first_pi_align].i;
    x[i] = (float)(parts_i[idx].x[0] - total_ci_shift[0]);
    y[i] = (float)(parts_i[idx].x[1] - total_ci_shift[1]);
    z[i] = (float)(parts_i[idx].x[2] - total_ci_shift[2]);
    h[i] = parts_i[idx].h;
    for(int j = 0; j < NUM_OF_FORCE_CACHE_FIELDS; j++) {
      fields[j][i] = *(float *)&(props[j].field[idx*sizePart]);
    }
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
    for(int j = 0; j < NUM_OF_FORCE_CACHE_FIELDS; j++) fields[j][i] = 1.f;
  }

  /* Construct a list of particle fields to read into the cache. */
  cache_read_particle_fields_force(parts_j, props, cj_cache);
 
  /* Construct a list of pointers to the arrays inside the cache. */
  for(int i=0; i<NUM_OF_FORCE_CACHE_FIELDS; i++) {
    fields[i] = props[i].cache_addr;
    swift_align_information(float, fields[i], SWIFT_CACHE_ALIGNMENT);
  }
  
  /* Let the compiler know that the data is aligned. */
  swift_declare_aligned_ptr(float, xj, cj_cache->x, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, yj, cj_cache->y, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, zj, cj_cache->z, SWIFT_CACHE_ALIGNMENT);
  swift_declare_aligned_ptr(float, hj, cj_cache->h, SWIFT_CACHE_ALIGNMENT);
  swift_align_information_loop(fields, MAX_NUM_OF_CACHE_FIELDS);

PRAGMA_IVDEP
  for (int i = 0; i <= last_pj_align; i++) {
    const int idx = sort_j[i].i;
    xj[i] = (float)(parts_j[idx].x[0] - total_cj_shift[0]);
    yj[i] = (float)(parts_j[idx].x[1] - total_cj_shift[1]);
    zj[i] = (float)(parts_j[idx].x[2] - total_cj_shift[2]);
    hj[i] = parts_j[idx].h;
    for(int j = 0; j < NUM_OF_FORCE_CACHE_FIELDS; j++) {
      fields[j][i] = *(float *)&(props[j].field[idx*sizePart]);
    }
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
    for(int j = 0; j < NUM_OF_FORCE_CACHE_FIELDS; j++) fields[j][i] = 1.f;
  }
}

#endif /* WITH_VECTORIZATION */

#endif /* SWIFT_CACHE_H */
