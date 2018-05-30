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
#ifndef SWIFT_MINIMAL_HYDRO_CACHE_H
#define SWIFT_MINIMAL_HYDRO_CACHE_H

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

#define NUM_VEC_PROC 2
#define C2_CACHE_SIZE (NUM_VEC_PROC * VEC_SIZE * 6) + (NUM_VEC_PROC * VEC_SIZE)
#define MAX_NUM_OF_CACHE_FIELDS 20
#define NUM_OF_DENSITY_CACHE_FIELDS 1
#define NUM_OF_DENSITY_UPDATE_CACHE_FIELDS 4
#define NUM_OF_FORCE_CACHE_FIELDS 8
#define NUM_OF_FORCE_UPDATE_CACHE_FIELDS 6

#ifdef __ICC
#define PRAGMA_IVDEP _Pragma("ivdep")
#define PRAGMA_NOUNROLL _Pragma("nounroll")
#define PRAGMA_UNROLL _Pragma("unroll")
#define PRAGMA_OMP_SIMD
#else
#define PRAGMA_IVDEP
#define PRAGMA_NOUNROLL
#define PRAGMA_UNROLL
#define PRAGMA_OMP_SIMD _Pragma("omp simd")
#endif

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

/**
 * @brief Specifies which particle fields to read from a dataset (density).
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param ci_cache Particle cache to populate.
 */
INLINE void cache_read_particle_fields_density(const struct part *restrict parts, struct cache_props* list,
                          struct cache *restrict const ci_cache) {

  /* List what we want to read */
  list[0] = cache_make_input_field("mass", parts, mass, ci_cache->m);
}

/**
 * @brief Specifies which particle fields to read from a dataset (force).
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param ci_cache Particle cache to populate.
 */
INLINE void cache_read_particle_fields_force(const struct part *restrict parts, struct cache_props* list,
                          struct cache *restrict const ci_cache) {

  /* List what we want to read */
  list[0] = cache_make_input_field("mass", parts, mass, ci_cache->m);
  list[1] = cache_make_input_field("vx", parts, v[0], ci_cache->vx);
  list[2] = cache_make_input_field("vy", parts, v[1], ci_cache->vy);
  list[3] = cache_make_input_field("vz", parts, v[2], ci_cache->vz);
  list[4] = cache_make_input_field("rho", parts, rho, ci_cache->rho);
  list[5] = cache_make_input_field("pressure", parts, force.pressure, ci_cache->pressure);
  list[6] = cache_make_input_field("grad_h", parts, force.f, ci_cache->grad_h);
  list[7] = cache_make_input_field("soundspeed", parts, force.soundspeed, ci_cache->soundspeed);
}

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

/* Cache to hold a list of vectors used to update particle properties after a density interaction. */
struct update_cache_density {
  vector v_rhoSum;
  vector v_rho_dhSum;
  vector v_wcountSum;
  vector v_wcount_dhSum;
};

/* Cache to hold a list of vectors used to update particle properties after a force interaction. */
struct update_cache_force {
  vector v_a_hydro_xSum;
  vector v_a_hydro_ySum;
  vector v_a_hydro_zSum;
  vector v_u_dtSum;
  vector v_h_dtSum;
  vector v_sigSum;
};

/* Input parameters needed for computing the density interaction. */
struct input_params_density {
  vector input[input_params_density_length];
};

/* Input parameters needed for computing the force interaction. */
struct input_params_force {
  vector input[input_params_force_length];
};

INLINE static void cache_read_particle_update_fields_density(const struct part *restrict parts, struct cache_props* list,
    struct update_cache_density *restrict const update_cache) {

  /* List what we want to read */
  list[0] = cache_make_output_field("rho", parts, rho, &update_cache->v_rhoSum.f[0], reduction_add);
  list[1] = cache_make_output_field("rho_dh", parts, density.rho_dh, &update_cache->v_rho_dhSum.f[0], reduction_add);
  list[2] = cache_make_output_field("wcount", parts, density.wcount, &update_cache->v_wcountSum.f[0], reduction_add);
  list[3] = cache_make_output_field("wcount_dh", parts, density.wcount_dh, &update_cache->v_wcount_dhSum.f[0], reduction_add);

}

INLINE static void cache_read_particle_update_fields_force(const struct part *restrict parts, struct cache_props* list,
    struct update_cache_force *restrict const update_cache) {

  /* List what we want to read */
  list[0] = cache_make_output_field("a_hydro_x", parts, a_hydro[0], &update_cache->v_a_hydro_xSum.f[0], reduction_add);
  list[1] = cache_make_output_field("a_hydro_y", parts, a_hydro[1], &update_cache->v_a_hydro_ySum.f[0], reduction_add);
  list[2] = cache_make_output_field("a_hydro_z", parts, a_hydro[2], &update_cache->v_a_hydro_zSum.f[0], reduction_add);
  list[3] = cache_make_output_field("u_dt", parts, u_dt, &update_cache->v_u_dtSum.f[0], reduction_add);
  list[4] = cache_make_output_field("h_dt", parts, force.h_dt, &update_cache->v_h_dtSum.f[0], reduction_add);
  list[5] = cache_make_output_field("v_sig", parts, force.v_sig, &update_cache->v_sigSum.f[0], reduction_max);

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

#endif /* SWIFT_MINIMAL_HYDRO_CACHE_H */
