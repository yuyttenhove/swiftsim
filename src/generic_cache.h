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
#ifndef SWIFT_GENERIC_CACHE_H
#define SWIFT_GENERIC_CACHE_H

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

#define CACHE_FIELD_BUFFER_SIZE 200

#ifdef WITH_VECTORIZATION

typedef void (*reduction_func)(vector,float*);

/**
 * @brief The properties of a given particle cache for vectorisation
 */
struct cache_props {

  /* Name */
  char name[CACHE_FIELD_BUFFER_SIZE];

  /* Pointer to the field of the first particle in the array */
  char* field;

  /* Pointer to cache field */
  float* cache_addr;

  /* Function pointer to reduction function */
  reduction_func reduction_f;
  
};

/**
 * @brief Constructs an #cache_props from its parameters which will be used to populate the particle cache.
 */
#define cache_make_input_field(name, part, field, cache_addr) \
  cache_make_input_field_(name, (char*)(&(part[0]).field), cache_addr)

/**
 * @brief Construct an #cache_props from its parameters
 *
 * @param name Name of the field to read
 * @param field Pointer to the field of the first particle
 * @param cache_addr Pointer to the field in the particle cache
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE struct cache_props cache_make_input_field_(
    char name[CACHE_FIELD_BUFFER_SIZE], char* field, float* cache_addr) {
  struct cache_props r;
  strcpy(r.name, name);
  r.field = field;
  r.cache_addr = cache_addr;

  return r;
}

/**
 * @brief Constructs an #cache_props from its parameters which will be used to hold particle updates.
 */
#define cache_make_output_field(name, part, field, cache_addr, reduction_op) \
  cache_make_output_field_(name, (char*)(&(part[0]).field), cache_addr, reduction_op)

/**
 * @brief Construct an #cache_props from its parameters
 * 
 * @param name Name of the field to read
 * @param field Pointer to the field of the first particle
 * @param cache_addr Pointer to the field in the particle cache
 * @param functPtr Function pointer to which reduction operation should be used on the particle sum vector
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE struct cache_props cache_make_output_field_(
    char name[CACHE_FIELD_BUFFER_SIZE], char* field, float* cache_addr, reduction_func funcPtr) {
  struct cache_props r;
  strcpy(r.name, name);
  r.field = field;
  r.cache_addr = cache_addr;
  r.reduction_f = funcPtr;

  return r;
}

/**
 * @brief Perform horizontal add on vector and update particle field.
 * 
 * @param field Vector to perform the reduction operation on
 * @param pi_update Pointer to particle field
 */
INLINE static void reduction_add(vector field, float *pi_update) {
  VEC_HADD(field, *pi_update);
}

/**
 * @brief Perform horizontal max on vector and update particle field.
 * 
 * @param field Vector to perform the reduction operation on
 * @param pi_update Pointer to particle field
 */
INLINE static void reduction_max(vector field, float *pi_update) {
  float hmax = 0.f;
  VEC_HMAX(field, hmax);
  *pi_update = max(*pi_update, hmax);
}

#endif /* WITH_VECTORIZATION */

#endif /* SWIFT_GENERIC_CACHE_H */
