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

#ifndef SWIFT_MESH_GRAVITY_PATCH_H
#define SWIFT_MESH_GRAVITY_PATCH_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "inline.h"
#include "hashmap.h"
#include "error.h"

/* Forward declarations */
struct cell;

/**
 * @brief Data structure for a patch of mesh covering a cell
 */
struct pm_mesh_patch {
  
  /*! Size of the full mesh */
  int N;

  /*! Inverse of the mesh cell size */
  double fac;
  
  /*! Minimum of coordinate range into which particles should be wrapped */
  double wrap_min[3];

  /*! Maximum of coordinate range into which particles should be wrapped */
  double wrap_max[3];

  /*! Size of the mesh in each dimension */
  int mesh_size[3];

  /*! Minimum integer coordinate of the mesh in each dimension */
  int mesh_min[3];

  /*! Maximum integer coordinate of the mesh in each dimension */
  int mesh_max[3];

  /*! Pointer to the mesh data */
  double *mesh;
};

void pm_mesh_patch_init(struct pm_mesh_patch *patch, const struct cell *cell, 
                        const int N, const double fac, const double dim[3],
                        const int boundary_size);

void pm_mesh_patch_zero(struct pm_mesh_patch *patch);

void pm_mesh_patch_set_values_from_hashmap(struct pm_mesh_patch *patch, hashmap_t *map);

void pm_mesh_patch_add_values_to_hashmap(struct pm_mesh_patch *patch, hashmap_t *map);

void pm_mesh_patch_clean(struct pm_mesh_patch *patch);

void pm_mesh_patch_evaluate(struct pm_mesh_patch *patch, double pos[3], double *accel, double *pot);

/**
 * @brief Return the array index in the patch corresponding to
 * coordinates in the full mesh
 *
 * @param patch Pointer to the patch
 * @param i Integer x coordinate in the full mesh
 * @param j Integer y coordinate in the full mesh
 * @param k Integer z coordinate in the full mesh
 *
 */
__attribute__((always_inline)) INLINE static int pm_mesh_patch_index(
      const struct pm_mesh_patch *patch, const int i, const int j, const int k) {

  /* Get coordinates in the mesh patch */
  const int ilocal = i - patch->mesh_min[0];
  const int jlocal = j - patch->mesh_min[1];
  const int klocal = k - patch->mesh_min[2];

  /* #ifdef SWIFT_DEBUG_CHECKS */
  if(ilocal < 0 || ilocal >= patch->mesh_size[0])error("Coordinate in local mesh out of range!");
  if(jlocal < 0 || jlocal >= patch->mesh_size[1])error("Coordinate in local mesh out of range!");
  if(klocal < 0 || klocal >= patch->mesh_size[2])error("Coordinate in local mesh out of range!");
  /* #endif */

  return (ilocal*patch->mesh_size[1]*patch->mesh_size[2]) + (jlocal*patch->mesh_size[2]) + klocal;
}

/**
 * @brief Increment a value in the specified mesh patch
 *
 * @param patch Pointer to the patch
 * @param i Integer x coordinate in the full mesh
 * @param j Integer y coordinate in the full mesh
 * @param k Integer z coordinate in the full mesh
 * @param mass Amount by which to increment the element
 *
 */
__attribute__((always_inline)) INLINE static void pm_mesh_patch_accumulate(
      struct pm_mesh_patch *patch, const int i, const int j, const int k, const double mass) {
  
  /* Calculate array index*/
  const int index = pm_mesh_patch_index(patch, i, j, k);

  /* Add the value */
  patch->mesh[index] += mass;
}

#endif
