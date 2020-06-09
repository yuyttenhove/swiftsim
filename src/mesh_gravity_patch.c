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

/* Config parameters. */
#include "../config.h"

#include <math.h>

/* This object's header. */
#include "mesh_gravity_patch.h"

/* Local includes. */
#include "hashmap.h"
#include "cell.h"
#include "error.h"
#include "row_major_id.h"

/**
 * @brief Initialize a mesh patch to cover a cell
 *
 * @param patch A pointer to the mesh patch
 * @param cell The cell which the mesh should cover
 * @param fac Inverse of the FFT mesh size
 * @param dim Size of the full volume in each dimension
 * @param boundary_size Size of the boundary layer to include
 */
void pm_mesh_patch_init(struct pm_mesh_patch *patch, const struct cell *cell,
                        const int N, const double fac, const double dim[3],
                        const int boundary_size) {

  const int gcount = cell->grav.count;
  const struct gpart *gp = cell->grav.parts;

  patch->N = N;
  patch->fac = fac;

  /* Will need to wrap particles to position nearest the cell centre */
  for(int i=0; i<3; i+=1) {
    patch->wrap_min[i] = cell->loc[i] + 0.5*cell->width[i] - 0.5*dim[i];
    patch->wrap_max[i] = cell->loc[i] + 0.5*cell->width[i] + 0.5*dim[i];
  }

  /* Find the extent of the particle distribution in the cell */
  double pos_min[3];
  double pos_max[3];
  for(int i=0; i<3; i+=1) {
    pos_min[i] = patch->wrap_max[i];
    pos_max[i] = patch->wrap_min[i];
  }
  for (int ipart = 0; ipart < gcount; ipart += 1) {
    for(int i=0; i<3; i+=1) {
      const double pos_wrap = box_wrap(gp[ipart].x[i], patch->wrap_min[i], patch->wrap_max[i]);
      if(pos_wrap < pos_min[i])pos_min[i] = pos_wrap;
      if(pos_wrap > pos_max[i])pos_max[i] = pos_wrap;
    }
  }

  /* Determine the integer size and coordinates of the mesh */
  int num_cells = 1;
  for(int i=0; i<3; i+=1) {
    patch->mesh_min[i] = floor(pos_min[i]*fac) - boundary_size;
    /* CIC interpolation requires one extra element in the positive direction */
    patch->mesh_max[i] = floor(pos_max[i]*fac) + boundary_size + 1;
    patch->mesh_size[i] = patch->mesh_max[i] - patch->mesh_min[i] + 1;
    num_cells *= patch->mesh_size[i];
  }

  /* Allocate the mesh */
  if (swift_memalign("mesh_patch", (void **) &patch->mesh, 32,
                     num_cells * sizeof(double)) != 0)
    error("Failed to allocate array for mesh patch!");

  return;
}

/**
 * @brief Set all values in a mesh patch to zero
 *
 * @param patch A pointer to the mesh patch
 */
void pm_mesh_patch_zero(struct pm_mesh_patch *patch) {

  int num = patch->mesh_size[0]*patch->mesh_size[1]*patch->mesh_size[2];
  for(int i=0; i<num; i+=1)
    patch->mesh[i] = 0.0;
}

/**
 * @brief Initialize mesh values using a hashmap
 *
 * @param patch Pointer to the pm_mesh_patch
 * @param map Pointer to the hashmap
 */
void pm_mesh_patch_set_values_from_hashmap(struct pm_mesh_patch *patch, hashmap_t *map) {

  /* Loop over all cells in the patch */
  for(int i=0; i<patch->mesh_size[0];i+=1) {
    for(int j=0; j<patch->mesh_size[1];j+=1) {
      for(int k=0; k<patch->mesh_size[2];k+=1) {

        /* Find array index in the mesh patch */
	const int local_index = pm_mesh_patch_index(patch, i, j, k);

        /* Find index in the full mesh */
	const size_t global_index = row_major_id_periodic_size_t_padded(
            i+patch->mesh_min[0], j+patch->mesh_min[1], k+patch->mesh_min[2], patch->N);

        /* Look up the value in the hashmap and store it in the mesh patch */
        hashmap_value_t *value = hashmap_lookup(map, global_index);
        if(!value) {
          error("Required cell is not present in potential hashmap");
        } else {
          patch->mesh[local_index] = value->value_dbl;
        }
      }
    }
  }
}

/**
 * @brief Accumulate values from the patch to a hashmap
 *
 * @param patch Pointer to the pm_mesh_patch
 * @param map Pointer to the hashmap
 */
void pm_mesh_patch_add_values_to_hashmap(struct pm_mesh_patch *patch, hashmap_t *map) {

  /* Loop over all cells in the patch */
  for(int i=0; i<patch->mesh_size[0];i+=1) {
    for(int j=0; j<patch->mesh_size[1];j+=1) {
      for(int k=0; k<patch->mesh_size[2];k+=1) {

        /* Find array index in the mesh patch */        
	const int local_index = pm_mesh_patch_index(patch, i, j, k);

        /* Find index in the full mesh */
	const size_t global_index = row_major_id_periodic_size_t_padded(
            i+patch->mesh_min[0], j+patch->mesh_min[1], k+patch->mesh_min[2], patch->N);

        /* Increment the hashmap entry, creating it if necessary */
        int created = 0;
        hashmap_value_t *value = hashmap_get_new(map, global_index, &created);
        if (created) {
          /* Key was not present, so this is a new element */
          value->value_dbl = patch->mesh[local_index];
        } else {
          /* Key was present, so add m_add to previous value */
          value->value_dbl += patch->mesh[local_index];
        }
      }
    }
  }
}

/**
 * @brief Free the memory associated with a mesh patch.
 *
 * @param patch A pointer to the mesh patch
 */
void pm_mesh_patch_clean(struct pm_mesh_patch *patch) {
  swift_free("mesh_patch", patch->mesh);
}


__attribute__((always_inline, const)) INLINE static double mesh_patch_CIC_get(
    const struct pm_mesh_patch *patch, const int i, const int j, const int k,
    const double tx, const double ty, const double tz, const double dx,
    const double dy, const double dz) {
  double temp;
  temp  = patch->mesh[pm_mesh_patch_index(patch, i + 0, j + 0, k + 0)] * tx * ty * tz;
  temp += patch->mesh[pm_mesh_patch_index(patch, i + 0, j + 0, k + 1)] * tx * ty * dz;
  temp += patch->mesh[pm_mesh_patch_index(patch, i + 0, j + 1, k + 0)] * tx * dy * tz;
  temp += patch->mesh[pm_mesh_patch_index(patch, i + 0, j + 1, k + 1)] * tx * dy * dz;
  temp += patch->mesh[pm_mesh_patch_index(patch, i + 1, j + 0, k + 0)] * dx * ty * tz;
  temp += patch->mesh[pm_mesh_patch_index(patch, i + 1, j + 0, k + 1)] * dx * ty * dz;
  temp += patch->mesh[pm_mesh_patch_index(patch, i + 1, j + 1, k + 0)] * dx * dy * tz;
  temp += patch->mesh[pm_mesh_patch_index(patch, i + 1, j + 1, k + 1)] * dx * dy * dz;
  return temp;
}

/**
 * @brief Evaluate acceleration and potential at a position using CIC
 *
 * @param patch Pointer to the pm_mesh_patch
 * @param pos Position at which to evaluate
 * @param accel Output array for acceleration vector
 * @param pot Output double for potential
 */

void pm_mesh_patch_evaluate(struct pm_mesh_patch *patch, double pos[3], double *accel, double *pot) {
  
  const double fac = patch->fac;

  /* Box wrap the position to nearest the cell centre */
  const double pos_x = box_wrap(pos[0], patch->wrap_min[0], patch->wrap_max[0]);
  const double pos_y = box_wrap(pos[1], patch->wrap_min[1], patch->wrap_max[1]);
  const double pos_z = box_wrap(pos[2], patch->wrap_min[2], patch->wrap_max[2]);

  /* Workout the CIC coefficients */
  int i = (int) floor(fac * pos_x);
  const double dx = fac * pos_x - i;
  const double tx = 1. - dx;

  int j = (int) floor(fac * pos_y);
  const double dy = fac * pos_y - j;
  const double ty = 1. - dy;
    
  int k = (int) floor(fac * pos_z);
  const double dz = fac * pos_z - k;
  const double tz = 1. - dz;

  /* Indices of (i,j,k) in the local copy of the mesh */
  const int ii = i - patch->mesh_min[0];
  const int jj = j - patch->mesh_min[1];
  const int kk = k - patch->mesh_min[2];

  /* Simple CIC for the potential itself */
  *pot = mesh_patch_CIC_get(patch, ii, jj, kk, tx, ty, tz, dx, dy, dz);

  accel[0] = 0.0;
  accel[1] = 0.0;
  accel[2] = 0.0;

  /* 5-point stencil along each axis for the accelerations */
  accel[0] += (1. / 12.) * mesh_patch_CIC_get(patch, ii + 2, jj, kk, tx, ty, tz, dx, dy, dz);
  accel[0] -= (2. / 3.)  * mesh_patch_CIC_get(patch, ii + 1, jj, kk, tx, ty, tz, dx, dy, dz);
  accel[0] += (2. / 3.)  * mesh_patch_CIC_get(patch, ii - 1, jj, kk, tx, ty, tz, dx, dy, dz);
  accel[0] -= (1. / 12.) * mesh_patch_CIC_get(patch, ii - 2, jj, kk, tx, ty, tz, dx, dy, dz);
  
  accel[1] += (1. / 12.) * mesh_patch_CIC_get(patch, ii, jj + 2, kk, tx, ty, tz, dx, dy, dz);
  accel[1] -= (2. / 3.)  * mesh_patch_CIC_get(patch, ii, jj + 1, kk, tx, ty, tz, dx, dy, dz);
  accel[1] += (2. / 3.)  * mesh_patch_CIC_get(patch, ii, jj - 1, kk, tx, ty, tz, dx, dy, dz);
  accel[1] -= (1. / 12.) * mesh_patch_CIC_get(patch, ii, jj - 2, kk, tx, ty, tz, dx, dy, dz);
  
  accel[2] += (1. / 12.) * mesh_patch_CIC_get(patch, ii, jj, kk + 2, tx, ty, tz, dx, dy, dz);
  accel[2] -= (2. / 3.)  * mesh_patch_CIC_get(patch, ii, jj, kk + 1, tx, ty, tz, dx, dy, dz);
  accel[2] += (2. / 3.)  * mesh_patch_CIC_get(patch, ii, jj, kk - 1, tx, ty, tz, dx, dy, dz);
  accel[2] -= (1. / 12.) * mesh_patch_CIC_get(patch, ii, jj, kk - 2, tx, ty, tz, dx, dy, dz);
}
