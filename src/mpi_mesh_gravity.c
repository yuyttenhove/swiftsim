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

/* This object's header. */
#include "mpi_mesh_gravity.h"

/* Local includes. */
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "part.h"
#include "space.h"

/**
 * @brief Data structure for communicating contributions to the mesh
 */
struct mesh_hashmap_data {
  hashmap_key_t key;
  double value;
}


/**
 * @brief Returns 1D index of a 3D NxNxN array using row-major style.
 *
 * Wraps around in the corresponding dimension if any of the 3 indices is >= N
 * or < 0.
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param N Size of the array along one axis.
 */
__attribute__((always_inline, const)) INLINE static int row_major_id_periodic(
    const int i, const int j, const int k, const int N) {

  return (((i + N) % N) * N * N + ((j + N) % N) * N + ((k + N) % N));
}


/**
 * @brief Increment the value associated with a hashmap key
 *
 * The hashmap entry specified by key is incremented by m_add
 * if it already exists, otherwise it is created and set equal
 * to m_add.
 *
 * @param map Pointer to the hash map.
 * @param key Which hash key to update.
 * @param m_add Amount by which to increment the value in the hash map.
 */
__attribute__((always_inline, const)) INLINE static void add_to_hashmap(hashmap_t *map, hashmap_key_t key, double m_add) {

  int created = 0;
  hashmap_value_t *value = hashmap_get_new(map, key, &created);
  if(created) {
    /* Key was not present, so this is a new element */
    value->value_dbl = m_add;
  } else {
    /* Key was present, so add m_add to previous value */
    value->value_dbl += m_add;
  }
}


/**
 * @brief Accumulate local contributions to the density field
 *
 * Creates a hashmap with the contributions to the density
 * mesh from local particles.
 *
 * TODO: parallelize. E.g. one hashmap per thread and combine
 * when done.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the cell size
 * @param s The #space containing the particles.
 * @param map The hashmap in which to store the results
 *
 */
void accumulate_local_gparts_to_hashmap(const int N, const double fac, 
                                        const struct space* s, hashmap_t *map) {

  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int* local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  /* Loop over local cells */
  for(int icell=0; icell < nr_local_cells; icell+=1) {
    /* Get a pointer to this cell */
    struct cell *cell = &(s->cells_top[local_cells[icell]]);
    /* Loop over particles in this cell */
    for(int ipart=0; ipart < cell->grav.count; ipart+=1) {

      struct gpart *gp = &(cell->grav.parts[ipart]);

      /* Box wrap the multipole's position */
      const double pos_x = box_wrap(gp->x[0], 0., dim[0]);
      const double pos_y = box_wrap(gp->x[1], 0., dim[1]);
      const double pos_z = box_wrap(gp->x[2], 0., dim[2]);

      /* Workout the CIC coefficients */
      int i = (int)(fac * pos_x);
      if (i >= N) i = N - 1;
      const double dx = fac * pos_x - i;
      const double tx = 1. - dx;

      int j = (int)(fac * pos_y);
      if (j >= N) j = N - 1;
      const double dy = fac * pos_y - j;
      const double ty = 1. - dy;

      int k = (int)(fac * pos_z);
      if (k >= N) k = N - 1;
      const double dz = fac * pos_z - k;
      const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
      if (i < 0 || i >= N) error("Invalid gpart position in x");
      if (j < 0 || j >= N) error("Invalid gpart position in y");
      if (k < 0 || k >= N) error("Invalid gpart position in z");
#endif

      /* Accumulate contributions to the hashmap */
      const double mass = gp->mass;
      add_to_hashmap(map, (hashmap_key_t) row_major_id_periodic(i + 0, j + 0, k + 0, N), mass * tx * ty * tz);
      add_to_hashmap(map, (hashmap_key_t) row_major_id_periodic(i + 0, j + 0, k + 1, N), mass * tx * ty * dz);
      add_to_hashmap(map, (hashmap_key_t) row_major_id_periodic(i + 0, j + 1, k + 0, N), mass * tx * dy * tz);
      add_to_hashmap(map, (hashmap_key_t) row_major_id_periodic(i + 0, j + 1, k + 1, N), mass * tx * dy * dz);
      add_to_hashmap(map, (hashmap_key_t) row_major_id_periodic(i + 1, j + 0, k + 0, N), mass * dx * ty * tz);
      add_to_hashmap(map, (hashmap_key_t) row_major_id_periodic(i + 1, j + 0, k + 1, N), mass * dx * ty * dz);
      add_to_hashmap(map, (hashmap_key_t) row_major_id_periodic(i + 1, j + 1, k + 0, N), mass * dx * dy * tz);
      add_to_hashmap(map, (hashmap_key_t) row_major_id_periodic(i + 1, j + 1, k + 1, N), mass * dx * dy * dz);

    } /* Next particle */
  } /* Next cell */
  
  return;
}
