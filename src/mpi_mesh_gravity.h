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

#ifndef SWIFT_MPI_MESH_GRAVITY_H
#define SWIFT_MPI_MESH_GRAVITY_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "hashmap.h"
#include "space.h"

/**
 * @brief Returns 1D index of a 3D NxNxN array using row-major style.
 *
 * Wraps around in the corresponding dimension if any of the 3 indices is >= N
 * or < 0. Note that indexes are of type size_t in this version!
 *
 * @param i Index along x.
 * @param j Index along y.
 * @param k Index along z.
 * @param N Size of the array along one axis.
 */
__attribute__((always_inline, const)) INLINE static size_t
row_major_id_periodic_size_t(const size_t i, const size_t j, const size_t k,
                             const size_t N) {

  return (((i + N) % N) * N * N + ((j + N) % N) * N + ((k + N) % N));
}

/**
 * @brief Accumulate local contributions to the density field
 *
 * Creates a hashmap with the contributions to the density
 * mesh from local particles.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the cell size
 * @param s The #space containing the particles.
 * @param map The hashmap in which to store the results
 *
 */
void accumulate_local_gparts_to_hashmap(const int N, const double fac,
                                        const struct space *s, hashmap_t *map);

/**
 * @brief Convert hashmaps to a slab-distributed 3D mesh
 *
 * For FFTW each rank needs to hold a slice of the full mesh.
 * This routine does the necessary communication to convert
 * the per-rank hashmaps into a slab-distributed mesh.
 *
 * @param e Pointer to the engine struct
 * @param N The size of the mesh
 * @param Nslice The thickness of the slice to store on this rank
 * @param map The hashmap with the local part of the mesh
 * @param mesh Pointer to the output data buffer
 *
 */
void hashmaps_to_slices(const int N, const int Nslice, hashmap_t *map,
                        double *mesh);

#endif
