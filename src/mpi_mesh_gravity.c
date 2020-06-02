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

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "mpi_mesh_gravity.h"

/* Local includes. */
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "part.h"
#include "space.h"

#if defined(WITH_MPI) && defined(HAVE_MPI_FFTW)

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
__attribute__((always_inline)) INLINE static void add_to_hashmap(
    hashmap_t *map, hashmap_key_t key, double m_add) {

  int created = 0;
  hashmap_value_t *value = hashmap_get_new(map, key, &created);
  if (created) {
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
 * mesh from local particles. Here we require that hash_key_t
 * can store values up to at least N*N*N.
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
void mpi_mesh_accumulate_gparts_to_hashmap(const int N, const double fac,
                                           const struct space *s, hashmap_t *map) {

  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  /* Loop over local cells */
  for (int icell = 0; icell < nr_local_cells; icell += 1) {
    /* Get a pointer to this cell */
    struct cell *cell = &(s->cells_top[local_cells[icell]]);
    /* Loop over particles in this cell */
    for (int ipart = 0; ipart < cell->grav.count; ipart += 1) {

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
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t_padded(i + 0, j + 0, k + 0, N),
          mass * tx * ty * tz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t_padded(i + 0, j + 0, k + 1, N),
          mass * tx * ty * dz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t_padded(i + 0, j + 1, k + 0, N),
          mass * tx * dy * tz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t_padded(i + 0, j + 1, k + 1, N),
          mass * tx * dy * dz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t_padded(i + 1, j + 0, k + 0, N),
          mass * dx * ty * tz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t_padded(i + 1, j + 0, k + 1, N),
          mass * dx * ty * dz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t_padded(i + 1, j + 1, k + 0, N),
          mass * dx * dy * tz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t_padded(i + 1, j + 1, k + 1, N),
          mass * dx * dy * dz);

    } /* Next particle */
  }   /* Next cell */

  return;
}

/**
 * @brief Store contributions to the mesh as (index, mass) pairs
 */
struct mesh_key_value {
  hashmap_key_t key;
  double value;
};

/**
 * @brief Comparison function to sort mesh_key_value by key
 *
 * @param a The first #mesh_key_value object.
 * @param b The second #mesh_key_value object.
 * @return 1 if a's key field is greater than b's key field,
 * -1 if a's key field is less than b's key field, and zero otherwise
 *
 */
int cmp_func_mesh_key_value(const void *a, const void *b) {
  struct mesh_key_value *a_key_value = (struct mesh_key_value *)a;
  struct mesh_key_value *b_key_value = (struct mesh_key_value *)b;
  if (a_key_value->key > b_key_value->key)
    return 1;
  else if (a_key_value->key < b_key_value->key)
    return -1;
  else
    return 0;
}

/**
 * @brief Data needed to iterate over a hashmap copying key-value pairs
 */
struct hashmap_mapper_data {
  size_t offset;
  struct mesh_key_value *buf;
};

/**
 * @brief Mapper function to copy elements from the hashmap
 *
 * @param key The key associated with this hashmap entry
 * @param value The value associated with this hashmap entry
 * @param data Contains pointer to output buffer and offset to
 * next element to write
 *
 */
void hashmap_copy_elements_mapper(hashmap_key_t key, hashmap_value_t *value,
                                  void *data) {

  struct hashmap_mapper_data *mapper_data = (struct hashmap_mapper_data *)data;
  struct mesh_key_value *element = &(mapper_data->buf[mapper_data->offset]);
  element->key = key;
  element->value = value->value_dbl;
  mapper_data->offset += 1;
}

/**
 * @brief Copy keys and values from a hashmap into an array of struct
 * mesh_key_value, with the elements sorted by key
 *
 * @param map The hashmap to copy into the array
 * @param array The output array
 * @param n The size of the array
 * @return The number of elements copied into the array
 *
 */
size_t hashmap_to_sorted_array(hashmap_t *map, struct mesh_key_value *array, size_t n) {

  /* Find how many elements are in the hashmap */
  size_t num_in_map = hashmap_size(map);
  if(num_in_map > n) {
    error("Array is to small to contain hash map elements!");
  }

  /* Copy the key-value pairs to the new array */
  struct hashmap_mapper_data mapper_data = {(size_t)0, array};
  hashmap_iterate(map, hashmap_copy_elements_mapper, &mapper_data);

  /* And sort the pairs by key */
  qsort(array, num_in_map, sizeof(struct mesh_key_value),
        cmp_func_mesh_key_value);

  return num_in_map;
}


/**
 * @brief Given an array of structs of size element_size, send 
 * nr_send[i] elements to each node i. Allocates the receive
 * buffer recvbuf to the appropriate size and returns its size
 * in nr_recv_tot.
 *
 * TODO: can/should we replace this with a call to engine_do_redistribute()?
 *
 * @param nr_send Number of elements to send to each node
 * @param nr_recv Number of elements to receive from each node
 * @param sendbuf The elements to send
 * @param recvbuf The output buffer
 *
 */
void exchange_structs(size_t *nr_send, char *sendbuf,
                      size_t *nr_recv, char *recvbuf,
                      size_t element_size) {

  /* Determine rank, number of ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Compute send offsets */
  size_t *send_offset = malloc(nr_nodes * sizeof(size_t));
  send_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    send_offset[i] = send_offset[i - 1] + nr_send[i - 1];
  }

  /* Compute receive offsets */
  size_t *recv_offset = malloc(nr_nodes * sizeof(size_t));
  recv_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    recv_offset[i] = recv_offset[i - 1] + nr_recv[i - 1];
  }

  /* Allocate request objects (one send and receive per node) */
  MPI_Request *request = malloc(2 * sizeof(MPI_Request) * nr_nodes);

  /* Make type to communicate mesh_key_value struct */
  MPI_Datatype mesh_key_value_mpi_type;
  if (MPI_Type_contiguous(element_size, MPI_BYTE, &mesh_key_value_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&mesh_key_value_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for mesh_key_value struct.");
  }

  /*
   * Post the send operations. This is an alltoallv really but
   * we want to avoid the limits imposed by int counts and offsets
   * in MPI_Alltoallv.
   */
  for (int i = 0; i < nr_nodes; i += 1) {
    if (nr_send[i] > 0) {

      /* TODO: handle very large messages */
      if(nr_send[i] > INT_MAX) error("exchange_structs() fails if nr_send > INT_MAX!");

      MPI_Isend(&(sendbuf[send_offset[i]*element_size]), (int)nr_send[i],
                mesh_key_value_mpi_type, i, 0, MPI_COMM_WORLD, &(request[i]));
    } else {
      request[i] = MPI_REQUEST_NULL;
    }
  }

  /* Post the receives */
  for (int i = 0; i < nr_nodes; i += 1) {
    if (nr_recv[i] > 0) {

      /* TODO: handle very large messages */
      if(nr_recv[i] > INT_MAX) error("exchange_structs() fails if nr_recv > INT_MAX!");

      MPI_Irecv(&(recvbuf[recv_offset[i]*element_size]), (int)nr_recv[i],
                mesh_key_value_mpi_type, i, 0, MPI_COMM_WORLD,
                &(request[i + nr_nodes]));
    } else {
      request[i + nr_nodes] = MPI_REQUEST_NULL;
    }
  }

  /* Wait for everything to complete */
  MPI_Waitall(2 * nr_nodes, request, MPI_STATUSES_IGNORE);

  /* Done with the MPI type */
  MPI_Type_free(&mesh_key_value_mpi_type);

  /* Tidy up */
  free(recv_offset);
  free(send_offset);
  free(request);
}

/**
 * @brief Convert hashmaps to a slab-distributed 3D mesh
 *
 * For FFTW each rank needs to hold a slice of the full mesh.
 * This routine does the necessary communication to convert
 * the per-rank hashmaps into a slab-distributed mesh.
 *
 * @param e Pointer to the engine struct
 * @param N The size of the mesh
 * @param local_n0 The thickness of the slice to store on this rank
 * @param map The hashmap with the local part of the mesh
 * @param mesh Pointer to the output data buffer
 *
 */
void mpi_mesh_hashmaps_to_slices(const int N, const int local_n0, hashmap_t *map,
                                 double *mesh) {

  /* Determine rank, number of ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Make an array with the (key, value) pairs from the hashmap.
   * The elements are sorted by key, which means they're sorted 
   * by x coordinate, then y coordinate, then z coordinate.
   * We're going to distribute them between ranks according to their
   * x coordinate, so this puts them in order of destination rank.
   */
  size_t map_size = hashmap_size(map);
  struct mesh_key_value *mesh_sendbuf;
  if (swift_memalign("mesh_sendbuf", (void **) &mesh_sendbuf, 32,
                     map_size * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate array for mesh send buffer!");
  size_t nr_send_tot = hashmap_to_sorted_array(map, mesh_sendbuf, map_size);

  /* Get width of the slice on each rank */
  int *slice_width = malloc(sizeof(int) * nr_nodes);
  MPI_Allgather(&local_n0, 1, MPI_INT, slice_width, 1, MPI_INT, MPI_COMM_WORLD);

  /* Determine offset to the slice on each rank */
  int *slice_offset = malloc(sizeof(int) * nr_nodes);
  slice_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    slice_offset[i] = slice_offset[i - 1] + slice_width[i - 1];
  }

  /* Compute how many elements are to be sent to each rank */
  size_t *nr_send = malloc(nr_nodes * sizeof(size_t));
  for (int i = 0; i < nr_nodes; i += 1) {
    nr_send[i] = 0;
  }
  int dest_node = 0;
  for (size_t i = 0; i < nr_send_tot; i += 1) {
    /* Get the x coordinate of this mesh cell in the global mesh */
    int mesh_x = get_xcoord_from_padded_row_major_id((size_t) mesh_sendbuf[i].key, N);
    /* Advance to the destination node that is to contain this x coordinate */
    while ((mesh_x >= slice_offset[dest_node] + slice_width[dest_node]) ||
           (slice_width[dest_node] == 0)) {
      dest_node += 1;
    }
    nr_send[dest_node] += 1;
  }

  /* Determine how many requests we'll receive from each MPI rank */
  size_t *nr_recv = malloc(sizeof(size_t)*nr_nodes);
  MPI_Alltoall(nr_send, sizeof(size_t), MPI_BYTE, nr_recv, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);
  size_t nr_recv_tot = 0;
  for(int i=0; i<nr_nodes; i+=1) {
    nr_recv_tot += nr_recv[i];
  }

  /* Allocate the receive buffer */
  struct mesh_key_value *mesh_recvbuf;
  if (swift_memalign("mesh_recvbuf", (void **) &mesh_recvbuf, 32,
                     nr_recv_tot * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate receive buffer for constructing MPI FFT mesh");

  /* Carry out the communication */
  exchange_structs(nr_send, (char *) mesh_sendbuf,
                   nr_recv, (char *) mesh_recvbuf, 
                   sizeof(struct mesh_key_value));

  /* Copy received data to the output buffer */
  for (size_t i = 0; i < nr_recv_tot; i += 1) {
#ifdef SWIFT_DEBUG_CHECKS
    const int xcoord = get_xcoord_from_padded_row_major_id(mesh_recvbuf[i].key, N);
    if(xcoord < slice_offset[nodeID])
      error("Received mesh cell is not in the local slice (xcoord too small)");
    if(xcoord >= slice_offset[nodeID] + slice_width[nodeID])
      error("Received mesh cell is not in the local slice (xcoord too large)");
#endif
    mesh[get_index_in_local_slice((size_t) mesh_recvbuf[i].key, N, slice_offset[nodeID])] += mesh_recvbuf[i].value;
  }

  /* Tidy up */
  free(slice_width);
  free(slice_offset);
  free(nr_send);
  free(nr_recv);
  swift_free("mesh_recvbuf", mesh_recvbuf);
  swift_free("mesh_sendbuf", mesh_sendbuf);
}


/**
 * @brief Retrieve the potential in the mesh cells we need to
 * compute the force on particles on this MPI rank. Result is
 * returned in the supplied hashmap, which should be initially
 * empty.
 *
 * We need all cells containing points two mesh cell widths
 * away from each particle along each axis to compute the
 * potential gradient. We also need to allow for the movement
 * of the particles between calculation of the mesh at rebuild time
 * and evaluation of the forces on each timestep. Particles can move
 * up to half of a top level cell size between updates of the mesh.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the FFT mesh cell size
 * @param s The #space containing the particles.
 * @param local_0_start Offset to the first mesh x coordinate on this rank
 * @param local_n0 Width of the mesh slab on this rank
 * @param potential_slice Array with the potential on the local slice of the mesh
 * @param potential_map A hashmap in which to store the potential data
 *
 */
void mpi_mesh_fetch_potential(const int N, const double fac,
                              const struct space *s,
                              int local_0_start, int local_n0,
                              double *potential_slice,
                              hashmap_t *potential_map) {

  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  /* Determine rank, number of MPI ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Create a hashmap to store the indexes of the cells we need */
  hashmap_t map;
  hashmap_init(&map);

  /* Loop over our local top level cells */
  for (int icell = 0; icell < nr_local_cells; icell += 1) {
    struct cell *cell = &(s->cells_top[local_cells[icell]]);
    if(cell->grav.count > 0) {

      /* Determine range of mesh cells we need for particles in this top level cell */
      int ixmin[3];
      int ixmax[3];
      for(int idim=0;idim<3;idim+=1) {
        const double xmin = cell->loc[idim] - 0.5*cell->width[idim] - 2.0/fac;
        const double xmax = cell->loc[idim] + 1.5*cell->width[idim] + 2.0/fac;
        ixmin[idim] = (int) floor(xmin*fac);
        ixmax[idim] = (int) floor(xmax*fac);
      }

      /* Add the required cells to the map */
      for(int i=ixmin[0]; i<=ixmax[0]; i+=1) {
        for(int j=ixmin[1]; j<=ixmax[1]; j+=1) {
          for(int k=ixmin[2]; k<=ixmax[2]; k+=1) {
            const size_t index = row_major_id_periodic_size_t_padded(i, j, k, N);
            /* We don't have a value associated with the entry yet */
	    hashmap_value_t *value = hashmap_get(&map, (hashmap_key_t) index);
	    value->value_dbl = 0.0;
          }
        }
      }
    }
  }

  /* Make an array with the cell IDs we need to request from other ranks */
  size_t nr_send_tot = hashmap_size(&map);
  struct mesh_key_value *send_cells;
  if (swift_memalign("send_cells", (void **) &send_cells, 32,
                     nr_send_tot * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate array for cells to request!");
  nr_send_tot = hashmap_to_sorted_array(&map, send_cells, nr_send_tot);
  hashmap_free(&map);

  /* Get width of the mesh slice on each rank */
  int *slice_width = malloc(sizeof(int) * nr_nodes);
  MPI_Allgather(&local_n0, 1, MPI_INT, slice_width, 1, MPI_INT, MPI_COMM_WORLD);

  /* Determine first mesh x coordinate stored on each rank */
  int *slice_offset = malloc(sizeof(int) * nr_nodes);
  slice_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    slice_offset[i] = slice_offset[i - 1] + slice_width[i - 1];
  }

  /* Count how many mesh cells we need to request from each MPI rank */
  int dest_rank = 0;
  size_t *nr_send = malloc(sizeof(size_t)*nr_nodes);
  for(int i=0; i<nr_nodes; i+=1) {
    nr_send[i] = 0;
  }
  for(size_t i=0; i<nr_send_tot; i+=1) {
    while(get_xcoord_from_padded_row_major_id(send_cells[i].key, N) >= (slice_offset[dest_rank]+slice_width[dest_rank]) || slice_width[dest_rank] == 0) {
      dest_rank += 1;
    }
#ifdef SWIFT_DEBUG_CHECKS
    if(dest_rank >= nr_nodes || dest_rank < 0)error("Destination rank out of range");
#endif
    nr_send[dest_rank] += 1;
  }

  /* Determine how many requests we'll receive from each MPI rank */
  size_t *nr_recv = malloc(sizeof(size_t)*nr_nodes);
  MPI_Alltoall(nr_send, sizeof(size_t), MPI_BYTE, nr_recv, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);
  size_t nr_recv_tot = 0;
  for(int i=0; i<nr_nodes; i+=1) {
    nr_recv_tot += nr_recv[i];
  }

  /* Allocate buffer to receive requests */
  struct mesh_key_value *recv_cells;
  if (swift_memalign("recv_cells", (void **) &recv_cells, 32,
                     nr_recv_tot * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate array for mesh receive buffer!");

  /* Send requests for cells to other ranks */
  exchange_structs(nr_send, (char *) send_cells,
                   nr_recv, (char *) recv_cells,
                   sizeof(struct mesh_key_value));
  
  /* Look up potential in the requested cells */
  for(size_t i=0; i<nr_recv_tot; i+=1) {
#ifdef SWIFT_DEBUG_CHECKS
    const size_t cells_in_slab  = ((size_t) N)*(2*(N/2+1));
    const size_t first_local_id = local_0_start*cells_in_slab;
    const size_t num_local_ids  = local_n0*cells_in_slab;
    if(recv_cells[i].key < first_local_id || recv_cells[i].key >= first_local_id+num_local_ids) {
      error("Requested potential mesh cell ID is out of range");
    }
#endif
    size_t local_id = get_index_in_local_slice(recv_cells[i].key, N, local_0_start);
#ifdef SWIFT_DEBUG_CHECKS
    const size_t Ns = N;
    if(local_id >= Ns*(2*(Ns/2+1))*local_n0)error("Local potential mesh cell ID is out of range");
#endif
    recv_cells[i].value = potential_slice[local_id];
  }

  /* Return the results */
  exchange_structs(nr_recv, (char *) recv_cells,
                   nr_send, (char *) send_cells,
                   sizeof(struct mesh_key_value));

  /* Store the results in the hashmap */
  for(size_t i=0; i<nr_send_tot; i+=1) {
    int created = 0;
    hashmap_value_t *value = hashmap_get_new(potential_map, send_cells[i].key,
					     &created);
    value->value_dbl = send_cells[i].value;
#ifdef SWIFT_DEBUG_CHECKS
    if(!created)error("Received duplicate potential hash map value");
    const size_t Ns = N;
    if(send_cells[i].key >= Ns*Ns*(2*(Ns/2+1)))error("Received potential mesh cell ID out of range");
#endif
  }

  /* Tidy up */
  swift_free("send_cells", send_cells);
  swift_free("recv_cells", recv_cells);
  free(slice_width);
  free(slice_offset);
  free(nr_send);
  free(nr_recv);
}

#endif
