//
// Created by yuyttenh on 02/07/2021.
//

#ifndef SWIFTSIM_VORONOI_FUNCTIONS_3D_H
#define SWIFTSIM_VORONOI_FUNCTIONS_3D_H

/* Forward declarations */
inline static int voronoi_new_face(struct voronoi *v, const struct delaunay *d,
                                   int left_part_idx_in_d,
                                   int right_part_idx_in_d, double *vertices,
                                   int n_vertices);
inline static void voronoi_check_grid(struct voronoi *restrict v);
inline static void voronoi_destroy(struct voronoi *restrict v);

/**
 * @brief Initializes a voronoi pair (i.e. a face). The vertices of the face are
 * only stored if VORONOI_STORE_FACES is defined.
 *
 * @param pair Voronoi pair to initialize
 * @param c Pointer to the SWIFT cell in which the right particle lives
 * (NULL if this is the same cell as the left particle)
 * @param left_part_pointer Pointer to the left particle of this pair
 * @param right_part_pointer Pointer to the right particle of this pair
 * @param vertices Vertices making up the face
 * @param n_vertices Number of vertices in the vertices array.
 */
inline static void voronoi_pair_init(struct voronoi_pair *pair,
                                     struct cell *restrict c,
                                     struct part *restrict left_part_pointer,
                                     struct part *restrict right_part_pointer,
                                     int left_part_idx_in_d,
                                     int right_part_idx_in_d, double *vertices,
                                     int n_vertices) {
  pair->right_cell = c;
  pair->left = left_part_pointer;
  pair->left_idx = left_part_idx_in_d;
  pair->right = right_part_pointer;
  pair->right_idx = right_part_idx_in_d;

  pair->surface_area =
      geometry3d_compute_centroid_area(vertices, n_vertices, pair->midpoint);

#ifdef VORONOI_STORE_FACES
  pair->vertices = (double *)malloc(3 * n_vertices * sizeof(double));
  pair->n_vertices = n_vertices;
#ifdef SWIFT_DEBUG_CHECKS
  assert(pair->n_vertices > 0);
#endif
  for (int i = 0; i < n_vertices; i++) {
    pair->vertices[3 * i] = vertices[3 * i];
    pair->vertices[3 * i + 1] = vertices[3 * i + 1];
    pair->vertices[3 * i + 2] = vertices[3 * i + 2];
  }
#endif
}

inline static void voronoi_pair_destroy(struct voronoi_pair *pair) {
#ifdef VORONOI_STORE_FACES
  free(pair->vertices);
  pair->vertices = NULL;
#endif
}

inline static void voronoi_malloc(struct voronoi *restrict v,
                                  int number_of_cells, double dmin,
                                  struct cell *restrict swift_cell) {
  v->number_of_cells = number_of_cells;
  /* allocate memory for the voronoi cells */
  v->cells = (struct voronoi_cell_new *)swift_malloc(
      "c.h.v.cells", v->number_of_cells * sizeof(struct voronoi_cell_new));
  v->cells_size = v->number_of_cells;

  /* Allocate memory for the voronoi pairs (faces). */
  for (int i = 0; i < 28; ++i) {
    v->pairs[i] = (struct voronoi_pair *)swift_malloc(
        "c.h.v.pairs", 10 * sizeof(struct voronoi_pair));
    v->pair_index[i] = 0;
    v->pair_size[i] = 10;
  }

  /* Allocate memory for the cell_pair connections */
  int2_lifo_queue_init(&v->cell_pair_connections, 6 * number_of_cells);

  v->swift_cell = swift_cell;
  v->min_surface_area = min_rel_voronoi_face_size * dmin;
  v->active = 1;
}

inline static void voronoi_init(struct voronoi *restrict v, int number_of_cells,
                                double dmin) {
  voronoi_assert(v->active);

  v->number_of_cells = number_of_cells;
  if (v->cells_size < v->number_of_cells) {
    /* allocate memory for the voronoi cells */
    v->cells = (struct voronoi_cell_new *)swift_realloc(
        "c.h.v.cells", v->cells,
        v->number_of_cells * sizeof(struct voronoi_cell_new));
    v->cells_size = v->number_of_cells;
  }

  /* reset indices for the voronoi pairs (faces). */
  for (int i = 0; i < 28; ++i) {
    /* Free any allocations made for the old pairs */
    for (int j = 0; j < v->pair_index[i]; j++) {
      voronoi_pair_destroy(&v->pairs[i][j]);
    }
    v->pair_index[i] = 0;
  }

  /* Reset the cell_pair connections */
  int2_lifo_queue_reset(&v->cell_pair_connections);

  v->min_surface_area = min_rel_voronoi_face_size * dmin;
}

/**
 * @brief Build the Voronoi grid based on the given Delaunay tessellation.
 *
 * This function allocates the memory for the Voronoi grid arrays and creates
 * the grid in linear time by
 *  1. Computing the grid vertices as the midpoints of the circumcircles of the
 *     Delaunay tetrahedra.
 *  2. Looping over all vertices and for each generator looping over all
 *     tetrahedra that link to that vertex. This is done by looping around all
 *     the Delaunay edges connected to that generator. While looping around an
 *     edge, for each tetrahedron, we add the edges of that tetrahedron which
 *     are connected to the current generator to a queue of next edges to loop
 *     around (if we did not already do so).
 *
 * During the second step, the geometrical properties (cell centroid, volume
 * and face midpoint, area) are computed as well.
 *
 * @param v Voronoi grid.
 * @param d Delaunay tessellation (read-only).
 */
inline static void voronoi_build(struct voronoi *v, struct delaunay *d) {

  /* the number of cells equals the number of non-ghost and non-dummy
     vertex_indices in the Delaunay tessellation */
  voronoi_assert(d->vertex_end > 0);
  voronoi_assert(v->number_of_cells == d->vertex_end - d->vertex_start);

  /* Allocate memory to store voronoi vertices (will be freed at end of this
   * function!) */
  double *voronoi_vertices =
      (double *)malloc(3 * (d->tetrahedra_index - 4) * sizeof(double));

  /* loop over the tetrahedra in the Delaunay tessellation and compute the
     midpoints of their circumspheres. These happen to be the vertices of
     the Voronoi grid (because they are the points of equal distance to 3
     generators, while the Voronoi edges are the lines of equal distance to 2
     generators) */
  for (int i = 0; i < d->tetrahedra_index - 4; i++) {
    struct tetrahedron *t = &d->tetrahedra[i + 4];
    /* Get the indices of the vertices of the tetrahedron */
    int v0 = t->vertices[0];
    int v1 = t->vertices[1];
    int v2 = t->vertices[2];
    int v3 = t->vertices[3];

    /* if the tetrahedron is inactive or not linked to a non-ghost, non-dummy
     * vertex, it is not a grid vertex and we can skip it. */
    if (!t->active || (v0 >= d->vertex_end && v1 >= d->vertex_end &&
                       v2 >= d->vertex_end && v3 >= d->vertex_end)) {
      voronoi_vertices[3 * i] = NAN;
      voronoi_vertices[3 * i + 1] = NAN;
      voronoi_vertices[3 * i + 2] = NAN;
      continue;
    }
    /* Check that the vertices are valid */
    voronoi_assert(v0 >= 0 && v1 >= 0 && v2 >= 0 && v3 >= 0);

    /* Extract coordinates from the Delaunay vertices (generators)
     * FUTURE NOTE: In swift we should read this from the particles themselves!
     * */
    double *v0d, *v1d, *v2d, *v3d;
    unsigned long *v0ul, *v1ul, *v2ul, *v3ul;
    if (v0 < d->vertex_end || v0 >= d->ngb_offset) {
      v0d = &d->rescaled_vertices[3 * v0];
      v0ul = &d->integer_vertices[3 * v0];
    } else {
      /* This could mean that a neighbouring cell of this grids cell is empty!
       * Or that we did not add all the necessary ghost vertex_indices to the
       * delaunay tesselation. */
      voronoi_error(
          "Vertex is part of tetrahedron with Dummy vertex! This could mean "
          "that one of the neighbouring cells is empty.");
    }
    if (v1 < d->vertex_end || v1 >= d->ngb_offset) {
      v1d = &d->rescaled_vertices[3 * v1];
      v1ul = &d->integer_vertices[3 * v1];
    } else {
      voronoi_error(
          "Vertex is part of tetrahedron with Dummy vertex! This could mean "
          "that one of the neighbouring cells is empty.");
    }
    if (v2 < d->vertex_end || v2 >= d->ngb_offset) {
      v2d = &d->rescaled_vertices[3 * v2];
      v2ul = &d->integer_vertices[3 * v2];
    } else {
      voronoi_error(
          "Vertex is part of tetrahedron with Dummy vertex! This could mean "
          "that "
          "one of the neighbouring cells is empty.");
    }
    if (v3 < d->vertex_end || v3 >= d->ngb_offset) {
      v3d = &d->rescaled_vertices[3 * v3];
      v3ul = &d->integer_vertices[3 * v3];
    } else {
      voronoi_error(
          "Vertex is part of tetrahedron with Dummy vertex! This could mean "
          "that one of the neighbouring cells is empty.");
    }

    geometry3d_compute_circumcenter_adaptive(
        &d->geometry, v0d, v1d, v2d, v3d, v0ul, v1ul, v2ul, v3ul,
        &voronoi_vertices[3 * i], d->side, d->anchor);

#ifdef VORONOI_CHECKS
    const double cx = voronoi_vertices[3 * i];
    const double cy = voronoi_vertices[3 * i + 1];
    const double cz = voronoi_vertices[3 * i + 2];

    const double r0 =
        sqrt((cx - d->vertices[3 * v0]) * (cx - d->vertices[3 * v0]) +
             (cy - d->vertices[3 * v0 + 1]) * (cy - d->vertices[3 * v0 + 1]) +
             (cz - d->vertices[3 * v0 + 2]) * (cz - d->vertices[3 * v0 + 2]));
    const double r1 =
        sqrt((cx - d->vertices[3 * v1]) * (cx - d->vertices[3 * v1]) +
             (cy - d->vertices[3 * v1 + 1]) * (cy - d->vertices[3 * v1 + 1]) +
             (cz - d->vertices[3 * v1 + 2]) * (cz - d->vertices[3 * v1 + 2]));
    const double r2 =
        sqrt((cx - d->vertices[3 * v2]) * (cx - d->vertices[3 * v2]) +
             (cy - d->vertices[3 * v2 + 1]) * (cy - d->vertices[3 * v2 + 1]) +
             (cz - d->vertices[3 * v2 + 2]) * (cz - d->vertices[3 * v2 + 2]));
    const double r3 =
        sqrt((cx - d->vertices[3 * v3]) * (cx - d->vertices[3 * v3]) +
             (cy - d->vertices[3 * v3 + 1]) * (cy - d->vertices[3 * v3 + 1]) +
             (cz - d->vertices[3 * v3 + 2]) * (cz - d->vertices[3 * v3 + 2]));
    voronoi_assert(double_cmp(r0, r1, 1e5) && double_cmp(r0, r2, 1e5) &&
                   double_cmp(r0, r3, 1e5));
#endif
  } /* loop over the Delaunay tetrahedra and compute the circumcenters */

  /* Allocate memory for the neighbour flags and initialize them to 0 (will be
   * freed at the end of this function!) */
  int *neighbour_flags = (int *)malloc(d->vertex_index * sizeof(int));
  for (int i = 0; i < d->vertex_index; i++) {
    neighbour_flags[i] = 0;
  }

  /* Allocate a tetrahedron_vertex_queue (will be freed at the end of this
   * function!) */
  struct int3_fifo_queue neighbour_info_q;
  int3_fifo_queue_init(&neighbour_info_q, 10);

  /* The size of the array used to temporarily store the vertices of the voronoi
   * faces in */
  int face_vertices_size = 10;
  /* Temporary array to store face vertices in (will be freed at the end of this
   * function!) */
  double *face_vertices =
      (double *)malloc(3 * face_vertices_size * sizeof(double));

  /* loop over all cell generators, and hence over all non-ghost, non-dummy
     Delaunay vertex_indices */
  for (int gen_idx_in_d = d->vertex_start; gen_idx_in_d < d->vertex_end;
       gen_idx_in_d++) {
    /* First reset the tetrahedron_vertex_queue */
    int3_fifo_queue_reset(&neighbour_info_q);
    /* Set the flag of the central generator so that we never pick it as
     * possible neighbour */
    neighbour_flags[gen_idx_in_d] = 1;

    /* Create a new voronoi cell for this generator */
    struct voronoi_cell_new *this_cell =
        &v->cells[gen_idx_in_d - d->vertex_start];
    this_cell->volume = 0.;
    this_cell->centroid[0] = 0.;
    this_cell->centroid[1] = 0.;
    this_cell->centroid[2] = 0.;
    int nface = 0;
    this_cell->pair_connections_offset = v->cell_pair_connections.index;

    /* get the generator position, we use it during centroid/volume
       calculations */
    voronoi_assert(gen_idx_in_d < d->vertex_end);
    double ax = d->vertices[3 * gen_idx_in_d];
    double ay = d->vertices[3 * gen_idx_in_d + 1];
    double az = d->vertices[3 * gen_idx_in_d + 2];

#ifdef VORONOI_STORE_GENERATORS
    this_cell->generator = d->part_pointers[gen_idx_in_d];
#endif

    /* Get a tetrahedron containing the central generator */
    int t_idx = d->vertex_tetrahedron_links[gen_idx_in_d];
    int gen_idx_in_t = d->vertex_tetrahedron_index[gen_idx_in_d];

    /* Pick another vertex (generator) from this tetrahedron and add it to the
     * queue */
    int other_v_idx_in_t = (gen_idx_in_t + 1) % 4;
    struct tetrahedron *t = &d->tetrahedra[t_idx];
    int other_v_idx_in_d = t->vertices[other_v_idx_in_t];
    int3 info = {._0 = t_idx, ._1 = other_v_idx_in_d, ._2 = other_v_idx_in_t};
    int3_fifo_queue_push(&neighbour_info_q, info);
    /* update flag of the other vertex */
    neighbour_flags[other_v_idx_in_d] = 1;

    while (!int3_fifo_queue_is_empty(&neighbour_info_q)) {
      /* Pop the next axis vertex and corresponding tetrahedron from the queue
       */
      info = int3_fifo_queue_pop(&neighbour_info_q);
      int first_t_idx = info._0;
      int axis_idx_in_d = info._1;
      int axis_idx_in_t = info._2;
      voronoi_assert(axis_idx_in_d >= 0 && (axis_idx_in_d < d->vertex_end ||
                                            axis_idx_in_d >= d->ngb_offset));
      struct tetrahedron *first_t = &d->tetrahedra[first_t_idx];

      /* Get a non axis vertex from first_t */
      int non_axis_idx_in_first_t = (axis_idx_in_t + 1) % 4;
      if (first_t->vertices[non_axis_idx_in_first_t] == gen_idx_in_d) {
        non_axis_idx_in_first_t = (non_axis_idx_in_first_t + 1) % 4;
      }
      int non_axis_idx_in_d = first_t->vertices[non_axis_idx_in_first_t];

      if (!neighbour_flags[non_axis_idx_in_d]) {
        /* Add this vertex and tetrahedron to the queue and update its flag */
        int3 new_info = {._0 = first_t_idx,
                         ._1 = non_axis_idx_in_d,
                         ._2 = non_axis_idx_in_first_t};
        int3_fifo_queue_push(&neighbour_info_q, new_info);
        neighbour_flags[non_axis_idx_in_d] |= 1;
      }

      /* Get a neighbouring tetrahedron of first_t sharing the axis */
      int cur_t_idx = first_t->neighbours[non_axis_idx_in_first_t];
      struct tetrahedron *cur_t = &d->tetrahedra[cur_t_idx];
      int prev_t_idx_in_cur_t =
          first_t->index_in_neighbour[non_axis_idx_in_first_t];

      /* Get a neighbouring tetrahedron of cur_t that is not first_t, sharing
       * the same axis */
      int next_t_idx_in_cur_t = (prev_t_idx_in_cur_t + 1) % 4;
      while (cur_t->vertices[next_t_idx_in_cur_t] == gen_idx_in_d ||
             cur_t->vertices[next_t_idx_in_cur_t] == axis_idx_in_d) {
        next_t_idx_in_cur_t = (next_t_idx_in_cur_t + 1) % 4;
      }
      int next_t_idx = cur_t->neighbours[next_t_idx_in_cur_t];

      /* Get the next non axis vertex and add it to the queue if necessary */
      int next_non_axis_idx_in_d = cur_t->vertices[next_t_idx_in_cur_t];
      if (!neighbour_flags[next_non_axis_idx_in_d]) {
        int3 new_info = {._0 = cur_t_idx,
                         ._1 = next_non_axis_idx_in_d,
                         ._2 = next_t_idx_in_cur_t};
        int3_fifo_queue_push(&neighbour_info_q, new_info);
        neighbour_flags[next_non_axis_idx_in_d] |= 1;
      }

      /* Get the coordinates of the voronoi vertex of the new face */
      int vor_vertex0_idx = first_t_idx - 4;
      face_vertices[0] = voronoi_vertices[3 * vor_vertex0_idx];
      face_vertices[1] = voronoi_vertices[3 * vor_vertex0_idx + 1];
      face_vertices[2] = voronoi_vertices[3 * vor_vertex0_idx + 2];
      int face_vertices_index = 1;

      /* Loop around the axis */
      while (next_t_idx != first_t_idx) {
        /* Get the coordinates of the voronoi vertex corresponding to cur_t and
         * next_t */
        if (face_vertices_index + 6 > face_vertices_size) {
          face_vertices_size <<= 1;
          face_vertices = (double *)realloc(
              face_vertices, 3 * face_vertices_size * sizeof(double));
        }
        const int vor_vertex1_idx = cur_t_idx - 4;
        face_vertices[3 * face_vertices_index] =
            voronoi_vertices[3 * vor_vertex1_idx];
        face_vertices[3 * face_vertices_index + 1] =
            voronoi_vertices[3 * vor_vertex1_idx + 1];
        face_vertices[3 * face_vertices_index + 2] =
            voronoi_vertices[3 * vor_vertex1_idx + 2];
        const int vor_vertex2_idx = next_t_idx - 4;
        face_vertices[3 * face_vertices_index + 3] =
            voronoi_vertices[3 * vor_vertex2_idx];
        face_vertices[3 * face_vertices_index + 4] =
            voronoi_vertices[3 * vor_vertex2_idx + 1];
        face_vertices[3 * face_vertices_index + 5] =
            voronoi_vertices[3 * vor_vertex2_idx + 2];
        face_vertices_index += 2;

        /* Update cell volume and tetrahedron_centroid */
        double tetrahedron_centroid[3] = {0., 0., 0.};
        const double V = geometry3d_compute_centroid_volume_tetrahedron(
            ax, ay, az, face_vertices[0], face_vertices[1], face_vertices[2],
            face_vertices[3 * face_vertices_index - 6],
            face_vertices[3 * face_vertices_index - 5],
            face_vertices[3 * face_vertices_index - 4],
            face_vertices[3 * face_vertices_index - 3],
            face_vertices[3 * face_vertices_index - 2],
            face_vertices[3 * face_vertices_index - 1], tetrahedron_centroid);
        this_cell->volume += V;
        this_cell->centroid[0] += V * tetrahedron_centroid[0];
        this_cell->centroid[1] += V * tetrahedron_centroid[1];
        this_cell->centroid[2] += V * tetrahedron_centroid[2];
        voronoi_assert(V >= 0.);

        /* Update variables */
        prev_t_idx_in_cur_t = cur_t->index_in_neighbour[next_t_idx_in_cur_t];
        cur_t_idx = next_t_idx;
        cur_t = &d->tetrahedra[cur_t_idx];
        next_t_idx_in_cur_t = (prev_t_idx_in_cur_t + 1) % 4;
        while (cur_t->vertices[next_t_idx_in_cur_t] == gen_idx_in_d ||
               cur_t->vertices[next_t_idx_in_cur_t] == axis_idx_in_d) {
          next_t_idx_in_cur_t = (next_t_idx_in_cur_t + 1) % 4;
        }
        next_t_idx = cur_t->neighbours[next_t_idx_in_cur_t];
        /* Get the next non axis vertex and add it to the queue if necessary */
        next_non_axis_idx_in_d = cur_t->vertices[next_t_idx_in_cur_t];
        if (!neighbour_flags[next_non_axis_idx_in_d]) {
          int3 new_info = {._0 = cur_t_idx,
                           ._1 = next_non_axis_idx_in_d,
                           ._2 = next_t_idx_in_cur_t};
          int3_fifo_queue_push(&neighbour_info_q, new_info);
          neighbour_flags[next_non_axis_idx_in_d] |= 1;
        }
      }
      if (voronoi_new_face(v, d, gen_idx_in_d, axis_idx_in_d, face_vertices,
                           face_vertices_index)) {
        /* The face is not degenerate */
        nface++;
      }
    }
    voronoi_assert(this_cell->volume > 0.);
    this_cell->centroid[0] /= this_cell->volume;
    this_cell->centroid[1] /= this_cell->volume;
    this_cell->centroid[2] /= this_cell->volume;
    this_cell->nface = nface;

    /* reset flags for all neighbours of this cell */
    neighbour_flags[gen_idx_in_d] = 0;
    for (int i = 0; i < neighbour_info_q.end; i++) {
      voronoi_assert(neighbour_info_q.values[i]._1 < d->vertex_index);
      neighbour_flags[neighbour_info_q.values[i]._1] = 0;
    }
#ifdef VORONOI_CHECKS
    for (int i = 0; i < d->vertex_index; i++) {
      voronoi_assert(neighbour_flags[i] == 0);
    }
#endif
  }
  free(voronoi_vertices);
  free(neighbour_flags);
  int3_fifo_queue_destroy(&neighbour_info_q);
  free(face_vertices);
  voronoi_check_grid(v);
}

/**
 * @brief Free up all memory used by the Voronoi grid.
 *
 * @param v Voronoi grid.
 */
inline static void voronoi_destroy(struct voronoi *restrict v) {
  /* Anything allocated? */
  if (!v->active) {
    return;
  }

  swift_free("c.h.v.cells", v->cells);
  v->cells = NULL;
  v->number_of_cells = 0;
  v->cells_size = 0;

  for (int i = 0; i < 28; ++i) {
    for (int j = 0; j < v->pair_index[i]; j++) {
      voronoi_pair_destroy(&v->pairs[i][j]);
      v->pair_index[i] = -1;
      v->pair_size[i] = 0;
    }
    swift_free("c.h.v.pairs", v->pairs[i]);
    v->pairs[i] = NULL;
  }
  int2_lifo_queue_destroy(&v->cell_pair_connections);
  v->swift_cell = NULL;
  v->active = 0;
}

/**
 * @brief Add a face (two particle pair) to the mesh.
 *
 * The grid connectivity is stored per cell sid: sid=13 corresponds to particle
 * pairs encountered during a self task (both particles are within the local
 * cell), while sid=0-12 and 14-26 correspond to particle interactions for which
 * the right neighbour is part of one of the 26 neighbouring cells.
 *
 * For each pair, we compute and store all the quantities required to compute
 * fluxes between the Voronoi cells: the surface area and midpoint of the
 * interface.
 *
 * @param v Voronoi grid.
 * @param sid 0 for pairs entirely in this cell, 1 for pairs between this cell
 * and a neighbouring cell (in SWIFT we use the convention from the
 * description).
 * @param cell Pointer to the cell of the right particle (NULL if the right
 * particle lives in the same cell as the left particle). For SWIFT only.
 * @param left_part_pointer Index of left particle in cell (particle in the
 * cell linked to this grid). FUTURE NOTE: For SWIFT, replace this with direct
 * pointer to the left particle.
 * @param right_part_pointer Index of right particle in cell (particle in the
 * cell linked to this grid), or -1 for ghost vertices. FUTURE NOTE: For SWIFT,
 * replace this with direct pointer to the right particle.
 * @param vertices Vertices of the interface.
 * @param n_vertices Number of vertices in the vertices array.
 * @returns 1 if a non-degenerate face was added or found, else 0
 */
inline static int voronoi_new_face(struct voronoi *v, const struct delaunay *d,
                                   int left_part_idx_in_d,
                                   int right_part_idx_in_d, double *vertices,
                                   int n_vertices) {
  int sid;
  struct cell *c;
  if (right_part_idx_in_d < d->ngb_offset) {
    if (right_part_idx_in_d < left_part_idx_in_d) {
      /* Pair was already added. Find it and add it to the cell_pair_connections
       * if necessary. If no pair is found, the face must have been degenerate.
       * Return early. */
      struct voronoi_cell_new *ngb_cell =
          &v->cells[right_part_idx_in_d - d->vertex_start];
      struct part *left_part = d->part_pointers[left_part_idx_in_d];
      for (int i = 0; i < ngb_cell->nface; i++) {
        int2 connection = v->cell_pair_connections
                              .values[ngb_cell->pair_connections_offset + i];
        if (v->pairs[connection._1][connection._0].right == left_part) {
          int2_lifo_queue_push(&v->cell_pair_connections, connection);
          return 1;
        }
      }
      return 0;
    }
    sid = 13;
    c = NULL;
  } else {
    sid = d->ngb_cell_sids[right_part_idx_in_d - d->ngb_offset];
    c = d->ngb_cell_ptrs[right_part_idx_in_d - d->ngb_offset];
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Sanity check */
  assert(c == NULL || !c->split);
#endif

  if (v->pair_index[sid] == v->pair_size[sid]) {
    v->pair_size[sid] <<= 1;
    v->pairs[sid] = (struct voronoi_pair *)swift_realloc(
        "c.h.v.pairs", v->pairs[sid],
        v->pair_size[sid] * sizeof(struct voronoi_pair));
  }

  /* Initialize pair */
  struct voronoi_pair *this_pair = &v->pairs[sid][v->pair_index[sid]];
  voronoi_pair_init(this_pair, c, d->part_pointers[left_part_idx_in_d],
                    d->part_pointers[right_part_idx_in_d], left_part_idx_in_d,
                    right_part_idx_in_d, vertices, n_vertices);

  /* is the face degenerate? */
  if (this_pair->surface_area < v->min_surface_area) {
    return 0;
  }
  /* increase index (i.e. add the face) */
  v->pair_index[sid]++;

  /* Add cell_pair_connection */
  int2 connection = {._0 = v->pair_index[sid], ._1 = sid};
  int2_lifo_queue_push(&v->cell_pair_connections, connection);

  return 1;
}

static inline double voronoi_compute_volume(const struct voronoi *restrict v) {
  double total_volume = 0.;
  for (int i = 0; i < v->number_of_cells; ++i) {
    total_volume += v->cells[i].volume;
  }
  return total_volume;
}

/**
 * @brief Sanity checks on the grid.
 *
 * Right now, this only checks the total volume of the cells.
 */
inline static void voronoi_check_grid(struct voronoi *restrict v) {
#ifdef VORONOI_CHECKS
  /* Check total volume */
  double total_volume = 0.;
  for (int i = 0; i < v->number_of_cells; i++) {
    total_volume += v->cells[i].volume;
  }
  //  fprintf(stderr, "Total volume: %g\n", total_volume);

  /* For each cell check that the total surface area is not bigger than the
   * surface area of a sphere with the same volume */
  double *surface_areas = malloc(v->number_of_cells * sizeof(double));
  for (int i = 0; i < v->number_of_cells; i++) {
    surface_areas[i] = 0;
  }
  for (int sid = 0; sid < 28; sid++) {
    for (int i = 0; i < v->pair_index[sid]; i++) {
      struct voronoi_pair *pair = &v->pairs[sid][i];
      surface_areas[pair->left_idx] += pair->surface_area;
      if (sid == 13) {
        surface_areas[pair->right_idx] += pair->surface_area;
      }
    }
  }
  for (int i = 0; i < v->number_of_cells; i++) {
    double sphere_surface_area =
        4. * M_PI * pow(3. * v->cells[i].volume / (4. * M_PI), 2. / 3.);
    voronoi_assert(sphere_surface_area < surface_areas[i]);
  }
  free(surface_areas);

#if defined(VORONOI_STORE_GENERATORS)
  /* Check connectivity */
  for (int i = 0; i < v->number_of_cells; i++) {
    struct voronoi_cell_new *this_cell = &v->cells[i];
    for (int j = this_cell->pair_connections_offset;
         j < this_cell->pair_connections_offset + this_cell->nface; j++) {
      int2 connection = v->cell_pair_connections.values[j];
      int pair_idx = connection._0;
      int sid = connection._1;
      struct voronoi_pair *pair = &v->pairs[sid][pair_idx];
      assert(i == pair->left_idx || (sid == 13 && i == pair->right_idx));
    }
  }
#endif

#endif
}

/**
 * @brief Write the Voronoi grid information to the given file.
 *
 * The output depends on the configuration. The maximal output contains 3
 * different types of output lines:
 *  - "G\tgx\tgx: x and y position of a single grid generator (optional).
 *  - "C\tcx\tcy\tV\tnface": centroid position, volume and (optionally) number
 *    of faces for a single Voronoi cell.
 *  - "F\tsid\tarea\tcx\tcy\tcz\t(v0x, v0y, v0z)\t...\t(vNx, vNy, vNz)": sid,
 *    area, coordinates of centroid, coordinates of vertices of face (optional).
 *
 * @param v Voronoi grid.
 * @param file File to write to.
 */
inline static void voronoi_write_grid(const struct voronoi *restrict v,
                                      FILE *file) {
  /* first write the cells (and generators, if those are stored) */
  for (int i = 0; i < v->number_of_cells; ++i) {
    struct voronoi_cell_new *this_cell = &v->cells[i];
#ifdef VORONOI_STORE_GENERATORS
    fprintf(file, "G\t%g\t%g\t%g\n", this_cell->generator->x[0],
            this_cell->generator->x[1], this_cell->generator->x[2]);
#endif
    fprintf(file, "C\t%g\t%g\t%g\t%g", this_cell->centroid[0],
            this_cell->centroid[1], this_cell->centroid[2], this_cell->volume);
    fprintf(file, "\t%i", this_cell->nface);
    fprintf(file, "\n");
  }

  /* now write the faces */
  for (int ngb = 0; ngb < 28; ++ngb) {
    for (int i = 0; i < v->pair_index[ngb]; ++i) {
      struct voronoi_pair *pair = &v->pairs[ngb][i];
      fprintf(file, "F\t%i\t%g\t%g\t%g\t%g", ngb, pair->surface_area,
              pair->midpoint[0], pair->midpoint[1], pair->midpoint[2]);
#ifdef VORONOI_STORE_FACES
      for (int j = 0; j < pair->n_vertices; j++) {
        fprintf(file, "\t(%g, %g, %g)", pair->vertices[3 * j],
                pair->vertices[3 * j + 1], pair->vertices[3 * j + 2]);
      }
#endif
      fprintf(file, "\n");
    }
  }
}

/**
 * @brief Print the Voronoi grid to a file with the given name.
 *
 * @param v Voronoi grid (read-only).
 * @param file_name Name of the output file.
 */
static inline void voronoi_print_grid(const struct voronoi *restrict v,
                                      const char *file_name) {

  FILE *file = fopen(file_name, "w");

  voronoi_write_grid(v, file);

  fclose(file);
}

#endif  // SWIFTSIM_VORONOI_FUNCTIONS_3D_H
