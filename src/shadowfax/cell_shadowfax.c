#include "cell_shadowfax.h"

#include "../space_getsid.h"

void cell_malloc_delaunay_tessellation_recursive(struct cell *c) {
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_malloc_delaunay_tessellation_recursive(c->progeny[k]);
      }
    }
  } else {
    cell_malloc_delaunay_tessellation(c);
  }
}

void cell_shadowfax_do_pair1_density_recursive(const struct engine *e,
                                               struct cell *ci, struct cell *cj,
                                               int sid, const double *shift,
                                               unsigned long sub_cell_mask_i,
                                               unsigned long sub_cell_mask_j) {
  int k;
  /* recurse? */
  if (ci->split) {
    sub_cell_mask_i *= sub_cell_mask_i;
    // TODO check if possible with max smoothing length
    for (k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        cell_shadowfax_do_pair1_density_recursive(e, ci->progeny[k], cj, sid,
                                                  shift, sub_cell_mask_i << k,
                                                  sub_cell_mask_j);
      }
    }
  } else if (cj->split) {
    sub_cell_mask_j *= sub_cell_mask_j;
    for (k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        cell_shadowfax_do_pair1_density_recursive(e, ci, cj->progeny[k], sid,
                                                  shift, sub_cell_mask_i,
                                                  sub_cell_mask_j << k);
      }
    }
  } else {
    cell_shadowfax_do_pair1_density(e, ci, cj, sid, shift, sub_cell_mask_i,
                                    sub_cell_mask_j);
  }
}

void cell_shadowfax_do_pair2_force_recursive(const struct engine *e,
                                             struct cell *restrict ci,
                                             struct cell *restrict cj, int sid,
                                             const double *shift) {
  int k;
  /* recurse */
  if (ci->split) {
    for (k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        cell_shadowfax_do_pair2_force_recursive(e, ci->progeny[k], cj, sid,
                                                shift);
      }
    }
  } else if (cj->split) {
    for (k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        cell_shadowfax_do_pair2_force_recursive(e, ci, cj->progeny[k], sid,
                                                shift);
      }
    }
  } else {
    cell_shadowfax_do_pair2_force(e, ci, cj, sid, shift);
  }
}

void cell_shadowfax_do_pair_subset_density_recursive(
    const struct engine *e, struct cell *ci, struct part *parts_i,
    const int *ind, int count, struct cell *cj, int sid, int flipped,
    const double *shift, unsigned long sub_cell_mask_i,
    unsigned long sub_cell_mask_j) {
  int k;
  /* recurse? */
  if (ci->split) {
    sub_cell_mask_i *= sub_cell_mask_i;
    int sub_start_ind = 0;
    for (k = 0; k < 8 && count > 0; k++) {
      // TODO check if possible with max smoothing length
      if (ci->progeny[k] != NULL) {
        struct cell *sub = ci->progeny[k];
        /* Does this sub contain some of the remaining particles?
         * (Particles are stored in order of sub cells) */
        if (&parts_i[ind[sub_start_ind]] >= &sub->hydro.parts[0] &&
            &parts_i[ind[sub_start_ind]] <
                &sub->hydro.parts[sub->hydro.count]) {
          int sub_count;
          /* Does this sub contain all of the remaining particles? */
          if (&parts_i[ind[count - 1]] < &sub->hydro.parts[sub->hydro.count]) {
            sub_count = count;
          } else {
            /* Some, but not all remaining particles are in this sub
             * Find the number of particles in this sub. */
            sub_count = 1;
            while (&parts_i[ind[sub_start_ind + sub_count]] <
                   &sub->hydro.parts[sub->hydro.count]) {
              sub_count++;
            }
          } /* Does this sub contain all of the remaining particles? */
          cell_shadowfax_do_pair_subset_density_recursive(
              e, sub, parts_i, &ind[sub_start_ind], sub_count, cj, sid, flipped,
              shift, sub_cell_mask_i << k, sub_cell_mask_j);
          sub_start_ind += sub_count;
          count -= sub_count;
        } /* Does this sub contain some of the remaining particles? */
      }   /* Progeny not NULL? */
    }
  } else if (cj->split) {
    sub_cell_mask_j *= sub_cell_mask_j;
    for (k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        // TODO check if possible with max smoothing length
        cell_shadowfax_do_pair_subset_density_recursive(
            e, ci, parts_i, ind, count, cj->progeny[k], sid, flipped, shift,
            sub_cell_mask_i, sub_cell_mask_j << k);
      }
    }
  } else {
    cell_shadowfax_do_pair_subset_density(e, ci, parts_i, ind, count, cj, sid,
                                          flipped, shift, sub_cell_mask_i,
                                          sub_cell_mask_j);
  }
}

void cell_shadowfax_do_self1_density_recursive(const struct engine *e,
                                               struct cell *restrict c) {
  double shift[3] = {0., 0., 0.};
  int sid;
  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      struct cell *ck = c->progeny[k];
      if (ck != NULL) {
        cell_shadowfax_do_self1_density_recursive(e, c->progeny[k]);
        for (int l = k + 1; l < 8; l++) {
          struct cell *cl = c->progeny[l];
          if (cl != NULL) {
            sid = space_getsid(e->s, &ck, &cl, shift);
            cell_shadowfax_do_pair1_density_recursive(e, ck, cl, sid, shift, 1,
                                                      1);
          }
        }
      }
    }
  } else {
    cell_shadowfax_do_self1_density(e, c);
  }
}

void cell_shadowfax_do_self2_force_recursive(const struct engine *e,
                                             struct cell *restrict c) {
  double shift[3] = {0., 0., 0.};
  int sid;
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      struct cell *ck = c->progeny[k];
      if (ck != NULL) {
        cell_shadowfax_do_self2_force_recursive(e, c->progeny[k]);
        for (int l = k + 1; l < 8; l++) {
          struct cell *cl = c->progeny[l];
          if (cl != NULL) {
            sid = space_getsid(e->s, &ck, &cl, shift);
            cell_shadowfax_do_pair2_force_recursive(e, ck, cl, sid, shift);
          }
        }
      }
    }
  } else {
    cell_shadowfax_do_self2_force(e, c);
  }
}

void cell_shadowfax_end_density_recursive(struct cell *restrict c) {
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_shadowfax_end_density_recursive(c->progeny[k]);
  } else {
    cell_shadowfax_end_density(c);
  }
}

void cell_shadowfax_write_tesselations(const struct cell *c, FILE *dfile,
                                       FILE *vfile, size_t *offset) {
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_shadowfax_write_tesselations(c->progeny[k], dfile, vfile, offset);
      }
    }
  } else {
    //    delaunay_write_tessellation(&c->hydro.deltess, dfile, offset);
    voronoi_write_grid(&c->hydro.vortess, vfile);
  }
}
