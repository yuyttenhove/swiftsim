#include "cell_shadowfax.h"

#include "../space_getsid.h"

void cell_malloc_delaunay_tessellation_recursive(struct cell *restrict c, const struct engine *restrict e) {
  /* anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) return;

  /* recurse? */
  if (c->split) {
    if (c->hydro.shadowfax_enabled) {
      /* Tessellations were allocated here in a previous iteration. We don't
       * need them any more, so free them now. */
      cell_destroy_tessellations(c);
    }
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_malloc_delaunay_tessellation_recursive(c->progeny[k], e);
      }
    }
  } else {
    /* Base case */
    cell_malloc_delaunay_tessellation(c);
  }
}

void cell_shadowfax_do_pair1_density_recursive(const struct engine *e,
                                               struct cell *restrict ci,
                                               struct cell *restrict cj,
                                               int sid, const double *shift) {
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;

  int k;
  /* recurse? */
  if (ci->split) {
    for (k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        cell_shadowfax_do_pair1_density_recursive(e, ci->progeny[k], cj, sid,
                                                  shift);
      }
    }
  } else if (cj->split) {
    for (k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        cell_shadowfax_do_pair1_density_recursive(e, ci, cj->progeny[k], sid,
                                                  shift);
      }
    }
  } else {
    cell_shadowfax_do_pair1_density(e, ci, cj, sid, shift);
  }
}

void cell_shadowfax_do_pair1_gradient_recursive(const struct engine *e,
                                                struct cell *restrict ci,
                                                struct cell *restrict cj,
                                                int sid, const double *shift) {
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;

  int k;
  /* recurse? */
  if (ci->split) {
    for (k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        cell_shadowfax_do_pair1_gradient_recursive(e, ci->progeny[k], cj, sid,
                                                   shift);
      }
    }
  } else if (cj->split) {
    for (k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        cell_shadowfax_do_pair1_gradient_recursive(e, ci, cj->progeny[k], sid,
                                                   shift);
      }
    }
  } else {
    cell_shadowfax_do_pair1_gradient(e, ci, cj, sid, shift);
  }
}

void cell_shadowfax_do_pair2_force_recursive(const struct engine *e,
                                             struct cell *restrict ci,
                                             struct cell *restrict cj, int sid,
                                             const double *shift) {
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;

  int k;
  /* recurse? */
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
    const struct engine *e, struct cell *restrict ci,
    struct part *restrict parts_i, const int *restrict ind, int count,
    struct cell *restrict cj, const int sid, const int flipped,
    const double *shift) {
  if (!cell_is_active_hydro(ci, e)) return;

  int k;
  /* recurse? */
  if (ci->split) {
    // TODO check if possible with max smoothing length
    int sub_start_ind = 0;
    for (k = 0; k < 8 && count > 0; k++) {
      if (ci->progeny[k] != NULL) {
        struct cell *sub = ci->progeny[k];
        /* Does this sub-cell contain some of the remaining particles?
         * (Particles are stored in order of sub-cells) */
        if (&parts_i[ind[sub_start_ind]] >= &sub->hydro.parts[0] &&
            &parts_i[ind[sub_start_ind]] <
                &sub->hydro.parts[sub->hydro.count]) {
          int sub_count;
          /* Does this sub-cell contain all of the remaining particles? */
          if (&parts_i[ind[count - 1]] < &sub->hydro.parts[sub->hydro.count]) {
            sub_count = count;
          } else {
            /* Some, but not all remaining particles are in this sub-cell
             * Find the number of particles in this sub-cell. */
            sub_count = 1;
            while (&parts_i[ind[sub_start_ind + sub_count]] <
                   &sub->hydro.parts[sub->hydro.count]) {
              sub_count++;
            }
          } /* Does this sub contain all of the remaining particles? */
          cell_shadowfax_do_pair_subset_density_recursive(
              e, sub, parts_i, &ind[sub_start_ind], sub_count, cj, sid, flipped,
              shift);
          sub_start_ind += sub_count;
          count -= sub_count;
        } /* Does this sub contain some of the remaining particles? */
      }   /* Progeny not NULL? */
    }
  } else if (cj->split) {
    for (k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        cell_shadowfax_do_pair_subset_density_recursive(
            e, ci, parts_i, ind, count, cj->progeny[k], sid, flipped, shift);
      }
    }
  } else {
    cell_shadowfax_do_pair_subset_density(e, ci, parts_i, ind, count, cj, sid,
                                          flipped, shift);
  }
}

void cell_shadowfax_do_self1_density_recursive(const struct engine *e,
                                               struct cell *restrict c) {
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    double shift[3] = {0., 0., 0.};
    int sid;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_shadowfax_do_self1_density_recursive(e, c->progeny[k]);
        for (int l = k + 1; l < 8; l++) {
          if (c->progeny[l] != NULL) {
            struct cell *ck = c->progeny[k];
            struct cell *cl = c->progeny[l];
            /* this might swap ck and cl! */
            sid = space_getsid(e->s, &ck, &cl, shift);
            cell_shadowfax_do_pair1_density_recursive(e, ck, cl, sid, shift);
          }
        }
      }
    }
  } else {
    cell_shadowfax_do_self1_density(e, c);
  }
}

void cell_shadowfax_do_self1_gradient_recursive(const struct engine *e,
                                                struct cell *restrict c) {
  if (!cell_is_active_hydro(c, e)) return;

  double shift[3] = {0., 0., 0.};
  int sid;
  /* recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_shadowfax_do_self1_gradient_recursive(e, c->progeny[k]);
        for (int l = k + 1; l < 8; l++) {
          if (c->progeny[l] != NULL) {
            struct cell *ck = c->progeny[k];
            struct cell *cl = c->progeny[l];
            /* this might swap ck and cl! */
            sid = space_getsid(e->s, &ck, &cl, shift);
            cell_shadowfax_do_pair1_gradient_recursive(e, ck, cl, sid, shift);
          }
        }
      }
    }
  } else {
    cell_shadowfax_do_self1_gradient(e, c);
  }
}

void cell_shadowfax_do_self2_force_recursive(const struct engine *e,
                                             struct cell *restrict c) {
  if (!cell_is_active_hydro(c, e)) return;

  double shift[3] = {0., 0., 0.};
  int sid;
  /* recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_shadowfax_do_self2_force_recursive(e, c->progeny[k]);
        for (int l = k + 1; l < 8; l++) {
          if (c->progeny[l] != NULL) {
            struct cell *ck = c->progeny[k];
            struct cell *cl = c->progeny[l];
            /* this might swap ck and cl! */
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

void cell_shadowfax_do_self_subset_density_recursive(
    const struct engine *e, struct cell *restrict c,
    struct part *restrict parts, const int *restrict ind, int count) {
  if (!cell_is_active_hydro(c, e)) return;

  int k, l, sid, flipped;
  double shift[3] = {0., 0., 0.};
  /* recurse? */
  if (c->split) {
    int sub_start_ind = 0;
    for (k = 0; k < 8 && count > 0; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *sub = c->progeny[k];
        /* Does this sub-cell contain some of the remaining particles?
         * (Particles are stored in order of sub-cells) */
        if (&parts[ind[sub_start_ind]] >= &sub->hydro.parts[0] &&
            &parts[ind[sub_start_ind]] < &sub->hydro.parts[sub->hydro.count]) {
          int sub_count;
          /* Does this sub-cell contain all of the remaining particles? */
          if (&parts[ind[count - 1]] < &sub->hydro.parts[sub->hydro.count]) {
            sub_count = count;
          } else {
            /* Some, but not all remaining particles are in this sub-cell
             * Find the number of particles in this sub-cell. */
            sub_count = 1;
            while (&parts[ind[sub_start_ind + sub_count]] <
                   &sub->hydro.parts[sub->hydro.count]) {
              sub_count++;
            }
          } /* Does this sub contain all of the remaining particles? */

          /* recursive self interaction of this sub cells */
          cell_shadowfax_do_self_subset_density_recursive(
              e, sub, parts, &ind[sub_start_ind], sub_count);

          /* pair interactions of this sub-cell with the other sub-cells of c*/
          for (l = 0; l < 8; l++) {
            if (l != k && c->progeny[l] != NULL) {
              struct cell *ck = c->progeny[k];
              struct cell *cj = c->progeny[l];
              sid = space_getsid(e->s, &ck, &cj, shift);
              flipped = (ck != sub);
              cell_shadowfax_do_pair_subset_density_recursive(
                  e, sub, parts, &ind[sub_start_ind], sub_count, c->progeny[l],
                  sid, flipped, shift);
            }
          }
          /* Update indices */
          sub_start_ind += sub_count;
          count -= sub_count;
        } /* Does this sub contain some of the remaining particles? */
      }
    }
  } else {
    cell_shadowfax_do_self_subset_density(e, c, parts, ind, count);
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
  } else if (c->hydro.shadowfax_enabled){
        delaunay_write_tessellation(&c->hydro.deltess, dfile, offset);
#ifdef SWIFT_DEBUG_CHECKS
    assert(c->hydro.vortess.active == c->hydro.deltess.active);
    assert(c->hydro.vortess.active);
#endif
    voronoi_write_grid(&c->hydro.vortess, vfile);
  }
}

double cell_shadowfax_voronoi_volume(const struct cell *c) {
  double total_volume = 0.;
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        total_volume += cell_shadowfax_voronoi_volume(c->progeny[k]);
      }
    }
  } else if (c->hydro.shadowfax_enabled) {
    total_volume = voronoi_compute_volume(&c->hydro.vortess);
  }
  return total_volume;
 }
