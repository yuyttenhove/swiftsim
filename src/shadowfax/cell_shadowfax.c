#include "cell_shadowfax.h"

#include "../space_getsid.h"

void cell_shadowfax_do_pair1_density_recursive(const struct engine *e,
                                               struct cell *restrict ci,
                                               struct cell *restrict cj,
                                               int sid, const double *shift) {
  int k;
  /* recurse? */
  if (ci->split) {
    // TODO check if possible with max smoothing length
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

void cell_shadowfax_do_pair_subset_density_recursive(
    const struct engine *e, struct cell *restrict ci,
    struct part *restrict parts_i, const int *restrict ind, int count,
    struct cell *restrict cj, const int sid, const int flipped,
    const double *shift) {
  int k;
  /* recurse? */
  if (ci->split) {
    // TODO check if possible with max smoothing length
    int sub_start_ind = 0;
    for (k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        struct cell *sub = ci->progeny[k];
        /* Does this sub contain some of the remaining particles? */
        if (&parts_i[ind[sub_start_ind]] >= &sub->hydro.parts[0] &&
            &parts_i[ind[sub_start_ind]] < &sub->hydro.parts[sub->hydro.count]) {
          /* Does this sub contain all of the remaining particles? */
          if (&parts_i[ind[count - 1]] < &sub->hydro.parts[sub->hydro.count]) {
            cell_shadowfax_do_pair_subset_density_recursive(
                e, sub, parts_i, ind, count, cj, sid, flipped, shift);
          } else {
            /* Some, but not all remaining particles are in this sub
             * Find the number of particles in this sub. */
            int sub_count = 1;
            while (&parts_i[ind[sub_start_ind + sub_count]] <
                   &sub->hydro.parts[sub->hydro.count]) {
              sub_count++;
            }
            cell_shadowfax_do_pair_subset_density_recursive(
                e, sub, parts_i, &ind[sub_start_ind], sub_count, cj, sid,
                flipped, shift);
            sub_start_ind += sub_count;
          } /* Does this sub contain all of the remaining particles? */
        }   /* Does this sub contain some of the remaining particles? */
      }     /* Progeny not NULL? */
    }
  } else if (cj->split) {
    for (k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        // TODO check if possible with max smoothing length
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
            cell_shadowfax_do_pair1_density_recursive(e, ck, cl, sid, shift);
          }
        }
      }
    }
  } else {
    cell_shadowfax_do_self1_density(e, c);
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
