#include "cell_shadowfax.h"

#include "../space_getsid.h"

void cell_shadowfax_do_pair1_density_recursive(const struct engine *e,
                                               struct cell *restrict ci,
                                               struct cell *restrict cj,
                                               int sid, const double *shift) {
  // TODO
}

void cell_shadowfax_do_pair_subset_density_recursive(
    const struct engine *e, struct cell *restrict ci,
    struct part *restrict parts_i, const int *restrict ind, int count,
    struct cell *restrict cj, const int sid, const int flipped,
    const double *shift) {
  // TODO
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
