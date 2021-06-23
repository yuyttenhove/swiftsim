#ifndef SWIFT_CELL_SHADOWFAX_H
#define SWIFT_CELL_SHADOWFAX_H

#include "active.h"
#include "cell.h"
#include "delaunay_functions.h"
#include "hydro/Shadowfax/hydro_gradients.h"
#include "hydro_shadowfax.h"
#include "voronoi_functions.h"

#ifdef SHADOWFAX_CHECK_DELAUNAY_FLAG
__attribute__((always_inline)) INLINE static void shadowfax_flag_particle_added(
    struct part *restrict p, int sid) {
  p->voronoi.flag |= (1 << sid);
}
#endif

#ifdef SHADOWFAX_CHECK_DELAUNAY_FLAG
__attribute__((always_inline)) INLINE static int shadowfax_particle_was_added(
    const struct part *restrict p, int sid) {
  return (p->voronoi.flag & (1 << sid)) > 0;
}
#endif

__attribute__((always_inline)) INLINE static void get_shift(
    const double *restrict p1, const double *restrict p2,
    struct space *restrict s, /*return*/ double *restrict shift) {
  const int periodic = s->periodic;
  double dx;
  for (int k = 0; k < 3; k++) {
    dx = p2[k] - p1[k];
    if (periodic && dx < -s->dim[k] / 2)
      shift[k] = s->dim[k];
    else if (periodic && dx > s->dim[k] / 2)
      shift[k] = -s->dim[k];
    else
      shift[k] = 0.0;
  }
}

__attribute__((always_inline)) INLINE static void
cell_malloc_delaunay_tessellation(struct cell *c) {

  const int count = c->hydro.count;
  double *loc, *width;

  if (c->hydro.super != NULL) {
    loc = c->hydro.super->loc;
    width = c->hydro.super->width;
  } else {
    error(
        "Trying to allocate a delaunay tesselation for a cell above the super "
        "level!");
  }

  if (c->hydro.shadowfax_enabled == 1) {
    delaunay_reset(&c->hydro.deltess, loc, width, count);
  } else {
    delaunay_init(&c->hydro.deltess, loc, width, count, 10 * count);
  }

  //  struct part *restrict parts = c->hydro.parts;
  //
  //  for (int pd = 0; pd < count; pd++) {
  //    /* Get a pointer to the ith particle. */
  //    struct part *restrict p = &parts[pd];
  //    p->voronoi.flag = 0;
  //    p->voronoi.nface = 0;
  //    p->voronoi.volume = 0;
  //  }

  c->hydro.shadowfax_enabled = 1;
}

void cell_malloc_delaunay_tessellation_recursive(struct cell *c);

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair1_density(const struct engine *e, struct cell *ci,
                                struct cell *cj, int sid, const double *shift) {

  if (ci == cj) error("Interacting cell with itself!");

  /*  cell_shadowfax_do_pair_naive(e, ci, cj, sid, shift);*/
  /*  return;*/

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

  /* Get some other useful values. */
  const double hi_max = ci->hydro.h_max * kernel_gamma - rshift;
  const double hj_max = cj->hydro.h_max * kernel_gamma;
  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;
  const double di_max = sort_i[count_i - 1].d - rshift;
  const double dj_min = sort_j[0].d;
  const float dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  if (cell_is_active_hydro(ci, e)) {

    /* Loop over the parts in ci. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[sort_i[pid].i];
      const float hi = pi->h;

      /* Skip inactive particles */
      if (!part_is_active(pi, e)) {
        /* TODO what should we do here?
         * For the moment, also build grid for inactive particles... */
        // continue;
      }

      /* Is there anything we need to interact with ? */
      const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
      if (di < dj_min) continue;

      /* Get some additional information about pi */
      const float hig2 = hi * hi * kernel_gamma2;
      const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
      const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
      const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float pjx = pj->x[0] - cj->loc[0];
        const float pjy = pj->x[1] - cj->loc[1];
        const float pjz = pj->x[2] - cj->loc[2];

        /* Compute the pairwise distance. */
        float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < hig2) {
          delaunay_add_new_vertex(&ci->hydro.deltess, pj->x[0] + shift[0],
                                  pj->x[1] + shift[1], 26 - sid, cj, pj);
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }     /* Cell ci is active */

  if (cell_is_active_hydro(cj, e)) {

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth part in cj. */
      struct part *pj = &parts_j[sort_j[pjd].i];
      const float hj = pj->h;

      /* Skip inactive particles */
      if (!part_is_active(pj, e)) {
        /* TODO what should we do here?
         * For the moment, also build grid for inactive particles... */
        // continue;
      }

      /* Is there anything we need to interact with ? */
      const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max + rshift;
      if (dj - rshift > di_max) continue;

      /* Get some additional information about pj */
      const float hjg2 = hj * hj * kernel_gamma2;
      const float pjx = pj->x[0] - cj->loc[0];
      const float pjy = pj->x[1] - cj->loc[1];
      const float pjz = pj->x[2] - cj->loc[2];

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
        const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
        const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

        /* Compute the pairwise distance. */
        float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < hjg2) {
          delaunay_add_new_vertex(&cj->hydro.deltess, pi->x[0] - shift[0],
                                  pi->x[1] - shift[1], sid, ci, pi);
        }
      } /* loop over the parts in ci. */
    }   /* loop over the parts in cj. */
  }     /* Cell cj is active */
}

void cell_shadowfax_do_pair1_density_recursive(const struct engine *e,
                                               struct cell *restrict ci,
                                               struct cell *restrict cj,
                                               int sid, const double *shift);

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair2_density(const struct engine *e,
                                struct cell *restrict ci,
                                struct cell *restrict cj, int sid,
                                const double *shift) {
  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair1_gradient(const struct engine *e,
                                 struct cell *restrict ci,
                                 struct cell *restrict cj, int sid,
                                 const double *shift) {
  if (ci == cj) error("Interacting cell with itself!");

  /* anything to do here? */
  if (!(cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e))) return;

  struct voronoi *vortess = &ci->hydro.vortess;
  /* loop over voronoi faces between ci and cj */
  for (int i = 0; i < vortess->pair_index[26 - sid]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[26 - sid][i];
    struct part *part_left = pair->left;
    struct part *part_right = pair->right;
    /* check if right particle in cj */
    if (pair->right_cell != cj) {
      continue;
    }
    if (part_left->force.active == 1 && part_right->force.active == 1) {
      hydro_shadowfax_gradients_collect(part_left, part_right, pair->midpoint,
                                        pair->surface_area, shift);
    } else if (part_left->force.active) {
      hydro_shadowfax_gradients_nonsym_collect(
          part_left, part_right, pair->midpoint, pair->surface_area, shift);
    } else if (part_right->force.active) {
      double inverse_shift[3];
      for (int k = 0; k < 3; k++) {
        inverse_shift[k] = -shift[k];
      }
      hydro_shadowfax_gradients_nonsym_collect(
          part_right, part_left, pair->midpoint, pair->surface_area,
          inverse_shift);
    }
  } /* loop over voronoi faces between ci and cj */
}

void cell_shadowfax_do_pair1_gradient_recursive(const struct engine *e,
                                                struct cell *restrict ci,
                                                struct cell *restrict cj,
                                                int sid, const double *shift);

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair2_gradient(const struct engine *e,
                                 struct cell *restrict ci,
                                 struct cell *restrict cj, int sid,
                                 const double *shift) {
  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_do_pair1_force(
    const struct engine *e, struct cell *restrict ci, struct cell *restrict cj,
    int sid, const double *shift) {
  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_do_pair2_force(
    const struct engine *e, struct cell *restrict ci, struct cell *restrict cj,
    int sid, const double *shift) {
#ifdef SWIFT_DEBUG_CHECKS
  int n_pairs_ci = 0;
  for (int i = 0; i < ci->hydro.vortess.pair_index[26 - sid]; i++) {
    if (ci->hydro.vortess.pairs[26 - sid][i].right_cell == cj) {
      n_pairs_ci++;
    }
  }
  int n_pairs_cj = 0;
  for (int i = 0; i < cj->hydro.vortess.pair_index[sid]; i++) {
    if (cj->hydro.vortess.pairs[sid][i].right_cell == ci) {
      n_pairs_cj++;
    }
  }
  assert(n_pairs_ci == n_pairs_cj);

  for (int i = 0; i < ci->hydro.vortess.pair_index[26 - sid]; i++) {
    struct voronoi_pair *pair_i = &ci->hydro.vortess.pairs[26 - sid][i];
    if (pair_i->right_cell != cj) continue;

    int found = 0;
    struct voronoi_pair *pair_j;
    for (int j = 0; !found && j < cj->hydro.vortess.pair_index[sid]; j++) {
      pair_j = &cj->hydro.vortess.pairs[sid][j];
      if (pair_i->left == pair_j->right && pair_j->left == pair_i->right) {
        found = 1;
      }
    }
    assert(found || pair_i->surface_area < 1e-15);
  }
  for (int j = 0; j < cj->hydro.vortess.pair_index[sid]; j++) {
    struct voronoi_pair *pair_j = &cj->hydro.vortess.pairs[sid][j];
    if (pair_j->right_cell != ci) continue;

    int found = 0;
    struct voronoi_pair *pair_i;
    for (int i = 0; !found && i < ci->hydro.vortess.pair_index[26-sid]; i++) {
      pair_i = &ci->hydro.vortess.pairs[26-sid][i];
      if (pair_j->left == pair_i->right && pair_i->left == pair_j->right) {
        found = 1;
      }
    }
    assert(found || pair_j->surface_area < 1e-15);
  }
#endif

  /* anything to do here? */
  if (!cell_is_active_hydro(ci, e)) return;

  struct voronoi *vortess = &ci->hydro.vortess;
  /* loop over voronoi faces between ci and cj */
  for (int i = 0; i < vortess->pair_index[26 - sid]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[26 - sid][i];
    struct part *part_left = pair->left;
    struct part *part_right = pair->right;
    if ((part_left->x[0] == 0.3833333333333333 &&
         part_left->x[1] == 0.68333333333333335) ||
        (part_right->x[0] == 0.3833333333333333 &&
         part_right->x[1] == 0.68333333333333335)) {
      int n_pairs = 0;
      for (int k = 0; k < 27; k++) {
        for (int j = 0; j < vortess->pair_index[k]; j++) {
          struct voronoi_pair *p = &vortess->pairs[k][j];
          if ((p->left->x[0] == 0.3833333333333333 &&
              p->left->x[1] == 0.68333333333333335) ||
              (p->right->x[0] == 0.3833333333333333 &&
                  p->right->x[1] == 0.68333333333333335)) {
            n_pairs++;
          }
        }
      }
      part_left->x[0] = part_left->x[0];
    }
    /* check if right particle in cj */
    if (pair->right_cell != cj) {
      continue;
    }
    if (part_left->force.active == 1 && part_right->force.active == 1) {
      hydro_shadowfax_flux_exchange(part_left, part_right, pair->midpoint,
                                    pair->surface_area, shift, 1);
    } else if (part_left->force.active) {
      hydro_shadowfax_flux_exchange(part_left, part_right, pair->midpoint,
                                    pair->surface_area, shift, 0);
    } else if (part_right->force.active) {
      double inverse_shift[3];
      for (int k = 0; k < 3; k++) {
        inverse_shift[k] = -shift[k];
      }
      hydro_shadowfax_flux_exchange(part_right, part_left, pair->midpoint,
                                    pair->surface_area, inverse_shift, 0);
    }
  } /* loop over voronoi faces between ci and cj */
}

void cell_shadowfax_do_pair2_force_recursive(const struct engine *e,
                                             struct cell *restrict ci,
                                             struct cell *restrict cj, int sid,
                                             const double *shift);

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair_subset_density(const struct engine *e,
                                      struct cell *restrict ci,
                                      struct part *restrict parts_i,
                                      const int *restrict ind, int count,
                                      struct cell *restrict cj, const int sid,
                                      const int flipped, const double *shift) {

  if (ci == cj) error("Interacting cell with itself!");

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Pick-out the sorted lists. */
  const struct sort_entry *sort_j = cell_get_hydro_sorts(cj, sid);
  const float dxj = cj->hydro.dx_max_sort;

  /* Parts are on the left? */
  if (!flipped) {

    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[ind[pid]];
      const double pix = pi->x[0] - (shift[0]);
      const double piy = pi->x[1] - (shift[1]);
      const double piz = pi->x[2] - (shift[2]);
      const float hi = pi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = hi * kernel_gamma + dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < hig2) {
          delaunay_add_new_vertex(&ci->hydro.deltess, pj->x[0] + shift[0],
                                  pj->x[1] + shift[1], 26 - sid, cj, pj);
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }

  /* Parts are on the right. */
  else {

    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[ind[pid]];
      const double pix = pi->x[0] - (shift[0]);
      const double piy = pi->x[1] - (shift[1]);
      const double piz = pi->x[2] - (shift[2]);
      const float hi = pi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = -hi * kernel_gamma - dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < hig2) {
          delaunay_add_new_vertex(&ci->hydro.deltess, pj->x[0] + shift[0],
                                  pj->x[1] + shift[1], sid, cj, pj);
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }
}

void cell_shadowfax_do_pair_subset_density_recursive(
    const struct engine *e, struct cell *restrict ci,
    struct part *restrict parts_i, const int *restrict ind, int count,
    struct cell *restrict cj, int sid, int flipped, const double *shift);

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_self1_density(const struct engine *e,
                                struct cell *restrict c) {

  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;

  /* Loop over the parts in c. */
  for (int i = 0; i < count; i++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict p = &parts[i];
    delaunay_add_local_vertex(&c->hydro.deltess, i, p->x[0], p->x[1], p);
  }
}

void cell_shadowfax_do_self1_density_recursive(const struct engine *e,
                                               struct cell *restrict c);

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_self2_density(const struct engine *e,
                                struct cell *restrict c) {
  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_self1_gradient(const struct engine *e,
                                 struct cell *restrict c) {
  double shift[3] = {0., 0., 0.};

  struct voronoi *vortess = &c->hydro.vortess;
  for (int i = 0; i < vortess->pair_index[13]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[13][i];
    struct part *part_left = pair->left;
    struct part *part_right = pair->right;
    if (part_left->force.active == 1 && part_right->force.active == 1) {
      hydro_shadowfax_gradients_collect(part_left, part_right, pair->midpoint,
                                        pair->surface_area, shift);
    } else if (part_left->force.active) {
      hydro_shadowfax_gradients_nonsym_collect(
          part_left, part_right, pair->midpoint, pair->surface_area, shift);
    } else if (part_right->force.active) {
      double inverse_shift[3];
      for (int k = 0; k < 3; k++) {
        inverse_shift[k] = -shift[k];
      }
      hydro_shadowfax_gradients_nonsym_collect(
          part_right, part_left, pair->midpoint, pair->surface_area,
          inverse_shift);
    }
  }
}

void cell_shadowfax_do_self1_gradient_recursive(const struct engine *e,
                                                struct cell *restrict c);

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_self2_gradient(const struct engine *e,
                                 struct cell *restrict c) {
  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_do_self1_force(
    const struct engine *e, struct cell *restrict c) {
  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_do_self2_force(
    const struct engine *e, struct cell *restrict c) {

  double shift[3] = {0., 0., 0.};

  struct voronoi *vortess = &c->hydro.vortess;
  for (int i = 0; i < vortess->pair_index[13]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[13][i];
    struct part *part_left = pair->left;
    struct part *part_right = pair->right;

    if ((part_left->x[0] == 0.3833333333333333 &&
         part_left->x[1] == 0.68333333333333335) ||
        (part_right->x[0] == 0.3833333333333333 &&
         part_right->x[1] == 0.68333333333333335)) {
      int n_pairs = 0;
      for (int k = 0; k < 27; k++) {
        for (int j = 0; j < vortess->pair_index[k]; j++) {
          struct voronoi_pair *p = &vortess->pairs[k][j];
          if ((p->left->x[0] == 0.3833333333333333 &&
               p->left->x[1] == 0.68333333333333335) ||
              (p->right->x[0] == 0.3833333333333333 &&
               p->right->x[1] == 0.68333333333333335)) {
            n_pairs++;
          }
        }
      }
      part_left->x[0] = part_left->x[0];
    }

    if (part_left->force.active == 1 && part_right->force.active == 1) {
      hydro_shadowfax_flux_exchange(part_left, part_right, pair->midpoint,
                                    pair->surface_area, shift, 1);
    } else if (part_left->force.active) {
      hydro_shadowfax_flux_exchange(part_left, part_right, pair->midpoint,
                                    pair->surface_area, shift, 0);
    } else if (part_right->force.active) {
      double inverse_shift[3];
      for (int k = 0; k < 3; k++) {
        inverse_shift[k] = -shift[k];
      }
      hydro_shadowfax_flux_exchange(part_right, part_left, pair->midpoint,
                                    pair->surface_area, inverse_shift, 0);
    }
  }
}

void cell_shadowfax_do_self2_force_recursive(const struct engine *e,
                                             struct cell *restrict c);

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_self_subset_density(const struct engine *e,
                                      struct cell *restrict c,
                                      struct part *restrict parts,
                                      const int *restrict ind, int count) {
  /* Nothing should happen here at the moment */
}

void cell_shadowfax_do_self_subset_density_recursive(const struct engine *e,
                                                     struct cell *c,
                                                     struct part *parts,
                                                     const int *ind, int count);

__attribute__((always_inline)) INLINE static void cell_shadowfax_end_density(
    struct cell *restrict c) {
  voronoi_init(&c->hydro.vortess, &c->hydro.deltess, c->hydro.parts);

  struct part *p;
  for (int i = 0; i < c->hydro.vortess.number_of_cells; i++) {
    p = &c->hydro.parts[i];
    p->density.wcount = 1.0f;
    p->voronoi.volume = c->hydro.vortess.cells[i].volume;
#ifdef VORONOI_STORE_CELL_STATS
    p->voronoi.nface = c->hydro.vortess.cells[i].nface;
#endif
    hydro_gradients_init(p);
    hydro_shadowfax_convert_conserved_to_primitive(p);
  }
}

void cell_shadowfax_end_density_recursive(struct cell *restrict c);

/* Naive functions */
__attribute__((always_inline)) INLINE static void cell_shadowfax_do_pair_naive(
    const struct engine *e, struct cell *restrict ci, struct cell *restrict cj,
    int sid, const double *shift) {

  if (sid < 0) {
    error("Doing a naive interaction!");
  }

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count_i; pid++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict pi = &parts_i[pid];
    delaunay_add_new_vertex(&cj->hydro.deltess, pi->x[0] - shift[0],
                            pi->x[1] - shift[1], sid, ci, pi);
  }

  /* Loop over the parts in cj. */
  for (int pjd = 0; pjd < count_j; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts_j[pjd];
    delaunay_add_new_vertex(&ci->hydro.deltess, pj->x[0] + shift[0],
                            pj->x[1] + shift[1], 26 - sid, cj, pj);
  }
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_do_pair1_naive(
    const struct engine *e, struct cell *restrict ci, struct cell *restrict cj,
    int sid, const double *shift) {

  if (ci == cj) error("Interacting cell with itself!");

  cell_shadowfax_do_pair_naive(e, ci, cj, sid, shift);
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_do_pair2_naive(
    const struct engine *e, struct cell *restrict ci, struct cell *restrict cj,
    int sid, const double *shift) {

  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair_subset_naive(const struct engine *e,
                                    struct cell *restrict ci,
                                    struct cell *restrict cj, int sid,
                                    const double *shift) {

  if (ci == cj) error("Interacting cell with itself!");

  cell_shadowfax_do_pair_naive(e, ci, cj, sid, shift);
}

void cell_shadowfax_write_tesselations(const struct cell *c, FILE *dfile,
                                       FILE *vfile, size_t *offset);

#endif /* SWIFT_CELL_SHADOWFAX_H */
