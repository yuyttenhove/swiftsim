#ifndef SWIFT_CELL_SHADOWFAX_H
#define SWIFT_CELL_SHADOWFAX_H

#include "active.h"
#include "cell.h"
#include "delaunay.h"
#include "hydro_shadowfax.h"
#include "voronoi.h"

__attribute__((always_inline)) INLINE static void shadowfax_flag_particle_added(
    struct part *restrict p, int sid) {
  p->voronoi.flag |= (1 << sid);
}

__attribute__((always_inline)) INLINE static int shadowfax_particle_was_added(
    const struct part *restrict p, int sid) {
  return (p->voronoi.flag & (1 << sid)) > 0;
}

__attribute__((always_inline)) INLINE static void
cell_malloc_delaunay_tessellation(struct cell *c,
                                  const struct hydro_space *hs) {

  if (c->hydro.shadowfax_enabled == 1) {
    delaunay_reset(&c->hydro.deltess, c->hydro.count);
  } else {
    delaunay_init(&c->hydro.deltess, c->loc, c->width, c->hydro.count,
                  10 * c->hydro.count);
  }

  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;

  for (int pd = 0; pd < count; pd++) {
    /* Get a pointer to the ith particle. */
    struct part *restrict p = &parts[pd];
    p->voronoi.flag = 0;
    p->voronoi.nface = 0;
    p->voronoi.volume = 0;
  }

  c->hydro.shadowfax_enabled = 1;
}

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_self1_density(const struct engine *e,
                                struct cell *restrict c) {

  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;

  /* Loop over the parts in c. */
  for (int pd = 0; pd < count; pd++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict p = &parts[pd];

    if (!shadowfax_particle_was_added(p, 13)) {
      delaunay_add_local_vertex(&c->hydro.deltess, pd, p->x[0], p->x[1]);
      shadowfax_flag_particle_added(p, 13);
    }
  }
}

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_self2_density(const struct engine *e,
                                struct cell *restrict c) {
  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_do_self1_force(
    const struct engine *e, struct cell *restrict c) {
  error("Shouldn't be using this function!");
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_do_self2_force(
    const struct engine *e, struct cell *restrict c) {

  struct part *restrict parts = c->hydro.parts;
  double shift[3] = {0., 0., 0.};

  struct voronoi *vortess = &c->hydro.vortess;
  for (int i = 0; i < vortess->pair_index[13]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[13][i];
    if (parts[pair->left].force.active == 1
        || parts[pair->right].force.active == 1) {
      hydro_shadowfax_flux_exchange(&parts[pair->left], &parts[pair->right],
                                    pair->midpoint, pair->surface_area, shift);
    }
  }
}

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair1_density(const struct engine *e,
                                struct cell *restrict ci,
                                struct cell *restrict cj, int sid,
                                const double *shift) {

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

          if (!shadowfax_particle_was_added(pj, 27 - sid)) {
            delaunay_add_new_vertex(&ci->hydro.deltess, pj->x[0] + shift[0],
                                    pj->x[1] + shift[1], 27 - sid,
                                    sort_j[pjd].i);
            shadowfax_flag_particle_added(pj, 27 - sid);
          }
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

          if (!shadowfax_particle_was_added(pi, sid)) {
            delaunay_add_new_vertex(&cj->hydro.deltess, pi->x[0] - shift[0],
                                    pi->x[1] - shift[1], sid, sort_i[pid].i);
            shadowfax_flag_particle_added(pi, sid);
          }
        }
      } /* loop over the parts in ci. */
    }   /* loop over the parts in cj. */
  }     /* Cell cj is active */
}

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair2_density(const struct engine *e,
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

  if (ci == cj) error("Interacting cell with itself!");

  /* anything to do here? */
  if (!(cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e))) return;

  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  struct voronoi *vortess = &ci->hydro.vortess;
  /* loop over voronoi faces between ci and cj */
  for (int i = 0; i < vortess->pair_index[27 - sid]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[27 - sid][i];
    /* at least one of the parts active? TODO check this later, for the moment, always continue! */
    if (1 || parts_i[pair->left].force.active == 1
        || parts_j[pair->right].force.active == 1) {
      hydro_shadowfax_flux_exchange(&parts_i[pair->left], &parts_j[pair->right],
                                    pair->midpoint, pair->surface_area, shift);
    } /* at least one of the parts active? */
  } /* loop over voronoi faces between ci and cj */
}

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
          if (!shadowfax_particle_was_added(pj, 27 - sid)) {
            delaunay_add_new_vertex(&ci->hydro.deltess, pj->x[0] + shift[0],
                                    pj->x[1] + shift[1], 27 - sid,
                                    sort_j[pjd].i);
            shadowfax_flag_particle_added(pj, 27 - sid);
          }
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
          if (!shadowfax_particle_was_added(pj, sid)) {
            delaunay_add_new_vertex(&ci->hydro.deltess, pj->x[0] + shift[0],
                                    pj->x[1] + shift[1], sid, sort_j[pjd].i);
            shadowfax_flag_particle_added(pj, sid);
          }
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }
}

__attribute__((always_inline)) INLINE static void cell_shadowfax_end_density(
    struct cell *restrict c) {
  voronoi_init(&c->hydro.vortess, &c->hydro.deltess, c->hydro.parts);

  struct part *p;
  for (int i = 0; i < c->hydro.vortess.number_of_cells; i++) {
    p = &c->hydro.parts[i];
    hydro_gradients_init(p);
    p->density.wcount = 1.0f;
    p->voronoi.volume = c->hydro.vortess.cells[i].volume;
    hydro_shadowfax_convert_conserved_to_primitive(p);
  }
}

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

    if (!shadowfax_particle_was_added(pi, sid)) {
      delaunay_add_new_vertex(&cj->hydro.deltess, pi->x[0] - shift[0],
                              pi->x[1] - shift[1], sid, pid);
      shadowfax_flag_particle_added(pi, sid);
    }
  }

  /* Loop over the parts in cj. */
  for (int pjd = 0; pjd < count_j; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts_j[pjd];

    if (!shadowfax_particle_was_added(pj, 27 - sid)) {
      delaunay_add_new_vertex(&ci->hydro.deltess, pj->x[0] + shift[0],
                              pj->x[1] + shift[1], 27 - sid, pjd);
      shadowfax_flag_particle_added(pj, 27 - sid);
    }
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
  if (ci == cj) error("Interacting cell with itself!");

  cell_shadowfax_do_pair_naive(e, ci, cj, sid, shift);
}

__attribute__((always_inline)) INLINE static void
cell_shadowfax_do_pair_subset_naive(const struct engine *e,
                                    struct cell *restrict ci,
                                    struct cell *restrict cj, int sid,
                                    const double *shift) {

  if (ci == cj) error("Interacting cell with itself!");

  cell_shadowfax_do_pair_naive(e, ci, cj, sid, shift);
}

#endif /* SWIFT_CELL_SHADOWFAX_H */
