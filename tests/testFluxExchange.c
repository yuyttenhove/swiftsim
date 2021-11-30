//
// Created by yuyttenh on 07/07/2021.
//

/* Config parameters. */
#include "../config.h"
/* Local headers. */
#include "swift.h"
#include "shadowfax/cell_shadowfax.h"
/* Some standard headers. */
#include <fenv.h>

#ifdef SHADOWFAX_NEW_SPH

#define N_GENERATORS_1D 2
#define CELL_SIZE_1D 1.
#define PERT 0.
#define H_PERT 0.
#define NODE_ID 0

struct cell *make_cell(double *offset, double h,
                       double density, long long *partId) {
  const size_t count = N_GENERATORS_1D * N_GENERATORS_1D * N_GENERATORS_1D;
  const double volume = CELL_SIZE_1D * CELL_SIZE_1D * CELL_SIZE_1D;
  float h_max = 0.f;
  struct cell *cell = (struct cell *)malloc(sizeof(struct cell));
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->hydro.parts, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->hydro.parts, count * sizeof(struct part));

  /* Construct the parts */
  struct part *part = cell->hydro.parts;
  for (size_t x = 0; x < N_GENERATORS_1D; ++x) {
    for (size_t y = 0; y < N_GENERATORS_1D; ++y) {
      for (size_t z = 0; z < N_GENERATORS_1D; ++z) {
        part->x[0] =
            offset[0] +
            CELL_SIZE_1D * (x + 0.5 + random_uniform(-0.5, 0.5) * PERT) / (float)N_GENERATORS_1D;
        part->x[1] =
            offset[1] +
            CELL_SIZE_1D * (y + 0.5 + random_uniform(-0.5, 0.5) * PERT) / (float)N_GENERATORS_1D;
        part->x[2] =
            offset[2] +
            CELL_SIZE_1D * (z + 0.5 + random_uniform(-0.5, 0.5) * PERT) / (float)N_GENERATORS_1D;

        /* Set velocities to zero */
        part->v[0] = 0.f;
        part->v[1] = 0.f;
        part->v[2] = 0.f;

        if (H_PERT)
          part->h = CELL_SIZE_1D * h * random_uniform(1.f, H_PERT) / (float)N_GENERATORS_1D;
        else
          part->h = CELL_SIZE_1D * h / (float)N_GENERATORS_1D;
        h_max = fmaxf(h_max, part->h);
        part->id = ++(*partId);

        part->conserved.mass = density * volume / count;

        part->time_bin = 1;

#ifdef SWIFT_DEBUG_CHECKS
        part->ti_drift = 8;
        part->ti_kick = 8;
#endif
        ++part;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->hydro.h_max = h_max;
  cell->hydro.h_max_active = h_max;
  cell->hydro.count = count;
  cell->hydro.dx_max_part = 0.;
  cell->hydro.dx_max_sort = 0.;
  cell->width[0] = CELL_SIZE_1D;
  cell->width[1] = CELL_SIZE_1D;
  cell->width[2] = CELL_SIZE_1D;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->hydro.ti_old_part = 8;
  cell->hydro.ti_end_min = 8;
  cell->nodeID = NODE_ID;

  shuffle_particles(cell->hydro.parts, cell->hydro.count);

  cell->hydro.sorted = 0;
  cell->hydro.sort = NULL;
  cell->hydro.super = cell;

  cell_malloc_tesselations(cell);

  return cell;
}

void clean_up(struct cell *ci) {
  free(ci->hydro.parts);
  free(ci->hydro.sort);
#ifdef SHADOWFAX_NEW_SPH
  cell_destroy_tessellations(ci);
#endif
  free(ci);
}

/* Just a forward declaration... */
void runner_dopair1_branch_density(struct runner *r, struct cell *ci,
                                   struct cell *cj);
void runner_doself1_branch_density(struct runner *r, struct cell *c);

/* And go... */
int main(int argc, char *argv[]) {

#ifdef HAVE_SETAFFINITY
  engine_pin();
#endif

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  srand(0);

  double h = 2, rho = 1.;

  /* Build the infrastructure */
  struct space space;
  space.periodic = 1;
  space.dim[0] = 3.;
  space.dim[1] = 3.;
  space.dim[2] = 3.;

  struct hydro_props hp;
  hydro_props_init_no_hydro(&hp);
  hp.eta_neighbours = h;
  hp.h_tolerance = 1e0;
  hp.h_max = FLT_MAX;
  hp.max_smoothing_iterations = 1;
  hp.CFL_condition = 0.1;

  struct engine engine;
  engine.s = &space;
  engine.time = 0.1f;
  engine.ti_current = 8;
  engine.max_active_bin = num_time_bins;
  engine.hydro_properties = &hp;
  engine.nodeID = NODE_ID;

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  engine.cosmology = &cosmo;

  struct runner runner;
  runner.e = &engine;

  /* Construct some cells */
  struct cell *cells[27];
  static long long partId = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        double offset[3] = {i * CELL_SIZE_1D, j * CELL_SIZE_1D, k * CELL_SIZE_1D};
        cells[i * 9 + j * 3 + k] =
            make_cell(offset, h, rho, &partId);

        runner_do_drift_part(&runner, cells[i * 9 + j * 3 + k], 0);

        runner_do_hydro_sort(&runner, cells[i * 9 + j * 3 + k], 0x1FFF, 0, 0);
      }
    }
  }

  /* Store the main cell for future use */
  struct cell *main_cell = cells[13];
  int code = 8 == main_cell->hydro.count;

  /* Build delaunay tesselation of main cell. */
  /* Run all the pairs */
  for (int j = 0; j < 27; ++j) {
    if (cells[j] != main_cell) {
      runner_dopair1_branch_density(&runner, main_cell, cells[j]);
    }
  }
  /* And now the self-interaction */
  runner_doself1_branch_density(&runner, main_cell);

  /* Clean things to make the sanitizer happy ... */
  for (int i = 0; i < 27; ++i) clean_up(cells[i]);

  return code;
}

#else

int main(int argc, char *argv[]) {
  return 0;
}

#endif
