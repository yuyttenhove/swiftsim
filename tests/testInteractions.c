/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "swift.h"
#include "cache.h"

/* Other schemes need to be added here if they are not vectorized, otherwise
 * this test will simply not compile. */

#if defined(WITH_VECTORIZATION) && (defined(GADGET2_SPH) || defined(MINIMAL_SPH))

#define array_align sizeof(float) * VEC_SIZE
#define ACC_THRESHOLD 1e-5

#ifndef IACT
#define IACT runner_iact_nonsym_density
#define IACT_VEC runner_iact_nonsym_1_vec_density
#define IACT_NAME "test_nonsym_density"
#define NUM_VEC_PROC_INT 1
#endif

/**
 * @brief Constructs an array of particles in a valid state prior to
 * a IACT_NONSYM and IACT_NONSYM_VEC call.
 *
 * @param count No. of particles to create
 * @param offset The position of the particle offset from (0,0,0).
 * @param spacing Particle spacing.
 * @param h The smoothing length of the particles in units of the inter-particle
 *separation.
 * @param partId The running counter of IDs.
 */
struct part *make_particles(size_t count, double *offset, double spacing,
                            double h, long long *partId) {

  struct part *particles;
  if (posix_memalign((void **)&particles, part_align,
                     count * sizeof(struct part)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(particles, count * sizeof(struct part));

  /* Construct the particles */
  struct part *p;

  /* Set test particle at centre of unit sphere. */
  p = &particles[0];

  /* Place the test particle at the centre of a unit sphere. */
  p->x[0] = 0.0f;
  p->x[1] = 0.0f;
  p->x[2] = 0.0f;

  p->h = h;
  p->id = ++(*partId);

#if !defined(GIZMO_MFV_SPH) && !defined(SHADOWFAX_SPH)
  p->mass = 1.0f;
#endif

  /* Place rest of particles around the test particle
   * with random position within a unit sphere. */
  for (size_t i = 1; i < count; ++i) {
    p = &particles[i];

    /* Randomise positions within a unit sphere. */
    p->x[0] = random_uniform(-1.0, 1.0);
    p->x[1] = random_uniform(-1.0, 1.0);
    p->x[2] = random_uniform(-1.0, 1.0);

    /* Randomise velocities. */
    p->v[0] = random_uniform(-0.05, 0.05);
    p->v[1] = random_uniform(-0.05, 0.05);
    p->v[2] = random_uniform(-0.05, 0.05);

    p->h = h;
    p->id = ++(*partId);
#if !defined(GIZMO_SPH) && !defined(SHADOWFAX_SPH)
    p->mass = 1.0f;
#endif
  }
  return particles;
}

/**
 * @brief Populates particle properties needed for the force calculation.
 */
void prepare_force(struct part *parts, size_t count) {

#if !defined(GIZMO_MFV_SPH) && !defined(SHADOWFAX_SPH) &&       \
    !defined(MINIMAL_SPH) && !defined(MINIMAL_MULTI_MAT_SPH) && \
    !defined(HOPKINS_PU_SPH)
  struct part *p;
  for (size_t i = 0; i < count; ++i) {
    p = &parts[i];
    p->rho = i + 1;
    p->force.soundspeed = random_uniform(2.0, 3.0);
    p->force.v_sig = 0.0f;
    p->force.h_dt = 0.0f;
#if !defined(MINIMAL_SPH)
    p->force.balsara = random_uniform(0.0, 1.0);
    p->force.P_over_rho2 = i + 1;
#endif
  }
#endif
}

/**
 * @brief Dumps all particle information to a file
 */
void dump_indv_particle_fields(char *fileName, struct part *p) {

  FILE *file = fopen(fileName, "a");

  fprintf(file,
          "%6llu %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f "
          "%8.5f "
          "%8.5f %8.5f %13e %13e %13e %13e %13e %8.5f %8.5f\n",
          p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], p->h,
          hydro_get_comoving_density(p),
#if defined(MINIMAL_SPH) || defined(MINIMAL_MULTI_MAT_SPH) || \
    defined(SHADOWFAX_SPH)
          0.f,
#else
          p->density.div_v,
#endif
          hydro_get_comoving_entropy(p), hydro_get_comoving_internal_energy(p),
          hydro_get_comoving_pressure(p), hydro_get_comoving_soundspeed(p),
          p->a_hydro[0], p->a_hydro[1], p->a_hydro[2], p->force.h_dt,
#if defined(GADGET2_SPH)
          p->force.v_sig, p->entropy_dt, 0.f
#elif defined(DEFAULT_SPH)
          p->force.v_sig, 0.f, p->force.u_dt
#elif defined(MINIMAL_SPH) || defined(HOPKINS_PU_SPH)
          p->force.v_sig, 0.f, p->u_dt
#else
          0.f, 0.f, 0.f
#endif
          );

  fclose(file);
}

/**
 * @brief Creates a header for the output file
 */
void write_header(char *fileName) {

  FILE *file = fopen(fileName, "w");
  /* Write header */
  fprintf(file,
          "# %4s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %13s %13s "
          "%13s %13s %13s %8s %8s\n",
          "ID", "pos_x", "pos_y", "pos_z", "v_x", "v_y", "v_z", "h", "rho",
          "div_v", "S", "u", "P", "c", "a_x", "a_y", "a_z", "h_dt", "v_sig",
          "dS/dt", "du/dt");

  fclose(file);
}

/**
 * @brief Compares the vectorised result against
 * the serial result of the interaction.
 *
 * @param serial_test_part Particle that has been updated serially
 * @param serial_parts Particle array that has been interacted serially
 * @param vec_test_part Particle that has been updated using vectors
 * @param vec_parts Particle array to be interacted using vectors
 * @param count No. of particles that have been interacted
 *
 * @return Non-zero value if difference found, 0 otherwise
 */
int check_results(struct part serial_test_part, struct part *serial_parts,
                  struct part vec_test_part, struct part *vec_parts,
                  int count) {
  int result = 0;
  result += compare_particles(serial_test_part, vec_test_part, ACC_THRESHOLD);

  for (int i = 0; i < count; i++)
    result += compare_particles(serial_parts[i], vec_parts[i], ACC_THRESHOLD);

  return result;
}

/*
 * @brief Calls the serial and vectorised version of the non-symmetrical density
 * interaction.
 *
 * @param test_part Particle that will be updated
 * @param parts Particle array to be interacted
 * @param count No. of particles to be interacted
 * @param serial_inter_func Serial interaction function to be called
 * @param vec_inter_func Vectorised interaction function to be called
 * @param runs No. of times to call interactions
 * @param num_vec_proc No. of vectors to use to process interaction
 *
 */
void test_interactions(struct part test_part, struct part *parts, size_t count,
                       char *filePrefix, int runs, int num_vec_proc) {

  ticks serial_time = 0;
  ticks vec_time = 0;

  const float a = 1.f;
  const float H = 0.f;

  char serial_filename[200] = "";
  char vec_filename[200] = "";

  strcpy(serial_filename, filePrefix);
  strcpy(vec_filename, filePrefix);
  sprintf(serial_filename + strlen(serial_filename), "_serial.dat");
  sprintf(vec_filename + strlen(vec_filename), "_%d_vec.dat", num_vec_proc);

  write_header(serial_filename);
  write_header(vec_filename);

  struct part pi_serial, pi_vec;
  struct part pj_serial[count], pj_vec[count];

  float r2[count] __attribute__((aligned(array_align)));
  float dx[3 * count] __attribute__((aligned(array_align)));

  struct part *piq[count], *pjq[count];
  for (size_t k = 0; k < count; k++) {
    piq[k] = NULL;
    pjq[k] = NULL;
  }

  /* Call serial interaction a set number of times. */
  for (int r = 0; r < runs; r++) {
    /* Reset particle to initial setup */
    pi_serial = test_part;
    for (size_t i = 0; i < count; i++) pj_serial[i] = parts[i];

    /* Perform serial interaction */
    for (size_t i = 0; i < count; i++) {
      /* Compute the pairwise distance. */
      r2[i] = 0.0f;
      for (int k = 0; k < 3; k++) {
        int ind = (3 * i) + k;
        dx[ind] = pi_serial.x[k] - pj_serial[i].x[k];
        r2[i] += dx[ind] * dx[ind];
      }
    }

    const ticks tic = getticks();
/* Perform serial interaction */
#ifdef __ICC
#pragma novector
#endif
    for (size_t i = 0; i < count; i++) {
      IACT(r2[i], &(dx[3 * i]), pi_serial.h, pj_serial[i].h, &pi_serial,
           &pj_serial[i], a, H);
    }
    serial_time += getticks() - tic;
  }

  /* Dump result of serial interaction. */
  dump_indv_particle_fields(serial_filename, &pi_serial);
  for (size_t i = 0; i < count; i++)
    dump_indv_particle_fields(serial_filename, &pj_serial[i]);

  /* Call vector interaction a set number of times. */
  for (int r = 0; r < runs; r++) {
    /* Reset particle to initial setup */
    pi_vec = test_part;
    for (size_t i = 0; i < count; i++) pj_vec[i] = parts[i];

    struct c2_cache int_cache;
    struct cache cell_cache;
    cell_cache.count = 0;
    cache_init(&cell_cache, count);

    /* Setup arrays for vector interaction. */
    for (size_t i = 0; i < count; i++) {
      /* Compute the pairwise distance. */
      float my_r2 = 0.0f;
      float my_dx[3];
      for (int k = 0; k < 3; k++) {
        my_dx[k] = pi_vec.x[k] - pj_vec[i].x[k];
        my_r2 += my_dx[k] * my_dx[k];
      }

      int_cache.r2q[i] = my_r2;
      int_cache.dxq[i] = my_dx[0];
      int_cache.dyq[i] = my_dx[1];
      int_cache.dzq[i] = my_dx[2];

      piq[i] = &pi_vec;
      pjq[i] = &pj_vec[i];

      int_cache.mq[i] = pj_vec[i].mass;
      cell_cache.m[i] = pj_vec[i].mass;
#ifdef GADGET2_SPH
      int_cache.vxq[i] = pj_vec[i].v[0];
      int_cache.vyq[i] = pj_vec[i].v[1];
      int_cache.vzq[i] = pj_vec[i].v[2];
      cell_cache.vx[i] = pj_vec[i].v[0];
      cell_cache.vy[i] = pj_vec[i].v[1];
      cell_cache.vz[i] = pj_vec[i].v[2];
#endif
    }

    struct input_params_density params;
    populate_input_params_density(piq[0], &params);

    /* Reset cumulative sums of update vectors. */
    struct update_cache_density sum_cache;
    struct cache_props update_props[NUM_OF_DENSITY_UPDATE_CACHE_FIELDS];
    cache_read_particle_update_fields_density(piq[0], update_props, &sum_cache);
    update_cache_init(NUM_OF_DENSITY_UPDATE_CACHE_FIELDS, update_props);

    mask_t mask, mask2;
    vec_init_mask_true(mask);
    vec_init_mask_true(mask2);

    const ticks vec_tic = getticks();

    for (size_t i = 0; i < count; i += num_vec_proc * VEC_SIZE) {

      /* Interleave two vectors for interaction. */
      if (num_vec_proc == 2) {
        runner_iact_nonsym_2_vec_density(
            &int_cache, i, &params, &sum_cache, mask, mask2, 0);
      } else { /* Only use one vector for interaction. */

        vector my_r2, my_dx, my_dy, my_dz;
        my_r2.v = vec_load(&(int_cache.r2q[i]));
        my_dx.v = vec_load(&(int_cache.dxq[i]));
        my_dy.v = vec_load(&(int_cache.dyq[i]));
        my_dz.v = vec_load(&(int_cache.dzq[i]));

        runner_iact_nonsym_1_vec_density(
            &my_r2, &my_dx, &my_dy, &my_dz, &params, &cell_cache, i, &sum_cache, mask);
      }
    }

    /* Perform horizontal adds on vector sums and store result in particle pi.*/
    update_particle(update_props, 0, NUM_OF_DENSITY_UPDATE_CACHE_FIELDS);

    vec_time += getticks() - vec_tic;
  }

  /* Dump result of serial interaction. */
  dump_indv_particle_fields(vec_filename, piq[0]);
  for (size_t i = 0; i < count; i++)
    dump_indv_particle_fields(vec_filename, pjq[i]);

#ifdef GADGET2_SPH
  /* Check serial results against the vectorised results. */
  if (check_results(pi_serial, pj_serial, pi_vec, pj_vec, count))
    message("Differences found...");
#endif

  message("The serial interactions took     : %15lli ticks.",
          serial_time / runs);
  message("The vectorised interactions took : %15lli ticks.", vec_time / runs);
  message("Speed up: %15fx.", (double)(serial_time) / vec_time);
}

/*
 * @brief Calls the serial and vectorised version of the non-symmetrical force
 * interaction.
 *
 * @param test_part Particle that will be updated
 * @param parts Particle array to be interacted
 * @param count No. of particles to be interacted
 * @param serial_inter_func Serial interaction function to be called
 * @param vec_inter_func Vectorised interaction function to be called
 * @param runs No. of times to call interactions
 *
 */
void test_force_interactions(struct part test_part, struct part *parts,
                             size_t count, char *filePrefix, int runs,
                             int num_vec_proc) {

  ticks serial_time = 0;
  ticks vec_time = 0;

  FILE *file;
  char serial_filename[200] = "";
  char vec_filename[200] = "";

  const float a = 1.f;
  const float H = 0.f;

  strcpy(serial_filename, filePrefix);
  strcpy(vec_filename, filePrefix);
  sprintf(serial_filename + strlen(serial_filename), "_serial.dat");
  sprintf(vec_filename + strlen(vec_filename), "_%d_vec.dat", num_vec_proc);

  write_header(serial_filename);
  write_header(vec_filename);

  struct part pi_serial, pi_vec;
  struct part pj_serial[count], pj_vec[count];

  float r2[count] __attribute__((aligned(array_align)));
  float dx[3 * count] __attribute__((aligned(array_align)));

  struct part *piq[count], *pjq[count];
  for (size_t k = 0; k < count; k++) {
    piq[k] = NULL;
    pjq[k] = NULL;
  }

  float r2q[count] __attribute__((aligned(array_align)));
  float dxq[count] __attribute__((aligned(array_align)));
  float dyq[count] __attribute__((aligned(array_align)));
  float dzq[count] __attribute__((aligned(array_align)));
  float hj_invq[count] __attribute__((aligned(array_align)));
  
  /* Call serial interaction a set number of times. */
  for (int r = 0; r < runs; r++) {
    /* Reset particle to initial setup */
    pi_serial = test_part;
    for (size_t i = 0; i < count; i++) pj_serial[i] = parts[i];

    /* Only dump data on first run. */
    if (r == 0) {
      /* Dump state of particles before serial interaction. */
      dump_indv_particle_fields(serial_filename, &pi_serial);
      for (size_t i = 0; i < count; i++)
        dump_indv_particle_fields(serial_filename, &pj_serial[i]);
    }

    /* Perform serial interaction */
    for (size_t i = 0; i < count; i++) {
      /* Compute the pairwise distance. */
      r2[i] = 0.0f;
      for (int k = 0; k < 3; k++) {
        int ind = (3 * i) + k;
        dx[ind] = pi_serial.x[k] - pj_serial[i].x[k];
        r2[i] += dx[ind] * dx[ind];
      }
    }

    const ticks tic = getticks();
/* Perform serial interaction */
#ifdef __ICC
#pragma novector
#endif
    for (size_t i = 0; i < count; i++) {
      runner_iact_nonsym_force(r2[i], &(dx[3 * i]), pi_serial.h, pj_serial[i].h,
                               &pi_serial, &pj_serial[i], a, H);
    }
    serial_time += getticks() - tic;
  }

  file = fopen(serial_filename, "a");
  fprintf(file, "\n# PARTICLES AFTER INTERACTION:\n");
  fclose(file);

  /* Dump result of serial interaction. */
  dump_indv_particle_fields(serial_filename, &pi_serial);
  for (size_t i = 0; i < count; i++)
    dump_indv_particle_fields(serial_filename, &pj_serial[i]);

  /* Call vector interaction a set number of times. */
  for (int r = 0; r < runs; r++) {
    /* Reset particle to initial setup */
    pi_vec = test_part;
    for (size_t i = 0; i < count; i++) pj_vec[i] = parts[i];

    struct cache cell_cache, cj_cache;
    cell_cache.count = 0;
    cj_cache.count = 0;
    cache_init(&cell_cache, count);
    cache_init(&cj_cache, count);

    /* Setup arrays for vector interaction. */
    for (size_t i = 0; i < count; i++) {
      /* Compute the pairwise distance. */
      float my_r2 = 0.0f;
      float my_dx[3];
      for (int k = 0; k < 3; k++) {
        my_dx[k] = pi_vec.x[k] - pj_vec[i].x[k];
        my_r2 += my_dx[k] * my_dx[k];
      }

      piq[i] = &pi_vec;
      pjq[i] = &pj_vec[i];

      r2q[i] = my_r2;
      dxq[i] = my_dx[0];
      dyq[i] = my_dx[1];
      dzq[i] = my_dx[2];

      cell_cache.h[i] = pi_vec.h;
      cell_cache.vx[i] = pi_vec.v[0];
      cell_cache.vy[i] = pi_vec.v[1];
      cell_cache.vz[i] = pi_vec.v[2];
      cell_cache.rho[i] = pi_vec.rho;
      cell_cache.grad_h[i] = pi_vec.force.f;
      cell_cache.soundspeed[i] = pi_vec.force.soundspeed;
#ifdef GADGET2_SPH
      cell_cache.pOrho2[i] = pi_vec.force.P_over_rho2;
      cell_cache.balsara[i] = pi_vec.force.balsara;
#endif

      hj_invq[i] = 1.f / pj_vec[i].h;
      cj_cache.m[i] = pj_vec[i].mass;
      cj_cache.vx[i] = pj_vec[i].v[0];
      cj_cache.vy[i] = pj_vec[i].v[1];
      cj_cache.vz[i] = pj_vec[i].v[2];
      cj_cache.rho[i] = pj_vec[i].rho;
      cj_cache.grad_h[i] = pj_vec[i].force.f;
      cj_cache.soundspeed[i] = pj_vec[i].force.soundspeed;
#ifdef GADGET2_SPH
      cj_cache.pOrho2[i] = pj_vec[i].force.P_over_rho2;
      cj_cache.balsara[i] = pj_vec[i].force.balsara;
#endif
    
    }

    /* Only dump data on first run. */
    if (r == 0) {
      /* Dump state of particles before vector interaction. */
      dump_indv_particle_fields(vec_filename, piq[0]);
      for (size_t i = 0; i < count; i++)
        dump_indv_particle_fields(vec_filename, pjq[i]);
    }

    struct input_params_force params;
    populate_input_params_force_cache(&cell_cache, 0, &params);
    
    /* Reset cumulative sums of update vectors. */
    struct update_cache_force sum_cache;
    struct cache_props update_props[NUM_OF_FORCE_UPDATE_CACHE_FIELDS];
    cache_read_particle_update_fields_force(piq[0], update_props, &sum_cache);
    update_cache_init(NUM_OF_FORCE_UPDATE_CACHE_FIELDS, update_props);

    /* Perform vector interaction. */
    mask_t mask;
    vec_init_mask_true(mask);

    const ticks vec_tic = getticks();

    for (size_t i = 0; i < count; i += num_vec_proc * VEC_SIZE) {

        vector my_r2, my_dx, my_dy, my_dz, hj, hj_inv;
        my_r2.v = vec_load(&(r2q[i]));
        my_dx.v = vec_load(&(dxq[i]));
        my_dy.v = vec_load(&(dyq[i]));
        my_dz.v = vec_load(&(dzq[i]));
        hj.v = vec_load(&hj_invq[i]);
        hj_inv = vec_reciprocal(hj);

        runner_iact_nonsym_1_vec_force(
            &my_r2, &my_dx, &my_dy, &my_dz, &params, &cj_cache, i, hj_inv, a, H, &sum_cache, mask);
    }

    /* Perform horizontal adds on vector sums and store result in particle pi.*/
    update_particle(update_props, 0, NUM_OF_FORCE_UPDATE_CACHE_FIELDS);

    vec_time += getticks() - vec_tic;
  }

  file = fopen(vec_filename, "a");
  fprintf(file, "\n# PARTICLES AFTER INTERACTION:\n");
  fclose(file);

  /* Dump result of serial interaction. */
  dump_indv_particle_fields(vec_filename, piq[0]);
  for (size_t i = 0; i < count; i++)
    dump_indv_particle_fields(vec_filename, pjq[i]);

#ifdef GADGET2_SPH
  /* Check serial results against the vectorised results. */
  if (check_results(pi_serial, pj_serial, pi_vec, pj_vec, count))
    message("Differences found...");
#endif

  message("The serial interactions took     : %15lli ticks.",
          serial_time / runs);
  message("The vectorised interactions took : %15lli ticks.", vec_time / runs);
  message("Speed up: %15fx.", (double)(serial_time) / vec_time);
}

/* And go... */
int main(int argc, char *argv[]) {
  size_t runs = 10000;
  double h = 1.0, spacing = 0.5;
  double offset[3] = {0.0, 0.0, 0.0};
  size_t count = 112;

  /* Get some randomness going */
  srand(0);

  char c;
  while ((c = getopt(argc, argv, "h:s:n:r:")) != -1) {
    switch (c) {
      case 'h':
        sscanf(optarg, "%lf", &h);
        break;
      case 's':
        sscanf(optarg, "%lf", &spacing);
        break;
      case 'n':
        sscanf(optarg, "%zu", &count);
        break;
      case 'r':
        sscanf(optarg, "%zu", &runs);
        break;
      case '?':
        error("Unknown option.");
        break;
    }
  }

  if (h < 0 || spacing < 0) {
    printf(
        "\nUsage: %s [OPTIONS...]\n"
        "\nGenerates a particle array with equal particle separation."
        "\nThese are then interacted using runner_iact_density and "
        "runner_iact_vec_density."
        "\n\nOptions:"
        "\n-h DISTANCE=1.2348 - Smoothing length in units of <x>"
        "\n-s SPACING=0.5     - Spacing between particles"
        "\n-n NUMBER=9        - No. of particles",
        argv[0]);
    exit(1);
  }

  /* Correct count so that VEC_SIZE of particles interact with the test
   * particle. */
  count = count - (count % VEC_SIZE) + 1;

  /* Build the infrastructure */
  static long long partId = 0;
  struct part test_particle;
  struct part *particles = make_particles(count, offset, spacing, h, &partId);

  test_particle = particles[0];
  /* Call the non-sym density test. */
  message("Testing %s interaction...", IACT_NAME);
  test_interactions(test_particle, &particles[1], count - 1, IACT_NAME, runs,
                    1);
  test_interactions(test_particle, &particles[1], count - 1, IACT_NAME, runs,
                    2);

  prepare_force(particles, count);

  test_force_interactions(test_particle, &particles[1], count - 1,
                          "test_nonsym_force", runs, 1);

  return 0;
}

#else

int main(int argc, char *argv[]) { return 1; }

#endif
