/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_DEBUGGING_GEAR_H
#define SWIFT_RT_DEBUGGING_GEAR_H

#ifdef SWIFT_RT_DEBUG_CHECKS

#include "rt_properties.h"

/**
 * @file src/rt/GEAR/rt_debugging.h
 * @brief Main header file for the GEAR radiative transfer scheme
 * extra debugging functions.
 */

/**
 * @brief Debugging checks loop over all star particles after each time step
 */
static void rt_debugging_end_of_step_stars_mapper(void *restrict map_data,
                                                  int scount,
                                                  void *restrict extra_data) {

  struct spart *restrict sparts = (struct spart *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  int emission_sum = 0;
  unsigned long long emission_sum_tot = 0;

  for (int k = 0; k < scount; k++) {
    struct spart *restrict sp = &sparts[k];
    emission_sum += sp->rt_data.debug_iact_hydro_inject;
    emission_sum_tot += sp->rt_data.debug_radiation_emitted_tot;

    for (int g = 0; g < RT_NGROUPS; g++) {
      /* also check now that we actually injected the correct
       * amount of energy
       * sp->rt_data.emission_this_step: energy we should distribute
       *                                 this step
       * sp->rt_data.debug_injected_energy: energy we actually did
       *                                    distribute this step    */
      if (sp->rt_data.debug_injected_energy[g] != 0.f) {
        float diff = 1.f - sp->rt_data.emission_this_step[g] /
                               sp->rt_data.debug_injected_energy[g];

        if (fabs(diff) > 1e-3) {
          message(
              "Incorrect injection ID %lld: "
              "group %d expected %.3g got %.3g diff %.3g",
              sp->id, g, sp->rt_data.emission_this_step[g],
              sp->rt_data.debug_injected_energy[g], diff);
        }
      }
    }

    /* reset stuff such that it is set even if the star particle
     * isn't active in the next step */
    sp->rt_data.debug_iact_hydro_inject = 0;
    for (int g = 0; g < RT_NGROUPS; g++) {
      sp->rt_data.debug_injected_energy[g] = 0.f;
    }
  }
  atomic_add(&e->rt_props->debug_radiation_emitted_this_step, emission_sum);
  atomic_add(&e->rt_props->debug_radiation_emitted_tot, emission_sum_tot);
}

/**
 * @brief Debugging checks loop over all hydro particles after each time step
 */
static void rt_debugging_end_of_step_hydro_mapper(void *restrict map_data,
                                                  int count,
                                                  void *restrict extra_data) {

  struct part *restrict parts = (struct part *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  int absorption_sum = 0;
  unsigned long long absorption_sum_tot = 0;
  for (int k = 0; k < count; k++) {
    struct part *restrict p = &parts[k];
    absorption_sum += p->rt_data.debug_iact_stars_inject;
    absorption_sum_tot += p->rt_data.debug_radiation_absorbed_tot;
    /* Reset all values here in case particles won't be active next step */
    p->rt_data.debug_iact_stars_inject = 0;
  }
  atomic_add(&e->rt_props->debug_radiation_absorbed_this_step, absorption_sum);
  atomic_add(&e->rt_props->debug_radiation_absorbed_tot, absorption_sum_tot);
}

/**
 * @brief At the end of each time step, loop over both hydro and star
 * particles and do whatever checks for this particular time step you
 * want done.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_checks_end_of_step(struct engine *e, int verbose) {

  const ticks tic = getticks();
  struct space *s = e->s;
  if (!(e->policy & engine_policy_rt)) return;

  /* reset values before the particle loops */
  e->rt_props->debug_radiation_emitted_this_step = 0ULL;
  e->rt_props->debug_radiation_absorbed_this_step = 0ULL;

  /* hydro particle loop */
  if (s->nr_parts > 0)
    threadpool_map(&e->threadpool, rt_debugging_end_of_step_hydro_mapper,
                   s->parts, s->nr_parts, sizeof(struct part),
                   threadpool_auto_chunk_size, /*extra_data=*/e);

  /* star particle loop */
  if (s->nr_sparts > 0)
    threadpool_map(&e->threadpool, rt_debugging_end_of_step_stars_mapper,
                   s->sparts, s->nr_sparts, sizeof(struct spart),
                   threadpool_auto_chunk_size, /*extra_data=*/e);

  /* Have we accidentally invented or deleted some radiation somewhere? */
  if ((e->rt_props->debug_radiation_emitted_this_step !=
       e->rt_props->debug_radiation_absorbed_this_step) ||
      (e->rt_props->debug_radiation_emitted_tot !=
       e->rt_props->debug_radiation_absorbed_tot))
    error(
        "Emitted and absorbed radiation vary.\n"
        "  This step: star emission %12d; gas absorption %12d\n"
        "Since start: star emission %12lld; gas absorption %12lld",
        e->rt_props->debug_radiation_emitted_this_step,
        e->rt_props->debug_radiation_absorbed_this_step,
        e->rt_props->debug_radiation_emitted_tot,
        e->rt_props->debug_radiation_absorbed_tot);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief This function is intended for debugging purposes only. It is called
 * during the self injection tasks, (regardless whether the particle actually
 * has neighbours to interact with) and intended to mark star or gas particles
 * to have been called during the step so further checks can be performed
 * further down the task system.
 *
 * @param p Hydro particle.
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_check_injection_part(struct part *restrict p,
                                  struct rt_props *props) {

  /* skip this for GEAR */
  /* if (props->do_all_parts_have_stars_checks) */
  /*   p->rt_data.debug_injection_check += 1; */
}

/**
 * @brief This function is intended for debugging purposes only. It is called
 * during the self injection tasks, (regardless whether the particle actually
 * has neighbours to interact with) and intended to mark star or gas particles
 * to have been called during the step so further checks can be performed
 * further down the task system.
 *
 * @param s Star particle.
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_check_injection_spart(struct spart *restrict s,
                                   struct rt_props *props) {

  /* skip this for GEAR */
  /* if (props->do_all_parts_have_stars_checks) */
  /*   s->rt_data.debug_injection_check += 1; */
}

/**
 * @brief Mark that a particle has been called during the gradient calls
 *
 * @param p Particle
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_count_gradient_call(struct part *restrict p) {

  p->rt_data.debug_calls_iact_gradient += 1;
}

/**
 * @brief Mark that a particle has been called during the transport calls
 *
 * @param p Particle
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_count_transport_call(struct part *restrict p) {

  p->rt_data.debug_calls_iact_transport += 1;
}
#endif /* SWIFT_RT_DEBUG_CHECKS */
#endif /* SWIFT_RT_DEBUGGING_GEAR_H */
