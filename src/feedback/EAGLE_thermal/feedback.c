/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* This file's header */
#include "feedback.h"

/* Local includes. */
#include "../EAGLE/enrichment.h"
#include "../EAGLE/yield_tables.h"
#include "hydro_properties.h"
#include "inline.h"
#include "random.h"
#include "timers.h"

/**
 * @brief Return the change in temperature (in internal units) to apply to a
 * gas particle affected by SNe feedback.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
double eagle_feedback_temperature_change(const struct spart* sp,
                                         const struct feedback_props* props) {

  /* In the EAGLE REF model, the change of temperature is constant */
  return props->SNe_deltaT_desired;
}

/**
 * @brief Computes the fraction of the available super-novae energy to
 * inject for a given event.
 *
 * Note that the fraction can be > 1.
 *
 * This allows a choice of different f_th scaling functions.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 * @param ngb_nH_cgs Hydrogen number density of the gas surrounding the star
 * (physical cgs units).
 * @param ngb_Z Metallicity (metal mass fraction) of the gas surrounding the
 * star.
 */
double eagle_feedback_energy_fraction(const struct spart* sp,
                                      const struct feedback_props* props,
                                      const double ngb_nH_cgs,
                                      const double ngb_Z) {

  /* Model parameters */
  const double f_E_max = props->f_E_max;
  const double f_E_min = props->f_E_min;
  const double Z_0 = props->Z_0;
  const double n_0 = props->n_0_cgs;
  const double n_Z = props->n_Z;
  const double n_n = props->n_n;

  /* Metallicity (metal mass fraction) at birth time of the star */
  const double Z_birth =
      chemistry_get_star_total_metal_mass_fraction_for_feedback(sp);

  /* Physical density of the gas at the star's birth time */
  const double rho_birth = sp->birth_density;
  const double n_birth_cgs = rho_birth * props->rho_to_n_cgs;

  /* Choose either the birth properties or current properties */
  const double nH =
      props->use_birth_density_for_f_th ? n_birth_cgs : ngb_nH_cgs;
  const double Z = props->use_birth_Z_for_f_th ? Z_birth : ngb_Z;

  /* Calculate f_E */
  const double Z_term = pow(max(Z, 1e-6) / Z_0, -n_Z);
  const double n_term = pow(nH / n_0, n_n);

  /* Different behaviour for different scaling functions ahead */
  double f_th;
  if (props->SNII_energy_scaling == SNII_scaling_independent) {
    /* Independent scaling of f_th with Z and n. Here, we have a second
     * parameter for the max f_th increase due to density, delta_E_n. */
    const double delta_E_n = props->SNII_delta_E_n;
    f_th = (f_E_max - (f_E_max - f_E_min) / (1. + Z_term)) *
           (delta_E_n - (delta_E_n - 1.) / (1. + n_term));

  } else if (props->SNII_energy_scaling == SNII_scaling_separable) {
    /* Separable scaling between fixed fE_min and fE_max */
    f_th = f_E_max - (f_E_max - f_E_min) / ((1. + Z_term) * (1. + n_term));

  } else if (props->SNII_energy_scaling == SNII_scaling_EAGLE) {
    /* Mixed scaling as described in Schaye et al. (2015) for EAGLE */
    f_th = f_E_max - (f_E_max - f_E_min) / (1. + Z_term * n_term);

  } else {
    error("Invalid SNII energy scaling model!");
    f_th = -1.;
  }

  return f_th;
}

/**
 * @brief Compute the properties of the SNII stochastic feedback energy
 * injection.
 *
 * Only does something if the particle reached the SNII age during this time
 * step.
 *
 * @param sp The star particle.
 * @param star_age Age of star at the beginning of the step in internal units.
 * @param dt Length of time-step in internal units.
 * @param ngb_gas_mass Total un-weighted mass in the star's kernel (internal
 * units)
 * @param num_gas_ngbs Total (integer) number of gas neighbours within the
 * star's kernel.
 * @param ngb_nH_cgs Hydrogen number density of the gas surrounding the star
 * (physical cgs units).
 * @param ngb_Z Metallicity (metal mass fraction) of the gas surrounding the
 * star.
 * @param feedback_props The properties of the feedback model.
 * @param min_dying_mass_Msun Minimal star mass dying this step (in solar
 * masses).
 * @param max_dying_mass_Msun Maximal star mass dying this step (in solar
 * masses).
 */
INLINE static void compute_SNII_feedback(
    struct spart* sp, const double star_age, const double dt,
    const int ngb_gas_N, const float ngb_gas_mass, const double ngb_nH_cgs,
    const double ngb_Z, const struct feedback_props* feedback_props,
    const double min_dying_mass_Msun, const double max_dying_mass_Msun,
    const integertime_t ti_begin) {

  /* Are we sampling the delay function or using a fixed delay? */
  const int SNII_sampled_delay = feedback_props->SNII_sampled_delay;

  /* Time after birth considered for SNII feedback (internal units)
   * when using a fixed delay */
  const double SNII_wind_delay = feedback_props->SNII_wind_delay;

  /* Are we doing feedback this step?
   * Note that since the ages are calculated using an interpolation table we
   * must allow some tolerance here*/
  if ((SNII_sampled_delay) || (star_age <= SNII_wind_delay &&
                               (star_age + 1.001 * dt) > SNII_wind_delay)) {

    /* Make sure a star does not do feedback twice
     * when using a fixed delay! */
    if (!SNII_sampled_delay && sp->f_E != -1.f) {
#ifdef SWIFT_DEBUG_CHECKS
      message("Star has already done feedback! sp->id=%lld age=%e d=%e", sp->id,
              star_age, dt);
#endif
      return;
    }

    /* Properties of the model (all in internal units) */
    const double delta_T =
        eagle_feedback_temperature_change(sp, feedback_props);
    const double E_SNe = feedback_props->E_SNII;
    const double f_E =
        eagle_feedback_energy_fraction(sp, feedback_props, ngb_nH_cgs, ngb_Z);

    /* Number of SNe at this time-step */
    double N_SNe;
    if (SNII_sampled_delay) {
      N_SNe = eagle_feedback_number_of_sampled_SNII(
          sp, feedback_props, min_dying_mass_Msun, max_dying_mass_Msun);
    } else {
      N_SNe = eagle_feedback_number_of_SNII(sp, feedback_props);
    }

    /* Abort if there are no SNe exploding this step */
    if (N_SNe <= 0.) return;

    /* Conversion factor from T to internal energy */
    const double conv_factor = feedback_props->temp_to_u_factor;

    /* Calculate the default heating probability (accounting for round-off) */
    double prob = f_E * E_SNe * N_SNe / (conv_factor * delta_T * ngb_gas_mass);
    prob = max(prob, 0.0);

    /* Calculate the change in internal energy of the gas particles that get
     * heated */
    double delta_u;

    /* Number of SNII events for this stellar particle */
    int number_of_SN_events = 0;

    if (prob <= 1.) {

      /* Normal case */
      delta_u = delta_T * conv_factor;

      for (int i = 0; i < ngb_gas_N; i++) {
        const double rand_thermal = random_unit_interval_part_ID_and_ray_idx(
            sp->id, i, ti_begin, random_number_stellar_feedback_3);
        if (rand_thermal < prob) number_of_SN_events++;
      }

    } else {

      /* Special case: we need to adjust the energy irrespective of the
         desired deltaT to ensure we inject all the available energy. */
      delta_u = f_E * E_SNe * N_SNe / ngb_gas_mass;

      /* Number of SNIa events is equal to the number of Ngbs */
      number_of_SN_events = ngb_gas_N;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (f_E < feedback_props->f_E_min || f_E > feedback_props->f_E_max)
      error("f_E is not in the valid range! f_E=%f sp->id=%lld", f_E, sp->id);
#endif

    /* If we have more heating events than the maximum number of
     * rays (eagle_feedback_number_of_rays), then we cannot
     * distribute all of the heating events (since 1 event = 1 ray), so we need
     * to increase the thermal energy per ray and make the number of events
     * equal to the number of rays */
    if (number_of_SN_events > eagle_SNII_feedback_num_of_rays) {
      const double alpha_thermal =
          (double)number_of_SN_events / (double)eagle_SNII_feedback_num_of_rays;
      delta_u *= alpha_thermal;
      number_of_SN_events = eagle_SNII_feedback_num_of_rays;
    }

    /* Current total f_E for this star */
    double star_f_E = sp->f_E * sp->number_of_SNII_events;

    /* New total */
    star_f_E = (star_f_E + f_E) / (sp->number_of_SNII_events + 1.);

    /* Store all of this in the star for delivery onto the gas and recording */
    sp->f_E = star_f_E;
    sp->number_of_SNII_events++;
    sp->feedback_data.to_distribute.SNII_delta_u = delta_u;
    sp->feedback_data.to_distribute.SNII_num_of_thermal_energy_inj =
        number_of_SN_events;
  }
}

/**
 * @brief calculates stellar mass in spart that died over the timestep, calls
 * functions to calculate feedback due to SNIa, SNII and AGB
 *
 * @param feedback_props feedback_props data structure
 * @param phys_const The physical constants in internal units.
 * @param cosmo The cosmological model.
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 * @param ti_begin The current integer time (for random number hashing).
 */
void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct phys_const* phys_const,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt, const integertime_t ti_begin) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (age < 0.f) error("Negative age for a star.");

  if (sp->count_since_last_enrichment != 0 && engine_current_step > 0)
    error("Computing feedback on a star that should not");
#endif

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_inv = 1. / (phys_const->const_year * 1e9);
  const double dt_Gyr = dt * Gyr_inv;
  const double star_age_Gyr = age * Gyr_inv;

  /* Get the birth mass of the star */
  const double M_init = sp->mass_init;

  /* Get the total metallicity (metal mass fraction) at birth time and impose a
   * minimum */
  const double Z =
      max(chemistry_get_star_total_metal_mass_fraction_for_feedback(sp),
          exp10(log10_min_metallicity));

  /* Get the individual abundances (mass fractions at birth time) */
  const float* const abundances =
      chemistry_get_star_metal_mass_fraction_for_feedback(sp);

  /* Properties collected in the stellar density loop. */
  const int ngb_Number = sp->feedback_data.to_collect.ngb_N;
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;
  const float ngb_gas_Z = sp->feedback_data.to_collect.ngb_Z;
  const float ngb_gas_rho = sp->feedback_data.to_collect.ngb_rho;
  const float ngb_gas_phys_nH_cgs =
      ngb_gas_rho * cosmo->a3_inv * feedback_props->rho_to_n_cgs;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Update the enrichment weights */
  const float enrichment_weight_inv =
      sp->feedback_data.to_collect.enrichment_weight_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_collect.enrichment_weight_inv < 0.)
    error("Negative inverse weight!");
#endif

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.to_distribute.enrichment_weight = enrichment_weight;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.enrichment_weight < 0.)
    error("Negative weight!");
#endif

  /* Calculate mass of stars that has died from the star's birth up to the
   * beginning and end of timestep */
  const double max_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr, Z, feedback_props);
  const double min_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr + dt_Gyr, Z, feedback_props);

#ifdef SWIFT_DEBUG_CHECKS
  /* Sanity check. Worth investigating if necessary as functions for evaluating
   * mass of stars dying might be strictly decreasing.  */
  if (min_dying_mass_Msun > max_dying_mass_Msun)
    error("min dying mass is greater than max dying mass");
#endif

  /* Compute properties of the stochastic SNII feedback model. */
  if (feedback_props->with_SNII_feedback) {
    compute_SNII_feedback(sp, age, dt, ngb_Number, ngb_gas_mass,
                          ngb_gas_phys_nH_cgs, ngb_gas_Z, feedback_props,
                          min_dying_mass_Msun, max_dying_mass_Msun, ti_begin);
  }

  /* Integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_Msun. Return without doing any
   * enrichment. */
  if (min_dying_mass_Msun == max_dying_mass_Msun) return;

  /* Life is better in log */
  const double log10_max_dying_mass_Msun = log10(max_dying_mass_Msun);
  const double log10_min_dying_mass_Msun = log10(min_dying_mass_Msun);

  /* Compute elements, energy and momentum to distribute from the
   *  three channels SNIa, SNII, AGB */
  if (feedback_props->with_SNIa_enrichment) {
    evolve_SNIa(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                feedback_props, star_age_Gyr, dt_Gyr, &sp->feedback_data);
  }
  if (feedback_props->with_SNII_enrichment) {
    evolve_SNII(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                abundances, feedback_props, &sp->feedback_data);
  }
  if (feedback_props->with_AGB_enrichment) {
    evolve_AGB(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
               abundances, feedback_props, &sp->feedback_data);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.mass != 0.f)
    error("Injected mass will be lost");
#endif

  /* Compute the total mass to distribute (H + He  metals) */
  sp->feedback_data.to_distribute.mass =
      sp->feedback_data.to_distribute.total_metal_mass +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_H] +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_He];

  /* Compute energy change due to kinetic energy of ejectas */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass *
      feedback_props->AGB_ejecta_specific_kinetic_energy;

  /* Compute energy change due to kinetic energy of the star */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass * 0.5f *
      (sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2]) *
      cosmo->a2_inv;

  TIMER_TOC(timer_do_star_evol);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void feedback_props_init(struct feedback_props* fp,
                         const struct phys_const* phys_const,
                         const struct unit_system* us,
                         struct swift_params* params,
                         const struct hydro_props* hydro_props,
                         const struct cosmology* cosmo) {

  /* Main operation modes ------------------------------------------------- */

  fp->with_SNII_feedback =
      parser_get_param_int(params, "EAGLEFeedback:use_SNII_feedback");

  fp->with_SNIa_feedback =
      parser_get_param_int(params, "EAGLEFeedback:use_SNIa_feedback");

  fp->with_AGB_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_AGB_enrichment");

  fp->with_SNII_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_SNII_enrichment");

  fp->with_SNIa_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_SNIa_enrichment");

  if (fp->with_SNIa_feedback && !fp->with_SNIa_enrichment) {
    error("Cannot run with SNIa feedback without SNIa enrichment.");
  }

  /* Properties of the IMF model ------------------------------------------ */

  /* Minimal and maximal mass considered */
  fp->imf_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:IMF_max_mass_Msun");
  fp->imf_min_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:IMF_min_mass_Msun");

  /* Check that it makes sense. */
  if (fp->imf_max_mass_msun < fp->imf_min_mass_msun) {
    error("Can't have the max IMF mass smaller than the min IMF mass!");
  }

  fp->log10_imf_max_mass_msun = log10(fp->imf_max_mass_msun);
  fp->log10_imf_min_mass_msun = log10(fp->imf_min_mass_msun);

  /* Properties of the SNII energy feedback model ------------------------- */

  char model[64];
  parser_get_param_string(params, "EAGLEFeedback:SNII_feedback_model", model);
  if (strcmp(model, "Random") == 0)
    fp->feedback_model = SNII_random_ngb_model;
  else if (strcmp(model, "Isotropic") == 0)
    fp->feedback_model = SNII_isotropic_model;
  else if (strcmp(model, "MinimumDistance") == 0)
    fp->feedback_model = SNII_minimum_distance_model;
  else if (strcmp(model, "MinimumDensity") == 0)
    fp->feedback_model = SNII_minimum_density_model;
  else
    error(
        "The SNII feedback model must be either 'Random', 'MinimumDistance', "
        "'MinimumDensity' or 'Isotropic', not %s",
        model);

  /* Are we sampling the SNII lifetimes for feedback or using a fixed delay? */
  fp->SNII_sampled_delay =
      parser_get_param_int(params, "EAGLEFeedback:SNII_sampled_delay");

  if (!fp->SNII_sampled_delay) {

    /* Set the delay time before SNII occur */
    fp->SNII_wind_delay =
        parser_get_param_double(params, "EAGLEFeedback:SNII_wind_delay_Gyr") *
        phys_const->const_year * 1e9;
  }

  /* Read the temperature change to use in stochastic heating */
  fp->SNe_deltaT_desired =
      parser_get_param_float(params, "EAGLEFeedback:SNII_delta_T_K");
  fp->SNe_deltaT_desired /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Energy released by supernova type II */
  fp->E_SNII_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_erg");
  fp->E_SNII =
      fp->E_SNII_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Stellar mass limits for SNII feedback */
  const double SNII_min_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_min_mass_Msun");
  const double SNII_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_max_mass_Msun");

  /* Check that it makes sense. */
  if (SNII_max_mass_msun < SNII_min_mass_msun) {
    error("Can't have the max SNII mass smaller than the min SNII mass!");
  }

  fp->SNII_min_mass_msun = SNII_min_mass_msun;
  fp->SNII_max_mass_msun = SNII_max_mass_msun;
  fp->log10_SNII_min_mass_msun = log10(SNII_min_mass_msun);
  fp->log10_SNII_max_mass_msun = log10(SNII_max_mass_msun);

  /* Properties of the energy fraction model */
  char energy_fraction[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "EAGLEFeedback:SNII_energy_fraction_function",
                          energy_fraction);

  if (strcmp(energy_fraction, "EAGLE") == 0) {
    fp->SNII_energy_scaling = SNII_scaling_EAGLE;
  } else if (strcmp(energy_fraction, "Separable") == 0) {
    fp->SNII_energy_scaling = SNII_scaling_separable;
  } else if (strcmp(energy_fraction, "Independent") == 0) {
    fp->SNII_energy_scaling = SNII_scaling_independent;
    fp->SNII_delta_E_n = parser_get_param_double(
        params, "EAGLEFeedback:SNII_energy_fraction_delta_E_n");
  } else {
    error(
        "Invalid value of "
        "EAGLEFeedback:SNII_energy_fraction_function: '%s'",
        energy_fraction);
  }

  fp->f_E_min =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_min");
  fp->f_E_max =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_max");
  fp->Z_0 =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_Z_0");
  fp->n_0_cgs = parser_get_param_double(
      params, "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3");

  /* Indicies can either be included in the parameter file as powers,
   * or as widths - but not both! The width and power are always related
   * by the relation power = 1.0 / (ln(10) * width) */

  int n_n_exists =
      parser_does_param_exist(params, "EAGLEFeedback:SNII_energy_fraction_n_n");
  int n_Z_exists =
      parser_does_param_exist(params, "EAGLEFeedback:SNII_energy_fraction_n_Z");
  int s_n_exists = parser_does_param_exist(
      params, "EAGLEFeedback:SNII_energy_fraction_sigma_n");
  int s_Z_exists = parser_does_param_exist(
      params, "EAGLEFeedback:SNII_energy_fraction_sigma_Z");

  if (n_n_exists && s_n_exists) {
    error("Cannot specify both n_n and sigma_n in SNII energy fraction.");
  }

  if (n_Z_exists && s_Z_exists) {
    error("Cannot specify both n_Z and sigma_Z in SNII energy fraction.");
  }

  if (n_n_exists) {
    fp->n_n = parser_get_param_double(params,
                                      "EAGLEFeedback:SNII_energy_fraction_n_n");
  } else {
    fp->n_n = 1.0 / (M_LN10 *
                     parser_get_param_double(
                         params, "EAGLEFeedback:SNII_energy_fraction_sigma_n"));
  }

  if (n_Z_exists) {
    fp->n_Z = parser_get_param_double(params,
                                      "EAGLEFeedback:SNII_energy_fraction_n_Z");
  } else {
    fp->n_Z = 1.0 / (M_LN10 *
                     parser_get_param_double(
                         params, "EAGLEFeedback:SNII_energy_fraction_sigma_Z"));
  }

  /* Check that it makes sense. */
  if (fp->f_E_max < fp->f_E_min) {
    error("Can't have the maximal energy fraction smaller than the minimal!");
  }

  /* Are we using the stars' birth properties or at feedback time? */
  fp->use_birth_density_for_f_th = parser_get_param_int(
      params, "EAGLEFeedback:SNII_energy_fraction_use_birth_density");
  fp->use_birth_Z_for_f_th = parser_get_param_int(
      params, "EAGLEFeedback:SNII_energy_fraction_use_birth_metallicity");

  /* Properties of the SNII enrichment model -------------------------------- */

  /* Set factors for each element adjusting SNII yield */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    char buffer[50];
    sprintf(buffer, "EAGLEFeedback:SNII_yield_factor_%s",
            chemistry_get_element_name((enum chemistry_element)elem));

    fp->SNII_yield_factor[elem] =
        parser_get_opt_param_float(params, buffer, 1.f);
  }

  /* Properties of the SNIa enrichment model -------------------------------- */

  fp->SNIa_DTD_delay_Gyr =
      parser_get_param_double(params, "EAGLEFeedback:SNIa_DTD_delay_Gyr");

  char temp[32] = {0};
  parser_get_param_string(params, "EAGLEFeedback:SNIa_DTD", temp);

  if (strcmp(temp, "Exponential") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_exponential;

    /* Read SNIa exponential DTD model parameters */
    fp->SNIa_DTD_exp_norm = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_exp_norm_p_Msun");
    fp->SNIa_DTD_exp_timescale_Gyr = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_exp_timescale_Gyr");
    fp->SNIa_DTD_exp_timescale_Gyr_inv = 1.f / fp->SNIa_DTD_exp_timescale_Gyr;

  } else if (strcmp(temp, "PowerLaw") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_power_law;

    /* Read SNIa power-law DTD model parameters */
    fp->SNIa_DTD_power_law_norm = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_power_law_norm_p_Msun");

    /* Renormalize everything such that the integral converges to
       'SNIa_DTD_power_law_norm' over 13.6 Gyr. */
    fp->SNIa_DTD_power_law_norm /= log(13.6 / fp->SNIa_DTD_delay_Gyr);

  } else {
    error("Invalid SNIa DTD model: '%s'", temp);
  }

  /* Energy released by supernova type Ia */
  fp->E_SNIa_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNIa_energy_erg");
  fp->E_SNIa =
      fp->E_SNIa_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Properties of the SNIa enrichment model -------------------------------- */

  /* Read AGB ejecta velocity */
  const float ejecta_velocity_km_p_s = parser_get_param_float(
      params, "EAGLEFeedback:AGB_ejecta_velocity_km_p_s");

  /* Convert to internal units */
  const float ejecta_velocity_cgs = ejecta_velocity_km_p_s * 1e5;
  const float ejecta_velocity =
      ejecta_velocity_cgs / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  /* Convert to specific thermal energy */
  fp->AGB_ejecta_specific_kinetic_energy =
      0.5f * ejecta_velocity * ejecta_velocity;

  /* Properties of the enrichment down-sampling ----------------------------- */

  fp->stellar_evolution_age_cut =
      parser_get_param_double(params,
                              "EAGLEFeedback:stellar_evolution_age_cut_Gyr") *
      phys_const->const_year * 1e9;

  fp->stellar_evolution_sampling_rate = parser_get_param_double(
      params, "EAGLEFeedback:stellar_evolution_sampling_rate");

  if (fp->stellar_evolution_sampling_rate < 1 ||
      fp->stellar_evolution_sampling_rate >= (1 << (8 * sizeof(char) - 1)))
    error("Stellar evolution sampling rate too large. Must be >0 and <%d",
          (1 << (8 * sizeof(char) - 1)));

  /* Check that we are not downsampling before reaching the SNII delay */
  if (!fp->SNII_sampled_delay &&
      fp->stellar_evolution_age_cut < fp->SNII_wind_delay)
    error(
        "Time at which the enrichment downsampling starts is lower than the "
        "SNII wind delay!");

  /* Gather common conversion factors --------------------------------------- */

  /* Calculate internal mass to solar mass conversion factor */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  fp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  fp->solar_mass_to_mass = 1. / fp->mass_to_solar_mass;

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  fp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  fp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Initialise the IMF ------------------------------------------------- */

  init_imf(fp);

  /* Calculate number of type II SN per unit solar mass based on our choice
   * of IMF and integration limits for type II SNe.
   * Note: No weighting by yields here. */
  fp->num_SNII_per_msun =
      integrate_imf(fp->log10_SNII_min_mass_msun, fp->log10_SNII_max_mass_msun,
                    eagle_imf_integration_no_weight,
                    /*(stellar_yields=)*/ NULL, fp);

  /* Initialise the yields ---------------------------------------------- */

  /* Read yield table filepath  */
  parser_get_param_string(params, "EAGLEFeedback:filename",
                          fp->yield_table_path);

  /* Allocate yield tables  */
  allocate_yield_tables(fp);

  /* Read the tables  */
  read_yield_tables(fp);

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (fp->log10_imf_max_mass_msun - fp->log10_imf_min_mass_msun) /
      (eagle_feedback_N_imf_bins - 1);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++)
    fp->yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + fp->log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(fp);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(fp);

  const double s_n = 1.0 / (M_LN10 * fp->n_n);
  const double s_Z = 1.0 / (M_LN10 * fp->n_Z);

  if (engine_rank == 0) {
    message("Feedback model is EAGLE (%s)", energy_fraction);
    message("Feedback energy fraction min=%f, max=%f", fp->f_E_min,
            fp->f_E_max);
    message("Feedback energy fraction powers: n_n=%f, n_Z=%f", fp->n_n,
            fp->n_Z);
    message("Feedback energy fraction widths: s_n=%f, s_Z=%f", s_n, s_Z);
    message("Feedback energy fraction pivots: Z_0=%f, n_0_cgs=%f", fp->Z_0,
            fp->n_0_cgs);
  }
}

/**
 * @brief Clean-up the memory allocated for the feedback routines
 *
 * We simply free all the arrays.
 *
 * @param fp the feedback data structure.
 */
void feedback_clean(struct feedback_props* fp) {

  swift_free("imf-tables", fp->imf);
  swift_free("imf-tables", fp->imf_mass_bin);
  swift_free("imf-tables", fp->imf_mass_bin_log10);
  swift_free("feedback-tables", fp->yields_SNIa);
  swift_free("feedback-tables", fp->yield_SNIa_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.mass);
  swift_free("feedback-tables", fp->yield_AGB.metallicity);
  swift_free("feedback-tables", fp->yield_AGB.yield);
  swift_free("feedback-tables", fp->yield_AGB.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.ejecta);
  swift_free("feedback-tables", fp->yield_AGB.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.total_metals);
  swift_free("feedback-tables", fp->yield_AGB.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.mass);
  swift_free("feedback-tables", fp->yield_SNII.metallicity);
  swift_free("feedback-tables", fp->yield_SNII.yield);
  swift_free("feedback-tables", fp->yield_SNII.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.ejecta);
  swift_free("feedback-tables", fp->yield_SNII.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.total_metals);
  swift_free("feedback-tables", fp->yield_SNII.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->lifetimes.mass);
  swift_free("feedback-tables", fp->lifetimes.metallicity);
  swift_free("feedback-tables", fp->yield_mass_bins);
  for (int i = 0; i < eagle_feedback_lifetime_N_metals; i++) {
    free(fp->lifetimes.dyingtime[i]);
  }
  free(fp->lifetimes.dyingtime);
  for (int i = 0; i < eagle_feedback_SNIa_N_elements; i++) {
    free(fp->SNIa_element_names[i]);
  }
  free(fp->SNIa_element_names);
  for (int i = 0; i < eagle_feedback_SNII_N_elements; i++) {
    free(fp->SNII_element_names[i]);
  }
  free(fp->SNII_element_names);
  for (int i = 0; i < eagle_feedback_AGB_N_elements; i++) {
    free(fp->AGB_element_names[i]);
  }
  free(fp->AGB_element_names);
}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the feedback routines. Helps debugging. */
  struct feedback_props feedback_copy = *feedback;

  /* zero AGB and SNII table pointers */
  zero_yield_table_pointers(&feedback_copy.yield_AGB);
  zero_yield_table_pointers(&feedback_copy.yield_SNII);

  /* zero SNIa table pointers */
  feedback_copy.yield_SNIa_IMF_resampled = NULL;
  feedback_copy.yields_SNIa = NULL;
  feedback_copy.yield_SNIa_total_metals_IMF_resampled = 0;

  /* zero element name tables */
  feedback_copy.SNIa_element_names = NULL;
  feedback_copy.SNII_element_names = NULL;
  feedback_copy.AGB_element_names = NULL;

  /* zero mass bins table */
  feedback_copy.yield_mass_bins = NULL;

  /* zero lifetime tracks */
  feedback_copy.lifetimes.mass = NULL;
  feedback_copy.lifetimes.metallicity = NULL;
  feedback_copy.lifetimes.dyingtime = NULL;

  /* zero IMF tables */
  feedback_copy.imf = NULL;
  feedback_copy.imf_mass_bin = NULL;
  feedback_copy.imf_mass_bin_log10 = NULL;

  restart_write_blocks((void*)&feedback_copy, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the feedback tables by
 * re-reading them.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {
  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");

  if (strlen(feedback->yield_table_path) != 0)
    feedback_restore_tables(feedback);
}
