//
// Created by yuyttenh on 25/10/2021.
//

#ifndef SWIFTSIM_HYDRO_GETTERS_H
#define SWIFTSIM_HYDRO_GETTERS_H

__attribute__((always_inline)) INLINE static void hydro_get_primitives(struct part *p, double *W) {
  W[0] = p->rho;
  W[1] = p->fluid_v[0];
  W[2] = p->fluid_v[1];
  W[3] = p->fluid_v[2];
  W[4] = p->P;
}

#endif  // SWIFTSIM_HYDRO_GETTERS_H
