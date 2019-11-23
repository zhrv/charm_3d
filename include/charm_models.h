//
// Created by zhrv on 30.07.19.
//

#ifndef CHARM_DG_CHARM_MODELS_H
#define CHARM_DG_CHARM_MODELS_H
#include "charm_globals.h"

#ifdef __cplusplus
extern "C" {
#endif


/** MODEL INVISCID */
void charm_model_euler_timestep_single(p4est_t *p4est, double *dt, p4est_ghost_t **_ghost, charm_data_t **_ghost_data);

double charm_model_euler_get_dt(p4est_t *p4est);


/** MODEL NAVIER-STOKES */
void charm_model_ns_timestep_single(p4est_t *p4est, double *dt, p4est_ghost_t **_ghost, charm_data_t **_ghost_data);

double charm_model_ns_get_dt(p4est_t *p4est);


#ifdef __cplusplus
}
#endif

#endif //CHARM_DG_CHARM_MODELS_H
