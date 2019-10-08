//
// Created by zhrv on 30.07.19.
//

#ifndef CHARM_DG_CHARM_MODELS_H
#define CHARM_DG_CHARM_MODELS_H
#include "charm_globals.h"
#include "charm_xml.h"


/** MODEL INVISCID */
void    charm_model_euler_timestep_single   (p4est_t * p4est, double *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data);
double  charm_model_euler_get_dt            (p4est_t * p4est);
void    charm_model_euler_init              (charm_ctx_t *ctx, mxml_node_t *node);

/** MODEL NAVIER-STOKES */
void    charm_model_ns_timestep_single  (p4est_t * p4est, double *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data);
double  charm_model_ns_get_dt           (p4est_t * p4est);
void    charm_model_ns_init             (charm_ctx_t *ctx, mxml_node_t *node);

/** MODEL NAVIER-STOKES LOW MACH */
void    charm_model_ns_low_mach_timestep_single  (p4est_t * p4est, double *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data);
double  charm_model_ns_low_mach_get_dt           (p4est_t * p4est);
void    charm_model_ns_low_mach_init             (charm_ctx_t *ctx, mxml_node_t *node);

#endif //CHARM_DG_CHARM_MODELS_H
