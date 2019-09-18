//
// Created by zhrv on 27.08.19.
//

#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"

void charm_model_ns_timestep_diff_grad(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_timestep_diff_integrals(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);

void charm_model_ns_timestep_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_model_ns_timestep_diff_grad(p4est, ghost, ghost_data);
    charm_model_ns_timestep_diff_integrals(p4est, ghost, ghost_data);
}
