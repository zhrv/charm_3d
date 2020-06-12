//
// Created by zhrv on 06.06.20.
//

#include <p8est_iterate.h>
#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"



void charm_model_ns_timestep_turb(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;

    if (ctx->model.ns.turb_fn) ctx->model.ns.turb_fn(p4est, ghost, ghost_data);
}