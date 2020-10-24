//
// Created by zhrv on 15.11.19.
//
#include "charm_bnd_cond.h"
#include "charm_fluxes.h"
#include "charm_globals.h"
#include "charm_eos.h"
#include "charm_limiter.h"
#include "charm_models.h"
#include "charm_amr.h"
#include "yaml-cpp/yaml.h"
#include <cstring>
#include <cstdlib>


void charm_model_euler_init(charm_ctx_t *ctx, YAML::Node model_node, YAML::Node yaml)
{
    ctx->get_dt_fn              = charm_model_euler_get_dt;
    ctx->timestep_single_fn     = charm_model_euler_timestep_single;
    ctx->amr_init_fn            = charm_adapt_init;
    ctx->amr_fn                 = charm_adapt;
    ctx->model_init_cond_fn     = nullptr;
}

