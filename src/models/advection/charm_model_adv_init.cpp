//
// Created by zhrv on 15.11.19.
//
#include "charm_bnd_cond.h"
#include "charm_fluxes.h"
#include "charm_globals.h"
#include "charm_eos.h"
#include "charm_limiter.h"
#include "charm_models.h"
#include "yaml-cpp/yaml.h"
#include <cstring>
#include <cstdlib>

extern "C" {
    charm_real_t charm_model_adv_get_dt (p4est_t * p4est);
    void charm_model_adv_timestep_single(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data);
}
void charm_model_adv_init(charm_ctx_t *ctx, YAML::Node model_node, YAML::Node yaml)
{
    ctx->get_dt_fn              = charm_model_adv_get_dt;
    ctx->timestep_single_fn     = charm_model_adv_timestep_single;
}

