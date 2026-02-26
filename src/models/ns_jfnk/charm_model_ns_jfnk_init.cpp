//
// Created by zhrv on 10.01.26.
//
#include "charm_globals.h"
#include "charm_models.h"
#include "charm_amr.h"
#include "yaml-cpp/yaml.h"
#include <cstring>

extern "C" {
    void charm_model_ns_jfnk_init_initial_condition(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);
}


void charm_model_ns_jfnk_init(charm_ctx_t *ctx, YAML::Node model_node, const YAML::Node &yaml)
{
    YAML::Node turb_node;

    ctx->get_dt_fn              = charm_model_ns_jfnk_get_dt;

    ctx->timestep_single_fn     = charm_model_ns_jfnk_timestep_single;
    ctx->model.ns_jfnk.use_visc = model_node["use_visc"].as<int>();
    ctx->model.ns_jfnk.use_diff = model_node["use_diffusion"].as<int>();
    ctx->model.ns_jfnk.t_ref    = model_node["t_ref"].as<charm_real_t>();

    ctx->amr_init_fn            = charm_adapt_init;
    ctx->amr_fn                 = charm_adapt;
    ctx->model_init_cond_fn     = charm_model_ns_jfnk_init_initial_condition;

}

