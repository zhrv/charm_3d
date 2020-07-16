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

extern "C" {
    void charm_model_ns_turb_sst(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
}

static void _charm_model_ns_chem_init_fetch_reaction(charm_ctx_t *ctx, YAML::Node node, charm_reaction_t *r)
{
    int i, comp;
    int c_count = ctx->comp->elem_count;

    YAML::Node left, right, n;
    r->a = node["a"].as<charm_real_t>();
    r->e = node["e"].as<charm_real_t>();
    r->n = node["n"].as<charm_real_t>();

    left = node["left"];
    i = 0;
    for (auto it : left) {
        if (i >= 3) {
            CHARM_LERROR("More than 3 components in reaction!\n");
            charm_abort(nullptr, 1);
        }
        comp = it.as<int>();
        if (!( -1 <= comp && comp < c_count )) {
            CHARM_LERRORF("Wrong component number in reaction: %d, but total components count: %d\n", comp, c_count);
            charm_abort(nullptr, 1);
        }
        r->left_comps[i] = comp;
        ++i;
    }

    right = node["right"];
    i = 0;
    for (auto it : right) {
        if (i >= 3) {
            CHARM_LERROR("More than 3 components in reaction!\n");
            charm_abort(nullptr, 1);
        }
        comp = it.as<int>();
        if (!( -1 <= comp && comp < c_count )) {
            CHARM_LERRORF("Wrong component number in reaction: %d, but total components count: %d\n", comp, c_count);
            charm_abort(nullptr, 1);
        }
        r->right_comps[i] = comp;
        ++i;
    }

}


static void _charm_model_ns_chem_init(charm_ctx_t *ctx, YAML::Node yaml)
{
    YAML::Node r_node;
    r_node = yaml["reactions"];
    if (!r_node.IsDefined()) {
        ctx->reactions = nullptr;
        return;
    }

    charm_reaction_t *r;
    ctx->reactions = sc_array_new(sizeof(charm_reaction_t));

    for (auto it : r_node) {
        r = (charm_reaction_t *) sc_array_push(ctx->reactions);
        _charm_model_ns_chem_init_fetch_reaction(ctx, it, r);
    }
}

static charm_turb_models_t _charm_turb_model_by_name(const char* name) {
    int i = 0;
    while (charm_turb_models[i] != NULL) {
        if (strcmp(charm_turb_models[i], name) == 0) {
            return (charm_turb_models_t)i;
        }
        i++;
    }
    return TURB_MODEL_UNKNOWN;
}

void _charm_model_ns_turb_sst_fetch_param(charm_ctx_t *ctx, YAML::Node par)
{
    ctx->model.ns.turb.param.sst.a1         = par["a1"].as<charm_real_t>();
    ctx->model.ns.turb.param.sst.sigma_k1   = par["sigma_k1"].as<charm_real_t>();
    ctx->model.ns.turb.param.sst.sigma_k2   = par["sigma_k2"].as<charm_real_t>();
    ctx->model.ns.turb.param.sst.sigma_w1   = par["sigma_w1"].as<charm_real_t>();
    ctx->model.ns.turb.param.sst.sigma_w2   = par["sigma_w2"].as<charm_real_t>();
    ctx->model.ns.turb.param.sst.beta_star  = par["beta_star"].as<charm_real_t>();
    ctx->model.ns.turb.param.sst.beta_1     = par["beta_1"].as<charm_real_t>();
    ctx->model.ns.turb.param.sst.beta_2     = par["beta_2"].as<charm_real_t>();
    ctx->model.ns.turb.param.sst.beta_2     = par["beta_2"].as<charm_real_t>();
}


static void _charm_model_ns_turb_init(charm_ctx_t *ctx, YAML::Node node)
{
    ctx->model.ns.turb.model_type = _charm_turb_model_by_name(node["model"].as<std::string>().c_str());
    switch (ctx->model.ns.turb.model_type) {
        case TURB_MODEL_SST:
            ctx->model.ns.turb.model_fn = charm_model_ns_turb_sst;
            _charm_model_ns_turb_sst_fetch_param(ctx, node["parameters"]);
            break;
    }
}


void charm_model_ns_init(charm_ctx_t *ctx, YAML::Node model_node, YAML::Node yaml)
{
    YAML::Node turb_node;

    ctx->get_dt_fn              = charm_model_ns_get_dt;
    ctx->timestep_single_fn     = charm_model_ns_timestep_single;
    ctx->model.ns.use_visc = model_node["use_visc"].as<int>();
    ctx->model.ns.use_diff = model_node["use_diffusion"].as<int>();
    ctx->model.ns.t_ref    = model_node["t_ref"].as<charm_real_t>();

    turb_node = model_node["turbulence"];
    if (turb_node) {
        _charm_model_ns_turb_init(ctx, turb_node);
    }

    ctx->amr_init_fn            = charm_adapt_init;
    ctx->amr_fn                 = charm_adapt;

    _charm_model_ns_chem_init(ctx, yaml);
}

