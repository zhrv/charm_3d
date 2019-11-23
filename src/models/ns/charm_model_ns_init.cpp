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


static void _charm_model_ns_chem_init_fetch_reaction(charm_ctx_t *ctx, YAML::Node node, charm_reaction_t *r)
{
    int i, comp;
    int c_count = ctx->comp->elem_count;

    YAML::Node left, right, n;
    r->a = node["a"].as<double>();
    r->e = node["e"].as<double>();
    r->n = node["n"].as<double>();

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


void charm_model_ns_init(charm_ctx_t *ctx, YAML::Node model_node, YAML::Node yaml)
{
    ctx->get_dt_fn              = charm_model_ns_get_dt;
    ctx->timestep_single_fn     = charm_model_ns_timestep_single;
    ctx->model.ns.use_visc = model_node["use_visc"].as<int>();
    ctx->model.ns.use_diff = model_node["use_diffusion"].as<int>();
    ctx->model.ns.t_ref    = model_node["t_ref"].as<double>();

    _charm_model_ns_chem_init(ctx, yaml);
}

