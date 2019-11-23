//
// Created by zhrv on 26.10.17.
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
#include <iostream>

#ifdef CHARM_CONFIG_YAML

void charm_model_ns_init(charm_ctx_t *ctx, YAML::Node model_node, YAML::Node yaml);
void charm_model_euler_init(charm_ctx_t *ctx, YAML::Node model_node, YAML::Node yaml);

static void _charm_init_fetch_bnd(YAML::Node node, charm_bnd_t *bnd)
{
    YAML::Node n2;
    
    strcpy(bnd->name, node["name"].as<std::string>().c_str());

    bnd->type = charm_bnd_type_by_name(node["type"].as<std::string>().c_str());
    bnd->params = nullptr;
    switch (bnd->type) {
        case BOUND_INLET:
            bnd->bnd_fn = charm_bnd_cond_fn_inlet;
            n2 = node["parameters"];
            bnd->params = CHARM_ALLOC(double, 5);
            bnd->params[0] = n2["V"]["x"].as<double>();
            bnd->params[1] = n2["V"]["y"].as<double>();
            bnd->params[2] = n2["V"]["z"].as<double>();
            bnd->params[3] = n2["T"].as<double>();
            bnd->params[4] = n2["P"].as<double>();
            break;
        case BOUND_OUTLET:
            bnd->bnd_fn = charm_bnd_cond_fn_outlet;
            break;
        case BOUND_WALL_SLIP:
            bnd->bnd_fn = charm_bnd_cond_fn_wall_slip;
            break;
        case BOUND_WALL_NO_SLIP: // @todo
            bnd->bnd_fn = charm_bnd_cond_fn_wall_no_slip;
            n2 = node["parameters"];
            bnd->params = CHARM_ALLOC(double, 1);
            bnd->params[0] = n2["T"].as<double>();
            break;
        default:
            CHARM_LERRORF("Unknown boundary type %d\n", bnd->type);

    }


}


static void _charm_init_bnd(charm_ctx_t *ctx, YAML::Node node)
{
    ctx->bnd = sc_array_new(sizeof(charm_bnd_t));
    for (auto it : node) {
        auto bnd = (charm_bnd_t *) sc_array_push(ctx->bnd);
        _charm_init_fetch_bnd(it, bnd);
    }

}


static void _charm_init_fetch_comp(YAML::Node node, charm_comp_t *comp)
{
    std::string str;
    comp->id = node["id"].as<int>();
    str = node["name"].as<std::string>();
    strcpy(comp->name, str.c_str());

    str = node["cp_type"].as<std::string>();
    if (str == "CONST") {
        comp->cp_type = COMP_CP_CONST;
    }
    else if (str == "POLYNOM") {
        comp->cp_type = COMP_CP_POLYNOM;
    }
    else {
        CHARM_LERRORF("Unknown Cp type '%s'. Use: CONST, POLYNOM.", str.c_str());
        charm_abort(nullptr, 1);
    }

    str = node["kp_type"].as<std::string>();
    if (str == "CONST") {
        comp->kp_type = COMP_KP_CONST;
    }
    else if (str == "SATHERLAND") {
        comp->kp_type = COMP_KP_SATHERLAND;
    }
    else {
        CHARM_LERRORF("Unknown KP type '%s'. Use: CONST, SATHERLAND.", str.c_str());
        charm_abort(nullptr, 1);
    }
    str = node["ml_type"].as<std::string>();
    if (str == "CONST") {
        comp->ml_type = COMP_ML_CONST;
    }
    else if (str == "SATHERLAND") {
        comp->ml_type = COMP_ML_SATHERLAND;
    }
    else {
        CHARM_LERRORF("Unknown ML type '%s'. Use: CONST, SATHERLAND.", str.c_str());
        charm_abort(nullptr, 1);
    }

    comp->m = node["M"].as<double>();
    comp->ml0 = node["ML0"].as<double>();
    comp->kp0 = node["KP0"].as<double>();
    comp->t0 = node["T0"].as<double>();
    comp->ts = node["TS"].as<double>();
    comp->sig = node["Sig"].as<double>();
    comp->ek = node["ek"].as<double>();
    comp->h0 = node["h0"].as<double>();

    comp->cp = sc_array_new(sizeof(double));
    YAML::Node cp = node["Cp"];
    for (auto it : cp) {
        auto tmp = (double*)sc_array_push(comp->cp);
        *tmp = it.as<double>();
    }
}


static void _charm_init_comps(charm_ctx_t *ctx, YAML::Node node)
{
    ctx->comp = sc_array_new(sizeof(charm_comp_t));
    for (auto c : node) {
        auto comp = (charm_comp_t *) sc_array_push(ctx->comp);
        _charm_init_fetch_comp(c, comp);
    }
}


static void _charm_init_fetch_mat(charm_ctx_t *ctx, YAML::Node node, charm_mat_t *mat)
{
    std::string str;
    mat->id = node["id"].as<int>();
    str = node["name"].as<std::string>();
    strcpy(mat->name, str.c_str());
    str = node["eof_type"].as<std::string>();
    if (str == "IDEAL") {
        mat->eos_fn = charm_mat_eos_ideal;
        if (ctx->comp->elem_count > 1) {
            CHARM_GLOBAL_ESSENTIAL("WARNING! There is more than one component in 'task.yaml'. First component's parameters is used for EOS. \n");
        }
    }
    else if (str == "MIX") {
        mat->eos_fn = charm_mat_eos_mix;
    }
    else if (str == "TABLE") {
        mat->eos_fn = charm_mat_eos_table;
    }
    else {
        CHARM_LERRORF("Unknown flux type '%s'. Use: LF, GODUNOV.", str.c_str());
        charm_abort(nullptr, 1);
    }
}


static void _charm_init_mat(charm_ctx_t *ctx, YAML::Node node)
{
    ctx->mat = sc_array_new(sizeof(charm_mat_t));
    for (auto it : node) {
        auto mat = (charm_mat_t *) sc_array_push(ctx->mat);
        _charm_init_fetch_mat(ctx, it, mat);
    }
}

static void _charm_init_fetch_reg(charm_ctx_t *ctx, YAML::Node node, charm_reg_t *reg)
{
    YAML::Node n1;
    int id, i;
    size_t idx;
    double c;

    std::string str;

    reg->id = node["id"].as<int>();
    str = node["name"].as<std::string>();
    strcpy(reg->name, str.c_str());
    reg->mat_id = node["material_id"].as<int>();

    n1 = node["parameters"];
    reg->v[0]   = n1["V"]["x"].as<double>();
    reg->v[1]   = n1["V"]["y"].as<double>();
    reg->v[2]   = n1["V"]["z"].as<double>();
    reg->t      = n1["T"].as<double>();
    reg->p      = n1["P"].as<double>();

    reg->grav[0]   = n1["G"]["x"].as<double>();
    reg->grav[1]   = n1["G"]["y"].as<double>();
    reg->grav[2]   = n1["G"]["z"].as<double>();

    memset(reg->c, 0, CHARM_MAX_COMPONETS_COUNT*sizeof(double));
    n1 = node["components"];
    for (auto it : n1) {
        id = it["id"].as<int>();
        c  = it["concentration"].as<double>();
        if (charm_comp_index_find_by_id(ctx, id, &idx)) {
            reg->c[idx] = c;
        }
        else {
            CHARM_LERRORF("Unknown component id %d for region '%s' in file 'task.yaml'\n", id, reg->name);
            charm_abort(nullptr, 1);
        }
    }

    c = 0;
    for (i = 0; i < CHARM_MAX_COMPONETS_COUNT; i++) {
        c += reg->c[i];
    }
    if (fabs(c)-1. > CHARM_EPS) {
        CHARM_LERRORF("Sum of concentrations for region '%s' is not equal to 1 in file 'task.yaml'\n", reg->name);
        charm_abort(nullptr, 1);
    }
}


static void _charm_init_reg(charm_ctx_t *ctx, YAML::Node node)
{
    ctx->reg = sc_array_new(sizeof(charm_reg_t));
    for (auto it : node) {
        auto reg = (charm_reg_t *) sc_array_push(ctx->reg);
        _charm_init_fetch_reg(ctx, it, reg);
    }
}


static void _charm_init_mesh_info(charm_ctx_t *ctx, YAML::Node node)
{
    charm_mesh_info_t *m;// = ctx->msh;
    std::string str;

    m = CHARM_ALLOC(charm_mesh_info_t, 1);

    str = node["files_type"].as<std::string>();
    m->type = charm_mesh_get_type_by_str((char*)str.c_str());
    str = node["name"].as<std::string>();
    strcpy(m->filename, str.c_str());

    ctx->msh = m;
}

extern "C" void charm_init_context_yaml(charm_ctx_t *ctx);

void charm_init_context_yaml(charm_ctx_t *ctx)
{
    std::string str;
    try {
        YAML::Node config = YAML::LoadFile("task.yaml");
        YAML::Node control = config["control"];

        str = control["FLUX_TYPE"].as<std::string>();
        if (str == "LF") {
            ctx->flux_fn = charm_calc_flux_lf;
        }
        else if (str == "GODUNOV") {
            ctx->flux_fn = charm_calc_flux_godunov;
        }
        else if (str == "HLLC") {
            ctx->flux_fn = charm_calc_flux_hllc;
        }
        else if (str == "CD") {
            ctx->flux_fn = charm_calc_flux_cd;
        }
        else {
            CHARM_LERRORF("Unknown flux type '%s'. Use: LF, GODUNOV.\n", str.c_str());
            charm_abort(nullptr, 1);
        }

        str = control["LIMITER"].as<std::string>();
        if (str == "NONE") {
            ctx->lim_fn = nullptr;
        }
        else if (str == "BJ") {
            ctx->lim_fn = charm_limiter_bj;
        }
        else {
            CHARM_LERRORF("Unknown limiter type '%s'. Use: NONE, BJ.\n", str.c_str());
            charm_abort(nullptr, 1);
        }


        ctx->max_err                = control["MAX_ERROR"].as<double>();
        ctx->refine_period          = control["REFINE_PERIOD"].as<int>();
        ctx->repartition_period     = control["REPARTITION_PERIOD"].as<int>();
        ctx->min_level              = control["MIN_LEVEL"].as<int>();
        ctx->max_level              = control["MAX_LEVEL"].as<int>();
        ctx->write_period           = control["FILE_OUTPUT_STEP"].as<int>();
        ctx->log_period             = control["LOG_OUTPUT_STEP"].as<int>();
        ctx->dt                     = control["TAU"].as<double>();
        ctx->CFL                    = control["CFL"].as<double>();
        ctx->time                   = control["TMAX"].as<double>();

        YAML::Node comps = config["components"];

        _charm_init_bnd(ctx, config["boundaries"]);




        _charm_init_bnd(       ctx, config["boundaries"]);
        _charm_init_comps(     ctx, config["components"]);
        _charm_init_mat(       ctx, config["materials"]);
        _charm_init_reg(       ctx, config["regions"]);
        _charm_init_mesh_info( ctx, config["mesh"]);

        ctx->model.ns.use_visc = 0;
        YAML::Node model = control["MODEL"];
        str = model["name"].as<std::string>();
        if (str == "EULER") {
            charm_model_euler_init(ctx, model, config);
        }
        else if (str == "NS") {
            charm_model_ns_init(ctx, model, config);
        }
        else {
            CHARM_LERRORF("Unknown model type '%s'. Use: EULER.\n", str.c_str());
            charm_abort(nullptr, 1);
        }

        ctx->timestep = 0;



    }
    catch(YAML::Exception &e) {
        std::cout << "Error: " << e.msg << std::endl;
        exit(1);
    }


}


#endif