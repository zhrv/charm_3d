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

static void _charm_init_fetch_bnd(charm_ctx_t *ctx, YAML::Node node, charm_bnd_t *bnd)
{
    YAML::Node n2, n3;
    charm_int_t i;
    charm_int_t c_count = ctx->comp->elem_count;
    strcpy(bnd->name, node["name"].as<std::string>().c_str());

    bnd->type = charm_bnd_type_by_name(node["type"].as<std::string>().c_str());
    bnd->params = nullptr;
    switch (bnd->type) {
        case BOUND_INLET:
            bnd->bnd_fn = charm_bnd_cond_fn_inlet;
            n2 = node["parameters"];
            bnd->params = CHARM_ALLOC(charm_real_t, 5);
            bnd->params[0] = n2["V"]["x"].as<charm_real_t>();
            bnd->params[1] = n2["V"]["y"].as<charm_real_t>();
            bnd->params[2] = n2["V"]["z"].as<charm_real_t>();
            bnd->params[3] = n2["T"].as<charm_real_t>();
            bnd->params[4] = n2["P"].as<charm_real_t>();
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
            bnd->params = CHARM_ALLOC(charm_real_t, 1);
            bnd->params[0] = n2["T"].as<charm_real_t>();
            break;
        case BOUND_MASS_FLOW: // @todo
            bnd->bnd_fn = charm_bnd_cond_fn_mass_flow;
            n2 = node["parameters"];
            bnd->params = CHARM_ALLOC(charm_real_t, 7+c_count);
            bnd->params[0] = n2["M"].as<charm_real_t>();
            bnd->params[1] = n2["P"].as<charm_real_t>();
            bnd->params[2] = n2["T"].as<charm_real_t>();
            bnd->params[3] = n2["CosX"].as<charm_real_t>();
            bnd->params[4] = n2["CosY"].as<charm_real_t>();
            bnd->params[5] = n2["CosZ"].as<charm_real_t>();
            bnd->params[6] = n2["P"].as<charm_real_t>();
            n3 = n2["components"];
            i = 7;
            for (auto it : n3) {
                if (i > 7+c_count) {
                    CHARM_LERRORF("BOUND_MASS_FLOW: Too many components specified. Must be %d\n", c_count);
                }
                bnd->params[i++] = it.as<charm_real_t>();
            }
            if (i < 7+c_count) {
                CHARM_LERRORF("BOUND_MASS_FLOW: Too few components specified. Must be %d\n", c_count);
            }

            break;
        case BOUND_UNKNOWN: // @todo
        default:
            CHARM_LERRORF("Unknown boundary type %s\n", bnd->name);

    }


}


static void _charm_init_bnd(charm_ctx_t *ctx, YAML::Node node)
{
    ctx->bnd = sc_array_new(sizeof(charm_bnd_t));
    for (auto it : node) {
        auto bnd = (charm_bnd_t *) sc_array_push(ctx->bnd);
        _charm_init_fetch_bnd(ctx, it, bnd);
    }

}


static void _charm_init_fetch_comp(charm_ctx_t *ctx, YAML::Node node, charm_comp_t *comp)
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

    comp->m = node["M"].as<charm_real_t>();
    comp->ml0 = node["ML0"].as<charm_real_t>();
    comp->kp0 = node["KP0"].as<charm_real_t>();
    comp->t0 = node["T0"].as<charm_real_t>();
    comp->ts = node["TS"].as<charm_real_t>();
    comp->sig = node["Sig"].as<charm_real_t>();
    comp->ek = node["ek"].as<charm_real_t>();
    comp->h0 = node["h0"].as<charm_real_t>();

    comp->cp = sc_array_new(sizeof(charm_real_t));
    YAML::Node cp = node["Cp"];
    for (auto it : cp) {
        auto tmp = (charm_real_t*)sc_array_push(comp->cp);
        *tmp = it.as<charm_real_t>();
    }
}


static void _charm_init_comps(charm_ctx_t *ctx, YAML::Node node)
{
    ctx->comp = sc_array_new(sizeof(charm_comp_t));
    for (auto c : node) {
        auto comp = (charm_comp_t *) sc_array_push(ctx->comp);
        _charm_init_fetch_comp(ctx, c, comp);
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
    charm_real_t c;

    std::string str;

    reg->id = node["id"].as<int>();
    str = node["name"].as<std::string>();
    strcpy(reg->name, str.c_str());
    reg->mat_id = node["material_id"].as<int>();

    n1 = node["parameters"];
    reg->v[0]   = n1["V"]["x"].as<charm_real_t>();
    reg->v[1]   = n1["V"]["y"].as<charm_real_t>();
    reg->v[2]   = n1["V"]["z"].as<charm_real_t>();
    reg->t      = n1["T"].as<charm_real_t>();
    reg->p      = n1["P"].as<charm_real_t>();

    reg->grav[0]   = n1["G"]["x"].as<charm_real_t>();
    reg->grav[1]   = n1["G"]["y"].as<charm_real_t>();
    reg->grav[2]   = n1["G"]["z"].as<charm_real_t>();

    memset(reg->c, 0, CHARM_MAX_COMPONETS_COUNT*sizeof(charm_real_t));
    n1 = node["components"];
    for (auto it : n1) {
        id = it["id"].as<int>();
        c  = it["concentration"].as<charm_real_t>();
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




bool _charm_init_yaml_check_version(std::string v)
{
    std::vector<int> ver, cver;
    std::string delimiter = ".";

    size_t pos = 0;
    std::string token;
    while ((pos = v.find(delimiter)) != std::string::npos) {
        token = v.substr(0, pos);
        ver.push_back(std::stoi(token));
        v.erase(0, pos + delimiter.length());
    }
    ver.push_back(std::stoi(v));

    pos = 0;
    std::string cv = YAML_VERSION;
    while ((pos = cv.find(delimiter)) != std::string::npos) {
        token = cv.substr(0, pos);
        cver.push_back(std::stoi(token));
        cv.erase(0, pos + delimiter.length());
    }
    cver.push_back(std::stoi(cv));

    if (ver.size() != cver.size()) {
        return false;
    }
    if (ver[0] > cver[0]) {
        return true;
    }
    else if (ver[0] == cver[0]) {
        if (ver[1] > cver[1]) {
            return true;
        }
        else if (ver[1] == cver[1]) {
            if (ver[2] >= cver[2]) {
                return true;
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

extern "C" void charm_init_context_yaml(charm_ctx_t *ctx);

void charm_init_context_yaml(charm_ctx_t *ctx)
{
    std::string str;
    try {
        YAML::Node config = YAML::LoadFile("task.yaml");
        str = config["method"].as<std::string>();
        if (str != "CHARM_3D") {
            throw YAML::Exception(YAML::Mark::null_mark(), "Method name must be 'CHARM_3D'");
        }
        
        str = config["version"].as<std::string>();

        if (!_charm_init_yaml_check_version(str)) {
            throw YAML::Exception(YAML::Mark::null_mark(), "Wrong YAML version. Must be "+std::string(YAML_VERSION)+" or higher.");
        }


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


        ctx->max_err                = control["MAX_ERROR"].as<charm_real_t>();
        ctx->refine_period          = control["REFINE_PERIOD"].as<int>();
        ctx->repartition_period     = control["REPARTITION_PERIOD"].as<int>();
        ctx->min_level              = control["MIN_LEVEL"].as<int>();
        ctx->max_level              = control["MAX_LEVEL"].as<int>();
        ctx->write_period           = control["FILE_OUTPUT_STEP"].as<int>();
        ctx->log_period             = control["LOG_OUTPUT_STEP"].as<int>();
        ctx->dt                     = control["TAU"].as<charm_real_t>();
        ctx->CFL                    = control["CFL"].as<charm_real_t>();
        ctx->time                   = control["TMAX"].as<charm_real_t>();


        _charm_init_comps(     ctx, config["components"]);
        _charm_init_bnd(       ctx, config["boundaries"]);
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
        std::cout << "Error (YAML): " << e.msg << std::endl;
        exit(1);
    }


}


#endif