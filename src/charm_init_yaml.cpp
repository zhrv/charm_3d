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
            bnd->params[4] = n2["P"].as<double>();
            break;
        default:
            CHARM_LERRORF("Unknown boundary type %d\n", bnd->type);

    }


}


static void _charm_init_bnd(charm_ctx_t *ctx, YAML::Node node)
{
    ctx->bnd = sc_array_new(sizeof(charm_bnd_t));
    for (auto it = node.begin(); it != node.end(); it++) {
        charm_bnd_t *bnd = (charm_bnd_t *) sc_array_push(ctx->bnd);
        _charm_init_fetch_bnd(*it, bnd);
    }

//    charm_bnd_t *bnd;
//    mxml_node_t *node1;
//    int i;
//
//    ctx->bnd = sc_array_new(sizeof(charm_bnd_t));
//
//    for (node1 = charm_xml_node_get_child(node, "boundCond");
//         node1 != NULL;
//         node1 = charm_xml_node_get_next_child(node1, node, "boundCond")) {
//        bnd = (charm_bnd_t *) sc_array_push(ctx->bnd);
//        _charm_init_fetch_bnd(node1, bnd);
//    }
}


static void charm_init_fetch_comp(mxml_node_t* node, charm_comp_t *comp)
{
//    mxml_node_t *n1, *n2;
//    char str[32];
//    double *tmp;
//
//
//    charm_xml_node_attr_int(node, "id", &(comp->id));
//    charm_xml_node_attr_str(node, "cp_type", str);
//    if (strcmp(str, "CONST") == 0) {
//        comp->cp_type = COMP_CP_CONST;
//    }
//    else if (strcmp(str, "POLYNOM") == 0) {
//        comp->cp_type = COMP_CP_POLYNOM;
//    }
//    else {
//        CHARM_LERRORF("Unknown Cp type '%s'. Use: CONST, POLYNOM.", str);
//        charm_abort(NULL, 1);
//    }
//    charm_xml_node_attr_str(node, "kp_type", str);
//    if (strcmp(str, "CONST") == 0) {
//        comp->kp_type = COMP_KP_CONST;
//    }
//    else if (strcmp(str, "SATHERLAND") == 0) {
//        comp->kp_type = COMP_KP_SATHERLAND;
//    }
//    else {
//        CHARM_LERRORF("Unknown KP type '%s'. Use: CONST, SATHERLAND.", str);
//        charm_abort(NULL, 1);
//    }
//    charm_xml_node_attr_str(node, "ml_type", str);
//    if (strcmp(str, "CONST") == 0) {
//        comp->ml_type = COMP_ML_CONST;
//    }
//    else if (strcmp(str, "SATHERLAND") == 0) {
//        comp->ml_type = COMP_ML_SATHERLAND;
//    }
//    else {
//        CHARM_LERRORF("Unknown ML type '%s'. Use: CONST, SATHERLAND.", str);
//        charm_abort(NULL, 1);
//    }
//    n1 = charm_xml_node_get_child(node, "name");
//    charm_xml_node_value_str(n1, comp->name);
//    //n1 = charm_xml_node_get_child(node, "parameters");
//    charm_xml_node_child_param_dbl(node, "M", &(comp->m));
//    charm_xml_node_child_param_dbl(node, "ML0", &(comp->ml0));
//    charm_xml_node_child_param_dbl(node, "KP0", &(comp->kp0));
//    charm_xml_node_child_param_dbl(node, "T0", &(comp->t0));
//    charm_xml_node_child_param_dbl(node, "TS", &(comp->ts));
//    charm_xml_node_child_param_dbl(node, "Sig", &(comp->sig));
//    charm_xml_node_child_param_dbl(node, "ek", &(comp->ek));
//    charm_xml_node_child_param_dbl(node, "h0", &(comp->h0));
//    comp->cp = sc_array_new(sizeof(double));
//    n1 = charm_xml_node_get_child(node, "Cp");
//    for (n2 = charm_xml_node_get_child(n1, "a");
//         n2 != NULL;
//         n2 = charm_xml_node_get_next_child(n2, n1, "a")) {
//        tmp = (double*)sc_array_push(comp->cp);
//        charm_xml_node_value_dbl(n2, tmp);
//    }
}


static void _charm_init_comps(charm_ctx_t *ctx, YAML::Node node) {
//    mxml_node_t *comp_node;
//    size_t mat_index;
//    charm_comp_t *comp;
//
//    ctx->comp = sc_array_new(sizeof(charm_comp_t));
//
//    for (comp_node = charm_xml_node_get_child(node, "component");
//         comp_node != NULL;
//         comp_node = charm_xml_node_get_next_child(comp_node, node, "component")) {
//        comp = (charm_comp_t *) sc_array_push(ctx->comp);
//        charm_init_fetch_comp(comp_node, comp);
//    }
}


static void charm_init_fetch_mat(charm_ctx_t *ctx, mxml_node_t* node, charm_mat_t *mat)
{
//    mxml_node_t *n1, *node1;
//    int id;
//    size_t *idx, idx_tmp;
//    char str[32];
//    charm_xml_node_attr_int(node, "id", &(mat->id));
//    charm_xml_node_child_param_str(node, "eof_type", str);
//    if (strcmp(str, "IDEAL") == 0) {
//        mat->eos_fn = charm_mat_eos_ideal;
//        if (ctx->comp->elem_count > 1) {
//            CHARM_GLOBAL_ESSENTIAL("WARNING! There is more than one component in 'task.xml'. First component's parameters is used for EOS. \n");
//        }
//    }
//    else if (strcmp(str, "MIX") == 0) {
//        mat->eos_fn = charm_mat_eos_mix;
//    }
//    else if (strcmp(str, "TABLE") == 0) {
//        mat->eos_fn = charm_mat_eos_table;
//    }
//    else {
//        CHARM_LERRORF("Unknown flux type '%s'. Use: LF, GODUNOV.", str);
//        charm_abort(NULL, 1);
//    }
//    n1 = charm_xml_node_get_child(node, "name");
//    charm_xml_node_value_str(n1, mat->name);
////    mat->comp_idx = sc_array_new(sizeof(size_t));
////    n1 = charm_xml_node_get_child(node, "components");
////    for (node1 = charm_xml_node_get_child(n1, "component");
////         node1 != NULL;
////         node1 = charm_xml_node_get_next_child(node1, n1, "component")) {
////
////        charm_xml_node_attr_int(node1, "id", &id);
////        if (charm_comp_index_find_by_id(ctx, id, &idx_tmp)) {
////            idx = (size_t *) sc_array_push(mat->comp_idx);
////            *idx = idx_tmp;
////        }
////        else {
////            CHARM_LERRORF("Unknown component id %d for material '%s' in file 'task.xml'\n", id, mat->name);
////            charm_abort(1);
////        }
////    }
}


static void _charm_init_mat(charm_ctx_t *ctx, YAML::Node node)
{
//    charm_mat_t *mat;
//    mxml_node_t *node1;
//    int i;
//
//    ctx->mat = sc_array_new(sizeof(charm_mat_t));
//
//    for (node1 = charm_xml_node_get_child(node, "material");
//         node1 != NULL;
//         node1 = charm_xml_node_get_next_child(node1, node, "material")) {
//        mat = (charm_mat_t *) sc_array_push(ctx->mat);
//        charm_init_fetch_mat(ctx, node1, mat);
//    }
}

static void charm_init_fetch_reg(charm_ctx_t *ctx, mxml_node_t* node, charm_reg_t *reg)
{
//    mxml_node_t * n1, *n2;
//    int tmp, id, i;
//    size_t idx;
//    double c;
//
//    charm_xml_node_attr_int(node, "id", &(reg->id));
//    n1 = charm_xml_node_get_child(node, "name");
//    charm_xml_node_value_str(n1, reg->name);
//    //charm_xml_node_child_param_int(node, "cell_type", &(reg->cell_type));
//    charm_xml_node_child_param_int(node, "material_id", &(reg->mat_id));
//
//    n1 = charm_xml_node_get_child(node, "parameters");
//    charm_xml_node_child_param_dbl(n1, "Vx", &(reg->v[0]));
//    charm_xml_node_child_param_dbl(n1, "Vy", &(reg->v[1]));
//    charm_xml_node_child_param_dbl(n1, "Vz", &(reg->v[2]));
//    charm_xml_node_child_param_dbl(n1, "T", &(reg->t));
//    charm_xml_node_child_param_dbl(n1, "P", &(reg->p));
//
//    charm_xml_node_child_param_dbl(n1, "Gx", &(reg->grav[0]));
//    charm_xml_node_child_param_dbl(n1, "Gy", &(reg->grav[1]));
//    charm_xml_node_child_param_dbl(n1, "Gz", &(reg->grav[2]));
//
//    memset(reg->c, 0, CHARM_MAX_COMPONETS_COUNT*sizeof(double));
//    n1 = charm_xml_node_get_child(node, "components");
//    for (n2 = charm_xml_node_get_child(n1, "component");
//         n2 != NULL;
//         n2 = charm_xml_node_get_next_child(n2, n1, "component")) {
//        charm_xml_node_attr_int(n2, "id", &id);
//        charm_xml_node_attr_dbl(n2, "concentration", &c);
//        if (charm_comp_index_find_by_id(ctx, id, &idx)) {
//            reg->c[idx] = c;
//        }
//        else {
//            CHARM_LERRORF("Unknown component id %d for region '%s' in file 'task.xml'\n", id, reg->name);
//            charm_abort(NULL, 1);
//        }
//    }
//
//    c = 0;
//    for (i = 0; i < CHARM_MAX_COMPONETS_COUNT; i++) {
//        c += reg->c[i];
//    }
//    if (fabs(c)-1. > CHARM_EPS) {
//        CHARM_LERRORF("Sum of concentrations for region '%s' is not equal to 1 in file 'task.xml'\n", reg->name);
//        charm_abort(NULL, 1);
//    }
}


static void _charm_init_reg(charm_ctx_t *ctx, YAML::Node node)
{
//    charm_reg_t *reg;
//    mxml_node_t *node1;
//    int i;
//
//    ctx->reg = sc_array_new(sizeof(charm_reg_t));
//
//    for (node1 = charm_xml_node_get_child(node, "region");
//         node1 != NULL;
//         node1 = charm_xml_node_get_next_child(node1, node, "region")) {
//        reg = (charm_reg_t *) sc_array_push(ctx->reg);
//        charm_init_fetch_reg(ctx, node1, reg);
//    }
}


static void _charm_init_mesh_info(charm_ctx_t *ctx, YAML::Node node)
{
//    charm_mesh_info_t *m;// = ctx->msh;
//    char str[128];
//
//    m = (charm_mesh_info_t*)malloc(sizeof(charm_mesh_info_t));
//
//    charm_xml_node_child_param_str(node, "name", m->filename);
//    charm_xml_node_child_param_str(node, "filesType", str);
//    m->type = charm_mesh_get_type_by_str(str);
//    ctx->msh = m;
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
            CHARM_LERRORF("Unknown limiter type '%s'. Use: NONE, BJ.\n", str);
            charm_abort(nullptr, 1);
        }


        ctx->max_err = control["MAX_ERROR"].as<double>();
        ctx->refine_period = control["REFINE_PERIOD"].as<int>();
        ctx->repartition_period = control["REPARTITION_PERIOD"].as<int>();
        ctx->min_level = control["MIN_LEVEL"].as<int>();
        ctx->max_level = control["MAX_LEVEL"].as<int>();
        ctx->write_period = control["FILE_OUTPUT_STEP"].as<int>();
        ctx->log_period = control["LOG_OUTPUT_STEP"].as<int>();
        ctx->dt = control["TAU"].as<double>();
        ctx->CFL = control["CFL"].as<double>();
        ctx->time = control["TMAX"].as<double>();

        YAML::Node comps = config["components"];

        _charm_init_bnd(ctx, config["boundaries"]);




        _charm_init_bnd(       ctx, config["boundaries"]);
        _charm_init_comps(     ctx, config["components"]);
        _charm_init_mat(       ctx, config["materials"]);
        _charm_init_reg(       ctx, config["regions"]);
        _charm_init_mesh_info( ctx, config["mesh"]);

//        ctx->use_visc = 0;
//        charm_xml_node_child_param_str(node, "MODEL", str);
//        if (strcmp(str, "EULER") == 0) {
//            charm_model_euler_init(ctx, charm_xml_node_get_child(node, "MODEL"), node_task);
//        }
//        else if (strcmp(str, "NS") == 0) {
//            charm_model_ns_init(ctx, charm_xml_node_get_child(node, "MODEL"), node_task);
//        }
//        else {
//            CHARM_LERRORF("Unknown model type '%s'. Use: EULER.\n", str);
//            charm_abort(NULL, 1);
//        }

        ctx->timestep = 0;



    }
    catch(YAML::Exception e) {
        std::cout << "Error: " << e.msg << std::endl;
        exit(1);
    }

    return;
}
