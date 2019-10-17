//
// Created by zhrv on 26.10.17.
//

#include "charm_xml.h"
#include "charm_bnd_cond.h"
#include "charm_fluxes.h"
#include "charm_globals.h"
#include "charm_eos.h"
#include "charm_limiter.h"
#include "charm_models.h"

void charm_init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *data    = (charm_data_t *) q->p.user_data;
    charm_param_t      *par     = &data->par;
    charm_prim_t        prim;
    charm_cons_t        cons;
    charm_tree_attr_t  *attr;
    charm_reg_t        *reg;
    int                 i;
    size_t              c_count = charm_get_comp_count(p4est);
    charm_mat_t        *mat;

    charm_geom_quad_calc(p4est, q, which_tree);

    attr = charm_get_tree_attr(p4est, which_tree);
    reg = attr->reg;
#ifdef POGGI
    double *x   = data->par.g.c;
    double pi   = 4.*atan(1.);
    double pi2  = pi*2.;
    double a0 = 0.5e-3;
    double lambda = 1.e-3;

    if (x[2] < -5.1e-3) {
        reg = charm_reg_find_by_id(ctx, 0);
    }
    else if ( x[2] > -a0*(1.-cos(pi2*x[0]/lambda))*(1.-cos(pi2*x[1]/lambda)) ) {
        reg = charm_reg_find_by_id(ctx, 2);
    }
    else {
        reg = charm_reg_find_by_id(ctx, 1);
    }
#endif
    prim.mat_id = reg->mat_id;
    prim.p   = reg->p;
    prim.t   = reg->t;

    prim.u   = reg->v[0];
    prim.v   = reg->v[1];
    prim.w   = reg->v[2];
    for (i = 0; i < c_count; i++) {
        prim.c[i] = reg->c[i];
    }
    mat = charm_mat_find_by_id(ctx, reg->mat_id);

    mat->eos_fn(p4est, &prim, 3); // (T,p) => (r, cz, e)
    prim.e_tot = prim.e + 0.5*(prim.u*prim.u+prim.v*prim.v+prim.w*prim.w);
    charm_param_prim_to_cons(p4est, &cons, &prim);
    memset(&(par->c), 0, sizeof(par->c));
    par->c.ru[0] = cons.ru;
    par->c.rv[0] = cons.rv;
    par->c.rw[0] = cons.rw;
    par->c.re[0] = cons.re;
    for (i = 0; i < c_count; i++) {
        par->c.rc[i][0] = cons.rc[i];
    }
    par->mat_id = reg->mat_id;
    for (i = 0; i < CHARM_DIM; i++) {
        par->grav[i] = reg->grav[i];
    }
}


void charm_init_fetch_bnd(mxml_node_t* node, charm_bnd_t *bnd)
{
    mxml_node_t * n1, *n2, *n3;
    char str[64];
    n1 = charm_xml_node_get_child(node, "name");
    charm_xml_node_value_str(n1, bnd->name);
    n1 = charm_xml_node_get_child(node, "type");
    charm_xml_node_value_str(n1, str);
    bnd->type = charm_bnd_type_by_name(str);
    bnd->params = NULL;
    switch (bnd->type) {
        case BOUND_INLET:
            bnd->bnd_fn = charm_bnd_cond_fn_inlet;
            n2 = charm_xml_node_get_child(node, "parameters");
            bnd->params = CHARM_ALLOC(double, 5);
            charm_xml_node_child_param_dbl(n2, "Vx", &(bnd->params[0]));
            charm_xml_node_child_param_dbl(n2, "Vy", &(bnd->params[1]));
            charm_xml_node_child_param_dbl(n2, "Vz", &(bnd->params[2]));
            charm_xml_node_child_param_dbl(n2, "T",  &(bnd->params[3]));
            charm_xml_node_child_param_dbl(n2, "P",  &(bnd->params[4]));
            break;
        case BOUND_OUTLET:
            bnd->bnd_fn = charm_bnd_cond_fn_outlet;
            break;
        case BOUND_WALL_SLIP:
            bnd->bnd_fn = charm_bnd_cond_fn_wall_slip;
            break;
        case BOUND_WALL_NO_SLIP: // @todo
            bnd->bnd_fn = charm_bnd_cond_fn_wall_no_slip;
            n2 = charm_xml_node_get_child(node, "parameters");
            bnd->params = CHARM_ALLOC(double, 1);
            charm_xml_node_child_param_dbl(n2, "T", &(bnd->params[0]));
            break;
        default:
            CHARM_LERRORF("Unknown boundary type %d\n", bnd->type);

    }


}


static void charm_init_bnd(charm_ctx_t *ctx, mxml_node_t *node)
{
    charm_bnd_t *bnd;
    mxml_node_t *node1;
    int i;

    ctx->bnd = sc_array_new(sizeof(charm_bnd_t));

    for (node1 = charm_xml_node_get_child(node, "boundCond");
         node1 != NULL;
         node1 = charm_xml_node_get_next_child(node1, node, "boundCond")) {
        bnd = (charm_bnd_t *) sc_array_push(ctx->bnd);
        charm_init_fetch_bnd(node1, bnd);
    }
}


static void charm_init_fetch_comp(mxml_node_t* node, charm_comp_t *comp)
{
    mxml_node_t * n1;
    char str[32];
    double *tmp;


    charm_xml_node_attr_int(node, "id", &(comp->id));
    charm_xml_node_attr_str(node, "cp_type", str);
    if (strcmp(str, "CONST") == 0) {
        comp->cp_type = COMP_CONST;
    }
    else if (strcmp(str, "POLYNOM") == 0) {
        comp->cp_type = COMP_POLYNOM;
    }
    else {
        CHARM_LERRORF("Unknown Cp type '%s'. Use: CONST, POLYNOM.", str);
        charm_abort(NULL, 1);
    }
    n1 = charm_xml_node_get_child(node, "name");
    charm_xml_node_value_str(n1, comp->name);
    //n1 = charm_xml_node_get_child(node, "parameters");
    charm_xml_node_child_param_dbl(node, "M", &(comp->m));
    charm_xml_node_child_param_dbl(node, "ML", &(comp->ml));
    charm_xml_node_child_param_dbl(node, "Lambda", &(comp->lambda));
    charm_xml_node_child_param_dbl(node, "K", &(comp->k));
    comp->cp = sc_array_new(sizeof(double));
    if (comp->cp_type == COMP_CONST) {
        tmp = sc_array_push(comp->cp);
        charm_xml_node_child_param_dbl(node, "Cp", tmp);
    }
    else {
        CHARM_LERROR("Cp type 'POLYNOM' is not released.\n");
        charm_abort(NULL, 1);
    }
}


static void charm_init_comps(charm_ctx_t* ctx, mxml_node_t* node) {
    mxml_node_t *comp_node;
    size_t mat_index;
    charm_comp_t *comp;

    ctx->comp = sc_array_new(sizeof(charm_comp_t));

    for (comp_node = charm_xml_node_get_child(node, "component");
         comp_node != NULL;
         comp_node = charm_xml_node_get_next_child(comp_node, node, "component")) {
        comp = (charm_comp_t *) sc_array_push(ctx->comp);
        charm_init_fetch_comp(comp_node, comp);
    }
}


static void charm_init_fetch_mat(charm_ctx_t *ctx, mxml_node_t* node, charm_mat_t *mat)
{
    mxml_node_t *n1, *node1;
    int id;
    size_t *idx, idx_tmp;
    char str[32];
    charm_xml_node_attr_int(node, "id", &(mat->id));
    charm_xml_node_child_param_str(node, "eof_type", str);
    if (strcmp(str, "IDEAL") == 0) {
        mat->eos_fn = charm_mat_eos_ideal;
        if (ctx->comp->elem_count > 1) {
            CHARM_GLOBAL_ESSENTIAL("WARNING! There is more than one component in 'task.xml'. First component's parameters is used for EOS. \n");
        }
    }
    else if (strcmp(str, "IDEAL_LOW_MACH") == 0) {
        mat->eos_fn = charm_mat_eos_ideal_low_mach;
        if (ctx->comp->elem_count > 1) {
            CHARM_GLOBAL_ESSENTIAL("WARNING! There is more than one component in 'task.xml'. First component's parameters is used for EOS. \n");
        }
    }
    else if (strcmp(str, "MIX") == 0) {
        mat->eos_fn = charm_mat_eos_mix;
    }
    else if (strcmp(str, "MIX_LOW_MACH") == 0) {
        mat->eos_fn = charm_mat_eos_mix_low_mach;
    }
    else if (strcmp(str, "TABLE") == 0) {
        mat->eos_fn = charm_mat_eos_table;
    }
    else {
        CHARM_LERRORF("Unknown flux type '%s'. Use: LF, GODUNOV.", str);
        charm_abort(NULL, 1);
    }
    n1 = charm_xml_node_get_child(node, "name");
    charm_xml_node_value_str(n1, mat->name);
//    mat->comp_idx = sc_array_new(sizeof(size_t));
//    n1 = charm_xml_node_get_child(node, "components");
//    for (node1 = charm_xml_node_get_child(n1, "component");
//         node1 != NULL;
//         node1 = charm_xml_node_get_next_child(node1, n1, "component")) {
//
//        charm_xml_node_attr_int(node1, "id", &id);
//        if (charm_comp_index_find_by_id(ctx, id, &idx_tmp)) {
//            idx = (size_t *) sc_array_push(mat->comp_idx);
//            *idx = idx_tmp;
//        }
//        else {
//            CHARM_LERRORF("Unknown component id %d for material '%s' in file 'task.xml'\n", id, mat->name);
//            charm_abort(1);
//        }
//    }
}


static void charm_init_mat(charm_ctx_t *ctx, mxml_node_t *node)
{
    charm_mat_t *mat;
    mxml_node_t *node1;
    int i;

    ctx->mat = sc_array_new(sizeof(charm_mat_t));

    for (node1 = charm_xml_node_get_child(node, "material");
         node1 != NULL;
         node1 = charm_xml_node_get_next_child(node1, node, "material")) {
        mat = (charm_mat_t *) sc_array_push(ctx->mat);
        charm_init_fetch_mat(ctx, node1, mat);
    }
}

static void charm_init_fetch_reg(charm_ctx_t *ctx, mxml_node_t* node, charm_reg_t *reg)
{
    mxml_node_t * n1, *n2;
    int tmp, id, i;
    size_t idx;
    double c;

    charm_xml_node_attr_int(node, "id", &(reg->id));
    n1 = charm_xml_node_get_child(node, "name");
    charm_xml_node_value_str(n1, reg->name);
    //charm_xml_node_child_param_int(node, "cell_type", &(reg->cell_type));
    charm_xml_node_child_param_int(node, "material_id", &(reg->mat_id));

    n1 = charm_xml_node_get_child(node, "parameters");
    charm_xml_node_child_param_dbl(n1, "Vx", &(reg->v[0]));
    charm_xml_node_child_param_dbl(n1, "Vy", &(reg->v[1]));
    charm_xml_node_child_param_dbl(n1, "Vz", &(reg->v[2]));
    charm_xml_node_child_param_dbl(n1, "T", &(reg->t));
    charm_xml_node_child_param_dbl(n1, "P", &(reg->p));

    charm_xml_node_child_param_dbl(n1, "Gx", &(reg->grav[0]));
    charm_xml_node_child_param_dbl(n1, "Gy", &(reg->grav[1]));
    charm_xml_node_child_param_dbl(n1, "Gz", &(reg->grav[2]));

    memset(reg->c, 0, CHARM_MAX_COMPONETS_COUNT*sizeof(double));
    n1 = charm_xml_node_get_child(node, "components");
    for (n2 = charm_xml_node_get_child(n1, "component");
         n2 != NULL;
         n2 = charm_xml_node_get_next_child(n2, n1, "component")) {
        charm_xml_node_attr_int(n2, "id", &id);
        charm_xml_node_attr_dbl(n2, "concentration", &c);
        if (charm_comp_index_find_by_id(ctx, id, &idx)) {
            reg->c[idx] = c;
        }
        else {
            CHARM_LERRORF("Unknown component id %d for region '%s' in file 'task.xml'\n", id, reg->name);
            charm_abort(NULL, 1);
        }
    }

    c = 0;
    for (i = 0; i < CHARM_MAX_COMPONETS_COUNT; i++) {
        c += reg->c[i];
    }
    if (fabs(c)-1. > CHARM_EPS) {
        CHARM_LERRORF("Sum of concentrations for region '%s' is not equal to 1 in file 'task.xml'\n", reg->name);
        charm_abort(NULL, 1);
    }
}


static void charm_init_reg(charm_ctx_t *ctx, mxml_node_t *node)
{
    charm_reg_t *reg;
    mxml_node_t *node1;
    int i;

    ctx->reg = sc_array_new(sizeof(charm_reg_t));

    for (node1 = charm_xml_node_get_child(node, "region");
         node1 != NULL;
         node1 = charm_xml_node_get_next_child(node1, node, "region")) {
        reg = (charm_reg_t *) sc_array_push(ctx->reg);
        charm_init_fetch_reg(ctx, node1, reg);
    }
}


static void charm_init_mesh_info(charm_ctx_t *ctx, mxml_node_t *node)
{
    charm_mesh_info_t *m;// = ctx->msh;
    char str[128];

    m = (charm_mesh_info_t*)malloc(sizeof(charm_mesh_info_t));

    charm_xml_node_child_param_str(node, "name", m->filename);
    charm_xml_node_child_param_str(node, "filesType", str);
    m->type = charm_mesh_get_type_by_str(str);
    ctx->msh = m;
}


void charm_init_context(charm_ctx_t *ctx)
{

    FILE *fp;
    mxml_node_t *tree;
    char str[64];

    fp = fopen("task.xml", "r");
    tree = mxmlLoadFile(NULL, fp,
                        MXML_TEXT_CALLBACK);
    fclose(fp);

    CHARM_GLOBAL_ESSENTIAL("Reading file 'task.xml'.\n");

    mxml_node_t *node_task, *node, *node1;
    double x;

    node_task = charm_xml_node_get_child(tree, "task");
    node = charm_xml_node_get_child(node_task, "control");

    charm_xml_node_child_param_str(node, "FLUX_TYPE", str);
    if (strcmp(str, "LF") == 0) {
        ctx->flux_fn = charm_calc_flux_lf;
    }
    else if (strcmp(str, "GODUNOV") == 0) {
        ctx->flux_fn = charm_calc_flux_godunov;
    }
    else if (strcmp(str, "HLLC") == 0) {
        ctx->flux_fn = charm_calc_flux_hllc;
    }
    else if (strcmp(str, "CD") == 0) {
        ctx->flux_fn = charm_calc_flux_cd;
    }
    else {
        CHARM_LERRORF("Unknown flux type '%s'. Use: LF, GODUNOV.\n", str);
        charm_abort(NULL, 1);
    }

    charm_xml_node_child_param_str(node, "LIMITER", str);
    if (strcmp(str, "NONE") == 0) {
        ctx->lim_fn = NULL;
    }
    else if (strcmp(str, "BJ") == 0) {
        ctx->lim_fn = charm_limiter_bj;
    }
    else {
        CHARM_LERRORF("Unknown limiter type '%s'. Use: NONE, BJ.\n", str);
        charm_abort(NULL, 1);
    }

    charm_xml_node_child_param_dbl(node, "MAX_ERROR", &(ctx->max_err));
    charm_xml_node_child_param_int(node, "REFINE_PERIOD", &(ctx->refine_period));
    charm_xml_node_child_param_int(node, "REPARTITION_PERIOD", &(ctx->repartition_period));
    charm_xml_node_child_param_int(node, "MIN_LEVEL", &(ctx->min_level));
    charm_xml_node_child_param_int(node, "MAX_LEVEL", &(ctx->max_level));
    charm_xml_node_child_param_int(node, "FILE_OUTPUT_STEP", &(ctx->write_period));
    charm_xml_node_child_param_int(node, "LOG_OUTPUT_STEP", &(ctx->log_period));

    charm_xml_node_child_param_dbl(node, "TAU", &(ctx->dt));
    charm_xml_node_child_param_dbl(node, "CFL", &(ctx->CFL));
    charm_xml_node_child_param_dbl(node, "TMAX", &(ctx->time));

    charm_xml_node_child_param_str(node, "MODEL", str);
    if (strcmp(str, "EULER") == 0) {
        ctx->model = CHARM_MODEL_EULER;
        charm_model_euler_init(ctx, charm_xml_node_get_child(node, "MODEL"));
    }
    else if (strcmp(str, "NS") == 0) {
        ctx->model = CHARM_MODEL_NS;
        charm_model_ns_init(ctx, charm_xml_node_get_child(node, "MODEL"));
    }
    else if (strcmp(str, "NS_LOW_MACH") == 0) {
        ctx->model = CHARM_MODEL_NS_LOW_MACH;
        charm_model_ns_low_mach_init(ctx, charm_xml_node_get_child(node, "MODEL"));
    }
    else {
        CHARM_LERRORF("Unknown model type '%s'. Use: EULER.\n", str);
        charm_abort(NULL, 1);
    }



    charm_init_bnd(       ctx, charm_xml_node_get_child(node_task, "boundaries"));
    charm_init_comps(     ctx, charm_xml_node_get_child(node_task, "components"));
    charm_init_mat(       ctx, charm_xml_node_get_child(node_task, "materials"));
    charm_init_reg(       ctx, charm_xml_node_get_child(node_task, "regions"));
    charm_init_mesh_info( ctx, charm_xml_node_get_child(node_task, "mesh"));

    ctx->timestep = 0;

}


