//
// Created by zhrv on 26.10.17.
//

#include "charm_init.h"
#include "charm_geom.h"
#include "charm_xml.h"
#include "charm_bnd_cond.h"


void charm_init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    charm_data_t       *data    = (charm_data_t *) q->p.user_data;
    charm_param_t      *par     = &data->par;
    charm_tree_attr_t  *attr;
    charm_reg_t        *reg;

    charm_geom_quad_calc(p4est, q, which_tree);

    attr = charm_get_tree_attr(p4est, which_tree);
    reg = attr->reg;

    par->p.p = reg->p;
    par->p.t = reg->t;
    par->p.u = reg->v[0];
    par->p.v = reg->v[1];
    par->p.w = reg->v[2];
    charm_mat_eos(reg->mat, par, 2);
    charm_mat_eos(reg->mat, par, 1);
    par->p.e_tot = par->p.e + 0.5*(par->p.u*par->p.u+par->p.v*par->p.v+par->p.w*par->p.w);
    charm_param_prim_to_cons(reg->mat, par);
}


void charm_init_fetch_bnd(mxml_node_t* node, charm_bnd_t *bnd)
{
    mxml_node_t * n1, *n2, *n3;
    char str[64];
    charm_xml_node_attr_int(node, "faceType", &(bnd->face_type));
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
            bnd->params = P4EST_ALLOC(double, 5);
            charm_xml_node_child_param_dbl(n2, "Vx", &(bnd->params[0]));
            charm_xml_node_child_param_dbl(n2, "Vy", &(bnd->params[1]));
            charm_xml_node_child_param_dbl(n2, "Vz", &(bnd->params[2]));
            charm_xml_node_child_param_dbl(n2, "T", &(bnd->params[3]));
            charm_xml_node_child_param_dbl(n2, "P", &(bnd->params[4]));
            break;
        case BOUND_OUTLET:
            bnd->bnd_fn = charm_bnd_cond_fn_outlet;
            break;
        case BOUND_WALL_SLIP: // @todo
            bnd->bnd_fn = charm_bnd_cond_fn_wall;
            break;
        case BOUND_WALL_NO_SLIP: // @todo
            bnd->bnd_fn = charm_bnd_cond_fn_wall;
            break;
        default:
            P4EST_LERRORF("Unknown boundary type %d\n", bnd->type);

    }


}


void charm_init_bnd(charm_ctx_t *ctx, mxml_node_t *node)
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


void charm_init_fetch_mat(mxml_node_t* node, charm_mat_t *mat)
{
    mxml_node_t * n1;
    charm_xml_node_attr_int(node, "id", &(mat->id));
    n1 = charm_xml_node_get_child(node, "name");
    charm_xml_node_value_str(n1, mat->name);
    n1 = charm_xml_node_get_child(node, "parameters");
    charm_xml_node_child_param_dbl(n1, "M", &(mat->m));
    charm_xml_node_child_param_dbl(n1, "Cp", &(mat->cp));
    charm_xml_node_child_param_dbl(n1, "ML", &(mat->ml));
    charm_xml_node_child_param_dbl(n1, "Lambda", &(mat->lambda));
    charm_xml_node_child_param_dbl(n1, "K", &(mat->k));
}


void charm_init_mat(charm_ctx_t *ctx, mxml_node_t *node)
{
    charm_mat_t *mat;
    mxml_node_t *node1;
    int i;

    ctx->mat = sc_array_new(sizeof(charm_mat_t));

    for (node1 = charm_xml_node_get_child(node, "material");
         node1 != NULL;
         node1 = charm_xml_node_get_next_child(node1, node, "material")) {
        mat = (charm_mat_t *) sc_array_push(ctx->mat);
        charm_init_fetch_mat(node1, mat);
    }
}


void charm_init_fetch_reg(charm_ctx_t* ctx, mxml_node_t* node, charm_reg_t *reg)
{
    mxml_node_t * n1;
    int tmp;
    charm_mat_t *mat;
    charm_xml_node_attr_int(node, "id", &(reg->id));
    n1 = charm_xml_node_get_child(node, "name");
    charm_xml_node_value_str(n1, reg->name);
    charm_xml_node_child_param_int(node, "cell_type", &(reg->cell_type));
    charm_xml_node_child_param_int(node, "material_id", &tmp);
    mat = charm_mat_find_by_id(ctx, tmp);
    reg->mat = mat;
    n1 = charm_xml_node_get_child(node, "parameters");
    charm_xml_node_child_param_dbl(n1, "Vx", &(reg->v[0]));
    charm_xml_node_child_param_dbl(n1, "Vy", &(reg->v[1]));
    charm_xml_node_child_param_dbl(n1, "Vz", &(reg->v[2]));
    charm_xml_node_child_param_dbl(n1, "T", &(reg->t));
    charm_xml_node_child_param_dbl(n1, "P", &(reg->p));
}


void charm_init_reg(charm_ctx_t *ctx, mxml_node_t *node)
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


void charm_init_mesh_info(charm_ctx_t *ctx, mxml_node_t *node)
{
    charm_mesh_info_t *m = ctx->msh;
    char str[128];

    charm_xml_node_child_param_str(node, "name", m->filename);
    charm_xml_node_child_param_str(node, "filesType", str);
    m->type = charm_mesh_get_type_by_str(str);
}


void charm_init_context(charm_ctx_t *ctx)
{

    FILE *fp;
    mxml_node_t *tree;

    fp = fopen("task.xml", "r");
    tree = mxmlLoadFile(NULL, fp,
                        MXML_TEXT_CALLBACK);
    fclose(fp);

    mxml_node_t *node_task, *node, *node1;
    double x;

    node_task = charm_xml_node_get_child(tree, "task");
    node = charm_xml_node_get_child(node_task, "control");

    charm_xml_node_child_param_dbl(node, "MAX_ERROR", &(ctx->max_err));
    charm_xml_node_child_param_int(node, "REFINE_PERIOD", &(ctx->refine_period));
    charm_xml_node_child_param_int(node, "REPARTITION_PERIOD", &(ctx->repartition_period));
    charm_xml_node_child_param_int(node, "MIN_LEVEL", &(ctx->min_level));
    charm_xml_node_child_param_int(node, "MAX_LEVEL", &(ctx->max_level));
    charm_xml_node_child_param_int(node, "FILE_OUTPUT_STEP", &(ctx->write_period));

    charm_xml_node_child_param_dbl(node, "TAU", &(ctx->dt));
    charm_xml_node_child_param_dbl(node, "TMAX", &(ctx->time));

    charm_init_bnd(ctx, charm_xml_node_get_child(node_task, "boundaries"));

    charm_init_mat(ctx, charm_xml_node_get_child(node_task, "materials"));

    charm_init_reg(ctx, charm_xml_node_get_child(node_task, "regions"));

    charm_init_mesh_info(ctx, charm_xml_node_get_child(node_task, "mesh"));
}


