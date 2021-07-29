//
// Created by zhrv on 14.10.2020.
//
#include "charm_globals.h"


void charm_model_ns_init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx     = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *data    = (charm_data_t *) q->p.user_data;
    charm_param_t      *par     = &data->par;
    charm_prim_t   prim;
    charm_cons_t        cons;
    charm_tree_attr_t  *attr;
    charm_reg_t        *reg;
    int                 i;
    size_t              c_count = charm_get_comp_count(p4est);
    charm_mat_t        *mat;

    charm_geom_quad_calc(p4est, q, which_tree);

    attr = charm_get_tree_attr(p4est, which_tree);
    reg = attr->reg;
    for (int i = 0; i < c_count; i++) {
        data->par.model.ns.d[i] = 0.;
    }

    if (ctx->model.ns.turb.init_cond_fn) ctx->model.ns.turb.init_cond_fn(p4est, which_tree, q);
}


