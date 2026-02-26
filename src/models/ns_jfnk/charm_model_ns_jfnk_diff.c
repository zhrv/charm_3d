//
// Created by zhrv on 10.01.26.
//

#include <p8est_iterate.h>
#include <charm_globals.h>
#include "charm_base_func.h"

void charm_model_ns_jfnk_dg_operator_diff_grad(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_jfnk_dg_operator_diff_integrals(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);

charm_real_t charm_model_ns_jfnk_get_visc_mu(p4est_t* p4est, charm_real_t *x, charm_data_t* data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    size_t c_count = charm_get_comp_count(p4est);
    charm_comp_t *comp;
    charm_cons_t cons;
    charm_prim_t prim;
    charm_real_t mu, cm, s;
    int i;

    charm_get_fields(data, x, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);
    s  = 0.;
    mu = 0.;
    for (i = 0; i < c_count; i++) {
        comp = charm_get_comp(p4est, i);
        cm = prim.c[i]/comp->m;
        s += cm;
        mu += cm*charm_comp_calc_ml(comp, prim.t);
    }
    mu /= s;

    return mu/s;
}

// charm_real_t charm_model_ns_jfnk_get_turb_mu(p4est_t* p4est, charm_real_t *x, charm_data_t* data)
// {
//     return data->par.model.turb.mu_t;
// }

charm_real_t charm_model_ns_jfnk_get_mu(p4est_t* p4est, charm_real_t *x, charm_data_t* data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    int mu = charm_model_ns_jfnk_get_visc_mu(p4est, x, data);
    // if (ctx->model.turb.model_type != TURB_MODEL_UNKNOWN) {
    //     mu += charm_model_ns_jfnk_get_turb_mu(p4est, x, data);;
    // }

    return mu;
}

charm_real_t charm_model_ns_jfnk_get_lambda(p4est_t* p4est, charm_data_t* data)
{
    return 0;
}




void charm_model_ns_jfnk_dg_operator_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    if (!ctx->model.ns_jfnk.use_visc) return;
    charm_model_ns_jfnk_dg_operator_diff_grad(p4est, ghost, ghost_data);
    charm_model_ns_jfnk_dg_operator_diff_integrals(p4est, ghost, ghost_data);
}
