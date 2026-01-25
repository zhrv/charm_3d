//
// Created by zhrv on 10.01.26.
//

#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"


// charm_real_t charm_model_ns_jfnk_newton_minres_step(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void _charm_model_ns_jfnk_copy_result_to_rhs_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_AXPY(data->par.model.ns_jfnk.c_result, data->par.model.ns_jfnk.c_rhs, -1.);
}


void charm_model_ns_jfnk_newton_rhs(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    int i;
    charm_real_t err;

    charm_model_ns_jfnk_dg_operator(p4est, dt, &ghost, &ghost_data);
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_jfnk_copy_result_to_rhs_quad_iter_fn, NULL, NULL, NULL);
}




static void _charm_model_ns_jfnk_copy_result_to_Jdelta_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_SET_SCALAR(data->par.model.ns_jfnk.c_Jdelta, 0.);
    FIELDS_AXPY(data->par.model.ns_jfnk.c_result, data->par.model.ns_jfnk.c_Jdelta, -1.);
}

static void _charm_model_ns_jfnk_add_delta_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_ctx_t    *ctx = (charm_ctx_t *) info->p4est->user_pointer;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_AXPY(data->par.model.ns_jfnk.c_delta, data->par.c, ctx->model.ns_jfnk.j_eps);
}

static void _charm_model_ns_jfnk_calc_Jdelta_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_ctx_t    *ctx = (charm_ctx_t *) info->p4est->user_pointer;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_AXPY(data->par.model.ns_jfnk.c_result, data->par.model.ns_jfnk.c_Jdelta, 1.);
    FIELDS_MULT(data->par.model.ns_jfnk.c_Jdelta, 1./ctx->model.ns_jfnk.j_eps)
}

void charm_model_ns_jfnk_newton_Jdelta(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_model_ns_jfnk_dg_operator(p4est, dt, &ghost, &ghost_data);
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_jfnk_copy_result_to_Jdelta_quad_iter_fn, NULL, NULL, NULL);
    charm_model_ns_jfnk_stash_push(p4est);
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_jfnk_add_delta_quad_iter_fn, NULL, NULL, NULL);
    charm_model_ns_jfnk_dg_operator(p4est, dt, &ghost, &ghost_data);
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_jfnk_calc_Jdelta_quad_iter_fn, NULL, NULL, NULL);
    charm_model_ns_jfnk_stash_pop(p4est);
}







static void _charm_model_ns_jfnk_copy_result_to_Jr_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_SET_SCALAR(data->par.model.ns_jfnk.c_Jr, 0.);
    FIELDS_AXPY(data->par.model.ns_jfnk.c_result, data->par.model.ns_jfnk.c_Jr, -1.);
}

static void _charm_model_ns_jfnk_add_r_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_ctx_t    *ctx = (charm_ctx_t *) info->p4est->user_pointer;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_AXPY(data->par.model.ns_jfnk.c_residual, data->par.c, ctx->model.ns_jfnk.j_eps);
}

static void _charm_model_ns_jfnk_calc_Jr_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_ctx_t    *ctx = (charm_ctx_t *) info->p4est->user_pointer;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_AXPY(data->par.model.ns_jfnk.c_result, data->par.model.ns_jfnk.c_Jr, 1.);
    FIELDS_MULT(data->par.model.ns_jfnk.c_Jr, 1./ctx->model.ns_jfnk.j_eps)
}

void charm_model_ns_jfnk_newton_Jr(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_model_ns_jfnk_dg_operator(p4est, dt, &ghost, &ghost_data);
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_jfnk_copy_result_to_Jr_quad_iter_fn, NULL, NULL, NULL);
    charm_model_ns_jfnk_stash_push(p4est);
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_jfnk_add_r_quad_iter_fn, NULL, NULL, NULL);
    charm_model_ns_jfnk_dg_operator(p4est, dt, &ghost, &ghost_data);
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_jfnk_calc_Jr_quad_iter_fn, NULL, NULL, NULL);
    charm_model_ns_jfnk_stash_pop(p4est);
}



static void _charm_model_ns_jfnk_residual_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_ctx_t    *ctx = (charm_ctx_t *) info->p4est->user_pointer;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_COPY(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_Jdelta);
    FIELDS_AXPY(data->par.model.ns_jfnk.c_rhs, data->par.model.ns_jfnk.c_residual, -1.);
}

void charm_model_ns_jfnk_newton_residual(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    p4est_iterate(p4est, NULL, NULL, _charm_model_ns_jfnk_residual_quad_iter_fn, NULL, NULL, NULL);
}



static void _charm_model_ns_jfnk_newton_zero_delta_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_ctx_t    *ctx = (charm_ctx_t *) info->p4est->user_pointer;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_SET_SCALAR(data->par.model.ns_jfnk.c_delta, 0.);
}

void charm_model_ns_jfnk_newton_zero_delta(p4est_t * p4est)
{
    p4est_iterate(p4est, NULL, NULL, _charm_model_ns_jfnk_newton_zero_delta_quad_iter_fn, NULL, NULL, NULL);
}

static void _charm_model_ns_jfnk_newton_minres_new_delta_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_ctx_t    *ctx = (charm_ctx_t *) info->p4est->user_pointer;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_real_t    tau = *((charm_real_t *) user_data);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_AXPY(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_delta, tau);
}



charm_real_t charm_model_ns_jfnk_newton_minres_step(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_real_t Jrr, JrJr, tau;
    // Ju = rhs

    charm_model_ns_jfnk_newton_Jdelta(p4est, dt, ghost, ghost_data);
    charm_model_ns_jfnk_newton_residual(p4est, dt, ghost, ghost_data);
    charm_model_ns_jfnk_newton_Jr(p4est, dt, ghost, ghost_data);

    Jrr = charm_model_ns_jfnk_newton_Jrr(p4est, dt, ghost, ghost_data);
    JrJr = charm_model_ns_jfnk_newton_JrJr(p4est, dt, ghost, ghost_data);

    tau = Jrr/JrJr;

    p4est_iterate(p4est, NULL, (void*)(&tau), _charm_model_ns_jfnk_newton_minres_new_delta_quad_iter_fn, NULL, NULL, NULL);
    // residual = Ju - rhs

    // Ju

    // 
}




static void _charm_model_ns_jfnk_newton_step_new_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_ctx_t    *ctx = (charm_ctx_t *) info->p4est->user_pointer;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    size_t          c_count = charm_get_comp_count(info->p4est);
    charm_int_t     i, j;

    FIELDS_AXPY(data->par.model.ns_jfnk.c_delta, data->par.c, 1.);
}

void charm_model_ns_jfnk_newton_step(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    int steps;
    charm_real_t err;
    steps = 0;
    charm_model_ns_jfnk_newton_zero_delta(p4est);
    while (1) {
        err = charm_model_ns_jfnk_newton_minres_step(p4est, dt, ghost, ghost_data);
        steps++;
        if (err < ctx->model.ns_jfnk.newton_tol || steps >= ctx->model.ns_jfnk.linsol_max_steps) break;
    }
    p4est_iterate(p4est, NULL, NULL, _charm_model_ns_jfnk_newton_step_new_quad_iter_fn, NULL, NULL, NULL);
}
