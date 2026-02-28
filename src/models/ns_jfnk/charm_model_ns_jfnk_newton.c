//
// Created by zhrv on 10.01.26.
//

#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"

void charm_model_ns_jfnk_dg_operator(p4est_t * p4est, p4est_ghost_t * _ghost, charm_data_t * _ghost_data);
void charm_model_ns_jfnk_dg_operator_stash_push(p4est_t * p4est);
void charm_model_ns_jfnk_dg_operator_stash_pop(p4est_t * p4est);


CHARM_DECL_QUAD_ITER_STATIC(_result_to_rhs, {
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_rhs, data->par.model.ns_jfnk.c_result, c_count);
    charm_fields_mult(data->par.model.ns_jfnk.c_rhs, -1., c_count);
    charm_fields_copy(data->par.model.ns_jfnk.c_dg_res, data->par.model.ns_jfnk.c_result, c_count);//???????
})


CHARM_DECL_FIELDS_DOT_FUNC(_rhs_norm2, data->par.model.ns_jfnk.c_rhs, data->par.model.ns_jfnk.c_rhs)

static void _calc_rhs(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = charm_get_ctx(p4est);
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    
    _result_to_rhs(p4est);
    ctx->tmp.ns_jfnk.rhs_norm = sqrt(_rhs_norm2(p4est));
}


CHARM_DECL_FIELDS_AXPY_FUNC(_calc_ju_add_arg, data->par.model.ns_jfnk.c_delta, data->par.c)

CHARM_DECL_QUAD_ITER_STATIC(_calc_ju_result, {
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_ju, data->par.model.ns_jfnk.c_result, c_count);
    charm_fields_axpy(data->par.model.ns_jfnk.c_dg_res, data->par.model.ns_jfnk.c_ju, -1., c_count);
    charm_fields_mult(data->par.model.ns_jfnk.c_ju, 1./ctx->model.ns_jfnk.newton.j_eps, c_count);
})


static inline void _calc_ju(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = charm_get_ctx(p4est);

    charm_model_ns_jfnk_dg_operator_stash_push(p4est);
    _calc_ju_add_arg(p4est, ctx->model.ns_jfnk.newton.j_eps);
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    _calc_ju_result(p4est);
    charm_model_ns_jfnk_dg_operator_stash_pop(p4est);        
}


CHARM_DECL_QUAD_ITER_STATIC(_calc_r, {
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_ju, c_count);
    charm_fields_sub(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_rhs, c_count);
})


CHARM_DECL_FIELDS_AXPY_FUNC(_calc_jr_add_arg, data->par.model.ns_jfnk.c_residual, data->par.c)

CHARM_DECL_QUAD_ITER_STATIC(_calc_jr_result, {
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_jr, data->par.model.ns_jfnk.c_result, c_count);
    charm_fields_axpy(data->par.model.ns_jfnk.c_dg_res, data->par.model.ns_jfnk.c_jr, -1., c_count);
    charm_fields_mult(data->par.model.ns_jfnk.c_jr, 1./ctx->model.ns_jfnk.newton.j_eps, c_count);
})


static inline void _calc_jr(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = charm_get_ctx(p4est);
    
    charm_model_ns_jfnk_dg_operator_stash_push(p4est);
    _calc_jr_add_arg(p4est, ctx->model.ns_jfnk.newton.j_eps);
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    _calc_jr_result(p4est);
    charm_model_ns_jfnk_dg_operator_stash_pop(p4est);        
}

static inline void _calc_ju_r_jr(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    // charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    // _copy_dg_result(p4est); // ???????????????????

    _calc_ju(p4est, ghost, ghost_data);
    _calc_r(p4est);
    _calc_jr(p4est, ghost, ghost_data);
}


CHARM_DECL_FIELDS_DOT_FUNC(_calc_jr_r,  data->par.model.ns_jfnk.c_jr, data->par.model.ns_jfnk.c_residual)
CHARM_DECL_FIELDS_DOT_FUNC(_calc_jr_jr, data->par.model.ns_jfnk.c_jr, data->par.model.ns_jfnk.c_jr)


static inline charm_real_t _calc_tau(p4est_t * p4est) 
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t        jr_r = _calc_jr_r(p4est);
    charm_real_t        jr_jr = _calc_jr_jr(p4est);
    
    return -jr_r/jr_jr;
}


CHARM_DECL_QUAD_ITER_STATIC(_zero_delta, {
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_zero(data->par.model.ns_jfnk.c_delta, c_count);
})

CHARM_DECL_QUAD_ITER_STATIC(_copy_delta_to_old, {
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_delta_old, data->par.model.ns_jfnk.c_delta, c_count);
})


CHARM_DECL_FIELDS_AXPY_FUNC(_update_delta, data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_delta)
CHARM_DECL_FIELDS_AXPY_FUNC(_delta_diff, data->par.model.ns_jfnk.c_delta, data->par.model.ns_jfnk.c_delta_old)
CHARM_DECL_FIELDS_DOT_FUNC(_delta_err2, data->par.model.ns_jfnk.c_delta_old, data->par.model.ns_jfnk.c_delta_old)


static void _calc_delta(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_real_t        tau_k, err;
    charm_int_t         km_step, km_stop;
    charm_ctx_t        *ctx = charm_get_ctx(p4est);

    km_stop = 0;
    km_step = 0;
    _zero_delta(p4est);
    while (!km_stop) {
        _copy_delta_to_old(p4est);
        _calc_ju_r_jr(p4est, ghost, ghost_data);

        tau_k = _calc_tau(p4est);
        _update_delta(p4est, tau_k);
        _delta_diff(p4est, -1.0);
        err = sqrt(_delta_err2(p4est));
        if (err < ctx->model.ns_jfnk.newton.solver_rtol*ctx->tmp.ns_jfnk.rhs_norm 
                                    || km_step >= ctx->model.ns_jfnk.newton.solver_max_step) {
            km_stop = 1;
        }
        km_step++;
    }
}

CHARM_DECL_FIELDS_AXPY_FUNC(_update_fld, data->par.model.ns_jfnk.c_delta, data->par.c)
CHARM_DECL_FIELDS_DOT_FUNC(_calc_err2, data->par.model.ns_jfnk.c_delta, data->par.model.ns_jfnk.c_delta)

charm_int_t charm_model_ns_jfnk_newton_step(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{    
    charm_ctx_t        *ctx = (charm_ctx_t*)p4est->user_pointer;
    charm_int_t         nm_stop = 0;
    charm_real_t        nm_err;

    _calc_rhs(p4est, ghost, ghost_data);
    _calc_delta(p4est, ghost, ghost_data);
    _update_fld(p4est, ctx->model.ns_jfnk.newton.relax);
    nm_err = sqrt(_calc_err2(p4est));
    if (nm_err < ctx->model.ns_jfnk.newton.rtol * ctx->tmp.ns_jfnk.fld_old_norm) {
        nm_stop = 1;
    }
    return nm_stop;
}
