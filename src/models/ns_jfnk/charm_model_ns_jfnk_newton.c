//
// Created by zhrv on 10.01.26.
//

#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"


void charm_model_ns_jfnk_dg_operator(p4est_t * p4est, charm_real_t dt, p4est_ghost_t * _ghost, charm_data_t * _ghost_data);


static void _result_to_rhs_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_rhs, data->par.model.ns_jfnk.c_result, c_count);
    charm_fields_mult(data->par.model.ns_jfnk.c_rhs, -1., c_count);
}


static void _calc_rhs(p4est_t * p4est, charm_real_t dt, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_model_ns_jfnk_dg_operator(p4est, dt, ghost, ghost_data);
    
    p4est_iterate (p4est, NULL, NULL, 
        _result_to_rhs_quad_iter_fn, 
        NULL, NULL, NULL);
}


static void _calc_delta(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    // TODO 
    
}


static void _update_fld_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_fields_axpy(data->par.model.ns_jfnk.c_delta, data->par.c, ctx->model.ns_jfnk.newton.relax, c_count);
}


static void _update_fld(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    p4est_iterate (p4est, NULL, NULL, 
        _update_fld_quad_iter_fn, 
        NULL, NULL, NULL);
        
}


static void _calc_err2_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_real_t   *err2 = (charm_real_t*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t          c_count = ctx->comp->elem_count;
    *err2 += charm_fields_get_norm2(data->par.model.ns_jfnk.c_delta, c_count);
}

static charm_real_t _calc_err2 (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t        loc_err2, glob_err2;
    int                 mpiret, i;

    loc_err2 = 0.0;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_err2,
                   _calc_err2_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_err2, &glob_err2, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return sqrt(glob_err2);
}


charm_int_t charm_model_ns_jfnk_newton_step(p4est_t * p4est, charm_real_t dt, charm_real_t fld_old_norm, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{    
    charm_ctx_t        *ctx = (charm_ctx_t*)p4est->user_pointer;
    charm_int_t         nm_stop = 0;
    charm_real_t        nm_err;

    _calc_rhs(p4est, dt, ghost, ghost_data);
    _calc_delta(p4est, ghost, ghost_data);
    _update_fld(p4est, ghost, ghost_data);
    nm_err = _calc_err2(p4est);
    if (nm_err < ctx->model.ns_jfnk.newton.rtol*fld_old_norm) {
        nm_stop = 1;
    }
    return nm_stop;
}
