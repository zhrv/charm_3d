//
// Created by zhrv on 27.02.26.
//

#include "charm_globals.h"





static void _result_to_rhs_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_rhs, data->par.model.ns_jfnk.c_result, c_count);
    charm_fields_mult(data->par.model.ns_jfnk.c_rhs, -1., c_count);
}


static void _rhs_norm2_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_real_t   *err2 = (charm_real_t*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t          c_count = ctx->comp->elem_count;
    *err2 += charm_fields_get_norm2(data->par.model.ns_jfnk.c_rhs, c_count);
}


static void _delta_err2_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_real_t   *err2 = (charm_real_t*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t          c_count = ctx->comp->elem_count;
    charm_fields_t  diff;

    charm_fields_copy(diff, data->par.model.ns_jfnk.c_delta, c_count);
    charm_fields_sub(diff, data->par.model.ns_jfnk.c_delta_old, c_count);

    *err2 += charm_fields_get_norm2(diff, c_count);
}


static void _copy_delta_to_old_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_delta_old, data->par.model.ns_jfnk.c_delta, c_count);
}


static void _calc_ju_copy_result_add_arg_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_jtmp1, data->par.model.ns_jfnk.c_result, c_count);
    charm_fields_axpy(data->par.model.ns_jfnk.c_delta, data->par.c, ctx->model.ns_jfnk.newton.j_eps, c_count);
}


static void _calc_ju_result_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_jtmp1, data->par.model.ns_jfnk.c_result, c_count);
    charm_fields_axpy(data->par.model.ns_jfnk.c_delta, data->par.c, ctx->model.ns_jfnk.newton.j_eps, c_count);
}


/*static void _calc_r_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_ju, c_count);
    charm_fields_sub(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_rhs, c_count);
}*/


static void _calc_jr_r_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;
    charm_real_t       *res = (charm_real_t*)user_data;
    
    *res += charm_fields_dot(data->par.model.ns_jfnk.c_jr, data->par.model.ns_jfnk.c_residual, c_count);
}


static void _calc_jr_jr_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;
    charm_real_t       *res = (charm_real_t*)user_data;
    
    *res += charm_fields_dot(data->par.model.ns_jfnk.c_jr, data->par.model.ns_jfnk.c_jr, c_count);
}


static void _update_delta_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;
    charm_real_t       *tau = (charm_real_t*)user_data;
    
    charm_fields_axpy(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_delta, *tau, c_count);
}


static void _update_fld_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_fields_axpy(data->par.model.ns_jfnk.c_delta, data->par.c, ctx->model.ns_jfnk.newton.relax, c_count);
}


static void _calc_err2_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_real_t   *err2 = (charm_real_t*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t          c_count = ctx->comp->elem_count;
    *err2 += charm_fields_get_norm2(data->par.model.ns_jfnk.c_delta, c_count);
}


