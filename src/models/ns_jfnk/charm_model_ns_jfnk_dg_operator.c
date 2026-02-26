//
// Created by zhrv on 10.01.26.
//

#include <p8est_iterate.h>
#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"


void charm_model_ns_jfnk_dg_operator_conv(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_jfnk_dg_operator_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void charm_model_ns_jfnk_dg_operator_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_fields_zero(data->integrals, c_count);
}


static void charm_model_ns_jfnk_dg_operator_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;
    charm_real_t        dt = *((charm_real_t *) user_data);
    charm_fields_t      rhs;

    charm_fields_zero(rhs, c_count);
    charm_fields_add(rhs, data->par.c, c_count);
    charm_fields_sub(rhs, data->par.c_old, c_count);
    charm_fields_mult(rhs, 1./dt, c_count);
    charm_matr_fields_mult(data->par.g.a, rhs, data->par.model.ns_jfnk.c_result, c_count);    
    charm_fields_add(data->par.model.ns_jfnk.c_result, data->integrals, c_count);
}


void charm_model_ns_jfnk_dg_operator(p4est_t * p4est, charm_real_t dt, p4est_ghost_t * ghost, charm_data_t * ghost_data) 
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    // p4est_ghost_t      *ghost       = *_ghost;
    // charm_data_t       *ghost_data  = *_ghost_data;

    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           
    p4est_iterate (p4est, ghost, (void *) ghost_data,               
            charm_model_ns_jfnk_dg_operator_zero_quad_iter_fn,         
            NULL, NULL, NULL);                                      
    charm_model_ns_jfnk_dg_operator_diff(p4est, ghost, ghost_data);    
    charm_model_ns_jfnk_dg_operator_conv(p4est, ghost, ghost_data);    
    p4est_iterate (p4est, NULL, (void *) &dt,                        
            charm_model_ns_jfnk_dg_operator_update_quad_iter_fn,       
            NULL, NULL, NULL);                                      
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           
    //charm_limiter(p4est, ghost, ghost_data);
}
