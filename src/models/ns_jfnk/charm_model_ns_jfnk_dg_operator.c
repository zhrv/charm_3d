//
// Created by zhrv on 10.01.26.
//

#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"


void charm_model_ns_jfnk_dg_operator_conv(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_jfnk_dg_operator_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void charm_model_ns_jfnk_dg_operator_result_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_real_t        dt = *((charm_real_t *) user_data);
    charm_real_t        _dt = 1./dt;
    charm_fields_t      rhs;
    size_t              c_count = charm_get_comp_count(info->p4est);
    int                 i, j;

    charm_matr_fields_mult(data->par.g.a_inv, data->int_r, rhs, c_count);

    FIELDS_COPY(rhs, data->par.c)
    FIELDS_AXPY(data->par.c_old, rhs, -1.);
    FIELDS_MULT(rhs, _dt);
    charm_matr_fields_mult(data->par.g.a, rhs, data->par.model.ns_jfnk.c_result, c_count);
    FIELDS_AXPY(rhs, data->par.model.ns_jfnk.c_result, 1.);
}


static void charm_model_ns_jfnk_dg_operator_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    size_t              c_count = charm_get_comp_count(info->p4est);
    charm_int_t         i, j;

    // FIELDS_COPY(data->par.model.ns_jfnk.c_stash, data->par.c, c_count);
    FIELDS_SET_SCALAR(data->int_r, 0.);
}


/**
 * 
 */
void charm_model_ns_jfnk_dg_operator(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t *ghost, charm_data_t *ghost_data) 
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;

    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           
    p4est_iterate (p4est, ghost, (void *) ghost_data,               
            charm_model_ns_jfnk_dg_operator_zero_quad_iter_fn,         
            NULL, NULL, NULL);                                      
    charm_model_ns_jfnk_dg_operator_diff(p4est, ghost, ghost_data);    
    charm_model_ns_jfnk_dg_operator_conv(p4est, ghost, ghost_data);    
    p4est_iterate (p4est, NULL, (void *) dt,                        
            charm_model_ns_jfnk_dg_operator_result_quad_iter_fn,       
            NULL, NULL, NULL);                                      
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           
    charm_limiter(p4est, ghost, ghost_data);
}
