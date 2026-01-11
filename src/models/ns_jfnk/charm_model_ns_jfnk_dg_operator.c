//
// Created by zhrv on 10.01.26.
//

#include <p8est_iterate.h>
#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"


void charm_model_ns_jfnk_dg_operator_conv(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_jfnk_dg_operator_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void charm_model_ns_jfnk_dg_operator_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_real_t              dt = *((charm_real_t *) user_data);
    charm_vect_t              rhs_ru;
    charm_vect_t              rhs_rv;
    charm_vect_t              rhs_rw;
    charm_vect_t              rhs_re;
    charm_vect_t              rhs_rc[CHARM_MAX_COMPONETS_COUNT];
    size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_matr_vect_mult(data->par.g.a_inv, data->int_ru, rhs_ru);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rv, rhs_rv);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rw, rhs_rw);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_re, rhs_re);

    for (j = 0; j < c_count; j++) {
        charm_matr_vect_mult(data->par.g.a_inv, data->int_rc[j], rhs_rc[j]);
    }

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ru[i] -= _NORM_(dt * rhs_ru[i]);
        data->par.c.rv[i] -= _NORM_(dt * rhs_rv[i]);
        data->par.c.rw[i] -= _NORM_(dt * rhs_rw[i]);
        data->par.c.re[i] -= _NORM_(dt * rhs_re[i]);
        for (j = 0; j < c_count; j++) {
            data->par.c.rc[j][i] -= _NORM_(dt * rhs_rc[j][i]);
        }
    }
}


static void charm_model_ns_jfnk_dg_operator_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i, j;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->int_ru[i] = 0.;
        data->int_rv[i] = 0.;
        data->int_rw[i] = 0.;
        data->int_re[i] = 0.;
        for (j = 0; j < CHARM_MAX_COMPONETS_COUNT; j++) {
            data->int_rc[j][i] = 0.;
        }
    }
}


void charm_model_ns_jfnk_dg_operator(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data) 
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    p4est_ghost_t      *ghost       = *_ghost;
    charm_data_t       *ghost_data  = *_ghost_data;

    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           
    p4est_iterate (p4est, ghost, (void *) ghost_data,               
            charm_model_ns_jfnk_dg_operator_zero_quad_iter_fn,         
            NULL, NULL, NULL);                                      
    charm_model_ns_jfnk_dg_operator_diff(p4est, ghost, ghost_data);    
    charm_model_ns_jfnk_dg_operator_conv(p4est, ghost, ghost_data);    
    p4est_iterate (p4est, NULL, (void *) dt,                        
            charm_model_ns_jfnk_dg_operator_update_quad_iter_fn,       
            NULL, NULL, NULL);                                      
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           
    charm_limiter(p4est, ghost, ghost_data);
}
