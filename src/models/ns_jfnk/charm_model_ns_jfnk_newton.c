//
// Created by zhrv on 10.01.26.
//

#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"

void charm_model_ns_jfnk_dg_operator(p4est_t * p4est, p4est_ghost_t * _ghost, charm_data_t * _ghost_data);
void charm_model_ns_jfnk_dg_operator_stash_push(p4est_t * p4est);
void charm_model_ns_jfnk_dg_operator_stash_pop(p4est_t * p4est);


static inline void _result_to_rhs_quad_iter_fn(p4est_iter_volume_info_t *info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(&(data->par.model.ns_jfnk.c_rhs), &(data->par.model.ns_jfnk.c_result), c_count);
    charm_fields_mult(&(data->par.model.ns_jfnk.c_rhs), -1., c_count);
    charm_fields_copy(&(data->par.model.ns_jfnk.c_dg_res), &(data->par.model.ns_jfnk.c_result), c_count);//???????
}

static inline void result_to_rhs(p4est_t * p4est) 
{
    p4est_iterate(p4est, NULL, NULL, _result_to_rhs_quad_iter_fn, NULL, NULL, NULL);
}




static inline void _rhs_norm2_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{                                                                                           
    charm_real_t   *res = (charm_real_t*) user_data;                                       
    charm_data_t   *data = charm_get_quad_data(info->quad);                                 
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;                          
    charm_size_t          c_count = ctx->comp->elem_count;                                        
    *res += charm_fields_dot(&(data->par.model.ns_jfnk.c_rhs), &(data->par.model.ns_jfnk.c_rhs), c_count);              
}                                                                                           
static inline charm_real_t rhs_norm2(p4est_t * p4est)                                                
{                                                                                           
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;                          
    charm_real_t        loc_res, glob_res;                                                  
    int                 mpiret;                                                             
                                                                                            
    loc_res = 0.0;                                                                          
    p4est_iterate (p4est, NULL,                                                             
                (void *) &loc_res,                                                          
                _rhs_norm2_quad_iter_fn,                                              
                NULL, NULL, NULL);                                                          
                                                                                            
    mpiret = sc_MPI_Allreduce (&loc_res, &glob_res, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);                                                                  
                                                                                            
    return glob_res;                                                                        
}



static void calc_rhs(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = charm_get_ctx(p4est);
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    
    result_to_rhs(p4est);
    ctx->tmp.ns_jfnk.rhs_norm = sqrt(rhs_norm2(p4est));
}


static void _calc_ju_add_arg_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t              c_count = ctx->comp->elem_count;
    charm_real_t       *a = (charm_real_t*)user_data;
    charm_fields_axpy(&(data->par.model.ns_jfnk.c_delta), &(data->par.c), *a, c_count);
}
static inline void calc_ju_add_arg(p4est_t * p4est, charm_real_t a) 
{
    p4est_iterate(p4est, NULL, (void*)&a, _calc_ju_add_arg_quad_iter_fn, NULL, NULL, NULL);
}



static inline void _calc_ju_result_quad_iter_fn(p4est_iter_volume_info_t *info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(&(data->par.model.ns_jfnk.c_ju), &(data->par.model.ns_jfnk.c_result), c_count);
    charm_fields_axpy(&(data->par.model.ns_jfnk.c_dg_res), &(data->par.model.ns_jfnk.c_ju), -1., c_count);
    charm_fields_mult(&(data->par.model.ns_jfnk.c_ju), 1./ctx->model.ns_jfnk.newton.j_eps, c_count);
}

static inline void calc_ju_result(p4est_t * p4est) 
{
    p4est_iterate(p4est, NULL, NULL, _calc_ju_result_quad_iter_fn, NULL, NULL, NULL);
}


static inline void calc_ju(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = charm_get_ctx(p4est);

    charm_model_ns_jfnk_dg_operator_stash_push(p4est);
    calc_ju_add_arg(p4est, ctx->model.ns_jfnk.newton.j_eps);
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    calc_ju_result(p4est);
    charm_model_ns_jfnk_dg_operator_stash_pop(p4est);        
}


static inline void _calc_r_quad_iter_fn(p4est_iter_volume_info_t *info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(&(data->par.model.ns_jfnk.c_residual), &(data->par.model.ns_jfnk.c_ju), c_count);
    charm_fields_sub(&(data->par.model.ns_jfnk.c_residual), &(data->par.model.ns_jfnk.c_rhs), c_count);
}

static inline void calc_r(p4est_t * p4est) 
{
    p4est_iterate(p4est, NULL, NULL, _calc_r_quad_iter_fn, NULL, NULL, NULL);
}




static void _calc_jr_add_arg_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data) 
{                                                                                           
    charm_data_t       *data = charm_get_quad_data(info->quad);                             
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;                      
    charm_size_t              c_count = ctx->comp->elem_count;                                    
    charm_real_t       *a = (charm_real_t*)user_data;                                       
    charm_fields_axpy(&(data->par.model.ns_jfnk.c_residual), &(data->par.c), *a, c_count);                                       \
}                                                                                           
static inline void calc_jr_add_arg(p4est_t * p4est, charm_real_t a)                                  
{                                                                                           
    p4est_iterate(p4est, NULL, (void*)&a, _calc_jr_add_arg_quad_iter_fn, NULL, NULL, NULL); 
}



static inline void _calc_jr_result_quad_iter_fn(p4est_iter_volume_info_t *info, void *user_data) 
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(&(data->par.model.ns_jfnk.c_jr), &(data->par.model.ns_jfnk.c_result), c_count);
    charm_fields_axpy(&(data->par.model.ns_jfnk.c_dg_res), &(data->par.model.ns_jfnk.c_jr), -1., c_count);
    charm_fields_mult(&(data->par.model.ns_jfnk.c_jr), 1./ctx->model.ns_jfnk.newton.j_eps, c_count);
}
static inline void calc_jr_result(p4est_t * p4est) 
{
    p4est_iterate(p4est, NULL, NULL, _calc_jr_result_quad_iter_fn, NULL, NULL, NULL);
}


static inline void calc_jr(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t        *ctx = charm_get_ctx(p4est);
    
    charm_model_ns_jfnk_dg_operator_stash_push(p4est);
    calc_jr_add_arg(p4est, ctx->model.ns_jfnk.newton.j_eps);
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    calc_jr_result(p4est);
    charm_model_ns_jfnk_dg_operator_stash_pop(p4est);        
}

static inline void calc_ju_r_jr(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    // charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    // _copy_dg_result(p4est); // ???????????????????

    calc_ju(p4est, ghost, ghost_data);
    calc_r(p4est);
    calc_jr(p4est, ghost, ghost_data);
}


static inline void _calc_jr_r_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{                                                                                           
    charm_real_t   *res = (charm_real_t*) user_data;                                       
    charm_data_t   *data = charm_get_quad_data(info->quad);                                 
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;                          
    charm_size_t          c_count = ctx->comp->elem_count;                                        
    *res += charm_fields_dot(&(data->par.model.ns_jfnk.c_jr), &(data->par.model.ns_jfnk.c_residual), c_count);
}                                                                                           
static inline charm_real_t calc_jr_r(p4est_t * p4est)                                                
{                                                                                           
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;                          
    charm_real_t        loc_res, glob_res;                                                  
    int                 mpiret;                                                             
                                                                                            
    loc_res = 0.0;                                                                          
    p4est_iterate (p4est, NULL,                                                             
                (void *) &loc_res,                                                          
                _calc_jr_r_quad_iter_fn,                                              
                NULL, NULL, NULL);                                                          
                                                                                            
    mpiret = sc_MPI_Allreduce (&loc_res, &glob_res, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);                                                                  
                                                                                            
    return glob_res;                                                                        
}


static inline void _calc_jr_jr_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{                                                                                           
    charm_real_t   *res = (charm_real_t*) user_data;                                       
    charm_data_t   *data = charm_get_quad_data(info->quad);                                 
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;                          
    charm_size_t          c_count = ctx->comp->elem_count;                                        
    *res += charm_fields_dot(&(data->par.model.ns_jfnk.c_jr), &(data->par.model.ns_jfnk.c_jr), c_count);
}                                                                                           
static inline charm_real_t calc_jr_jr(p4est_t * p4est)                                                
{                                                                                           
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;                          
    charm_real_t        loc_res, glob_res;                                                  
    int                 mpiret;                                                             
                                                                                            
    loc_res = 0.0;                                                                          
    p4est_iterate (p4est, NULL,                                                             
                (void *) &loc_res,                                                          
                _calc_jr_jr_quad_iter_fn,                                              
                NULL, NULL, NULL);                                                          
                                                                                            
    mpiret = sc_MPI_Allreduce (&loc_res, &glob_res, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);                                                                  
                                                                                            
    return glob_res;                                                                        
}


static inline charm_real_t calc_tau(p4est_t * p4est) 
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t        jr_r = calc_jr_r(p4est);
    charm_real_t        jr_jr = calc_jr_jr(p4est);
    
    return -jr_r/jr_jr;
}


static inline void _zero_delta_quad_iter_fn(p4est_iter_volume_info_t *info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_zero(&(data->par.model.ns_jfnk.c_delta), c_count);
}
static inline void zero_delta(p4est_t * p4est) 
{
    p4est_iterate(p4est, NULL, NULL, _zero_delta_quad_iter_fn, NULL, NULL, NULL);
}


static inline void _copy_delta_to_old_quad_iter_fn(p4est_iter_volume_info_t *info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t        c_count = ctx->comp->elem_count;

    charm_fields_copy(&(data->par.model.ns_jfnk.c_delta_old), &(data->par.model.ns_jfnk.c_delta), c_count);
}
static inline void copy_delta_to_old(p4est_t * p4est) 
{
    p4est_iterate(p4est, NULL, NULL, _copy_delta_to_old_quad_iter_fn, NULL, NULL, NULL);
}


static void _update_delta_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data) 
{                                                                                           
    charm_data_t       *data = charm_get_quad_data(info->quad);                             
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;                      
    charm_size_t              c_count = ctx->comp->elem_count;                                    
    charm_real_t       *a = (charm_real_t*)user_data;                                       
    charm_fields_axpy(&(data->par.model.ns_jfnk.c_residual), &(data->par.model.ns_jfnk.c_delta), *a, c_count);                                       
}                                                                                           
static inline void update_delta(p4est_t * p4est, charm_real_t a)                                  
{                                                                                           
    p4est_iterate(p4est, NULL, (void*)&a, _update_delta_quad_iter_fn, NULL, NULL, NULL); 
}


static void _delta_diff_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data) 
{                                                                                           
    charm_data_t       *data    = charm_get_quad_data(info->quad);                             
    charm_ctx_t        *ctx     = (charm_ctx_t*)info->p4est->user_pointer;                      
    charm_size_t        c_count = ctx->comp->elem_count;                                    
    charm_real_t       *a       = (charm_real_t*)user_data;                                       
    charm_fields_axpy(&(data->par.model.ns_jfnk.c_delta), &(data->par.model.ns_jfnk.c_delta_old), -1., c_count);                                       
}
                                                                                     
static inline void delta_diff(p4est_t * p4est)                                  
{                                                                                           
    p4est_iterate(p4est, NULL, NULL, _delta_diff_quad_iter_fn, NULL, NULL, NULL); 
}



static inline void _delta_err2_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{                                                                                           
    charm_real_t   *res = (charm_real_t*) user_data;                                       
    charm_data_t   *data = charm_get_quad_data(info->quad);                                 
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;                          
    charm_size_t    c_count = ctx->comp->elem_count;                                        
    *res += charm_fields_dot(&(data->par.model.ns_jfnk.c_delta_old), &(data->par.model.ns_jfnk.c_delta_old), c_count);              
}             

static inline charm_real_t delta_err2(p4est_t * p4est)                                                
{                                                                                           
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;                          
    charm_real_t        loc_res, glob_res;                                                  
    int                 mpiret;                                                             
                                                                                            
    loc_res = 0.0;                                                                          
    p4est_iterate (p4est, NULL,                                                             
                (void *) &loc_res,                                                          
                _delta_err2_quad_iter_fn,                                              
                NULL, NULL, NULL);                                                          
                                                                                            
    mpiret = sc_MPI_Allreduce (&loc_res, &glob_res, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);                                                                  
                                                                                            
    return glob_res;                                                                        
}


static void calc_delta(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_real_t        tau_k, err;
    charm_int_t         km_step, km_stop;
    charm_ctx_t        *ctx = charm_get_ctx(p4est);

    km_stop = 0;
    km_step = 0;
    zero_delta(p4est);
    while (!km_stop) {
        copy_delta_to_old(p4est);
        calc_ju_r_jr(p4est, ghost, ghost_data);

        tau_k = calc_tau(p4est);
        update_delta(p4est, tau_k);
        delta_diff(p4est);
        err = sqrt(delta_err2(p4est));
        if (err < ctx->model.ns_jfnk.newton.solver_rtol*ctx->tmp.ns_jfnk.rhs_norm 
                                    || km_step >= ctx->model.ns_jfnk.newton.solver_max_step) {
            km_stop = 1;
        }
        km_step++;
    }
}

static void _update_fld_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data) 
{                                                                                           
    charm_data_t       *data = charm_get_quad_data(info->quad);                             
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;                      
    charm_size_t              c_count = ctx->comp->elem_count;                                    
    charm_real_t       *a = (charm_real_t*)user_data;                                       
    charm_fields_axpy(&(data->par.model.ns_jfnk.c_delta), &(data->par.c), *a, c_count);                                       
}                                                                                           
static inline void update_fld(p4est_t * p4est, charm_real_t a)                                  
{                                                                                           
    p4est_iterate(p4est, NULL, (void*)&a, _update_fld_quad_iter_fn, NULL, NULL, NULL); 
}



static inline void _calc_err2_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{                                                                                           
    charm_real_t   *res = (charm_real_t*) user_data;                                       
    charm_data_t   *data = charm_get_quad_data(info->quad);                                 
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;                          
    charm_size_t          c_count = ctx->comp->elem_count;                                        
    *res += charm_fields_dot(&(data->par.model.ns_jfnk.c_delta), &(data->par.model.ns_jfnk.c_delta), c_count);              
}                                                                                           
static inline charm_real_t calc_err2(p4est_t * p4est)                                                
{                                                                                           
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;                          
    charm_real_t        loc_res, glob_res;                                                  
    int                 mpiret;                                                             
                                                                                            
    loc_res = 0.0;                                                                          
    p4est_iterate (p4est, NULL,                                                             
                (void *) &loc_res,                                                          
                _calc_err2_quad_iter_fn,                                              
                NULL, NULL, NULL);                                                          
                                                                                            
    mpiret = sc_MPI_Allreduce (&loc_res, &glob_res, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);                                                                  
                                                                                            
    return glob_res;                                                                        
}




charm_int_t charm_model_ns_jfnk_newton_step(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{    
    charm_ctx_t        *ctx = (charm_ctx_t*)p4est->user_pointer;
    charm_int_t         nm_stop = 0;
    charm_real_t        nm_err;

    calc_rhs(p4est, ghost, ghost_data);
    calc_delta(p4est, ghost, ghost_data);
    update_fld(p4est, ctx->model.ns_jfnk.newton.relax);
    nm_err = sqrt(calc_err2(p4est))*ctx->model.ns_jfnk.newton.relax;
    if (nm_err < ctx->model.ns_jfnk.newton.rtol * ctx->tmp.ns_jfnk.fld_old_norm) {
        nm_stop = 1;
    }
    return nm_stop;
}
