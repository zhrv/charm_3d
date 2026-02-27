//
// Created by zhrv on 10.01.26.
//

#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"

#include "charm_model_ns_jfnk_newton_iter_fn.h"


void charm_model_ns_jfnk_dg_operator(p4est_t * p4est, p4est_ghost_t * _ghost, charm_data_t * _ghost_data);
void charm_model_ns_jfnk_dg_operator_stash_push(p4est_t * p4est);
void charm_model_ns_jfnk_dg_operator_stash_pop(p4est_t * p4est);



static void _calc_rhs(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_real_t        loc_norm2, glob_norm2;
    charm_ctx_t        *ctx = charm_get_ctx(p4est);
    int                 mpiret;
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    
    p4est_iterate (p4est, NULL, NULL, 
        _result_to_rhs_quad_iter_fn, 
        NULL, NULL, NULL);

    loc_norm2 = 0.0;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_norm2,
                   _rhs_norm2_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_norm2, &glob_norm2, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    ctx->tmp.ns_jfnk.rhs_norm = sqrt(glob_norm2);
}






static charm_real_t _calc_delta_err(p4est_t * p4est)
{
    charm_real_t        loc_err2, glob_err2;
    charm_ctx_t        *ctx = charm_get_ctx(p4est);
    int                 mpiret;
    
    loc_err2 = 0.0;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_err2,
                   _delta_err2_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_err2, &glob_err2, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return sqrt(glob_err2);
}

/*
__host__ 
inline void calc_Ju(data_t *d, fields_t *fld_u, fields_t *fld_out) {
    real_t epsilon = EPS_J; // TODO epsilon
    d->copy_fld_to(d->d.fld_tmp1);
    calc_dg_operator(d);
    d->copy_dg_res_to(d->d.fld_tmp2);
    krnl_axpy<<<BLOCKS, THREADS>>>(fld_u, d->d.fld, epsilon); 
    calc_dg_operator(d);
    d->copy_dg_res_to(fld_out);
    krnl_axpy<<<BLOCKS, THREADS>>>(d->d.fld_tmp2, fld_out, -1.0); 
    krnl_mult_scalar<<<BLOCKS, THREADS>>>(fld_out, 1./epsilon); 
    d->copy_fld_from(d->d.fld_tmp1);
}

*/
static inline void _calc_ju(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_model_ns_jfnk_dg_operator_stash_push(p4est);
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    p4est_iterate (p4est, NULL, NULL,
                   _calc_ju_copy_result_add_arg_quad_iter_fn,
                   NULL, NULL, NULL);
    charm_model_ns_jfnk_dg_operator(p4est, ghost, ghost_data);
    p4est_iterate (p4est, NULL, NULL,
                   _calc_ju_result_quad_iter_fn,
                   NULL, NULL, NULL);
    charm_model_ns_jfnk_dg_operator_stash_pop(p4est);
        
}


/*static inline void _calc_r(p4est_t * p4est)
{
    p4est_iterate (p4est, NULL, NULL,
                   _calc_r_quad_iter_fn,
                   NULL, NULL, NULL);
}*/

CHARM_DECL_QUAD_ITER_STATIC(_calc_r, {
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    size_t              c_count = ctx->comp->elem_count;

    charm_fields_copy(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_ju, c_count);
    charm_fields_sub(data->par.model.ns_jfnk.c_residual, data->par.model.ns_jfnk.c_rhs, c_count);
})


static inline void _calc_jr(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{

}


static inline void _calc_ju_r_jr(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    _calc_ju(p4est, ghost, ghost_data);
    _calc_r(p4est);
    _calc_jr(p4est, ghost, ghost_data);
}



static charm_real_t _calc_tau(p4est_t * p4est) 
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t        loc_res;
    int                 mpiret;
    charm_real_t        jr_r, jr_jr;
    
    // (Jr,r)
    loc_res = 0.0;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_res,
                   _calc_jr_r_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_res, &jr_r, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    // (Jr,Jr)
    loc_res = 0.0;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_res,
                   _calc_jr_jr_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_res, &jr_jr, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return -jr_r/jr_jr;

}


static void _calc_delta(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_real_t        tau_k, err;
    charm_int_t         km_step, km_stop;
    charm_ctx_t        *ctx = charm_get_ctx(p4est);

    km_stop = 0;
    km_step = 0;
    while (!km_stop) {
        p4est_iterate (p4est, NULL, NULL, 
            _copy_delta_to_old_quad_iter_fn, 
            NULL, NULL, NULL);
        _calc_ju_r_jr(p4est, ghost, ghost_data);

        tau_k = _calc_tau(p4est);
        p4est_iterate (p4est, NULL, &tau_k, 
            _update_delta_quad_iter_fn, 
            NULL, NULL, NULL);
        err = _calc_delta_err(p4est);
        if (err < ctx->model.ns_jfnk.newton.solver_rtol*ctx->tmp.ns_jfnk.rhs_norm 
                                    || km_step >= ctx->model.ns_jfnk.newton.solver_max_step) {
            km_stop = 1;
        }
        km_step++;
    }
}


static void _update_fld(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    p4est_iterate (p4est, NULL, NULL, 
        _update_fld_quad_iter_fn, 
        NULL, NULL, NULL);
        
}


static charm_real_t _calc_err2 (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t        loc_err2, glob_err2;
    int                 mpiret;

    loc_err2 = 0.0;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_err2,
                   _calc_err2_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_err2, &glob_err2, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return sqrt(glob_err2);
}


charm_int_t charm_model_ns_jfnk_newton_step(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{    
    charm_ctx_t        *ctx = (charm_ctx_t*)p4est->user_pointer;
    charm_int_t         nm_stop = 0;
    charm_real_t        nm_err;

    _calc_rhs(p4est, ghost, ghost_data);
    _calc_delta(p4est, ghost, ghost_data);
    _update_fld(p4est, ghost, ghost_data);
    nm_err = _calc_err2(p4est);
    if (nm_err < ctx->model.ns_jfnk.newton.rtol*ctx->tmp.ns_jfnk.fld_old_norm) {
        nm_stop = 1;
    }
    return nm_stop;
}
