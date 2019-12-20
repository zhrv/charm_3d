//
// Created by zhrv on 27.08.19.
//

#include <p8est_iterate.h>
#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"


void charm_model_ns_li_timestep_conv(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_li_timestep_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_li_timestep_chem(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_li_timestep_diffusion(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void _charm_model_ns_li_timestep_min_dt_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    double         *dt = (double*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*) info->p4est->user_pointer;
    double          dt_loc;
    charm_cons_t    cons;
    charm_prim_t    prim;

    charm_get_fields(data, data->par.g.c, &cons);
    charm_param_cons_to_prim(info->p4est, &prim, &cons);

    dt_loc = ctx->CFL * data->par.g.volume / (sqrt(_MAG_(prim.u, prim.v, prim.w)) + prim.cz);

    *dt = SC_MIN(*dt, dt_loc);
}

/** Compute the timestep.
 *
 * \param [in] p4est the forest
 * \return the timestep.
 */
double charm_model_ns_li_get_dt (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    double              loc_dt, glob_dt;
    int                 mpiret, i;

    return ctx->dt;
    loc_dt = ctx->dt;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_dt,
                   _charm_model_ns_li_timestep_min_dt_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_dt, &glob_dt, 1, sc_MPI_DOUBLE, sc_MPI_MIN, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return glob_dt;
}







static void charm_model_ns_li_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    double             *dt  = (double *) user_data;

    double              rhs_ru[CHARM_BASE_FN_COUNT];
    double              rhs_rv[CHARM_BASE_FN_COUNT];
    double              rhs_rw[CHARM_BASE_FN_COUNT];
    double              rhs_re[CHARM_BASE_FN_COUNT];
    double              rhs_rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];

    double              lhs_ru[CHARM_BASE_FN_COUNT];
    double              lhs_rv[CHARM_BASE_FN_COUNT];
    double              lhs_rw[CHARM_BASE_FN_COUNT];
    double              lhs_re[CHARM_BASE_FN_COUNT];
    double              lhs_rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];

    size_t              c_count = ctx->comp->elem_count;
    int                 i, j;
    double              _tau;

    charm_vect_copy(rhs_ru, data->par.c.ru);
    charm_vect_copy(rhs_rv, data->par.c.rv);
    charm_vect_copy(rhs_rw, data->par.c.rw);
    charm_vect_copy(rhs_re, data->par.c.re);
    for (j = 0; j < c_count; j++) {
        charm_vect_copy(rhs_rc[j], data->par.c.rc[j]);
    }

    charm_vect_sub(rhs_ru, data->par.c_old.ru);
    charm_vect_sub(rhs_rv, data->par.c_old.rv);
    charm_vect_sub(rhs_rw, data->par.c_old.rw);
    charm_vect_sub(rhs_re, data->par.c_old.re);
    for (j = 0; j < c_count; j++) {
        charm_vect_sub(rhs_rc[j], data->par.c_old.rc[j]);
    }

    _tau = 1./dt[0];
    charm_vect_scalar_mult(rhs_ru, _tau);
    charm_vect_scalar_mult(rhs_rv, _tau);
    charm_vect_scalar_mult(rhs_rw, _tau);
    charm_vect_scalar_mult(rhs_re, _tau);
    for (j = 0; j < c_count; j++) {
        charm_vect_scalar_mult(rhs_rc[j], _tau);
    }

    charm_matr_vect_mult(data->par.g.a, rhs_ru, lhs_ru);
    charm_matr_vect_mult(data->par.g.a, rhs_rv, lhs_rv);
    charm_matr_vect_mult(data->par.g.a, rhs_rw, lhs_rw);
    charm_matr_vect_mult(data->par.g.a, rhs_re, lhs_re);
    for (j = 0; j < c_count; j++) {
        charm_matr_vect_mult(data->par.g.a, rhs_rc[j], lhs_rc[j]);
    }

    charm_vect_add(lhs_ru, data->int_ru);
    charm_vect_add(lhs_rv, data->int_rv);
    charm_vect_add(lhs_rw, data->int_rw);
    charm_vect_add(lhs_re, data->int_re);
    for (j = 0; j < c_count; j++) {
        charm_vect_add(lhs_rc[j], data->int_rc[j]);
    }




    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ru[i] -= _NORM_(dt[1] * lhs_ru[i]);
        data->par.c.rv[i] -= _NORM_(dt[1] * lhs_rv[i]);
        data->par.c.rw[i] -= _NORM_(dt[1] * lhs_rw[i]);
        data->par.c.re[i] -= _NORM_(dt[1] * lhs_re[i]);
        for (j = 0; j < c_count; j++) {
            data->par.c.rc[j][i] -= _NORM_(dt[1] * lhs_rc[j][i]);
        }
    }
}


static void charm_model_ns_li_timestep_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
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


static void _charm_model_ns_li_timestep_0(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i, j;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c_old.ru[i] = data->par.c.ru[i];
        data->par.c_old.rv[i] = data->par.c.rv[i];
        data->par.c_old.rw[i] = data->par.c.rw[i];
        data->par.c_old.re[i] = data->par.c.re[i];
        for (j = 0; j < CHARM_MAX_COMPONETS_COUNT; j++) {
            data->par.c_old.rc[j][i] = data->par.c.rc[j][i];
        }
    }
}


void charm_model_ns_li_timestep_single(p4est_t * p4est, double *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    int                 refine_period = ctx->refine_period;
    int                 repartition_period = ctx->repartition_period;
    int                 write_period = ctx->write_period;
    int                 allowcoarsening = 1;
    int                 i_sig;
    p4est_ghost_t      *ghost       = *_ghost;
    charm_data_t       *ghost_data  = *_ghost_data;
    double              dt_sig[2], *sig;

    if (refine_period) {
        if (!(ctx->timestep % refine_period)) {
            if (ctx->timestep) {
                charm_adapt(p4est, ghost, ghost_data); /* adapt */
                if (ghost) {
                    p4est_ghost_destroy(ghost);
                    CHARM_FREE (ghost_data);
                    ghost = NULL;
                    ghost_data = NULL;
                }
            }
            *dt = ctx->get_dt_fn(p4est);

        }
    }
    else {
        *dt = ctx->get_dt_fn(p4est);
    }

    /* repartition */
    if (repartition_period) {
        if (ctx->timestep && !(ctx->timestep % repartition_period)) {

            p4est_partition(p4est, allowcoarsening, NULL);

            if (ghost) {
                p4est_ghost_destroy(ghost);
                CHARM_FREE (ghost_data);
                ghost = NULL;
                ghost_data = NULL;
            }
        }
    }

    /* write out solution */
    if (!(ctx->timestep % write_period)) {
        charm_write_solution (p4est);
        CHARM_GLOBAL_ESSENTIALF (" File for step #%d is saved \n", ctx->timestep);
    }

    /* synchronize the ghost data */
    if (!ghost) {
        ghost = p4est_ghost_new (p4est, CHARM_CONNECT_FULL);
        ghost_data = CHARM_ALLOC (charm_data_t, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }


    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_li_timestep_0, NULL, NULL, NULL);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

    dt_sig[0] = *dt;
    for (i_sig = 0; i_sig < ctx->model.ns_li.li_sigma->elem_count; i_sig++) {
        sig = (double*) sc_array_index(ctx->model.ns_li.li_sigma, i_sig);
        dt_sig[1] = *sig;

        charm_model_ns_li_timestep_chem(p4est, ghost, ghost_data);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
        p4est_iterate (p4est, ghost, (void *) ghost_data,
                charm_model_ns_li_timestep_zero_quad_iter_fn,
                NULL, NULL, NULL);
        charm_model_ns_li_timestep_diff(p4est, ghost, ghost_data);
        charm_model_ns_li_timestep_conv(p4est, ghost, ghost_data);
        p4est_iterate (p4est, NULL, (void *) dt_sig,
                charm_model_ns_li_timestep_update_quad_iter_fn,
                NULL, NULL, NULL);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
        charm_limiter(p4est, ghost, ghost_data);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }


    *_ghost       = ghost;
    *_ghost_data  = ghost_data;

}


