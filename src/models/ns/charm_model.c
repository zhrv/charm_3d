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


void charm_model_ns_timestep_conv(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_timestep_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_timestep_chem(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_timestep_diffusion(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void _charm_model_ns_timestep_min_dt_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_real_t         *dt = (charm_real_t*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*) info->p4est->user_pointer;
    charm_real_t          dt_loc;
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
charm_real_t charm_model_ns_get_dt (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t              loc_dt, glob_dt;
    int                 mpiret, i;

    return ctx->dt;
    loc_dt = ctx->dt;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_dt,
                   _charm_model_ns_timestep_min_dt_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_dt, &glob_dt, 1, sc_MPI_DOUBLE, sc_MPI_MIN, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return glob_dt;
}







static void charm_model_ns_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_real_t              dt = *((charm_real_t *) user_data);
    charm_real_t              rhs_ru[CHARM_BASE_FN_COUNT];
    charm_real_t              rhs_rv[CHARM_BASE_FN_COUNT];
    charm_real_t              rhs_rw[CHARM_BASE_FN_COUNT];
    charm_real_t              rhs_re[CHARM_BASE_FN_COUNT];
    charm_real_t              rhs_rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];
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


static void charm_model_ns_timestep_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
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


static void _charm_model_ns_timestep_rk_0(p4est_iter_volume_info_t * info, void *user_data)
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


static void _charm_model_ns_timestep_rk_1(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i, j;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ru[i] *= 0.25;
        data->par.c.rv[i] *= 0.25;
        data->par.c.rw[i] *= 0.25;
        data->par.c.re[i] *= 0.25;

        data->par.c.ru[i] += 0.75*data->par.c_old.ru[i];
        data->par.c.rv[i] += 0.75*data->par.c_old.rv[i];
        data->par.c.rw[i] += 0.75*data->par.c_old.rw[i];
        data->par.c.re[i] += 0.75*data->par.c_old.re[i];

        for (j = 0; j < CHARM_MAX_COMPONETS_COUNT; j++) {
            data->par.c.rc[j][i] *= 0.25;
            data->par.c.rc[j][i] += 0.75*data->par.c_old.rc[j][i];
        }
    }
}


static void _charm_model_ns_timestep_rk_2(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i, j;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ru[i] *= 2.;
        data->par.c.rv[i] *= 2.;
        data->par.c.rw[i] *= 2.;
        data->par.c.re[i] *= 2.;
        data->par.c.ru[i] /= 3.;
        data->par.c.rv[i] /= 3.;
        data->par.c.rw[i] /= 3.;
        data->par.c.re[i] /= 3.;
        data->par.c.ru[i] += data->par.c_old.ru[i] / 3.;
        data->par.c.rv[i] += data->par.c_old.rv[i] / 3.;
        data->par.c.rw[i] += data->par.c_old.rw[i] / 3.;
        data->par.c.re[i] += data->par.c_old.re[i] / 3.;
        for (j = 0; j < CHARM_MAX_COMPONETS_COUNT; j++) {
            data->par.c.rc[j][i] *= 2.;
            data->par.c.rc[j][i] /= 3.;
            data->par.c.rc[j][i] += data->par.c_old.rc[j][i] / 3.;
        }
    }
}


void charm_model_ns_timestep_single(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    int                 refine_period = ctx->refine_period;
    int                 repartition_period = ctx->repartition_period;
    int                 write_period = ctx->write_period;
    int                 allowcoarsening = 1;
    p4est_ghost_t      *ghost       = *_ghost;
    charm_data_t       *ghost_data  = *_ghost_data;

    if (refine_period) {
        if (!(ctx->timestep % refine_period)) {
            if (ctx->timestep) {
                ctx->amr_fn(p4est, ghost, ghost_data); /* adapt */
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

#define CHARM_RUNGE_KUTTA_STEP()                                                    \
                    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           \
                    charm_model_ns_timestep_chem(p4est, ghost, ghost_data);         \
                    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           \
                    p4est_iterate (p4est, ghost, (void *) ghost_data,               \
                            charm_model_ns_timestep_zero_quad_iter_fn,              \
                            NULL, NULL, NULL);                                      \
                    charm_model_ns_timestep_diff(p4est, ghost, ghost_data);         \
                    charm_model_ns_timestep_conv(p4est, ghost, ghost_data);         \
                    p4est_iterate (p4est, NULL, (void *) dt,                        \
                            charm_model_ns_timestep_update_quad_iter_fn,            \
                            NULL, NULL, NULL);                                      \
                    p4est_ghost_exchange_data (p4est, ghost, ghost_data);           \
                    charm_limiter(p4est, ghost, ghost_data);

    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_timestep_rk_0, NULL, NULL, NULL);
    CHARM_RUNGE_KUTTA_STEP()
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    CHARM_RUNGE_KUTTA_STEP()
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_timestep_rk_1, NULL, NULL, NULL);
    CHARM_RUNGE_KUTTA_STEP()
    p4est_iterate (p4est, NULL, NULL, _charm_model_ns_timestep_rk_2, NULL, NULL, NULL);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

#undef CHARM_RUNGE_KUTTA_STEP

    *_ghost       = ghost;
    *_ghost_data  = ghost_data;

}


