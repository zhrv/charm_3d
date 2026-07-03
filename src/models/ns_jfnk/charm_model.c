//
// Created by zhrv on 10.01.26.
//

#include <p8est_iterate.h>
#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"


void charm_model_ns_jfnk_timestep_conv(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_jfnk_timestep_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_jfnk_geom_calc(p4est_t *p4est);

charm_int_t charm_model_ns_jfnk_newton_step(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void charm_model_ns_jfnk_timestep_min_dt_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
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
charm_real_t charm_model_ns_jfnk_get_dt (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t              loc_dt, glob_dt;
    int                 mpiret, i;

    return ctx->dt;
    loc_dt = ctx->dt;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_dt,
                   charm_model_ns_jfnk_timestep_min_dt_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_dt, &glob_dt, 1, sc_MPI_DOUBLE, sc_MPI_MIN, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return glob_dt;
}


static void charm_model_ns_jfnk_timestep_copy_to_old_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t              c_count = ctx->comp->elem_count;

    charm_fields_copy(&(data->par.c_old), &(data->par.c), c_count);
}


static void _charm_model_ns_jfnk_calc_old_norm2_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_real_t   *err2 = (charm_real_t*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t          c_count = ctx->comp->elem_count;
    *err2 += charm_fields_get_norm2(&(data->par.c_old), c_count);
}

static void charm_model_ns_jfnk_calc_old_norm2 (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t        loc_err2, glob_err2;
    int                 mpiret, i;

    loc_err2 = 0.0;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_err2,
                   _charm_model_ns_jfnk_calc_old_norm2_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_err2, &glob_err2, 1, sc_MPI_DOUBLE, sc_MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    ctx->tmp.ns_jfnk.fld_old_norm = sqrt(glob_err2);
}


void charm_model_ns_jfnk_timestep_single(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    int                 refine_period = ctx->refine_period;
    int                 repartition_period = ctx->repartition_period;
    int                 write_period = ctx->write_period;
    int                 allowcoarsening = 1;
    p4est_ghost_t      *ghost       = *_ghost;
    charm_data_t       *ghost_data  = *_ghost_data;
    charm_int_t         nm_stop;
    charm_int_t         nm_step;
    charm_real_t        fld_old_norm;

    if (!ctx->timestep) {
        charm_model_ns_jfnk_geom_calc(p4est);
    }
    if (refine_period) {
        if (!(ctx->timestep % refine_period)) {
            if (ctx->timestep) {
                ctx->amr_fn(p4est, ghost, ghost_data); /* adapt */
                charm_model_ns_jfnk_geom_calc(p4est); //@todo выяснить нужно или нет
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

    ctx->tmp.ns_jfnk.dt = *dt;

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
        //p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }

    p4est_iterate (p4est, NULL, NULL, 
        charm_model_ns_jfnk_timestep_copy_to_old_quad_iter_fn, 
        NULL, NULL, NULL);
    charm_model_ns_jfnk_calc_old_norm2 (p4est);

    nm_stop = 0;
    nm_step = 0;
    while (!nm_stop  && nm_step < ctx->model.ns_jfnk.newton.max_step) { // итерации метода Ньютона
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
        nm_stop = charm_model_ns_jfnk_newton_step(p4est, ghost, ghost_data);
        nm_step++;
    }


    *_ghost       = ghost;
    *_ghost_data  = ghost_data;
}

