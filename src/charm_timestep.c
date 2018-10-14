//
// Created by zhrv on 26.10.17.
//

#include <p8est_iterate.h>
#include "charm_timestep.h"
#include "charm_bnd_cond.h"
#include "charm_base_func.h"
#include "charm_timestep_conv.h"
#include "charm_limiter.h"
#include "charm_timestep_press.h"
#include "charm_timestep_correct_velosity.h"

static double _charm_get_timestep (p4est_t * p4est);
static void _charm_timestep_single(p4est_t * p4est, int step, double time, double *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data);


/** Timestep the advection problem.
 *
 * Update the state, refine, repartition, and write the solution to file.
 *
 * \param [in,out] p4est the forest, whose state is updated
 * \param [in] time      the end time
 */
void charm_timesteps(p4est_t * p4est, double time)
{
    double              t;
    double              dt;
    int                 i;
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *ghost_data;
    p4est_ghost_t      *ghost;
    double              calc_time;


    ghost = p4est_ghost_new (p4est, CHARM_CONNECT_FULL);
    ghost_data = CHARM_ALLOC (charm_data_t, ghost->ghosts.elem_count);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);


    dt = _charm_get_timestep(p4est);

    CHARM_GLOBAL_ESSENTIAL("Starting time steps...\n");
    for (t = 0., i = 0; t < time; t += dt, i++) {
        calc_time = sc_MPI_Wtime();
        _charm_timestep_single(p4est, i, t, &dt, &ghost, &ghost_data);
        calc_time = sc_MPI_Wtime() - calc_time;
        if (i % ctx->log_period == 0) {
            charm_log_statistics(p4est, i, t, dt, calc_time);
        }
    }

    CHARM_FREE (ghost_data);
    p4est_ghost_destroy (ghost);
}


static void _charm_timestep_single(p4est_t * p4est, int step, double time, double *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data)
{
    double              t = 0.;
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    int                 refine_period = ctx->refine_period;
    int                 repartition_period = ctx->repartition_period;
    int                 write_period = ctx->write_period;
    int                 allowcoarsening = 1;
    int                 mpiret;
    double              orig_max_err = ctx->max_err;
    double              umax, global_umax;
    int                 ref_flag = 0;
    p4est_ghost_t      *ghost       = *_ghost;
    charm_data_t       *ghost_data  = *_ghost_data;

    /* refine */
    ref_flag = 0;
    if (refine_period) {
        if (!(step % refine_period)) {
            if (step) {
                charm_adapt(p4est, ghost, ghost_data); /* adapt */
                if (ghost) {
                    p4est_ghost_destroy(ghost);
                    CHARM_FREE (ghost_data);
                    ghost = NULL;
                    ghost_data = NULL;
                }

                ref_flag = 1;
            }
            *dt = _charm_get_timestep(p4est);

        }
    }
    else {
        *dt = _charm_get_timestep(p4est);
    }

    /* repartition */
    if (repartition_period) {
        if (step && !(step % repartition_period)) {

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
    if (!(step % write_period)) {
        charm_write_solution (p4est, step);
    }

    /* synchronize the ghost data */
    if (!ghost) {
        ghost = p4est_ghost_new (p4est, CHARM_CONNECT_FULL);
        ghost_data = CHARM_ALLOC (charm_data_t, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }

    charm_timestep_conv(p4est, ghost, ghost_data, dt);

    charm_timestep_press(p4est, ghost, ghost_data, dt);
    charm_timestep_correct_velosity(p4est, ghost, ghost_data, dt);
    charm_limiter(p4est, ghost, ghost_data);

    *_ghost       = ghost;
    *_ghost_data  = ghost_data;

}


static void _charm_timestep_min_dt_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    double         *dt = (double*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*) info->p4est->user_pointer;
    double          dt_loc;
    charm_cons_t    cons;
    charm_prim_t    prim;

    charm_get_fields(data, data->par.g.c, &cons);
    charm_param_cons_to_prim(info->p4est, &prim, &cons);

    dt_loc = ctx->CFL * data->par.g.volume / (sqrt(_MAG_(prim.u, prim.v, prim.w)));

    *dt = SC_MIN(*dt, dt_loc);
}

/** Compute the timestep.
 *
 * Find the smallest quadrant and scale the timestep based on that length and
 * the advection velocity.
 *
 * \param [in] p4est the forest
 * \return the timestep.
 */
static double _charm_get_timestep (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    double              loc_dt, glob_dt;
    int                 mpiret, i;

    return ctx->dt;
    loc_dt = ctx->dt;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_dt,
                   _charm_timestep_min_dt_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_dt, &glob_dt, 1, sc_MPI_DOUBLE, sc_MPI_MIN, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return glob_dt;
}



