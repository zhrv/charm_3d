//
// Created by zhrv on 26.10.17.
//

#include <p8est_iterate.h>
#include "charm_timestep.h"
#include "charm_bnd_cond.h"
#include "charm_base_func.h"
#include "charm_timestep_conv.h"
#include "charm_limiter.h"

static double _charm_get_timestep (p4est_t * p4est);
static void charm_timestep_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data);
static void charm_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data);
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
        CHARM_GLOBAL_ESSENTIALF (" File for step #%d is saved \n", step);
    }

    /* synchronize the ghost data */
    if (!ghost) {
        ghost = p4est_ghost_new (p4est, CHARM_CONNECT_FULL);
        ghost_data = CHARM_ALLOC (charm_data_t, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   charm_timestep_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   charm_convect_volume_int_iter_fn,
                   charm_convect_surface_int_iter_fn,
                   NULL, NULL);


//    p4est_iterate (p4est,                                /* the forest */
//                   ghost,                                /* the ghost layer */
//                   (void *) ghost_data,                  /* the synchronized ghost data */
//                   NULL,                                 /* callback to compute each quad's interior contribution to du/dt */
//                   charm_diffusion_flux_face_iter_fn,    /* callback to compute each quads' faces' contributions to du/du */
//                   NULL,
//                   NULL);                                /* there is no callback for the
//                                                                corners between quadrants */

    /* update u */
    p4est_iterate (p4est, NULL,                              /* ghosts are not needed for this loop */
                   (void *) dt,                             /* pass in dt */
                   charm_timestep_update_quad_iter_fn,       /* update each cell */
                   NULL,                                     /* there is no callback for the faces between quadrants */
                   NULL,                                     /* there is no callback for the faces between quadrants */
                   NULL);                                    /* there is no callback for the corners between quadrants */

    /* synchronize the ghost data */
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

    charm_limiter(p4est, ghost, ghost_data);
    /* update grad */
//    charm_calc_grad(p4est, ghost, ghost_data);

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

    dt_loc = ctx->CFL * data->par.g.volume / (sqrt(_MAG_(prim.u, prim.v, prim.w)) + prim.cz);

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


static void charm_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    double              dt = *((double *) user_data);
    double              rhs_ro[CHARM_BASE_FN_COUNT];
    double              rhs_ru[CHARM_BASE_FN_COUNT];
    double              rhs_rv[CHARM_BASE_FN_COUNT];
    double              rhs_rw[CHARM_BASE_FN_COUNT];
    double              rhs_re[CHARM_BASE_FN_COUNT];
    double              rhs_rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];
    size_t              c_count = ctx->comp->elem_count;
    int                 i;

    charm_matr_vect_mult(data->par.g.a_inv, data->int_ro, rhs_ro);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_ru, rhs_ru);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rv, rhs_rv);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rw, rhs_rw);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_re, rhs_re);

    for (int j = 0; j < c_count; j++) {
        charm_matr_vect_mult(data->par.g.a_inv, data->int_rc[j], rhs_rc[j]);
    }

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ro[i] -= _NORM_(dt * rhs_ro[i]);
        data->par.c.ru[i] -= _NORM_(dt * rhs_ru[i]);
        data->par.c.rv[i] -= _NORM_(dt * rhs_rv[i]);
        data->par.c.rw[i] -= _NORM_(dt * rhs_rw[i]);
        data->par.c.re[i] -= _NORM_(dt * rhs_re[i]);
        for (int j = 0; j < c_count; j++) {
            data->int_rc[j][i] -= _NORM_(dt * rhs_rc[j][i]);;
        }
    }
}


static void charm_timestep_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->int_ro[i] = 0.;
        data->int_ru[i] = 0.;
        data->int_rv[i] = 0.;
        data->int_rw[i] = 0.;
        data->int_re[i] = 0.;
        for (int j = 0; j < CHARM_MAX_COMPONETS_COUNT; j++) {
            data->int_rc[j][i] = 0.;
        }
    }
}

