//
// Created by zhrv on 26.10.17.
//

#include <p8est_iterate.h>
#include "charm_timestep.h"
#include "charm_bnd_cond.h"
#include "charm_base_func.h"
#include "charm_timestep_conv.h"

static double _charm_get_timestep (p4est_t * p4est);
static void charm_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data);
static void _charm_timestep_single(p4est_t * p4est, int step, double time, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data);


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


    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);


    dt = _charm_get_timestep(p4est);

    CHARM_GLOBAL_ESSENTIAL("Starting time steps...\n");
    for (t = 0., i = 0; t < time; t += dt, i++) {
        _charm_timestep_single(p4est, i, t, &ghost, &ghost_data);
        if (i % ctx->log_period == 0) {
            charm_log_statistics(p4est, i, t);
        }
    }

    P4EST_FREE (ghost_data);
    p4est_ghost_destroy (ghost);
}


static void _charm_timestep_single(p4est_t * p4est, int step, double time, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data)
{
    double              t = 0.;
    double              dt;
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
                charm_adapt(p4est); /* adapt */
                if (ghost) {
                    p4est_ghost_destroy(ghost);
                    P4EST_FREE (ghost_data);
                    ghost = NULL;
                    ghost_data = NULL;
                }

                ref_flag = 1;
            }
            dt = _charm_get_timestep(p4est);

        }
    }

    /* repartition */
    if (step && !(step % repartition_period)) {

        p4est_partition (p4est, allowcoarsening, NULL);

        if (ghost) {
            p4est_ghost_destroy (ghost);
            P4EST_FREE (ghost_data);
            ghost = NULL;
            ghost_data = NULL;
        }
    }

    /* write out solution */
    if (!(step % write_period)) {
        charm_write_solution (p4est, step);
        CHARM_GLOBAL_ESSENTIALF ("**************** File for step #%d is saved ***************\n", step);
    }

    /* synchronize the ghost data */
    if (!ghost) {
        ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
        ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   charm_convect_volume_int_iter_fn,
                   NULL,
                   NULL,
                   NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   NULL,
                   charm_convect_surface_int_iter_fn,
                   NULL,
                   NULL);


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
                   (void *) &dt,                             /* pass in dt */
                   charm_timestep_update_quad_iter_fn,       /* update each cell */
                   NULL,                                     /* there is no callback for the faces between quadrants */
                   NULL,                                     /* there is no callback for the faces between quadrants */
                   NULL);                                    /* there is no callback for the corners between quadrants */

    /* synchronize the ghost data */
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

    /* update grad */
//    charm_calc_grad(p4est, ghost, ghost_data);

    *_ghost       = ghost;
    *_ghost_data  = ghost_data;

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
//    p4est_topidx_t      t, flt, llt;
//    p4est_tree_t       *tree;
//    int                 max_level, global_max_level;
//    int                 mpiret, i;
//    double              min_h, vnorm;
//    double              dt;
//
//    /* compute the timestep by finding the smallest quadrant */
//    flt = p4est->first_local_tree;
//    llt = p4est->last_local_tree;
//
//    max_level = 0;
//    for (t = flt; t <= llt; t++) {
//        tree = p4est_tree_array_index (p4est->trees, t);
//        max_level = SC_MAX (max_level, tree->maxlevel);
//
//    }
//    mpiret =
//            sc_MPI_Allreduce (&max_level, &global_max_level, 1, sc_MPI_INT,
//                              sc_MPI_MAX, p4est->mpicomm);
//    SC_CHECK_MPI (mpiret);
//
//    min_h = CHARM_GET_H(global_max_level);
//
//    vnorm = 0;
//    for (i = 0; i < P4EST_DIM; i++) {
////        vnorm += ctx->v[i] * ctx->v[i];
//    }
//    vnorm = sqrt (vnorm);
//
//    min_h = CHARM_GET_H(ctx->allowed_level);
//
//    dt = ctx->CFL* min_h  / ctx->v_ref;

    return ctx->dt;
}




static void charm_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    double              dt = *((double *) user_data);
    double              rhs_ro[CHARM_BASE_FN_COUNT];
    double              rhs_ru[CHARM_BASE_FN_COUNT];
    double              rhs_rv[CHARM_BASE_FN_COUNT];
    double              rhs_rw[CHARM_BASE_FN_COUNT];
    double              rhs_re[CHARM_BASE_FN_COUNT];
    int                 i;

    charm_matr_vect_mult(data->par.g.a_inv, data->int_ro, rhs_ro);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_ru, rhs_ru);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rv, rhs_rv);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rw, rhs_rw);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_re, rhs_re);

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ro[i] += dt * rhs_ro[i];
        data->par.c.ru[i] += dt * rhs_ru[i];
        data->par.c.rv[i] += dt * rhs_rv[i];
        data->par.c.rw[i] += dt * rhs_rw[i];
        data->par.c.re[i] += dt * rhs_re[i];
    }
}

