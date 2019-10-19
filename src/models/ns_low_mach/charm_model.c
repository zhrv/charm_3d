//
// Created by zhrv on 27.08.19.
//

#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"
#include "charm_xml.h"


void charm_model_ns_low_mach_timestep_conv(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_low_mach_timestep_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_low_mach_pressure(p4est_t * p4est, double *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void _charm_model_ns_low_mach_timestep_min_dt_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_t        *p4est = info->p4est;
    double         *dt = (double*) user_data;
    charm_data_t   *data = charm_get_quad_data(info->quad);
    charm_ctx_t    *ctx = (charm_ctx_t*) info->p4est->user_pointer;
    double          dt_loc;
    charm_cons_t    cons;
    charm_prim_t    prim;

    charm_get_fields(p4est, data, data->par.g.c, &cons);
    charm_param_cons_to_prim(info->p4est, &prim, &cons);

    dt_loc = ctx->CFL * data->par.g.volume / (sqrt(_MAG_(prim.u, prim.v, prim.w)) + prim.cz);

    *dt = SC_MIN(*dt, dt_loc);
}

/** Compute the timestep.
 *
 * \param [in] p4est the forest
 * \return the timestep.
 */
double charm_model_ns_low_mach_get_dt (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    double              loc_dt, glob_dt;
    int                 mpiret, i;

    return ctx->dt;
    loc_dt = ctx->dt;
    p4est_iterate (p4est, NULL,
                   (void *) &loc_dt,
                   _charm_model_ns_low_mach_timestep_min_dt_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_dt, &glob_dt, 1, sc_MPI_DOUBLE, sc_MPI_MIN, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return glob_dt;
}







static void charm_model_ns_low_mach_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    double              dt = *((double *) user_data);
    double              rhs_ru[CHARM_BASE_FN_COUNT];
    double              rhs_rv[CHARM_BASE_FN_COUNT];
    double              rhs_rw[CHARM_BASE_FN_COUNT];
    double              rhs_rh[CHARM_BASE_FN_COUNT];
    double              rhs_rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];
    size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_matr_vect_mult(data->par.g.a_inv, data->int_ru, rhs_ru);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rv, rhs_rv);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rw, rhs_rw);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rh, rhs_rh);

    for (j = 0; j < c_count; j++) {
        charm_matr_vect_mult(data->par.g.a_inv, data->int_rc[j], rhs_rc[j]);
    }

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ru[i] -= _NORM_(dt * rhs_ru[i]);
        data->par.c.rv[i] -= _NORM_(dt * rhs_rv[i]);
        data->par.c.rw[i] -= _NORM_(dt * rhs_rw[i]);
        data->par.c.rh[i] -= _NORM_(dt * rhs_rh[i]);
        for (j = 0; j < c_count; j++) {
            data->par.c.rc[j][i] -= _NORM_(dt * rhs_rc[j][i]);
        }
    }
}


static void charm_model_ns_low_mach_timestep_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i, j;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->int_ru[i] = 0.;
        data->int_rv[i] = 0.;
        data->int_rw[i] = 0.;
        data->int_rh[i] = 0.;
        for (j = 0; j < CHARM_MAX_COMPONETS_COUNT; j++) {
            data->int_rc[j][i] = 0.;
        }
    }
}



void charm_model_ns_low_mach_timestep_single(p4est_t * p4est, double *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data)
{
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
        if (!(ctx->timestep % refine_period)) {
            if (ctx->timestep) {
                charm_adapt(p4est, ghost, ghost_data); /* adapt */
                if (ghost) {
                    p4est_ghost_destroy(ghost);
                    CHARM_FREE (ghost_data);
                    ghost = NULL;
                    ghost_data = NULL;
                }

                ref_flag = 1;
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

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   charm_model_ns_low_mach_timestep_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    charm_model_ns_low_mach_timestep_conv(p4est, ghost, ghost_data);
    //charm_model_ns_low_mach_timestep_diff(p4est, ghost, ghost_data);

    p4est_iterate (p4est, NULL,
                   (void *) dt,
                   charm_model_ns_low_mach_timestep_update_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_ghost_exchange_data (p4est, ghost, ghost_data); /* synchronize the ghost data */

    charm_model_ns_low_mach_pressure(p4est, dt, ghost, ghost_data);

    charm_limiter(p4est, ghost, ghost_data);

    *_ghost       = ghost;
    *_ghost_data  = ghost_data;

}


void charm_model_ns_low_mach_init(charm_ctx_t *ctx, mxml_node_t *node)
{
    ctx->get_dt_fn              = charm_model_ns_low_mach_get_dt;
    ctx->timestep_single_fn     = charm_model_ns_low_mach_timestep_single;

    charm_xml_node_child_param_int(node, "PRESSURE_ITERATIONS", &(ctx->pressure_iterations));
    charm_xml_node_child_param_dbl(node, "PRESSURE_TAU",        &(ctx->pressure_tau));
}

