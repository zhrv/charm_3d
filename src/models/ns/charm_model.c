//
// Created by zhrv on 27.08.19.
//

#include <p8est_iterate.h>
#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_limiter.h"


void charm_model_ns_timestep_conv(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_timestep_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_timestep_chem(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_timestep_diffusion(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_geom_calc(p4est_t *p4est);


static void charm_model_ns_timestep_min_dt_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
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
                   charm_model_ns_timestep_min_dt_quad_iter_fn,
                   NULL, NULL, NULL);

    mpiret = sc_MPI_Allreduce (&loc_dt, &glob_dt, 1, sc_MPI_DOUBLE, sc_MPI_MIN, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    return glob_dt;
}




static void charm_model_ns_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t              c_count = ctx->comp->elem_count;
    charm_real_t        dt = *((charm_real_t *) user_data);
    charm_fields_t      rhs;
    int                 i, j;

    charm_matr_fields_mult(data->par.g.a_inv, data->integrals, rhs, c_count);

    charm_fields_axpy(rhs, data->par.c, dt, c_count);
}


static void charm_model_ns_timestep_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_fields_zero(data->integrals, c_count);
}


static void charm_model_ns_timestep_rk_0(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_fields_copy(data->par.c_old, data->par.c, c_count);
}


static void charm_model_ns_timestep_rk_1(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_fields_mult(data->par.c, 0.25, c_count);
    charm_fields_axpy(data->par.c_old, data->par.c, 0.75, c_count);
}


static void charm_model_ns_timestep_rk_2(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_size_t              c_count = ctx->comp->elem_count;
    int                 i, j;

    charm_fields_mult(data->par.c, 2./3., c_count);
    charm_fields_axpy(data->par.c_old, data->par.c, 1./3., c_count);
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

    if (!ctx->timestep) {
        charm_model_ns_geom_calc(p4est);
    }
    if (refine_period) {
        if (!(ctx->timestep % refine_period)) {
            if (ctx->timestep) {
                ctx->amr_fn(p4est, ghost, ghost_data); /* adapt */
                charm_model_ns_geom_calc(p4est); //@todo выяснить нужно или нет
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

    p4est_iterate (p4est, NULL, NULL, charm_model_ns_timestep_rk_0, NULL, NULL, NULL);
    CHARM_RUNGE_KUTTA_STEP()
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    CHARM_RUNGE_KUTTA_STEP()
    p4est_iterate (p4est, NULL, NULL, charm_model_ns_timestep_rk_1, NULL, NULL, NULL);
    CHARM_RUNGE_KUTTA_STEP()
    p4est_iterate (p4est, NULL, NULL, charm_model_ns_timestep_rk_2, NULL, NULL, NULL);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

#undef CHARM_RUNGE_KUTTA_STEP

    *_ghost       = ghost;
    *_ghost_data  = ghost_data;

}


