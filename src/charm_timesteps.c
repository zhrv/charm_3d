//
// Created by zhrv on 26.10.17.
//

#include <p8est_iterate.h>
#include "charm_bnd_cond.h"
#include "charm_base_func.h"
#include "charm_models.h"
#include "charm_limiter.h"
#include "charm_globals.h"
#include "charm_fluxes.h"
#include "charm_amr.h"



/** Timestep
 *
 * Update the state, refine, repartition, and write the solution to file.
 *
 * \param [in,out] p4est the forest, whose state is updated
 * \param [in] time      the end time
 */
void charm_timesteps(p4est_t * p4est)
{
    charm_real_t              dt;
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *ghost_data;
    p4est_ghost_t      *ghost;
    charm_real_t              calc_time;


    ghost = p4est_ghost_new (p4est, CHARM_CONNECT_FULL);
    ghost_data = CHARM_ALLOC (charm_data_t, ghost->ghosts.elem_count);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);


    dt = ctx->get_dt_fn(p4est);

    CHARM_GLOBAL_ESSENTIAL("Starting time steps...\n");
    for (ctx->t = 0., ctx->timestep = 0; ctx->t < ctx->time; ctx->t += dt, ctx->timestep++) {
        calc_time = sc_MPI_Wtime();
        ctx->timestep_single_fn(p4est, &dt, &ghost, &ghost_data);
        calc_time = sc_MPI_Wtime() - calc_time;
        if (ctx->timestep % ctx->log_period == 0) {
            charm_log_statistics(p4est, ctx->timestep+1, ctx->t+dt, dt, calc_time);
        }
    }

    CHARM_FREE (ghost_data);
    p4est_ghost_destroy (ghost);
}


