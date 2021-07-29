//
// Created by appmath on 29.09.18.
//

#include "charm_globals.h"



void charm_limiter(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    if (ctx->lim_fn != NULL) {
        ctx->lim_fn(p4est, ghost, ghost_data);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);   /* synchronize the ghost data */
    }
}