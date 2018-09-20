//
// Created by zhrv on 26.10.17.
//

#ifndef CHAMR_3D_CHARM_INIT_H
#define CHAMR_3D_CHARM_INIT_H

#include "charm_globals.h"

void charm_initial_condition (double x[], double u[FLD_COUNT], double du[FLD_COUNT][CHARM_DIM], charm_ctx_t * ctx);
void charm_init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree,
                                   p4est_quadrant_t * q);

void charm_init_context(charm_ctx_t *ctx);


#endif //CHAMR_3D_CHARM_INIT_H
