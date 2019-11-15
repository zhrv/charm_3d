//
// Created by appmath on 29.09.18.
//

#ifndef CHARM_DG_CHARM_LIMITER_H
#define CHARM_DG_CHARM_LIMITER_H

#include "charm_globals.h"
#ifdef __cplusplus
extern "C" {
#endif

void charm_limiter(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data);

void charm_limiter_bj(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data);

#ifdef __cplusplus
}
#endif

#endif //CHARM_DG_CHARM_LIMITER_H
