//
// Created by zhrv on 26.10.17.
//

#ifndef CHAMR_3D_CHARM_AMR_H
#define CHAMR_3D_CHARM_AMR_H
#include "charm_globals.h"
#include "charm_init.h"

void charm_adapt_init(p4est_t *p4est);
void charm_adapt(p4est_t *p4est, p4est_ghost_t **ghost, charm_data_t **ghost_data);

#endif //CHAMR_3D_CHARM_AMR_H
