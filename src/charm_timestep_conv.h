//
// Created by zhrv on 15.09.18.
//

#ifndef CHARM_DG_CHARM_TIMESTEP_CONV_H
#define CHARM_DG_CHARM_TIMESTEP_CONV_H
#include "charm_globals.h"

void charm_timestep_conv(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data, double *dt);

#endif //CHARM_DG_CHARM_TIMESTEP_CONV_H
