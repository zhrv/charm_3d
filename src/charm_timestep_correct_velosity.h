//
// Created by appmath on 14.10.18.
//

#ifndef CHARM_DG_CHARM_TIMESTEP_CORRECT_VELOSITY_H
#define CHARM_DG_CHARM_TIMESTEP_CORRECT_VELOSITY_H

#include "charm_globals.h"

void charm_timestep_correct_velosity(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data, double *dt);

#endif //CHARM_DG_CHARM_TIMESTEP_CORRECT_VELOSITY_H
