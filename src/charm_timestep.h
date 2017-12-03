//
// Created by zhrv on 26.10.17.
//

#ifndef CHAMR_3D_CHARM_TIMESTEP_H
#define CHAMR_3D_CHARM_TIMESTEP_H
#include "charm_globals.h"
#include "charm_fluxes.h"
#include "charm_amr.h"
#include "charm_output.h"


void charm_timesteps (p4est_t * p4est, double time);
void charm_calc_grad(p4est_t * p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data);


#endif //CHAMR_3D_CHARM_TIMESTEP_H
