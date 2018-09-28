//
// Created by zhrv on 26.10.17.
//

#ifndef CHAMR_3D_CHARM_OUTPUT_H
#define CHAMR_3D_CHARM_OUTPUT_H
#include "charm_globals.h"


void charm_write_solution (p4est_t * p4est, int timestep);
void charm_log_statistics(p4est_t * p4est, int timestep, double time, double dt);

#endif //CHAMR_3D_CHARM_OUTPUT_H
