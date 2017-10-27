//
// Created by zhrv on 27.10.17.
//

#ifndef CHAMR_3D_CHARM_GEOM_H
#define CHAMR_3D_CHARM_GEOM_H

#include "charm_globals.h"

void charm_geom_quad_calc(p4est_t * p4est, p4est_quadrant_t* q, p4est_topidx_t treeid);
void charm_geom_calc(p4est_t * p4est);

#endif //CHAMR_3D_CHARM_GEOM_H
