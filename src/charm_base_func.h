//
// Created by appmath on 29.08.18.
//

#ifndef CHARM_DG_CHARM_BASE_FUNC_H
#define CHARM_DG_CHARM_BASE_FUNC_H

#include "charm_globals.h"


double charm_base_func(double* x, int k, p4est_quadrant_t* q);
double charm_base_func_dx(double* x, int k, p4est_quadrant_t* q);
double charm_base_func_dy(double* x, int k, p4est_quadrant_t* q);
double charm_base_func_dz(double* x, int k, p4est_quadrant_t* q);




#endif //CHARM_DG_CHARM_BASE_FUNC_H
