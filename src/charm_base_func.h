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

double charm_get_field_ro(p4est_quadrant_t* q, double* x);
double charm_get_field_ru(p4est_quadrant_t* q, double* x);
double charm_get_field_rv(p4est_quadrant_t* q, double* x);
double charm_get_field_rw(p4est_quadrant_t* q, double* x);
double charm_get_field_re(p4est_quadrant_t* q, double* x);
double charm_get_field_rc(p4est_quadrant_t* q, double* x, int k);

void charm_get_fields(p4est_quadrant_t* q, double* x, charm_cons_t* c);



#endif //CHARM_DG_CHARM_BASE_FUNC_H
