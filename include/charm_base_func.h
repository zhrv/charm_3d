//
// Created by appmath on 29.08.18.
//

#ifndef CHARM_DG_CHARM_BASE_FUNC_H
#define CHARM_DG_CHARM_BASE_FUNC_H

#include "charm_globals.h"


double charm_base_func(double* x, int k, charm_data_t* p);
double charm_base_func_dx(double* x, int k, charm_data_t* p);
double charm_base_func_dy(double* x, int k, charm_data_t* p);
double charm_base_func_dz(double* x, int k, charm_data_t* p);

double charm_get_field_ro(charm_data_t* p, double* x);
double charm_get_field_ru(charm_data_t* p, double* x);
double charm_get_field_rv(charm_data_t* p, double* x);
double charm_get_field_rw(charm_data_t* p, double* x);
double charm_get_field_re(charm_data_t* p, double* x);
double charm_get_field_rh(charm_data_t* p, double* x);
double charm_get_field_p(charm_data_t* p, double* x);
double charm_get_field_rc(charm_data_t* p, double* x, int k);

void charm_get_fields(p4est_t *p4est, charm_data_t* p, double* x, charm_cons_t* c);
void charm_get_fields_avg(p4est_t *p4est, charm_data_t* p, charm_cons_t* c);
void charm_get_fields_arr(p4est_t *p4est, charm_data_t* p, double* fld[5]);

void charm_get_visc_tau(charm_data_t *p, double* x, charm_tensor_t *tau);
void charm_get_heat_q(charm_data_t *p, double* x, double *q);
void charm_get_grad_p(charm_data_t *p, double* x, double *grad_p);

#endif //CHARM_DG_CHARM_BASE_FUNC_H
