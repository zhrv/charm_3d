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
double charm_get_field_rc(charm_data_t* p, double* x, int k);

void charm_get_fields(charm_data_t* p, double* x, charm_cons_t* c);
void charm_get_fields_avg(charm_data_t* p, charm_cons_t* c);
void charm_get_fields_arr(charm_data_t* p, double* fld[5]);


#endif //CHARM_DG_CHARM_BASE_FUNC_H
