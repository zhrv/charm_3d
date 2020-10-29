//
// Created by appmath on 29.08.18.
//

#ifndef CHARM_DG_CHARM_BASE_FUNC_H
#define CHARM_DG_CHARM_BASE_FUNC_H

#include "charm_globals.h"


charm_real_t charm_base_func(charm_real_t* x, int k, charm_data_t* p);
charm_real_t charm_base_func_dx(charm_real_t* x, int k, charm_data_t* p);
charm_real_t charm_base_func_dy(charm_real_t* x, int k, charm_data_t* p);
charm_real_t charm_base_func_dz(charm_real_t* x, int k, charm_data_t* p);

charm_real_t charm_get_field_ro(charm_data_t* p, charm_real_t* x);
charm_real_t charm_get_field_ru(charm_data_t* p, charm_real_t* x);
charm_real_t charm_get_field_rv(charm_data_t* p, charm_real_t* x);
charm_real_t charm_get_field_rw(charm_data_t* p, charm_real_t* x);
charm_real_t charm_get_field_re(charm_data_t* p, charm_real_t* x);
charm_real_t charm_get_field_rc(charm_data_t* p, charm_real_t* x, int k);

void charm_get_fields(charm_data_t* p, charm_real_t* x, charm_cons_t* c);
void charm_get_fields_avg(charm_data_t* p, charm_cons_t* c);
void charm_get_fields_arr(charm_data_t* p, charm_real_t* fld[5]);

void charm_get_visc_tau(charm_data_t *p, charm_vec_t x, charm_tensor_t *tau);
void charm_get_heat_q(charm_data_t *p, charm_vec_t x, charm_vec_t q);

#endif //CHARM_DG_CHARM_BASE_FUNC_H
