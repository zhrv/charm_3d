//
// Created by zhrv on 26.10.17.
//

#ifndef CHAMR_3D_CHARM_BND_COND_H
#define CHAMR_3D_CHARM_BND_COND_H

#include "charm_globals.h"


void charm_bnd_cond_fn_inlet(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, double* param, double n[CHARM_DIM]);

void charm_bnd_cond_fn_outlet(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, double* param, double n[CHARM_DIM]);

void charm_bnd_cond_fn_wall_slip(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, double* param, double n[CHARM_DIM]);

void charm_bnd_cond_fn_wall_no_slip(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, double* param, double n[CHARM_DIM]);

void charm_bnd_cond(p4est_t* p4est, p4est_topidx_t treeid, int8_t face,
                    charm_prim_t *par_in, charm_prim_t *par_out, double n[CHARM_DIM]);


int charm_bnd_type_by_name(const char* name);

#endif //CHAMR_3D_CHARM_BND_COND_H
