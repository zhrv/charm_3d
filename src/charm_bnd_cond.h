//
// Created by zhrv on 26.10.17.
//

#ifndef CHAMR_3D_CHARM_BND_COND_H
#define CHAMR_3D_CHARM_BND_COND_H

#include "charm_globals.h"


void charm_bnd_cond_fn_inlet(double ro, double ru, double rv, double rw, double re,
                             double* ro_, double* ru_, double* rv_, double* rw_, double* re_,
                             double* n, double* param);

void charm_bnd_cond_fn_outlet(double ro, double ru, double rv, double rw, double re,
                              double* ro_, double* ru_, double* rv_, double* rw_, double* re_,
                              double* n, double* param);

void charm_bnd_cond_fn_wall(double ro, double ru, double rv, double rw, double re,
                            double* ro_, double* ru_, double* rv_, double* rw_, double* re_,
                            double* n, double* param);



void charm_bnd_cond(p4est_t* p4est, p4est_topidx_t treeid, int8_t face,
                    double ro, double ru, double rv, double rw, double re,
                    double* ro_, double* ru_, double* rv_, double* rw_, double* re_);


#endif //CHAMR_3D_CHARM_BND_COND_H
