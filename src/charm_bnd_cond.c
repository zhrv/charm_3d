//
// Created by zhrv on 26.10.17.
//

#include "charm_bnd_cond.h"

int charm_bnd_type_by_name(const char* name) {
    int i = 0;
    while (charm_bnd_types[i] != NULL) {
        if (strcmp(charm_bnd_types[i], name) == 0) {
            return i;
        }
        i++;
    }
    return -1;
}


void charm_bnd_cond(p4est_t* p4est, p4est_topidx_t treeid, int8_t face,
                    double  ro,  double  ru,  double  rv,  double  rw,  double  re,
                    double* ro_, double* ru_, double* rv_, double* rw_, double* re_)
{
    p4est_topidx_t f_type;
    *ro_ = ro;
    *ru_ = ru;
    *rv_ = rv;
    *rw_ = rw;
    *re_ = re;
    //f_type = charm_conn_get
}


void charm_bnd_cond_fn_inlet(double ro, double ru, double rv, double rw, double re,
                             double* ro_, double* ru_, double* rv_, double* rw_, double* re_,
                             double* n, double* param)
{

}

void charm_bnd_cond_fn_outlet(double ro, double ru, double rv, double rw, double re,
                              double* ro_, double* ru_, double* rv_, double* rw_, double* re_,
                              double* n, double* param)
{

}

void charm_bnd_cond_fn_wall(double ro, double ru, double rv, double rw, double re,
                            double* ro_, double* ru_, double* rv_, double* rw_, double* re_,
                            double* n, double* param)
{

}
