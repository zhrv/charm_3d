//
// Created by appmath on 29.08.18.
//

#include "charm_base_func.h"



double charm_base_func(double* x, int k, p4est_quadrant_t* q) {
    P4EST_ASSERT(k < CHARM_BASE_FN_COUNT);
    charm_data_t *p = (charm_data_t*)q->p.user_data;
    double * c = p->par.g.c;
    switch (k) {
        case 0:
            return 1.;
        case 1:
            return (x[0]-c[0])/p->par.g.dh[0];
        case 2:
            return (x[1]-c[1])/p->par.g.dh[1];
        case 3:
            return (x[2]-c[2])/p->par.g.dh[2];
        case 4:
            return _SQR_(x[0]-c[0])/_SQR_(p->par.g.dh[0]);
        case 5:
            return _SQR_(x[1]-c[1])/_SQR_(p->par.g.dh[1]);
        case 6:
            return _SQR_(x[2]-c[2])/_SQR_(p->par.g.dh[2]);
        case 7:
            return (x[0]-c[0])*(x[1]-c[1])/(p->par.g.dh[0]*p->par.g.dh[1]);
        case 8:
            return (x[0]-c[0])*(x[2]-c[2])/(p->par.g.dh[0]*p->par.g.dh[2]);
        case 9:
            return (x[1]-c[1])*(x[2]-c[2])/(p->par.g.dh[1]*p->par.g.dh[2]);
        default:
            return 0;
    }
}


double charm_base_func_dx(double* x, int k, p4est_quadrant_t* q) {
    P4EST_ASSERT(k < CHARM_BASE_FN_COUNT);
    charm_data_t *p = (charm_data_t*)q->p.user_data;
    double * c = p->par.g.c;
    switch (k) {
        case 0:
            return 0.;
        case 1:
            return 1./p->par.g.dh[0];
        case 2:
            return 0.;
        case 3:
            return 0.;
        case 4:
            return 2.*(x[0]-c[0])/_SQR_(p->par.g.dh[0]);
        case 5:
            return 0.;
        case 6:
            return 0.;
        case 7:
            return (x[1]-c[1])/(p->par.g.dh[0]*p->par.g.dh[1]);
        case 8:
            return (x[2]-c[2])/(p->par.g.dh[0]*p->par.g.dh[2]);
        case 9:
            return 0.;
        default:
            return 0;
    }
}


double charm_base_func_dy(double* x, int k, p4est_quadrant_t* q) {
    P4EST_ASSERT(k < CHARM_BASE_FN_COUNT);
    charm_data_t *p = (charm_data_t*)q->p.user_data;
    double * c = p->par.g.c;
    switch (k) {
        case 0:
            return 0.;
        case 1:
            return 0.;
        case 2:
            return 1./p->par.g.dh[1];
        case 3:
            return 0.;
        case 4:
            return 0.;
        case 5:
            return 2.*(x[1]-c[1])/_SQR_(p->par.g.dh[1]);
        case 6:
            return 0.;
        case 7:
            return (x[0]-c[0])/(p->par.g.dh[0]*p->par.g.dh[1]);
        case 8:
            return 0.;
        case 9:
            return (x[2]-c[2])/(p->par.g.dh[1]*p->par.g.dh[2]);
        default:
            return 0;
    }
}


double charm_base_func_dz(double* x, int k, p4est_quadrant_t* q) {
    P4EST_ASSERT(k < CHARM_BASE_FN_COUNT);
    charm_data_t *p = (charm_data_t*)q->p.user_data;
    double * c = p->par.g.c;
    switch (k) {
        case 0:
            return 0.;
        case 1:
            return 0.;
        case 2:
            return 0.;
        case 3:
            return 1./p->par.g.dh[2];
        case 4:
            return 0.;
        case 5:
            return 0.;
        case 6:
            return 2.*(x[2]-c[2])/_SQR_(p->par.g.dh[2]);
        case 7:
            return 0.;
        case 8:
            return (x[0]-c[0])/(p->par.g.dh[0]*p->par.g.dh[2]);
        case 9:
            return (x[1]-c[1])/(p->par.g.dh[1]*p->par.g.dh[2]);
        default:
            return 0;
    }
}
