//
// Created by appmath on 29.08.18.
//

#include "charm_base_func.h"



double charm_base_func(double* x, int k, p4est_quadrant_t* q) {
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
            return _SQR_(x[0]-c[0])/p->par.g.dh[0];
        case 5:
            return (x[1]-c[1])/p->par.g.dh[1];
        case 6:
            return (x[2]-c[2])/p->par.g.dh[2];
    }
}


double charm_base_func_dx(double x, int k, p4est_quadrant_t* q) {
    charm_data_t *p = (charm_data_t*)q->p.user_data;
    double * c = p->par.g.c;
    switch (k) {
        case 0:
            return 0.;
        case 1:
            return 1./p->par.g.dh[0];
        case 2:
            return 1./p->par.g.dh[1];
        case 3:
            return 1./p->par.g.dh[2];
    }
}