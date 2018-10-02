//
// Created by appmath on 29.08.18.
//

#include "charm_base_func.h"


double charm_base_func(double* x, int k, charm_data_t *p) {
    CHARM_ASSERT(k < CHARM_BASE_FN_COUNT);
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


double charm_base_func_dx(double* x, int k, charm_data_t *p) {
    CHARM_ASSERT(k < CHARM_BASE_FN_COUNT);
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


double charm_base_func_dy(double* x, int k, charm_data_t *p) {
    CHARM_ASSERT(k < CHARM_BASE_FN_COUNT);
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


double charm_base_func_dz(double* x, int k, charm_data_t *p) {
    CHARM_ASSERT(k < CHARM_BASE_FN_COUNT);
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

double charm_get_field_ro(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.ro[i]*charm_base_func(x, i, p);
    }
    return result;
}


double charm_get_field_ru(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.ru[i]*charm_base_func(x, i, p);
    }
    return result;
}


double charm_get_field_rv(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.rv[i]*charm_base_func(x, i, p);
    }
    return result;
}


double charm_get_field_rw(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.rw[i]*charm_base_func(x, i, p);
    }
    return result;
}


double charm_get_field_re(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.re[i]*charm_base_func(x, i, p);
    }
    return result;
}


double charm_get_field_rc(charm_data_t* p, double* x, int k)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.rc[k][i]*charm_base_func(x, i, p);
    }
    return result;
}

void charm_get_fields(charm_data_t* p, double* x, charm_cons_t* c){
    int k;
    c->ro = charm_get_field_ro(p, x);
    c->ru = charm_get_field_ru(p, x);
    c->rv = charm_get_field_rv(p, x);
    c->rw = charm_get_field_rw(p, x);
    c->re = charm_get_field_re(p, x);
    for (k = 0; k < CHARM_MAX_COMPONETS_COUNT; k++) {
        c->rc[k] = charm_get_field_rc(p, x, k);
    }
    c->mat_id = p->par.mat_id;
}


void charm_get_fields_arr(charm_data_t* p, double* fld[5])
{
    fld[0] = p->par.c.ro;
    fld[1] = p->par.c.ru;
    fld[2] = p->par.c.rv;
    fld[3] = p->par.c.rw;
    fld[4] = p->par.c.re;
}


double charm_get_avg_ro(charm_data_t* p)
{
    double result = 0.;
    int i;
    double vol = charm_quad_get_volume(p);
    double *x;

    for (i = 0; i < CHARM_QUAD_GP_COUNT; i++) {
        x = p->par.g.quad_gp[i];
        result += charm_get_field_ro(p, x)*p->par.g.quad_gw[i]*p->par.g.quad_gj[i];
    }
    return result/vol;
}


double charm_get_avg_ru(charm_data_t* p)
{
    double result = 0.;
    int i;
    double vol = charm_quad_get_volume(p);
    double *x;

    for (i = 0; i < CHARM_QUAD_GP_COUNT; i++) {
        x = p->par.g.quad_gp[i];
        result += charm_get_field_ru(p, x)*p->par.g.quad_gw[i]*p->par.g.quad_gj[i];
    }
    return result/vol;
}


double charm_get_avg_rv(charm_data_t* p)
{
    double result = 0.;
    int i;
    double vol = charm_quad_get_volume(p);
    double *x;

    for (i = 0; i < CHARM_QUAD_GP_COUNT; i++) {
        x = p->par.g.quad_gp[i];
        result += charm_get_field_rv(p, x)*p->par.g.quad_gw[i]*p->par.g.quad_gj[i];
    }
    return result/vol;
}


double charm_get_avg_rw(charm_data_t* p)
{
    double result = 0.;
    int i;
    double vol = charm_quad_get_volume(p);
    double *x;

    for (i = 0; i < CHARM_QUAD_GP_COUNT; i++) {
        x = p->par.g.quad_gp[i];
        result += charm_get_field_rw(p, x)*p->par.g.quad_gw[i]*p->par.g.quad_gj[i];
    }
    return result/vol;
}


double charm_get_avg_re(charm_data_t* p)
{
    double result = 0.;
    int i;
    double vol = charm_quad_get_volume(p);
    double *x;

    for (i = 0; i < CHARM_QUAD_GP_COUNT; i++) {
        x = p->par.g.quad_gp[i];
        result += charm_get_field_re(p, x)*p->par.g.quad_gw[i]*p->par.g.quad_gj[i];
    }
    return result/vol;
}

