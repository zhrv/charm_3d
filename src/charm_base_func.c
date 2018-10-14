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


double charm_get_field_rh(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.rh[i]*charm_base_func(x, i, p);
    }
    return result;
}


double charm_get_field_p(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.p[i]*charm_base_func(x, i, p);
    }
    return (result > CHARM_EPS ? result : 0.);
}


double charm_get_field_p0(charm_data_t* p)
{
    return p->par.p0;;
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


double charm_get_field_grad_p_x(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.grad_p[0][i]*charm_base_func(x, i, p);
    }
    return result;
}


double charm_get_field_grad_p_y(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.grad_p[1][i]*charm_base_func(x, i, p);
    }
    return result;
}


double charm_get_field_grad_p_z(charm_data_t* p, double* x)
{
    double result = 0.;
    int i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        result += p->par.c.grad_p[2][i]*charm_base_func(x, i, p);
    }
    return result;
}


void charm_get_fields(charm_data_t* p, double* x, charm_cons_t* c){
    int k;
    size_t c_count = CHARM_MAX_COMPONETS_COUNT; // @todo fix by real components count
    c->ru = charm_get_field_ru(p, x);
    c->rv = charm_get_field_rv(p, x);
    c->rw = charm_get_field_rw(p, x);
    c->rh = charm_get_field_rh(p, x);
    c->p  = charm_get_field_p(p, x);
    c->p0 = charm_get_field_p0(p);
    for (k = 0; k < c_count; k++) {
        c->rc[k] = charm_get_field_rc(p, x, k);
    }

    c->mat_id = p->par.mat_id;
}


void charm_get_fields_avg(charm_data_t* p, charm_cons_t* c)
{
    int k, igp;
    size_t c_count = CHARM_MAX_COMPONETS_COUNT; // @todo fix by real components count
    double *gx, gjw;

    c->ru = 0.;
    c->rv = 0.;
    c->rw = 0.;
    c->rh = 0.;
    c->p  = 0.;
    for (k = 0; k < c_count; k++) {
        c->rc[k] = 0.;
    }

    for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
        gx = p->par.g.quad_gp[igp];
        gjw = p->par.g.quad_gj[igp]*p->par.g.quad_gw[igp]/p->par.g.volume;
        c->ru += charm_get_field_ru(p, gx) * gjw;
        c->rv += charm_get_field_rv(p, gx) * gjw;
        c->rw += charm_get_field_rw(p, gx) * gjw;
        c->rh += charm_get_field_rh(p, gx) * gjw;
        c->p  += charm_get_field_p (p, gx) * gjw;
        for (k = 0; k < c_count; k++) {
            c->rc[k] += charm_get_field_rc(p, gx, k) * gjw;
        }
    }
    c->p0 = charm_get_field_p0(p);
    c->mat_id = p->par.mat_id;
}
