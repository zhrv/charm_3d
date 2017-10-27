//
// Created by zhrv on 26.10.17.
//

#include "charm_init.h"
#include "charm_geom.h"
#include "yaml.h"


void charm_initial_condition (double x[], double u[FLD_COUNT], double du[FLD_COUNT][P4EST_DIM], charm_ctx_t * ctx)
{
    int                 i;
    double r_,p_,u_,v_,w_,e_;

    double pi = 4.0*atan(1.0);
    double pi2 = pi*2.0;

    if (x[2] < 0.0) {
        r_ = 1.;
        u_ = 0.;
        v_ = 0.;
        w_ = 0.;
        p_ = 1.;
    }
    else {
        r_ = 0.125;
        u_ = 0.;
        v_ = 0.;
        w_ = 0.;
        p_ = 0.1;
    }
    u[0] = r_;
    u[1] = r_*u_;
    u[2] = r_*v_;
    u[3] = r_*w_;
    u[4] = p_/0.4+0.5*r_*(u_*u_+v_*v_+w_*w_);

    if (du) {
        for (i = 0; i < P4EST_DIM; i++) {
            du[0][i] = 0.0;
            du[1][i] = 0.0;
            du[2][i] = 0.0;
            du[3][i] = 0.0;
            du[4][i] = 0.0;
        }
    }
}

void charm_init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree,
                                          p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    double              midpoint[3];
    int i;

    double du[FLD_COUNT][P4EST_DIM], u[FLD_COUNT];

    charm_geom_quad_calc(p4est, q, which_tree);

    charm_quad_get_center (q, midpoint);
    charm_initial_condition (midpoint, u, du, ctx);

    data->par.c.ro = u[0];
    data->par.c.ru = u[1];
    data->par.c.rv = u[2];
    data->par.c.rw = u[3];
    data->par.c.re = u[4];

//    for (i = 0; i < P4EST_DIM; i++) {
//        data->dro[i] = du[0][i];
//        data->dru[i] = du[1][i];
//        data->drv[i] = du[2][i];
//        data->drw[i] = du[3][i];
//        data->dre[i] = du[4][i];
//    }
}


void charm_init_context(charm_ctx_t *ctx)
{
    ctx->max_err                = 1.e-2;

    ctx->refine_period          = 10;
    ctx->repartition_period     = 100;

    ctx->min_level              = 2;
    ctx->allowed_level          = 2;

    ctx->write_period           = 1;

    ctx->v_ref                  = 20.;
    ctx->CFL                    = 0.1;

    ctx->dt                     = 1.e-8;

    ctx->time                   = 3.5e-3;
}


