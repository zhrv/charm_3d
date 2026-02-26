//
// Created by zhrv on 26.10.17.
//

#include "charm_bnd_cond.h"
#include "charm_fluxes.h"
#include "charm_globals.h"
#include "charm_eos.h"
#include "charm_limiter.h"
#include "charm_models.h"


#ifdef TAYLOR_GREEN
void _charm_init_initial_condition_taylor_green (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);
#endif

void charm_init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *data    = (charm_data_t *) q->p.user_data;
    charm_param_t      *par     = &data->par;
    charm_prim_t        prim;
    charm_cons_t        cons;
    charm_tree_attr_t  *attr;
    charm_reg_t        *reg;
    int                 i;
    size_t              c_count = charm_get_comp_count(p4est);
    charm_mat_t        *mat;

    charm_geom_quad_calc(p4est, q, which_tree);

    attr = charm_get_tree_attr(p4est, which_tree);
    reg = attr->reg;
    prim.mat_id = reg->mat_id;

#ifdef POGGI
    charm_real_t *x   = data->par.g.c;
    charm_real_t pi   = 4.*atan(1.);
    charm_real_t pi2  = pi*2.;
    charm_real_t a0 = 0.5e-3;
    charm_real_t lambda = 1.e-3;

    if (x[2] < -5.1e-3) {
        reg = charm_reg_find_by_id(ctx, 0);
    }
    else if ( x[2] > -a0*(1.-cos(pi2*x[0]/lambda))*(1.-cos(pi2*x[1]/lambda)) ) {
        reg = charm_reg_find_by_id(ctx, 2);
    }
    else {
        reg = charm_reg_find_by_id(ctx, 1);
    }
#endif /* POGGI */

    prim.p   = reg->p;
    prim.t   = reg->t;

    prim.u   = reg->v[0];
    prim.v   = reg->v[1];
    prim.w   = reg->v[2];
    for (i = 0; i < c_count; i++) {
        prim.c[i] = reg->c[i];
    }
    mat = charm_mat_find_by_id(ctx, reg->mat_id);

    mat->eos_fn(p4est, &prim, 3); // (T,p) => (r, cz, e)
    prim.e_tot = prim.e + 0.5*(prim.u*prim.u+prim.v*prim.v+prim.w*prim.w);

    charm_param_prim_to_cons(p4est, &cons, &prim);
    memset(&(par->c), 0, sizeof(par->c));
    par->c.ru[0] = cons.ru;
    par->c.rv[0] = cons.rv;
    par->c.rw[0] = cons.rw;
    par->c.re[0] = cons.re;
    for (i = 0; i < c_count; i++) {
        par->c.rc[i][0] = cons.rc[i];
    }
    par->mat_id = reg->mat_id;
    for (i = 0; i < CHARM_DIM; i++) {
        par->grav[i] = reg->grav[i];
    }

#ifdef TAYLOR_GREEN
    _charm_init_initial_condition_taylor_green(p4est, which_tree, q);
#endif /* TAYLOR_GREEN */

    if (ctx->model_init_cond_fn) ctx->model_init_cond_fn(p4est, which_tree, q);
}


#ifdef TAYLOR_GREEN
void _charm_init_initial_condition_taylor_green (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx     = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *data    = (charm_data_t *) q->p.user_data;
    charm_param_t      *par     = &data->par;
    charm_prim_t        prim;
    charm_cons_t        cons;
    charm_tree_attr_t  *attr    = charm_get_tree_attr(p4est, which_tree);
    charm_reg_t        *reg     = attr->reg;
    int                 i;
    size_t              c_count = charm_get_comp_count(p4est);
    charm_mat_t        *mat;
    charm_comp_t       *comp    = sc_array_index(ctx->comp, 0);

    charm_geom_quad_calc(p4est, q, which_tree);

    prim.mat_id = reg->mat_id;
    mat = charm_mat_find_by_id(ctx, reg->mat_id);

    charm_real_t L    = 0.016/M_PI;
    charm_real_t Re   = 100; // 100, 280, 1600
    charm_real_t T    = 273.;
    charm_real_t mu0  = 1.67e-5;
    charm_real_t Ma   = 0.1;
    charm_real_t Cz   = 337.;
    charm_real_t U0   = Ma*Cz;
    charm_real_t rho  = Re/(U0*L/mu0);
    charm_real_t *x   = par->g.c;
    charm_real_t R    = 8.31446261815324; // Дж/(моль К)
    charm_real_t M    = 0.0280134;
    charm_real_t Cp   = 1038.810682;

    // charm_real_t *cp = sc_array_index(comp->cp, 0);
    // *cp = Cp;
    // comp->id = 0;
    // comp->m = M;
    // comp->ml0 = mu0;
    // comp->kp0 = 0.;
    // comp->t0 = T;
    // comp->ts = T;
    // comp->sig = 0.;
    // comp->ek = 0.;
    // comp->h0 = 0.;
    // comp->cp_type = COMP_CP_CONST;
    // comp->kp_type = COMP_KP_CONST;
    // comp->ml_type = COMP_ML_CONST;

    prim.r   = rho;
    prim.t   = T;
    memset(prim.c, 0, sizeof(prim.c));
    prim.c[0] = 1.;
    mat->eos_fn(p4est, &prim, EOS_R_T_TO_P_CZ_E); 

    prim.p  +=  (rho*U0*U0/16.)*(cos(2.*x[0]/L)+cos(2.*x[1]/L))*(cos(2.*x[2]/L)+2.);
    prim.u   =  U0*sin(x[0]/L)*cos(x[1]/L)*cos(x[2]/L);
    prim.v   = -U0*cos(x[0]/L)*sin(x[1]/L)*cos(x[2]/L);
    prim.w   =  0.;
    mat->eos_fn(p4est, &prim, 3); // (T,p) => (r, cz, e)
    prim.e_tot = prim.e + 0.5*(prim.u*prim.u+prim.v*prim.v+prim.w*prim.w);


    charm_param_prim_to_cons(p4est, &cons, &prim);
    memset(&(par->c), 0, sizeof(par->c));
    par->c.ru[0] = cons.ru;
    par->c.rv[0] = cons.rv;
    par->c.rw[0] = cons.rw;
    par->c.re[0] = cons.re;
    for (i = 0; i < c_count; i++) {
        par->c.rc[i][0] = cons.rc[i];
    }
    par->mat_id = reg->mat_id;
    for (i = 0; i < CHARM_DIM; i++) {
        par->grav[i] = reg->grav[i];
    }

}
#endif /* TAYLOR_GREEN */