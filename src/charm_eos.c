//
// Created by zhrv on 05.10.18.
//

#include "charm_eos.h"


double gR = 8.31446261815324; // Дж/(моль К)

double charm_eos_get_r() { return gR; }

static void _calc_t( p4est_t * p4est, charm_prim_t * p )
{
    charm_comp_t * comp;
    double tt = 1.0;		// начальное приближение для температуры
    double e  = p->e;		// энергия
    int c_count = charm_get_comp_count(p4est);
    int i, j;
    double rm, cp, cp_dt, ft, ft_dt, tt1;

    double m_ = 0.0;		//  формулы   M_ = 1 / M   где   M = 1 / SUM( Ci / Mi )


    for (i=0; i<c_count; i++)	{
        comp = charm_get_comp(p4est, i);
        m_ += p->c[i] / comp->m;
    }
    rm = gR * m_;


    for (j=0; j<100; j++ )
    {
        cp = 0.0;
        cp_dt = 0.0;
        for (i=0; i < c_count; i++)
        {
            comp = charm_get_comp(p4est, i);
            cp  += p->c[i] * charm_comp_calc_cp(comp, tt);
            cp_dt += p->c[i] * charm_comp_calc_cp_dt(comp, tt);
        }


        ft    = ( cp - rm ) * tt - e;
        ft_dt = ( cp - rm ) + cp_dt * tt;

        tt1 = tt - ft / ft_dt;
        if ( pow( tt1 - tt, 2 ) < CHARM_EPS ) {
            tt = tt1;
            break;
        }
        tt = tt1;
    }
    p->t = tt;
} // Newton_Method


static void _charm_mat_eos_switch(charm_prim_t * p, int flag) {
    switch (flag)
    {
        case 0:		// p=p(r,e)
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            p->p = p->r*p->e*(p->gam-1);
            p->cz = sqrt(p->gam*p->p/p->r);
            break;

        case 1:		// e=e(r,p)
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            if (p->p < CHARM_EPS) p->p = CHARM_EPS;
            p->e = p->p/(p->r*(p->gam-1));
            p->t = p->e/p->cv;
            break;

        case 2:		// r=r(T,p)
            if (p->p < CHARM_EPS) p->p = CHARM_EPS;
            p->r = p->p*p->m/(p->t*gR);
            p->cz = sqrt(p->gam*p->p/p->r);
            break;

        case 3: // (T,p) => (r, cz, e)
            if (p->p < CHARM_EPS) p->p = CHARM_EPS;
            p->r  = p->p*p->m/(p->t*gR);
            p->cz = sqrt(p->gam*p->p/p->r);
            p->e  = p->p/(p->r*(p->gam-1));
            break;

        case 4: // (r,e) => (p, cz, T)
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            p->p  = p->r*p->e*(p->gam-1);
            p->cz = sqrt(p->gam*p->p/p->r);
            p->t  = p->e/p->cv;
            break;

        default:
                    CHARM_ASSERT(flag < 5);
    }
}

void charm_mat_eos_ideal(p4est_t * p4est, charm_prim_t * p, int flag)
{
    int id;
    charm_ctx_t  *ctx = (charm_ctx_t *)p4est->user_pointer;
    charm_comp_t *comp = sc_array_index(ctx->comp, 0);
    double t;

    if (comp->cp_type == COMP_CP_POLYNOM && (flag == 1 || flag == 4)) {
        _calc_t(p4est, p);
        t = p->t;
    }
    else {
        t = p->t;
    }

    p->cp  = charm_comp_calc_cp(comp, t);
    p->m   = comp->m;
    p->cv  = p->cp-gR/p->m;
    p->gam = p->cp/p->cv;
    _charm_mat_eos_switch(p, flag);
}

void charm_mat_eos_mix(p4est_t * p4est, charm_prim_t * p, int flag)
{
    charm_ctx_t *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_mat_t *mat = charm_mat_find_by_id(ctx, p->mat_id);
    size_t c_count   = charm_get_comp_count(p4est);
    int i;
    charm_comp_t *comp = charm_get_comp(p4est, 0);
    double t;

    if (comp->cp_type == COMP_CP_POLYNOM && (flag == 1 || flag == 4)) {
        _calc_t(p4est, p);
        t = p->t;
    }
    else {
        t = p->t;
    }

    p->cp  = 0.;
    double M_  = 0.;
    for (i = 0; i < c_count; i++) {
        comp = charm_get_comp(p4est, i);
        M_ += p->c[i]/comp->m;
        p->cp += p->c[i]*charm_comp_calc_cp(comp, t);
    }
    p->m   = 1./M_;
    p->cv  = p->cp-gR/p->m;
    p->gam = p->cp/p->cv;
    _charm_mat_eos_switch(p, flag);
}

void charm_mat_eos_table(p4est_t * p4est, charm_prim_t * p, int flag)
{
    CHARM_LERROR("TABLE EOS IS NOT RELEASED\n");
    charm_abort(NULL, 1);
}


