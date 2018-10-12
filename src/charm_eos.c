//
// Created by zhrv on 05.10.18.
//

#include "charm_eos.h"


double gR = 8.314472;

double charm_comp_calc_cp(p4est_t *p4est, charm_comp_t * comp, double t)
{
    double *cp = sc_array_index(comp->cp, 0);
    return *cp;
}


void charm_mat_eos_ideal(p4est_t * p4est, charm_prim_t * p, int flag)
{
    int id;
    charm_ctx_t  *ctx = (charm_ctx_t *)p4est->user_pointer;
    charm_comp_t *comp = sc_array_index(ctx->comp, 0);

    double Cp = charm_comp_calc_cp(p4est, comp, 273.); // @todo temperature
    double M  = comp->m;
    double Cv = Cp-gR/M;
    double gam = Cp/Cv;
    p->gam = gam;
    p->cp = Cp;
    p->cv = Cv;
    switch (flag)
    {
        case 0:		// p=p(r,e)
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            CHARM_ASSERT(0); // p->p = p->r*p->e*(gam-1);
            p->cz = sqrt(gam*p->p0/p->r);
            break;

        case 1:		// e=e(r,p)
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            CHARM_ASSERT(p->p0 < CHARM_EPS); // p->p0 = CHARM_EPS;
            p->e = p->p0/(p->r*(gam-1));
            p->t = p->e/Cv;
            break;

        case 2:		// r=r(T,p)
            CHARM_ASSERT(p->p0 < CHARM_EPS); // p->p0 = CHARM_EPS;
            p->r = p->p0*M/(p->t*gR);
            p->cz = sqrt(gam*p->p0/p->r);
            break;

        case 3:
            CHARM_ASSERT(p->p0 < CHARM_EPS); // p->p0 = CHARM_EPS;
            p->r  = p->p0*M/(p->t*gR);
            p->cz = sqrt(gam*p->p0/p->r);
            p->e  = p->p0/(p->r*(gam-1));
            break;

        case 4:
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            CHARM_ASSERT(0); // p->p  = p->r*p->e*(gam-1);
            p->cz = sqrt(gam*p->p0/p->r);
            p->t  = p->e/Cv;
            break;

        default:
            CHARM_ASSERT(flag < 5);
    }

}

void charm_mat_eos_mix(p4est_t * p4est, charm_prim_t * p, int flag)
{
    charm_ctx_t *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_mat_t *mat = charm_mat_find_by_id(ctx, p->mat_id);
    size_t c_count   = charm_get_comp_count(p4est);
    int i;
    charm_comp_t *comp;

    double Cp  = 0.;
    double M_  = 0.;
    for (i = 0; i < c_count; i++) {
        comp = charm_get_comp(p4est, i);
        M_ += p->c[i]/comp->m;
        Cp += p->c[i]*charm_comp_calc_cp(p4est, comp, 273.); // @todo temperature
    }
    double M   = 1./M_;
    double Cv  = Cp-gR/M;
    double gam = Cp/Cv;
    p->gam     = gam;
    p->cp      = Cp;
    p->cv      = Cv;
    switch (flag)
    {
        case 0:		// p=p(r,e)
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            CHARM_ASSERT(0); // p->p = p->r*p->e*(gam-1);
            p->cz = sqrt(gam*p->p0/p->r);
            break;

        case 1:		// e=e(r,p)
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            CHARM_ASSERT(p->p0 < CHARM_EPS); // p->p = CHARM_EPS;
            p->e = p->p0/(p->r*(gam-1));
            p->t = p->e/Cv;
            break;

        case 2:		// r=r(T,p)
            CHARM_ASSERT(p->p0 < CHARM_EPS); // p->p = CHARM_EPS;
            p->r = p->p0*M/(p->t*gR);
            p->cz = sqrt(gam*p->p0/p->r);
            break;

        case 3:
            CHARM_ASSERT(p->p0 < CHARM_EPS); // p->p = CHARM_EPS;
            p->r  = p->p0*M/(p->t*gR);
            p->cz = sqrt(gam*p->p0/p->r);
            p->e  = p->p0/(p->r*(gam-1));
            break;

        case 4:
            if (p->r < CHARM_EPS) p->r = CHARM_EPS;
            CHARM_ASSERT(0); // p->p  = p->r*p->e*(gam-1);
            p->cz = sqrt(gam*p->p0/p->r);
            p->t  = p->e/Cv;
            break;

        default:
            CHARM_ASSERT(flag < 5);
    }
}

void charm_mat_eos_table(p4est_t * p4est, charm_prim_t * p, int flag)
{
    CHARM_LERROR("TABLE EOS IS NOT RELEASED\n");
    charm_abort(1);
}


