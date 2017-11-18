//
// Created by zhrv on 26.10.17.
//

#include "charm_bnd_cond.h"
#include "charm_grad.h"

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
                    charm_param_t *par_in, charm_param_t *par_out)
{
    charm_tree_attr_t *attr = charm_get_tree_attr(p4est, treeid);
    charm_bnd_t *bnd = attr->bnd[face];
    P4EST_ASSERT(bnd);
    charm_param_cons_to_prim(attr->reg->mat, par_in);
    bnd->bnd_fn(par_in, par_out, face, bnd->params);
    charm_mat_eos(attr->reg->mat, par_out, 2);
    charm_mat_eos(attr->reg->mat, par_out, 1);
    charm_param_prim_to_cons(attr->reg->mat, par_out);
}


void charm_bnd_cond_fn_inlet(charm_param_t *par_in, charm_param_t *par_out, int8_t face, double* param)
{
    P4EST_ASSERT(param);

    par_out->p.u = param[0];
    par_out->p.v = param[1];
    par_out->p.w = param[2];
    par_out->p.t = param[3];
    par_out->p.p = param[4];
}

void charm_bnd_cond_fn_outlet(charm_param_t *par_in, charm_param_t *par_out, int8_t face, double* param)
{
#define REAL double
    charm_prim_cpy(par_out, par_in);
    REAL fG = 1.4;
    REAL fQ = par_in->p.u*par_in->p.u+par_in->p.v*par_in->p.v+par_in->p.w*par_in->p.w;
    REAL fR = fQ + 2.0 * par_in->p.cz / ( fG - 1.0 );
    REAL fS = par_in->p.p / pow( par_in->p.r, fG );

    par_out->p.p = 46066.;
    par_out->p.r = pow( par_out->p.p / fS, 1.0 / fG );
    par_out->p.e = par_out->p.p / ( par_out->p.r * ( fG - 1.0 ) );

    REAL fBeta = ( fR - 2.0 / ( fG - 1.0 ) * sqrt( fG * par_out->p.p / par_out->p.r ) ) / fQ;

    par_out->p.u = fBeta * par_in->p.u;
    par_out->p.v = fBeta * par_in->p.v;
    par_out->p.w = fBeta * par_in->p.w;
#undef REAL

}

void charm_bnd_cond_fn_wall_slip(charm_param_t *par_in, charm_param_t *par_out, int8_t face, double* param)
{
    int i;
    double *n    = par_in->g.n[face];
    double  v[3] = {par_in->p.u, par_in->p.v, par_in->p.w};

    charm_prim_cpy(par_out, par_in);

    double   svn = scalar_prod( v, n );
    double   vv[3] = {n[0]*svn, n[1]*svn, n[3]*svn};
    for (i = 0; i < 3; i++) {
        v[i] -= vv[i];
    }
    par_out->p.u = v[0];
    par_out->p.v = v[1];
    par_out->p.w = v[2];
}


// @todo
void charm_bnd_cond_fn_wall_no_slip(charm_param_t *par_in, charm_param_t *par_out, int8_t face, double* param)
{
    int i;
    double *n    = par_in->g.n[face];
    double  v[3] = {par_in->p.u, par_in->p.v, par_in->p.w};

    charm_prim_cpy(par_out, par_in);
    par_out->p.u = 0.;
    par_out->p.v = 0.;
    par_out->p.w = 0.;
}
