//
// Created by zhrv on 26.10.17.
//

#include "charm_bnd_cond.h"

charm_bnd_types_t charm_bnd_type_by_name(const char* name) {
    int i = 0;
    while (charm_bnd_types[i] != NULL) {
        if (strcmp(charm_bnd_types[i], name) == 0) {
            return (charm_bnd_types_t)i;
        }
        i++;
    }
    return BOUND_UNKNOWN;
}


void charm_bnd_cond(p4est_t* p4est, p4est_topidx_t treeid, int8_t face,
                    charm_prim_t *par_in, charm_prim_t *par_out, charm_real_t n[CHARM_DIM])
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    charm_tree_attr_t *attr = charm_get_tree_attr(p4est, treeid);
    charm_bnd_t *bnd = attr->bnd[face];
    charm_mat_t *mat = charm_mat_find_by_id(ctx, par_in->mat_id);
    CHARM_ASSERT(bnd);
    bnd->bnd_fn(par_in, par_out, face, bnd->params, n);
    par_out->mat_id = par_in->mat_id;
    mat->eos_fn(p4est, par_out, 3); // (T,p) => (r, cz, e)
    par_out->e_tot = par_out->e+0.5*_MAG_(par_out->u, par_out->v, par_out->w);
}


void charm_bnd_cond_fn_inlet(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, charm_real_t* param, charm_real_t n[CHARM_DIM])
{
    CHARM_ASSERT(param);

    charm_prim_cpy(par_out, par_in); //@todo

    par_out->u = param[0];
    par_out->v = param[1];
    par_out->w = param[2];
    par_out->t = param[3];
    par_out->p = param[4];
}


void charm_bnd_cond_fn_outlet(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, charm_real_t* param, charm_real_t n[CHARM_DIM])
{
    charm_prim_cpy(par_out, par_in);
}


void charm_bnd_cond_fn_wall_slip(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, charm_real_t* param, charm_real_t n[CHARM_DIM])
{
    int i;
    charm_real_t  v[3] = {par_in->u, par_in->v, par_in->w};

    charm_prim_cpy(par_out, par_in);

    charm_real_t   svn = scalar_prod( v, n );
    charm_real_t   vv[3] = {n[0]*svn, n[1]*svn, n[2]*svn};
    for (i = 0; i < 3; i++) {
        v[i] -= vv[i];
        v[i] -= vv[i];
    }
    par_out->u = v[0];
    par_out->v = v[1];
    par_out->w = v[2];
}


// @todo
void charm_bnd_cond_fn_wall_no_slip(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, charm_real_t* param, charm_real_t n[CHARM_DIM])
{
    int i;
    charm_real_t  v[3] = {par_in->u, par_in->v, par_in->w};

    charm_prim_cpy(par_out, par_in);

    charm_real_t   svn = scalar_prod( v, n );
    charm_real_t   vv[3] = {n[0]*svn, n[1]*svn, n[2]*svn};
    for (i = 0; i < 3; i++) {
        v[i] -= vv[i];
        // v[i] -= vv[i];
    }
    par_out->u = v[0];
    par_out->v = v[1];
    par_out->w = v[2];
    par_out->t = param[0];
}


void charm_bnd_cond_fn_massflow(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, charm_real_t* param, charm_real_t n[CHARM_DIM])
{
}
