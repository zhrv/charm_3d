//
// Created by zhrv on 26.10.17.
//
#define GLOBAL_H_FILE
#include "charm_globals.h"

const char *charm_bnd_types[] ={
        "BOUND_INLET",
        "BOUND_OUTLET",
        "BOUND_WALL_SLIP",
        "BOUND_WALL_NO_SLIP",
        NULL
};


double scalar_prod(double v1[3], double v2[3])
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}


double vect_length(double v[3])
{
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

void vect_prod(double v1[3], double v2[3], double res[3])
{
    res[0] =  v1[1]*v2[2]-v1[2]*v2[1];
    res[1] = -v1[0]*v2[2]+v1[2]*v2[0];
    res[2] =  v1[0]*v2[1]-v1[1]*v2[0];
}


double charm_face_get_normal(p4est_quadrant_t* q, int8_t face, double* n)
{
    charm_data_t* d = (charm_data_t*) q->p.user_data;
    memcpy(n, d->par.g.n[face], 3*sizeof(double));
    return d->par.g.area[face];
}

void charm_quad_get_center(p4est_quadrant_t* q, double* c)
{
    charm_data_t* d = (charm_data_t*) q->p.user_data;
    memcpy(c, d->par.g.c, 3*sizeof(double));
}

void charm_face_get_center(p4est_quadrant_t* q, int8_t face, double* c)
{
    charm_data_t* d = (charm_data_t*) q->p.user_data;
    memcpy(c, d->par.g.fc[face], 3*sizeof(double));
}

double charm_quad_get_volume(p4est_quadrant_t* q)
{
    charm_data_t* d = (charm_data_t*) q->p.user_data;
    return d->par.g.volume;
}


charm_mat_t * charm_mat_find_by_id(charm_ctx_t *ctx, int id)
{
    int i;
    sc_array_t *arr = ctx->mat;
    charm_mat_t * mat;


    for (i = 0; i < arr->elem_count; i++) {
        mat = sc_array_index(arr, i);
        if (mat->id == id) {
            return mat;
        }
    }

    return NULL;
}


charm_bnd_t * charm_bnd_find_by_face_type(charm_ctx_t *ctx, int type)
{
    int i;
    sc_array_t *arr = ctx->bnd;
    charm_bnd_t * bnd;


    for (i = 0; i < arr->elem_count; i++) {
        bnd = sc_array_index(arr, i);
        if (bnd->face_type == type) {
            return bnd;
        }
    }

    return NULL;
}


charm_mesh_type_t charm_mesh_get_type_by_str(char *str)
{
    if (strcmp(str, "gmsh_msh") == 0) {
        return CHARM_MESH_GMSH_MSH;
    }
    else if (strcmp(str, "gmsh_unv") == 0) {
        return CHARM_MESH_GMSH_UNV;
    }
    else if (strcmp(str, "gmsh_inp") == 0) {
        return CHARM_MESH_GMSH_INP;
    }
    else {
        return CHARM_MESH_UNKNOWN;
    }
}

charm_tree_attr_t * charm_get_tree_attr(p4est_t * p4est, p4est_topidx_t which_tree)
{
    p4est_connectivity_t * conn = p4est->connectivity;
    return (charm_tree_attr_t *)&(conn->tree_to_attr[sizeof(charm_tree_attr_t)*which_tree]);
}

double gR = 8.314472;
void charm_mat_eos(charm_mat_t * mat, charm_param_t * p, int flag)
{
    double Cp = mat->cp;
    double M  = mat->m;
    double Cv = Cp-gR/M;
    double gam = Cp/Cv;
    switch (flag)
    {
        case 0:		// p=p(r,e)
            p->p.p = p->p.r*p->p.e*(gam-1);
            p->p.cz = sqrt(gam*p->p.p/p->p.r);
            break;

        case 1:		// e=e(r,p)
            p->p.e = p->p.p/(p->p.r*(gam-1));
            p->p.t = p->p.e/Cv;
            break;

        case 2:		// r=r(T,p)
            p->p.r = p->p.p*M/(p->p.t*gR);
            p->p.cz = sqrt(gam*p->p.p/p->p.r);
            break;
    }

}


void charm_param_cons_to_prim(charm_mat_t * mat, charm_param_t * p)
{
    p->p.r      = p->c.ro;
    p->p.u      = p->c.ru/p->c.ro;
    p->p.v      = p->c.rv/p->c.ro;
    p->p.e_tot  = p->c.re/p->c.ro;
    p->p.e      = p->p.e_tot-0.5*(p->p.u*p->p.u+p->p.v*p->p.v+p->p.w*p->p.w);

    charm_mat_eos(mat, p, 0);
    charm_mat_eos(mat, p, 1);
}


void charm_param_prim_to_cons(charm_mat_t * mat, charm_param_t * p)
{
    p->c.ro = p->p.r;
    p->c.ru = p->p.r*p->p.u;
    p->c.rv = p->p.r*p->p.v;
    p->c.rw = p->p.r*p->p.w;
    p->c.re = p->p.r*(p->p.e+0.5*(p->p.u*p->p.u+p->p.v*p->p.v+p->p.w*p->p.w));
}


void charm_prim_cpy(charm_param_t * dest, charm_param_t * src)
{
    dest->p.r = src->p.r;
    dest->p.p = src->p.p;
    dest->p.u = src->p.u;
    dest->p.v = src->p.v;
    dest->p.w = src->p.w;
    dest->p.t = src->p.t;
    dest->p.cz = src->p.cz;
    dest->p.e = src->p.e;
    dest->p.e_tot = src->p.e_tot;
}