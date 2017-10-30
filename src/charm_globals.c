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