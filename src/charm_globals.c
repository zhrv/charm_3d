//
// Created by zhrv on 26.10.17.
//

#include "charm_globals.h"

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
