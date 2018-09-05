//
// Created by zhrv on 27.10.17.
//

#include "charm_geom.h"

static double charm_face_calc_normal(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, double* n)
{
    const int ftv[6][4] =
            {{ 0, 2, 4, 6 },
             { 1, 3, 5, 7 },
             { 0, 1, 4, 5 },
             { 2, 3, 6, 7 },
             { 0, 1, 2, 3 },
             { 4, 5, 6, 7 }};

    int i,k;
    double x[4][3], v[2][3], nl;
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    for (i = 0; i < 4; i++) {
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l * (ftv[face][i] % 2), q->y + l * ((ftv[face][i] / 2) % 2), q->z + l * (ftv[face][i] / 4), x[i]);
    }

    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[0][i];
        v[1][i] = x[2][i]-x[0][i];
    }

    vect_prod(v[0], v[1], n);
    nl = vect_length(n);
    for (i = 0; i < 3; i++) {
        n[i] /= nl;
    }

    return nl;
}


static double charm_face_calc_area(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face)
{
    const int ftv[6][4] =
            {{ 0, 2, 4, 6 },
             { 1, 3, 5, 7 },
             { 0, 1, 4, 5 },
             { 2, 3, 6, 7 },
             { 0, 1, 2, 3 },
             { 4, 5, 6, 7 }};

    int i,k;
    double x[4][3], v[2][3], n[3], s1, s2;
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    for (i = 0; i < 4; i++) {
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l * (ftv[face][i] % 2), q->y + l * ((ftv[face][i] / 2) % 2), q->z + l * (ftv[face][i] / 4), x[i]);
    }

    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[0][i];
        v[1][i] = x[2][i]-x[0][i];
    }

    vect_prod(v[0], v[1], n);
    s1 = 0.5*vect_length(n);
    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[3][i];
        v[1][i] = x[2][i]-x[3][i];
    }

    vect_prod(v[0], v[1], n);
    s2 = 0.5*vect_length(n);

    return s1+s2;
}



static void charm_quad_calc_center(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, double* c)
{
    p4est_qcoord_t l2 = P4EST_QUADRANT_LEN(q->level) / 2;
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l2, q->y + l2, q->z + l2, c);


}

static double charm_quad_calc_j(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, double* xyz)
{
    return 0; // @todo
}

static void charm_quad_calc_gp(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, double** gp, double* gw)
{
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(q->level);
    charm_data_t *p = (charm_data_t*)q->p.user_data;
    int i, j;
    double a[3], b[3];
    double sqrt3 = 1./sqrt(3.);
    double t[8][3] = {
            {-sqrt3, -sqrt3, -sqrt3},
            { sqrt3, -sqrt3, -sqrt3},
            {-sqrt3,  sqrt3, -sqrt3},
            { sqrt3,  sqrt3, -sqrt3},
            {-sqrt3, -sqrt3,  sqrt3},
            { sqrt3, -sqrt3,  sqrt3},
            {-sqrt3,  sqrt3,  sqrt3},
            { sqrt3,  sqrt3,  sqrt3}
    };

    p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x, q->y, q->z, a);
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l, q->y + l, q->z + l, b);

    for (i = 0; i < 8; i++) {
        for (j = 0; j < 3; j++) {
            gp[i][j] = 0.5*(a[j]+b[j]+(b[j]-a[j])*t[i][j]);
        }
        gw[i] = charm_quad_calc_j(p4est, q, treeid, gp[i]);
    }

}


static void charm_face_calc_center(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, double* c)
{
    p4est_qcoord_t l  = P4EST_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    p4est_qcoord_t fc[6][3] = {
            {0, l2, l2},
            {l, l2, l2},
            {l2, 0, l2},
            {l2, l, l2},
            {l2, l2, 0},
            {l2, l2, l},
    };
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + fc[face][0], q->y + fc[face][1], q->z + fc[face][2], c);
}

static void charm_face_calc_gp(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, double** gp, double* gw)
{
    p4est_qcoord_t l  = P4EST_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    p4est_qcoord_t fc[6][3] = {
            {0, l2, l2},
            {l, l2, l2},
            {l2, 0, l2},
            {l2, l, l2},
            {l2, l2, 0},
            {l2, l2, l},
    };
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + fc[face][0], q->y + fc[face][1], q->z + fc[face][2], c);
}

double charm_tet_calc_volume(double v[4][3]) // @todo проверить корректность
{
    int8_t i;
    double v1[3], v2[3], v3[3], v12[3];
    for (i = 0; i < 3; i++) {
        v1[i] = v[1][i] - v[0][i];
        v2[i] = v[2][i] - v[0][i];
        v3[i] = v[3][i] - v[0][i];
    }
    vect_prod(v1, v2, v12);
    return fabs(scalar_prod(v12, v3));
}


double charm_quad_calc_volume(p4est_t * p4est, p4est_quadrant_t* q, p4est_topidx_t treeid)
{
    const int ttv[6][4] =
            {{ 0, 2, 4, 1 },
             { 5, 2, 4, 1 },
             { 5, 2, 4, 6 },
             { 3, 2, 7, 1 },
             { 5, 2, 7, 1 },
             { 5, 2, 7, 6 }};
    int8_t i, j, k;
    p4est_qcoord_t l  = P4EST_QUADRANT_LEN(q->level);
    double v[8][3], tet[4][3], vol;
    for (i = 0; i < P4EST_CHILDREN; i++) {
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l * (i % 2), q->y + l * ((i / 2) % 2), q->z + l * (i / 4), v[i]);
    }
    vol = 0.;
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < P4EST_DIM; k++) {
                tet[j][k] = v[ttv[i][j]][k];
            }
        }
        vol += charm_tet_calc_volume(tet);
    }
    return vol;
}


void charm_geom_quad_calc(p4est_t * p4est, p4est_quadrant_t* q, p4est_topidx_t treeid)
{
    int8_t i;
    charm_data_t *p = (charm_data_t*)q->p.user_data;

    for (i = 0; i < P4EST_FACES; i++) {
        charm_face_calc_center(p4est, q, treeid, i, p->par.g.fc[i]);
        charm_face_calc_normal(p4est, q, treeid, i, p->par.g.n[i]);
        charm_face_calc_gp(p4est, q, treeid, i, p->par.g.face_gp[i], p->par.g.face_gw[i]);
        p->par.g.area[i] = charm_face_calc_area(p4est, q, treeid, i);

    }
    p->par.g.volume = charm_quad_calc_volume(p4est, q, treeid);
    charm_quad_calc_center(p4est, q, treeid, p->par.g.c);
    charm_quad_calc_gp(p4est, q, treeid, p->par.g.quad_gp, p->par.g.quad_gw);

}


static void charm_geom_quad_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_geom_quad_calc(info->p4est, q, info->treeid);
}

void charm_geom_calc(p4est_t * p4est)
{
    p4est_iterate (p4est,
                   NULL,                      /* ghosts are not needed for this loop */
                   NULL,                      /*  NULL */
                   charm_geom_quad_fn,        /* update each cell */
                   NULL,                      /* there is no callback for the faces between quadrants */
                   NULL,                      /* there is no callback for the faces between quadrants */
                   NULL);                     /* there is no callback for the corners between quadrants */

}

