//
// Created by zhrv on 27.10.17.
//

#include "charm_base_func.h"


/*
 *  Cells
 */

void charm_quad_get_vertices(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, charm_vec_t v[8])
{
    p4est_qcoord_t l = CHARM_QUADRANT_LEN(q->level);
    p4est_qcoord_t qx, qy, qz;
    charm_data_t *p = (charm_data_t*)q->p.user_data;
    int i, iz, iy, ix;

    i = 0;
    for (iz = 0; iz < 2; iz++) {
        qz = q->z + iz*l;
        for (iy = 0; iy < 2; iy++) {
            qy = q->y + iy*l;
            for (ix = 0; ix < 2; ix++) {
                qx = q->x + ix*l;
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, qx, qy, qz, v[i]);
                i++;
            }
        }
    }
}

static void charm_quad_calc_center(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, charm_real_t* c)
{
    p4est_qcoord_t l2 = CHARM_QUADRANT_LEN(q->level) / 2;
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l2, q->y + l2, q->z + l2, c);
}


static void charm_quad_calc_gp_at_point(charm_real_t vertices[8][3], charm_vec_t ref_p, charm_vec_t gp, charm_real_t* gj)
{
    int                 ix, iy, iz, vindex;
    charm_real_t              wx[2], wy[2], wz[2];
    charm_real_t              xfactor, yfactor;
    charm_real_t              jacob[3][3], a;

    wx[0] = 1.-ref_p[0];
    wx[1] = 1.+ref_p[0];
    wy[0] = 1.-ref_p[1];
    wy[1] = 1.+ref_p[1];
    wz[0] = 1.-ref_p[2];
    wz[1] = 1.+ref_p[2];
    gp[0] = gp[1] = gp[2] = 0.;
    vindex = 0;
    for (iz = 0; iz < 2; iz++) {
        for (iy = 0; iy < 2; iy++) {
            yfactor = wz[iz]*wy[iy];
            for (ix = 0; ix < 2; ix++) {
                xfactor = yfactor*wx[ix]/8.;
                gp[0] += xfactor * vertices[vindex][0];
                gp[1] += xfactor * vertices[vindex][1];
                gp[2] += xfactor * vertices[vindex][2];
                vindex++;
            }
        }
    }

    // d/dksi
    wx[0] = -1.;
    wx[1] =  1.;
    jacob[0][0] = jacob[1][0] = jacob[2][0] = 0.;
    vindex = 0;
    for (iz = 0; iz < 2; iz++) {
        for (iy = 0; iy < 2; iy++) {
            yfactor = wz[iz]*wy[iy];
            for (ix = 0; ix < 2; ix++) {
                xfactor = yfactor*wx[ix]/8.;
                jacob[0][0] += xfactor * vertices[vindex][0];
                jacob[1][0] += xfactor * vertices[vindex][1];
                jacob[2][0] += xfactor * vertices[vindex][2];
                vindex++;
            }
        }
    }

    // d/deta
    wx[0] = 1.-ref_p[0];
    wx[1] = 1.+ref_p[0];
    wy[0] = -1.;
    wy[1] =  1.;
    jacob[0][1] = jacob[1][1] = jacob[2][1] = 0.;
    vindex = 0;
    for (iz = 0; iz < 2; iz++) {
        for (iy = 0; iy < 2; iy++) {
            yfactor = wz[iz]*wy[iy];
            for (ix = 0; ix < 2; ix++) {
                xfactor = yfactor*wx[ix]/8.;
                jacob[0][1] += xfactor * vertices[vindex][0];
                jacob[1][1] += xfactor * vertices[vindex][1];
                jacob[2][1] += xfactor * vertices[vindex][2];
                vindex++;
            }
        }
    }

    // d/dmu
    wy[0] = 1.-ref_p[1];
    wy[1] = 1.+ref_p[1];
    wz[0] = -1.;
    wz[1] =  1.;
    jacob[0][2] = jacob[1][2] = jacob[2][2] = 0.;
    vindex = 0;
    for (iz = 0; iz < 2; iz++) {
        for (iy = 0; iy < 2; iy++) {
            yfactor = wz[iz]*wy[iy];
            for (ix = 0; ix < 2; ix++) {
                xfactor = yfactor*wx[ix]/8.;
                jacob[0][2] += xfactor * vertices[vindex][0];
                jacob[1][2] += xfactor * vertices[vindex][1];
                jacob[2][2] += xfactor * vertices[vindex][2];
                vindex++;
            }
        }
    }

    *gj = charm_matr3_det(jacob);
    CHARM_ASSERT(*gj != 0.);
}

static void charm_quad_calc_gp(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, charm_vec_t gp[CHARM_QUAD_GP_COUNT], charm_real_t gw[CHARM_QUAD_GP_COUNT], charm_real_t gj[CHARM_QUAD_GP_COUNT])
{
    p4est_qcoord_t l = CHARM_QUADRANT_LEN(q->level);
    p4est_qcoord_t qx, qy, qz;
    charm_data_t *p = (charm_data_t*)q->p.user_data;
    int i, j, iz, iy, ix;
    charm_real_t v[8][3];
    charm_real_t vmin[3], vmax[3];
    charm_real_t sqrt3 = 1./sqrt(3.);
    charm_real_t t[8][3] = {
            {-sqrt3, -sqrt3, -sqrt3},
            { sqrt3, -sqrt3, -sqrt3},
            {-sqrt3,  sqrt3, -sqrt3},
            { sqrt3,  sqrt3, -sqrt3},
            {-sqrt3, -sqrt3,  sqrt3},
            { sqrt3, -sqrt3,  sqrt3},
            {-sqrt3,  sqrt3,  sqrt3},
            { sqrt3,  sqrt3,  sqrt3}
    };

    charm_quad_get_vertices(p4est, q, treeid, v);

    for (i = 0; i < 8; i++) {
        charm_quad_calc_gp_at_point(v, t[i], gp[i], &gj[i]);
        gw[i] = 1.;
    }

    for (j = 0; j < 3; j++) {
        vmin[j] = v[0][j];
        vmax[j] = v[0][j];
    }
    for (i = 1; i < 8; i++) {
        for (j = 0; j < 3; j++) {
            if (vmin[j] > v[i][j]) vmin[j] = v[i][j];
            if (vmax[j] < v[i][j]) vmax[j] = v[i][j];
        }
    }
    for (j = 0; j < 3; j++) {
        p->par.g.dh[j] = vmax[j] - vmin[j];
    }


}



/*
 *  Faces
 */


static charm_real_t charm_face_calc_normal(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, charm_real_t* n)
{
    const int ftv[6][4] =
            {{ 0, 2, 4, 6 },
             { 1, 3, 5, 7 },
             { 0, 1, 4, 5 },
             { 2, 3, 6, 7 },
             { 0, 1, 2, 3 },
             { 4, 5, 6, 7 }};

    int i,k;
    charm_real_t x[4][3], v[2][3], nl;
    p4est_qcoord_t l = CHARM_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    for (i = 0; i < 4; i++) {
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l * (ftv[face][i] % 2), q->y + l * ((ftv[face][i] / 2) % 2), q->z + l * (ftv[face][i] / 4), x[i]);
    }

    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[0][i];
        v[1][i] = x[2][i]-x[0][i];
    }

    vector_prod(v[0], v[1], n);
    nl = vector_length(n);
    for (i = 0; i < 3; i++) {
        n[i] /= nl;
    }

    return nl;
}


static charm_real_t charm_face_calc_area(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face)
{
    const int ftv[6][4] =
            {{ 0, 2, 4, 6 },
             { 1, 3, 5, 7 },
             { 0, 1, 4, 5 },
             { 2, 3, 6, 7 },
             { 0, 1, 2, 3 },
             { 4, 5, 6, 7 }};

    int i,k;
    charm_real_t x[4][3], v[2][3], n[3], s1, s2;
    p4est_qcoord_t l = CHARM_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    for (i = 0; i < 4; i++) {
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l * (ftv[face][i] % 2), q->y + l * ((ftv[face][i] / 2) % 2), q->z + l * (ftv[face][i] / 4), x[i]);
    }

    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[0][i];
        v[1][i] = x[2][i]-x[0][i];
    }

    vector_prod(v[0], v[1], n);
    s1 = 0.5*vector_length(n);
    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[3][i];
        v[1][i] = x[2][i]-x[3][i];
    }

    vector_prod(v[0], v[1], n);
    s2 = 0.5*vector_length(n);



    return s1+s2;


}


static void charm_face_calc_center(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, charm_real_t* c)
{
    p4est_qcoord_t l  = CHARM_QUADRANT_LEN(q->level);
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

static void charm_geom_face_get_v(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, charm_real_t v[4][3])
{
    const int ftv[6][4] =
            {{ 0, 2, 4, 6 },
             { 1, 3, 5, 7 },
             { 0, 1, 4, 5 },
             { 2, 3, 6, 7 },
             { 0, 1, 2, 3 },
             { 4, 5, 6, 7 }};
    p4est_qcoord_t l  = CHARM_QUADRANT_LEN(q->level);
    int i;

    for (i = 0; i < 4; i++) {
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l * (ftv[face][i] % 2), q->y + l * ((ftv[face][i] / 2) % 2), q->z + l * (ftv[face][i] / 4), v[i]);
    }
}

static void charm_face_calc_gp_at_point(charm_real_t vertices[8][3], int8_t face, int i_gp, charm_vec_t ref_p, charm_vec_t gp, charm_real_t* gj)
{
    int                 ix, iy, iz, vindex, i, j;
    charm_real_t              wx[2], wy[2], wz[2];
    charm_real_t              xfactor, yfactor;
    charm_real_t              jacob[2][2], a;
    charm_real_t              v[4][3];
    charm_real_t              sqrt3 = 1./sqrt(3.);
    int                 ind[2];
    const int ftv[6][4] =
            {{ 0, 2, 4, 6 },
             { 1, 3, 5, 7 },
             { 0, 1, 4, 5 },
             { 2, 3, 6, 7 },
             { 0, 1, 2, 3 },
             { 4, 5, 6, 7 }};
    charm_real_t ref_gp[4][3] =
            {{-sqrt3, -sqrt3},
             { sqrt3, -sqrt3},
             {-sqrt3,  sqrt3},
             { sqrt3,  sqrt3}};


    wx[0] = 1.-ref_p[0];
    wx[1] = 1.+ref_p[0];
    wy[0] = 1.-ref_p[1];
    wy[1] = 1.+ref_p[1];
    wz[0] = 1.-ref_p[2];
    wz[1] = 1.+ref_p[2];
    gp[0] = gp[1] = gp[2] = 0.;
    vindex = 0;
    for (iz = 0; iz < 2; iz++) {
        for (iy = 0; iy < 2; iy++) {
            yfactor = wz[iz]*wy[iy];
            for (ix = 0; ix < 2; ix++) {
                xfactor = yfactor*wx[ix]/8.;
                gp[0] += xfactor * vertices[vindex][0];
                gp[1] += xfactor * vertices[vindex][1];
                gp[2] += xfactor * vertices[vindex][2];
                vindex++;
            }
        }
    }

    for (i = 0; i < 4; i++) {
        for (j = 0; j < CHARM_DIM; j++) {
            v[i][j] = vertices[ftv[face][i]][j];
        }
    }

    if (face == 0 || face == 1) {
        ind[0] = 1;
        ind[1] = 2;
    }
    else if (face == 2 || face == 3) {
        ind[0] = 0;
        ind[1] = 2;
    }
    else if (face == 4 || face == 5) {
        ind[0] = 0;
        ind[1] = 1;
    }

    // d/dksi
    wx[0] = -1.;
    wx[1] =  1.;
    wy[0] = 1.-ref_gp[i_gp][1];
    wy[1] = 1.+ref_gp[i_gp][1];
    jacob[0][0] = jacob[1][0] = 0.;
    vindex = 0;
    for (iy = 0; iy < 2; iy++) {
        yfactor = wy[iy];
        for (ix = 0; ix < 2; ix++) {
            xfactor = yfactor*wx[ix]/4.;
            jacob[0][0] += xfactor * v[vindex][ind[0]];
            jacob[1][0] += xfactor * v[vindex][ind[1]];
            vindex++;
        }
    }

    // d/deta
    wx[0] = 1.-ref_gp[i_gp][0];
    wx[1] = 1.+ref_gp[i_gp][0];
    wy[0] = -1.;
    wy[1] =  1.;
    jacob[0][1] = jacob[1][1] = 0.;
    vindex = 0;
    for (iy = 0; iy < 2; iy++) {
        yfactor = wy[iy];
        for (ix = 0; ix < 2; ix++) {
            xfactor = yfactor*wx[ix]/4.;
            jacob[0][1] += xfactor * v[vindex][ind[0]];
            jacob[1][1] += xfactor * v[vindex][ind[1]];
            vindex++;
        }
    }


    *gj = jacob[0][0]*jacob[1][1]-jacob[0][1]*jacob[1][0];
    CHARM_ASSERT(*gj != 0.);
}

static void charm_face_calc_gp(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, charm_vec_t gp[CHARM_FACE_GP_COUNT], charm_real_t gw[CHARM_FACE_GP_COUNT], charm_real_t gj[CHARM_FACE_GP_COUNT])
{
    int i;
    charm_real_t sqrt3 = 1./sqrt(3.);
    charm_real_t v[8][3];
    charm_real_t ref_gp[6][4][3] = {
            {
                    {-1., -sqrt3, -sqrt3},
                    {-1.,  sqrt3, -sqrt3},
                    {-1., -sqrt3,  sqrt3},
                    {-1.,  sqrt3,  sqrt3}
            },
            {
                    { 1., -sqrt3, -sqrt3},
                    { 1.,  sqrt3, -sqrt3},
                    { 1., -sqrt3,  sqrt3},
                    { 1.,  sqrt3,  sqrt3}
            },
            {
                    {-sqrt3, -1., -sqrt3},
                    { sqrt3, -1., -sqrt3},
                    {-sqrt3, -1.,  sqrt3},
                    { sqrt3, -1.,  sqrt3}
            },
            {
                    {-sqrt3,  1., -sqrt3},
                    { sqrt3,  1., -sqrt3},
                    {-sqrt3,  1.,  sqrt3},
                    { sqrt3,  1.,  sqrt3}
            },
            {
                    {-sqrt3, -sqrt3, -1.},
                    { sqrt3, -sqrt3, -1.},
                    {-sqrt3,  sqrt3, -1.},
                    { sqrt3,  sqrt3, -1.}
            },
            {
                    {-sqrt3, -sqrt3,  1.},
                    { sqrt3, -sqrt3,  1.},
                    {-sqrt3,  sqrt3,  1.},
                    { sqrt3,  sqrt3,  1.}
            }
    };

    charm_quad_get_vertices(p4est, q, treeid, v);

    for (i = 0; i < CHARM_FACE_GP_COUNT; i++) {
        charm_face_calc_gp_at_point(v, face, i, ref_gp[face][i], gp[i], &gj[i]);
        gw[i] = 1.;
    }


}

static void charm_face_calc_gp_by_tri(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, charm_vec_t gp[CHARM_FACE_GP_COUNT], charm_real_t gw[CHARM_FACE_GP_COUNT], charm_real_t gj[CHARM_FACE_GP_COUNT])
{
    charm_real_t	pt[3][3];
    charm_real_t	pGP[3][3];
    charm_real_t  v[4][3];
    charm_real_t	a = 2.0/3.0;
    charm_real_t	b = 1.0/6.0;
    int i,j,k ;


    charm_geom_face_get_v(p4est, q, treeid, face, v);

    for (i = 0; i < CHARM_FACE_GP_COUNT; i++)
    {
        gw[i] = 1./6.;
    }

    for ( i = 0; i < 2; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++) {
                pt[j][k] = v[i + j][k];
            }
        }

        charm_real_t a1 = pt[0][0] - pt[2][0];	charm_real_t a2 = pt[1][0] - pt[2][0];	charm_real_t a3 = pt[2][0];
        charm_real_t b1 = pt[0][1] - pt[2][1];	charm_real_t b2 = pt[1][1] - pt[2][1];	charm_real_t b3 = pt[2][1];
        charm_real_t c1 = pt[0][2] - pt[2][2];	charm_real_t c2 = pt[1][2] - pt[2][2];	charm_real_t c3 = pt[2][2];

        charm_real_t E = a1*a1+b1*b1+c1*c1;
        charm_real_t F = a2*a2+b2*b2+c2*c2;
        charm_real_t G = a1*a2+b1*b2+c1*c2;

        gj[3*i+0] = gj[3*i+1] = gj[3*i+2] = sqrt(E*F-G*G);

        gp[3*i+0][0] = a * a1 + b * a2 + a3;
        gp[3*i+0][1] = a * b1 + b * b2 + b3;
        gp[3*i+0][2] = a * c1 + b * c2 + c3;

        gp[3*i+1][0] = b * a1 + a * a2 + a3;
        gp[3*i+1][1] = b * b1 + a * b2 + b3;
        gp[3*i+1][2] = b * c1 + a * c2 + c3;

        gp[3*i+2][0] = b * a1 + b * a2 + a3;
        gp[3*i+2][1] = b * b1 + b * b2 + b3;
        gp[3*i+2][2] = b * c1 + b * c2 + c3;
    }



}

void charm_geom_quad_calc(p4est_t * p4est, p4est_quadrant_t* q, p4est_topidx_t treeid)
{
    int8_t i, j, ig;
    charm_data_t *p = (charm_data_t*)q->p.user_data;

    for (i = 0; i < CHARM_FACES; i++) {
        charm_face_calc_center(p4est, q, treeid, i, p->par.g.fc[i]);
        charm_face_calc_normal(p4est, q, treeid, i, p->par.g.n[i]);
        charm_face_calc_gp_by_tri(p4est, q, treeid, i, p->par.g.face_gp[i], p->par.g.face_gw[i], p->par.g.face_gj[i]);
        p->par.g.area[i] = 0.;
        for (j = 0; j < CHARM_FACE_GP_COUNT; j++) {
            p->par.g.area[i] += p->par.g.face_gw[i][j]*p->par.g.face_gj[i][j];
        }
    }
    charm_quad_calc_center(p4est, q, treeid, p->par.g.c);
    charm_quad_calc_gp(p4est, q, treeid, p->par.g.quad_gp, p->par.g.quad_gw, p->par.g.quad_gj);
    p->par.g.volume = 0.;
    for (i = 0; i < CHARM_QUAD_GP_COUNT; i++) {
        p->par.g.volume += p->par.g.quad_gw[i]*p->par.g.quad_gj[i];
    }

    // mass matr
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        for (j = 0; j < CHARM_BASE_FN_COUNT; j++) {
            p->par.g.a[i][j] = 0.;
            for (ig = 0; ig < CHARM_QUAD_GP_COUNT; ig++) {
                p->par.g.a[i][j] += p->par.g.quad_gw[ig]*p->par.g.quad_gj[ig]
                                    * charm_base_func(p->par.g.quad_gp[ig], i, p)
                                    * charm_base_func(p->par.g.quad_gp[ig], j, p);
            }
        }
    }

    charm_matr_inv(p->par.g.a, p->par.g.a_inv);
}
