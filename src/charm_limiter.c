//
// Created by appmath on 29.09.18.
//

#include "charm_limiter.h"
#include "charm_globals.h"
#include "charm_base_func.h"
#include "charm_bnd_cond.h"
#include "charm_geom.h"

static void _charm_limiter_init_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t *p = charm_get_quad_data(info->quad);

    p->par.l.count = 1;
    p->par.l.ro[0] = p->par.c.ro[0];
    p->par.l.ru[0] = p->par.c.ru[0];
    p->par.l.rv[0] = p->par.c.rv[0];
    p->par.l.rw[0] = p->par.c.rw[0];
    p->par.l.re[0] = p->par.c.re[0];
}

static void _charm_limiter_neigh_iter_bnd(p4est_iter_face_info_t * info, void *user_data)
{
    int i, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    double n[3];
    double qr, qu, qv, qw, qe;
    double bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);


    int8_t face;
    double c[2][3], l[3];
    double r_[2], p_[2], u_[2], v_[2], w_[2], e_[2];
    charm_cons_t cons[2];
    charm_prim_t prim[2];
    double *x, gw, gj;
    double intg[2][5];


    CHARM_ASSERT(info->tree_boundary);

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    CHARM_ASSERT(!side[0]->is_hanging);

    if (side[0]->is.full.is_ghost) {
        CHARM_ASSERT(0);
        udata = &(ghost_data[side[0]->is.full.quadid]);
    } else {
        udata = charm_get_quad_data(side[0]->is.full.quad);//(charm_data_t *) side[0]->is.full.quad->p.user_data;
    }
    face = side[0]->face;
    charm_face_get_normal(udata, face, n);
    charm_quad_get_center(udata, c[0]);
    charm_face_get_center(udata, face, c[1]);

    for (i = 0; i < 3; i++) {
        l[i] = c[1][i] - c[0][i];
    }

    if (scalar_prod(n, l) < 0) {
        for (i = 0; i < 3; i++) {
            n[i] *= -1.0;
        }
    }

    charm_get_fields(udata, c[0], &(cons[0]));
    charm_param_cons_to_prim(p4est, &(prim[0]), &(cons[0]));
    charm_bnd_cond(p4est, side[0]->treeid, face, &(prim[0]), &(prim[1]), n);
    charm_param_prim_to_cons(p4est, &(cons[1]), &(prim[1]));
    i = udata->par.l.count;
    udata->par.l.ro[i] = cons[1].ro;
    udata->par.l.ru[i] = cons[1].ru;
    udata->par.l.rv[i] = cons[1].rv;
    udata->par.l.rw[i] = cons[1].rw;
    udata->par.l.re[i] = cons[1].re;
    udata->par.l.count++;
}


static void _charm_limiter_neigh_iter_inner(p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, k,  h_side;
    p4est_t                *p4est = info->p4est;
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    double                  c[3];
    int8_t                  face[2];

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    face[0] = side[0]->face;
    face[1] = side[1]->face;

    h_side = -1;
    if (side[0]->is_hanging || side[1]->is_hanging) { // @todo
                CHARM_ASSERT(0);
//        for (j = 0; j < CHARM_HALF; j++) {
//            for (i = 0; i < 2; i++) {
//                if (side[i]->is_hanging) {
//                    if (side[i]->is.hanging.is_ghost[j]) {
//                        udata[i] = &(ghost_data[side[i]->is.hanging.quadid[j]]);
//                    }
//                    else {
//                        udata[i] = (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
//                    }
//                    h_side = i;
//                }
//                else {
//                    if (side[i]->is.full.is_ghost) {
//                        udata[i] = &ghost_data[side[i]->is.full.quadid];
//                    }
//                    else {
//                        udata[i] = (charm_data_t *) side[i]->is.full.quad->p.user_data;
//                    }
//                }
//            }
//
//            CHARM_ASSERT(h_side != -1);
//
//            facearea = charm_face_get_area(udata[h_side], side[h_side]->face);
//            charm_face_get_normal(udata[0], face[0], n);
//            charm_quad_get_center(udata[0], c[0]);
//            charm_face_get_center(udata[0], face[0], c[1]);
//
//            for (i = 0; i < 3; i++) {
//                l[i] = c[1][i]-c[0][i];
//            }
//
//            if (scalar_prod(n, l) < 0) {
//                for (i = 0; i < 3; i++) {
//                    n[i] *= -1.0;
//                }
//            }
//
//            for (i = 0; i < 2; i++) {
////                charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[i]->treeid);
////                charm_param_cons_to_prim(attr->reg->mat, &(cons[i]), &());
////                r_[i] = udata[i]->par.p.r;
////                u_[i] = udata[i]->par.p.u;
////                v_[i] = udata[i]->par.p.v;
////                w_[i] = udata[i]->par.p.w;
////                p_[i] = udata[i]->par.p.p;
//            }
//
//            /* flux from side 0 to side 1 */
//            charm_calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, n);
//
//            for (i = 0; i < 2; i++) {
//                if (side[i]->is_hanging) {
//                    if (side[i]->is.hanging.is_ghost[j]) {
//                        continue;
//                    }
//
//                }
//                else {
//                    if (side[i]->is.full.is_ghost) {
//                        continue;
//                    }
//                }
////                udata[i]->drodt += qr * facearea * (i ? 1. : -1.);
////                udata[i]->drudt += qu * facearea * (i ? 1. : -1.);
////                udata[i]->drvdt += qv * facearea * (i ? 1. : -1.);
////                udata[i]->drwdt += qw * facearea * (i ? 1. : -1.);
////                udata[i]->dredt += qe * facearea * (i ? 1. : -1.);
//
//            }
//
//        }
    }
    else {

        for (i = 0; i < 2; i++) {
            if (side[i]->is.full.is_ghost) {
                udata[i] = &(ghost_data[side[i]->is.full.quadid]);
            }
            else {
                udata[i] = charm_get_quad_data(side[i]->is.full.quad);//(charm_data_t *) side[i]->is.full.quad->p.user_data;
            }
        }

        for (i = 0; i < 2; i++) {
            charm_quad_get_center(udata[i], c);
            charm_get_fields(udata[i], c, &(cons[i]));
            k = (i+1) % 2;
            if (!side[k]->is.full.is_ghost) {
                j = udata[k]->par.l.count;
                udata[k]->par.l.ro[j] = cons[i].ro;
                udata[k]->par.l.ru[j] = cons[i].ru;
                udata[k]->par.l.rv[j] = cons[i].rv;
                udata[k]->par.l.rw[j] = cons[i].rw;
                udata[k]->par.l.re[j] = cons[i].re;
                udata[k]->par.l.count++;
            }

        }
    }
}


static void _charm_limiter_neigh_iter_fn(p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_limiter_neigh_iter_bnd(info, user_data);
    }
    else {
        _charm_limiter_neigh_iter_inner(info, user_data);
    }


}


static void _charm_limiter_calc_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t *p = charm_get_quad_data(info->quad);
    double *u[5], f[5][8];
    double u_min[5], u_max[5], psi[5], psi_tmp;
    int i,j;
    double v[8][CHARM_DIM];
    charm_cons_t cons;

    CHARM_ASSERT(p->par.l.count == 7);

    u[0] = p->par.l.ro;
    u[1] = p->par.l.ru;
    u[2] = p->par.l.rv;
    u[3] = p->par.l.rw;
    u[4] = p->par.l.re;

    for (j = 0; j < 5; j++) {
        u_min[j] = u_max[j] = u[j][0];
    }


    for (j = 1; j < 5; j++) {
        for (i = 0; i < p->par.l.count; i++) {
            if (u_min[j] > u[j][i]) u_min[j] = u[j][i];
            if (u_max[j] < u[j][i]) u_max[j] = u[j][i];
        }
    }

    charm_quad_get_vertices(info->p4est, info->quad, info->treeid, v);
    for (j = 0; j < 5; j++) {
        psi[j] = 1.;
    }
    for (i = 0; i < 8; i++) {
        charm_get_fields(p, v[i], &cons);
        f[0][i] = cons.ro;
        f[1][i] = cons.ru;
        f[2][i] = cons.rv;
        f[3][i] = cons.rw;
        f[4][i] = cons.re;
        for (j = 0; j < 5; j++) {

            if (u[j][i]-u[j][0] > 0) {
                psi_tmp = _MIN_(1, (u_max[j]-u[j][0])/(u[j][i]-u[j][0]));
            }
            else if (u[j][i]-u[j][0] < 0) {
                psi_tmp = _MIN_(1, (u_min[j]-u[j][0])/(u[j][i]-u[j][0]));
            }
            else {
                psi_tmp = 1.;
            }
            if (psi_tmp < psi[j]) psi[j] = psi_tmp;
        }
    }
    for (j = 1; j < CHARM_BASE_FN_COUNT; j++) {
        p->par.c.ro[j] *= psi[0];
        p->par.c.ru[j] *= psi[1];
        p->par.c.rv[j] *= psi[2];
        p->par.c.rw[j] *= psi[3];
        p->par.c.re[j] *= psi[4];
    }
}


void charm_limiter(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data)
{
    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   _charm_limiter_init_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est, ghost, (void *) ghost_data, NULL,
                   _charm_limiter_neigh_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   _charm_limiter_calc_iter_fn,
                   NULL, NULL, NULL);

}