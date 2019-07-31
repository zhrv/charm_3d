//
// Created by zhrv on 31.07.19.
//

#include "charm_limiter.h"
#include "charm_globals.h"
#include "charm_base_func.h"
#include "charm_bnd_cond.h"


static void _charm_limiter_bj_init_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t *p = charm_get_quad_data(info->quad);
    size_t c_count = charm_get_comp_count(info->p4est);
    charm_cons_t cons;
    int j;

    charm_get_fields_avg(p, &cons);
    p->par.l.count = 1;
    p->par.l.ru[0] = cons.ru;
    p->par.l.rv[0] = cons.rv;
    p->par.l.rw[0] = cons.rw;
    p->par.l.re[0] = cons.re;
    for (j =0; j < c_count; j++) {
        p->par.l.rc[j][0] = cons.rc[j];
    }
}


static void _charm_limiter_bj_neigh_iter_bnd(p4est_iter_face_info_t * info, void *user_data)
{
    int i, j, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    size_t c_count = charm_get_comp_count(info->p4est);
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
        udata = charm_get_quad_data(side[0]->is.full.quad);
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

    charm_get_fields_avg(udata, &(cons[0]));
    charm_param_cons_to_prim(p4est, &(prim[0]), &(cons[0]));
    charm_bnd_cond(p4est, side[0]->treeid, face, &(prim[0]), &(prim[1]), n);
    charm_param_prim_to_cons(p4est, &(cons[1]), &(prim[1]));
    i = udata->par.l.count;
    udata->par.l.ru[i] = cons[1].ru;
    udata->par.l.rv[i] = cons[1].rv;
    udata->par.l.rw[i] = cons[1].rw;
    udata->par.l.re[i] = cons[1].re;
    for (j =0; j < c_count; j++) {
        udata->par.l.rc[j][i] = cons[1].rc[j];
    }
    udata->par.l.count++;
}


static void _charm_limiter_bj_neigh_iter_inner(p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, k,  h_side, f_side, cj;
    p4est_t                *p4est = info->p4est;
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2], cons_j;
    double                  c[3];
    int8_t                  face[2];
    double                  vol, svol;
    size_t                  c_count = charm_get_comp_count(info->p4est);

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    face[0] = side[0]->face;
    face[1] = side[1]->face;

    if (side[0]->is_hanging || side[1]->is_hanging) { // @todo
        if (side[0]->is_hanging) {
            h_side = 0;
            f_side = 1;
        }
        else {
            h_side = 1;
            f_side = 0;
        }

        for (i = 0; i < 2; i++) {
            if (side[i]->is_hanging) {
                svol = 0.;
                cons[i].ru = 0.;
                cons[i].rv = 0.;
                cons[i].rw = 0.;
                cons[i].re = 0.;
                for (cj = 0; cj < c_count; cj++) {
                    cons[i].rc[cj] = 0.;
                }
                for (j = 0; j < CHARM_HALF; j++) {
                    if (side[i]->is.hanging.is_ghost[j]) {
                        udata[i] = &(ghost_data[side[i]->is.hanging.quadid[j]]);
                    }
                    else {
                        udata[i] = (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
                    }
                    vol = charm_quad_get_volume(udata[i]);
                    svol += vol;
                    charm_quad_get_center(udata[i], c);
                    charm_get_fields(udata[i], c, &cons_j);
                    cons[i].ru += cons_j.ru*vol;
                    cons[i].rv += cons_j.rv*vol;
                    cons[i].rw += cons_j.rw*vol;
                    cons[i].re += cons_j.re*vol;
                    for (cj = 0; cj < c_count; cj++) {
                        cons[i].rc[cj] += cons_j.rc[cj]*vol;
                    }
                }
                cons[i].ru /= svol;
                cons[i].rv /= svol;
                cons[i].rw /= svol;
                cons[i].re /= svol;
                for (cj = 0; cj < c_count; cj++) {
                    cons[i].rc[cj] /= svol;
                }
            }
            else {
                if (side[i]->is.full.is_ghost) {
                    udata[i] = &ghost_data[side[i]->is.full.quadid];
                }
                else {
                    udata[i] = (charm_data_t *) side[i]->is.full.quad->p.user_data;
                }
                charm_quad_get_center(udata[i], c);
                charm_get_fields(udata[i], c, &(cons[i]));
            }
        }


        if (!side[f_side]->is.full.is_ghost) {
            udata[f_side] = (charm_data_t *) side[f_side]->is.full.quad->p.user_data;
            j = udata[f_side]->par.l.count;
            udata[f_side]->par.l.ru[j] = cons[h_side].ru;
            udata[f_side]->par.l.rv[j] = cons[h_side].rv;
            udata[f_side]->par.l.rw[j] = cons[h_side].rw;
            udata[f_side]->par.l.re[j] = cons[h_side].re;
            udata[f_side]->par.l.count++;
        }
        for (j = 0; j < CHARM_HALF; j++) {
            if (!side[h_side]->is.hanging.is_ghost[j]) {
                udata[h_side] = (charm_data_t *) side[h_side]->is.hanging.quad[j]->p.user_data;
                k = udata[h_side]->par.l.count;
                udata[h_side]->par.l.ru[k] = cons[f_side].ru;
                udata[h_side]->par.l.rv[k] = cons[f_side].rv;
                udata[h_side]->par.l.rw[k] = cons[f_side].rw;
                udata[h_side]->par.l.re[k] = cons[f_side].re;
                udata[h_side]->par.l.count++;
            }
        }
    }
    else {

        for (i = 0; i < 2; i++) {
            if (side[i]->is.full.is_ghost) {
                udata[i] = &(ghost_data[side[i]->is.full.quadid]);
            }
            else {
                udata[i] = charm_get_quad_data(side[i]->is.full.quad);
            }
        }

        for (i = 0; i < 2; i++) {
            charm_quad_get_center(udata[i], c);
            charm_get_fields_avg(udata[i], &(cons[i]));
            k = (i+1) % 2;
            if (!side[k]->is.full.is_ghost) {
                j = udata[k]->par.l.count;
                udata[k]->par.l.ru[j] = cons[i].ru;
                udata[k]->par.l.rv[j] = cons[i].rv;
                udata[k]->par.l.rw[j] = cons[i].rw;
                udata[k]->par.l.re[j] = cons[i].re;
                for (cj =0; cj < c_count; cj++) {
                    udata[k]->par.l.rc[cj][j] = cons[i].rc[cj];
                }
                udata[k]->par.l.count++;
            }
        }
    }
}


static void _charm_limiter_bj_neigh_iter_fn(p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_limiter_bj_neigh_iter_bnd(info, user_data);
    }
    else {
        _charm_limiter_bj_neigh_iter_inner(info, user_data);
    }


}


static void _charm_limiter_bj_calc_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t *p = charm_get_quad_data(info->quad);
    double **u, *f;
    double *u_min, *u_max, *psi, psi_tmp;
    int i,j;
    double v[8][CHARM_DIM];
    charm_cons_t cons;
    size_t                  c_count = charm_get_comp_count(info->p4est);
    size_t                  f_count = c_count+4;
    CHARM_ASSERT(p->par.l.count == 7);

    u_min = CHARM_ALLOC(double, f_count);
    u_max = CHARM_ALLOC(double, f_count);
    psi = CHARM_ALLOC(double, f_count);
    u = CHARM_ALLOC(double*, f_count);
    f = CHARM_ALLOC(double, f_count);

    u[0] = p->par.l.ru;
    u[1] = p->par.l.rv;
    u[2] = p->par.l.rw;
    u[3] = p->par.l.re;
    for (j = 4; j < f_count; j++) {
        u[j] = p->par.l.rc[j-4];
    }

    for (j = 0; j < f_count; j++) {
        u_min[j] = u_max[j] = u[j][0];
        for (i = 1; i < p->par.l.count; i++) {
            if (u_min[j] > u[j][i]) u_min[j] = u[j][i];
            if (u_max[j] < u[j][i]) u_max[j] = u[j][i];
        }
    }

    charm_quad_get_vertices(info->p4est, info->quad, info->treeid, v);
    for (j = 0; j < f_count; j++) {
        psi[j] = 1.;
    }
    for (i = 0; i < 8; i++) {
        charm_get_fields(p, v[i], &cons);
        f[0] = cons.ru;
        f[1] = cons.rv;
        f[2] = cons.rw;
        f[3] = cons.re;
        for (j = 4; j < f_count; j++) {
            f[j] = cons.rc[j-4];
        }
        for (j = 0; j < f_count; j++) {
            psi_tmp = 1.;
            if (fabs(f[j]-u[j][0]) > CHARM_EPS) {
                if (f[j] - u[j][0] > 0) {
                    psi_tmp = _MIN_(1, (u_max[j] - u[j][0]) / (f[j] - u[j][0]));
                } else if (f[j] - u[j][0] < 0) {
                    psi_tmp = _MIN_(1, (u_min[j] - u[j][0]) / (f[j] - u[j][0]));
                }
            }
            if (psi_tmp < psi[j]) psi[j] = psi_tmp;
        }
    }
    for (j = 1; j < CHARM_BASE_FN_COUNT; j++) {
        p->par.c.ru[j] *= psi[0];
        p->par.c.rv[j] *= psi[1];
        p->par.c.rw[j] *= psi[2];
        p->par.c.re[j] *= psi[3];
        for (i = 0; i < c_count; i++) {
            p->par.c.rc[i][j] *= psi[4+i];
        }
    }

    CHARM_FREE(u_min);
    CHARM_FREE(u_max);
    CHARM_FREE(psi);
    CHARM_FREE(u);
    CHARM_FREE(f);
}


void charm_limiter_bj(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data)
{
    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   _charm_limiter_bj_init_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est, ghost, (void *) ghost_data, NULL,
                   _charm_limiter_bj_neigh_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   _charm_limiter_bj_calc_iter_fn,
                   NULL, NULL, NULL);

}


