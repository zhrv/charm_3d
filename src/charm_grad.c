//
// Created by zhrv on 27.10.17.
//

#include "charm_grad.h"
#include "charm_bnd_cond.h"


static void charm_grad_reset_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q    = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    int i;

    for (i = 0; i < P4EST_DIM; i++) {
        data->par.grad.r[i] = 0.0;
        data->par.grad.u[i] = 0.0;
        data->par.grad.v[i] = 0.0;
        data->par.grad.w[i] = 0.0;
        data->par.grad.p[i] = 0.0;
    }
}


static void charm_grad_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    int                 i, j, k;
    p4est_t            *p4est = info->p4est;
    charm_data_t       *ghost_data = (charm_data_t *) user_data;
    charm_data_t       *udata, edata;
    p4est_quadrant_t   *quad;
    double              n[3];
    double              ro_avg[2], ru_avg[2], rv_avg[2], rw_avg[2], re_avg[2];
    double              qr, qu, qv, qw, qp;
    double              facearea;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);

    int8_t face[2];
    double c[2][3], l[3];

    if (sides->elem_count != 2) {

        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        P4EST_ASSERT(!side[0]->is_hanging);

        if (side[0]->is.full.is_ghost) {
            udata = &ghost_data[side[0]->is.full.quadid];
        }
        else {
            udata = (charm_data_t *) side[0]->is.full.quad->p.user_data;
        }
        face[0] = side[0]->face;
        facearea = charm_face_get_normal(udata, face[0], n);
        charm_quad_get_center(udata, c[0]);
        charm_face_get_center(udata, face[0], c[1]);

        for (i = 0; i < 3; i++) {
            l[i] = c[1][i]-c[0][i];
        }

        if (scalar_prod(n, l) < 0) {
            for (i = 0; i < 3; i++) {
                n[i] *= -1.0;
            }
        }

        ro_avg[0] = udata->par.c.ro;
        ru_avg[0] = udata->par.c.ru;
        rv_avg[0] = udata->par.c.rv;
        rw_avg[0] = udata->par.c.rw;
        re_avg[0] = udata->par.c.re;

        charm_bnd_cond(p4est, side[0]->treeid, face[0], &udata->par, &edata.par);

        ro_avg[1] = edata.par.c.ro;
        ru_avg[1] = edata.par.c.ru;
        rv_avg[1] = edata.par.c.rv;
        rw_avg[1] = edata.par.c.rw;
        re_avg[1] = edata.par.c.re;

    }
    else {
        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        side[1] = p4est_iter_fside_array_index_int(sides, 1);
        face[0] = side[0]->face;
        face[1] = side[1]->face;

        i = side[0]->is_hanging ? 1 : 0;

        if (side[i]->is.full.is_ghost) {
            udata = &ghost_data[side[i]->is.full.quadid];
        }
        else {
            udata = (charm_data_t *) side[i]->is.full.quad->p.user_data;
        }
        face[i] = side[i]->face;
        facearea = charm_face_get_normal(udata, face[i], n);
        charm_quad_get_center(udata, c[0]);
        charm_face_get_center(udata, face[i], c[1]);

        for (i = 0; i < 3; i++) {
            l[i] = c[1][i]-c[0][i];
        }

        if (scalar_prod(n, l) < 0) {
            for (i = 0; i < 3; i++) {
                n[i] *= (side[0]->is_hanging ? 1.0 : -1.0);
            }
        }


        for (i = 0; i < 2; i++) {
            ro_avg[i] = 0;
            ru_avg[i] = 0;
            rv_avg[i] = 0;
            rw_avg[i] = 0;
            re_avg[i] = 0;
            if (side[i]->is_hanging) {
                for (j = 0; j < P4EST_HALF; j++) {
                    if (side[i]->is.hanging.is_ghost[j]) {
                        udata = &ghost_data[side[i]->is.hanging.quadid[j]];
                    }
                    else {
                        udata = (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
                    }
                    ro_avg[i] += udata->par.c.ro;
                    ru_avg[i] += udata->par.c.ru;
                    rv_avg[i] += udata->par.c.rv;
                    rw_avg[i] += udata->par.c.rw;
                    re_avg[i] += udata->par.c.re;
                }
                ro_avg[i] /= P4EST_HALF;
                ru_avg[i] /= P4EST_HALF;
                rv_avg[i] /= P4EST_HALF;
                rw_avg[i] /= P4EST_HALF;
                re_avg[i] /= P4EST_HALF;
            }
            else {
                if (side[i]->is.full.is_ghost) {
                    udata = &ghost_data[side[i]->is.full.quadid];
                }
                else {
                    udata = (charm_data_t *) side[i]->is.full.quad->p.user_data;
                }
                ro_avg[i] = udata->par.c.ro;
                ru_avg[i] = udata->par.c.ru;
                rv_avg[i] = udata->par.c.rv;
                rw_avg[i] = udata->par.c.rw;
                re_avg[i] = udata->par.c.re;
            }
        }

    }

    double r_[2], p_[2], u_[2], v_[2], w_[2];
    for (i = 0; i < 2; i++) {
        r_[i] = ro_avg[i];
        u_[i] = ru_avg[i]/r_[i];
        v_[i] = rv_avg[i]/r_[i];
        w_[i] = rw_avg[i]/r_[i];
        p_[i] = (re_avg[i]-0.5*r_[i]*(u_[i]*u_[i]+v_[i]*v_[i]+w_[i]*w_[i]))*(GAM-1.0);
    }

//    /* flux from side 0 to side 1 */
//    calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, n, bnd);
//
    qr = 0.5*(r_[0]+r_[1]);
    qu = 0.5*(u_[0]+u_[1]);
    qv = 0.5*(v_[0]+v_[1]);
    qw = 0.5*(w_[0]+w_[1]);
    qp = 0.5*(p_[0]+p_[1]);
    for (i = 0; i < sides->elem_count; i++) {
        for (k = 0; k < CHARM_DIM; k++) {
            if (side[i]->is_hanging) {
                for (j = 0; j < P4EST_HALF; j++) {
                    quad = side[i]->is.hanging.quad[j];
                    if (!side[i]->is.hanging.is_ghost[j]) {
                        udata = (charm_data_t *) quad->p.user_data;
                        facearea = charm_face_get_area(udata, side[i]->face);
                        udata->par.grad.r[k] += qr * n[k] * facearea * (i ? 1. : -1.);
                        udata->par.grad.u[k] += qu * n[k] * facearea * (i ? 1. : -1.);
                        udata->par.grad.v[k] += qv * n[k] * facearea * (i ? 1. : -1.);
                        udata->par.grad.w[k] += qw * n[k] * facearea * (i ? 1. : -1.);
                        udata->par.grad.p[k] += qp * n[k] * facearea * (i ? 1. : -1.);
                    }
                }
            }
            else {
                quad = side[i]->is.full.quad;
                if (!side[i]->is.full.is_ghost) {
                    udata = (charm_data_t *) quad->p.user_data;
                    facearea = charm_face_get_area(udata, side[i]->face);
                    udata->par.grad.r[k] += qr * n[k] * facearea * (i ? 1. : -1.);
                    udata->par.grad.u[k] += qu * n[k] * facearea * (i ? 1. : -1.);
                    udata->par.grad.v[k] += qv * n[k] * facearea * (i ? 1. : -1.);
                    udata->par.grad.w[k] += qw * n[k] * facearea * (i ? 1. : -1.);
                    udata->par.grad.p[k] += qp * n[k] * facearea * (i ? 1. : -1.);
                }
            }
        }
    }
}

static void charm_grad_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q    = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    double              vol  = charm_quad_get_volume(data);
    int i;

    for (i = 0; i < P4EST_DIM; i++) {
        data->par.grad.r[i] /= vol;
        data->par.grad.u[i] /= vol;
        data->par.grad.v[i] /= vol;
        data->par.grad.w[i] /= vol;
        data->par.grad.p[i] /= vol;
    }
}


void charm_calc_grad(p4est_t * p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data)
{
    DBG_CH(p4est->mpirank);
    int my_ghost = 0;
    if (!ghost) {
        DBG_CH(p4est->mpirank);

        ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
        ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
        my_ghost = 1;
    }

    DBG_CH(p4est->mpirank);
    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   charm_grad_reset_quad_iter_fn,
                   charm_grad_face_iter_fn,
                   NULL,
                   NULL);


    DBG_CH(p4est->mpirank);
    p4est_iterate (p4est, ghost, (void *) ghost_data,     /* pass in ghost data that we just exchanged */
                   charm_grad_quad_iter_fn,       /* blank the previously calculated derivatives */
                   NULL,
                   NULL,
                   NULL);

    if (ghost && my_ghost) {
        p4est_ghost_destroy (ghost);
        P4EST_FREE (ghost_data);
        ghost      = NULL;
        ghost_data = NULL;
    }

    DBG_CH(p4est->mpirank);
}

