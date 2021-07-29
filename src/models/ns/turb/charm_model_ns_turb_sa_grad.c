//
// Created by zhrv on 27.09.2020.
//

#include "charm_globals.h"
#include "charm_base_func.h"
#include "charm_bnd_cond.h"

static void charm_model_ns_turb_sa_grad_zero_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t *data = (charm_data_t *) info->quad->p.user_data;
    int i;
    for (i = 0; i < CHARM_DIM; i++) {
        data->par.model.ns.turb.model.sa.grad_nu_[i] = 0.;
        memset(data->par.model.ns.turb.model.sa.grad_u[i], 0, sizeof(charm_vec_t));
    }
}


static void charm_model_ns_turb_sa_grad_update_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t *data = (charm_data_t *) info->quad->p.user_data;
    charm_real_t volume = data->par.g.volume;
    int i;
    for (i = 0; i < CHARM_DIM; i++) {
        data->par.model.ns.turb.model.sa.grad_nu_[i] /= volume;
        data->par.model.ns.turb.model.sa.grad_u[i][0] /= volume;
        data->par.model.ns.turb.model.sa.grad_u[i][1] /= volume;
        data->par.model.ns.turb.model.sa.grad_u[i][2] /= volume;
    }
}


static void charm_model_ns_turb_sa_grad_surface_int_iter_bnd(p4est_iter_face_info_t * info, void *user_data)
{
    int i, j, k;
    p4est_t                    *p4est = info->p4est;
    charm_ctx_t                *ctx = charm_get_ctx(p4est);
    charm_data_t               *ghost_data = (charm_data_t *) user_data;
    charm_data_t               *udata;
    charm_real_t                n[3];
    charm_real_t                qu;
    charm_vec_t                 gu;
    p4est_iter_face_side_t     *side[2];
    sc_array_t                 *sides = &(info->sides);
    size_t                      c_count = charm_get_comp_count(info->p4est);


    int8_t face;
    charm_real_t  c[2][3], l[3];
    charm_real_t  nu_[2], *int_nu;
    charm_cons_t  cons;
    charm_prim_t  prim[2];
    charm_real_t  s;
    charm_vec_t x;
    charm_real_t  intg[2][5];


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

    int_nu = &(udata->par.model.ns.turb.model.sa.int_nu_);

    for (i = 0; i < 3; i++) {
        l[i] = c[1][i] - c[0][i];
    }

    if (scalar_prod(n, l) < 0) {
        for (i = 0; i < 3; i++) {
            n[i] *= -1.0;
        }
    }

    charm_face_get_center(udata, face, x);
    s = charm_face_get_area(udata, face);
    nu_[0] = udata->par.model.ns.turb.model.sa.nu_;
    charm_get_fields(udata, x, &cons);
    charm_param_cons_to_prim(p4est, &(prim[0]), &cons);
    charm_bnd_cond(p4est, side[0]->treeid, face, &(prim[0]), &(prim[1]), n);
    qu = 0.;
    gu[0] = gu[1] = gu[2] = 0.;
    for (i = 0; i < 2; i++) {
        charm_get_fields(udata, x, &cons);
        charm_param_cons_to_prim(p4est, &(prim[i]), &cons);
        qu += udata->par.model.ns.turb.model.sa.nu_; // TODO оптимизировать
        gu[0] += prim[i].u;
        gu[1] += prim[i].v;
        gu[2] += prim[i].w;
    }
    qu *= 0.5;
    gu[0] *= 0.5;
    gu[1] *= 0.5;
    gu[2] *= 0.5;
    if (!side[0]->is.full.is_ghost) {
        for (k = 0; k < CHARM_DIM; k++) {
            udata->par.model.ns.turb.model.sa.grad_nu_[i] += qu * s * n[i];
            udata->par.model.ns.turb.model.sa.grad_u[0][i] += gu[0] * s * n[i];
            udata->par.model.ns.turb.model.sa.grad_u[1][i] += gu[1] * s * n[i];
            udata->par.model.ns.turb.model.sa.grad_u[2][i] += gu[2] * s * n[i];
        }
    }

}


static void charm_model_ns_turb_sa_grad_surface_int_iter_inner(p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, k, h_side,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    charm_real_t            n[3];
    charm_real_t            qu;
    charm_vec_t             gu;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_real_t            s;
    charm_real_t            c[2][3];
    charm_real_t            l[3];
    int8_t                  face[2];
    charm_prim_t            prim[2];
    charm_cons_t            cons[2];
    charm_vec_t           x;

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    face[0] = side[0]->face;
    face[1] = side[1]->face;

    h_side = -1;
    if (side[0]->is_hanging || side[1]->is_hanging) { // @todo
        for (j = 0; j < CHARM_HALF; j++) {
            for (i = 0; i < 2; i++) {
                if (side[i]->is_hanging) {
                    if (side[i]->is.hanging.is_ghost[j]) {
                        udata[i] = &(ghost_data[side[i]->is.hanging.quadid[j]]);
                    }
                    else {
                        udata[i] = (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
                    }
                    h_side = i;
                }
                else {
                    if (side[i]->is.full.is_ghost) {
                        udata[i] = &ghost_data[side[i]->is.full.quadid];
                    }
                    else {
                        udata[i] = (charm_data_t *) side[i]->is.full.quad->p.user_data;
                    }
                }
            }

            CHARM_ASSERT(h_side != -1);

            charm_face_get_normal(udata[0], face[0], n);
            charm_quad_get_center(udata[0], c[0]);
            charm_face_get_center(udata[0], face[0], c[1]);

            for (i = 0; i < 3; i++) {
                l[i] = c[1][i]-c[0][i];
            }

            if (scalar_prod(n, l) < 0) {
                for (i = 0; i < 3; i++) {
                    n[i] *= -1.0;
                }
            }

            charm_face_get_center(udata[h_side], face[h_side], x);
            s = charm_face_get_area(udata[h_side], face[h_side]);
            qu = 0.;
            gu[0] = gu[1] = gu[2] = 0.;
            for (i = 0; i < 2; i++) {
                charm_get_fields(udata[i], x, &(cons[i]));
                charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                qu += udata[i]->par.model.ns.turb.model.sa.nu_;
                gu[0] += prim[i].u;
                gu[1] += prim[i].v;
                gu[2] += prim[i].w;
            }
            qu *= 0.5;
            gu[0] *= 0.5;
            gu[1] *= 0.5;
            gu[2] *= 0.5;
            for (i = 0; i < 2; i++) {
                if (i == h_side) {
                    if (!side[i]->is.hanging.is_ghost[j]) {
                        for (k = 0; k < CHARM_DIM; k++) {
                            udata[i]->par.model.ns.turb.model.sa.grad_nu_[i] += qu * (i ? -1. : 1.) * s * n[i];
                            udata[i]->par.model.ns.turb.model.sa.grad_u[0][i] += gu[0] * (i ? -1. : 1.) * s * n[i];
                            udata[i]->par.model.ns.turb.model.sa.grad_u[1][i] += gu[1] * (i ? -1. : 1.) * s * n[i];
                            udata[i]->par.model.ns.turb.model.sa.grad_u[2][i] += gu[2] * (i ? -1. : 1.) * s * n[i];
                        }
                    }
                }
                else {
                    if (!side[i]->is.full.is_ghost) {
                        for (k = 0; k < CHARM_DIM; k++) {
                            udata[i]->par.model.ns.turb.model.sa.grad_nu_[i] += qu * (i ? -1. : 1.) * s * n[i];
                            udata[i]->par.model.ns.turb.model.sa.grad_u[0][i] += gu[0] * (i ? -1. : 1.) * s * n[i];
                            udata[i]->par.model.ns.turb.model.sa.grad_u[1][i] += gu[1] * (i ? -1. : 1.) * s * n[i];
                            udata[i]->par.model.ns.turb.model.sa.grad_u[2][i] += gu[2] * (i ? -1. : 1.) * s * n[i];
                        }
                    }
                }
            }
        }
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
        charm_face_get_normal(udata[0], face[0], n);
        charm_quad_get_center(udata[0], c[0]);
        charm_face_get_center(udata[0], face[0], c[1]);

        for (i = 0; i < 3; i++) {
            l[i] = c[1][i]-c[0][i];
        }

        if (scalar_prod(n, l) < 0) {
            for (i = 0; i < 3; i++) {
                n[i] *= -1.0;
            }
        }

        charm_face_get_center(udata[0], face[0], x);
        s = charm_face_get_area(udata[0], face[0]);
        qu = 0.;
        gu[0] = gu[1] = gu[2] = 0.;
        for (i = 0; i < 2; i++) {
            charm_get_fields(udata[i], x, &(cons[i]));
            charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
            qu    += udata[i]->par.model.ns.turb.model.sa.nu_;
            gu[0] += prim[i].u;
            gu[1] += prim[i].v;
            gu[2] += prim[i].w;
        }
        qu    *= 0.5;
        gu[0] *= 0.5;
        gu[1] *= 0.5;
        gu[2] *= 0.5;

        for (i = 0; i < 2; i++) {
            if (!side[i]->is.full.is_ghost) {
                for (k = 0; k < CHARM_DIM; k++) {
                    udata[i]->par.model.ns.turb.model.sa.grad_nu_[i] += qu * (i ? -1. : 1.) * s * n[i];
                    udata[i]->par.model.ns.turb.model.sa.grad_u[0][i] += gu[0] * (i ? -1. : 1.) * s * n[i];
                    udata[i]->par.model.ns.turb.model.sa.grad_u[1][i] += gu[1] * (i ? -1. : 1.) * s * n[i];
                    udata[i]->par.model.ns.turb.model.sa.grad_u[2][i] += gu[2] * (i ? -1. : 1.) * s * n[i];
                }
            }
        }
    }
}


static void charm_model_ns_turb_sa_grad_surface_int_iter_fn(p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        charm_model_ns_turb_sa_grad_surface_int_iter_bnd(info, user_data);
    }
    else {
        charm_model_ns_turb_sa_grad_surface_int_iter_inner(info, user_data);
    }

}

void charm_model_ns_turb_sa_grad(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   charm_model_ns_turb_sa_grad_zero_quad_iter_fn, NULL, NULL, NULL);

    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   NULL, charm_model_ns_turb_sa_grad_surface_int_iter_fn, NULL, NULL);

    p4est_iterate (p4est, NULL, NULL,
                   charm_model_ns_turb_sa_grad_update_quad_iter_fn, NULL, NULL, NULL);

    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

}