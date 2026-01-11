//
// Created by zhrv on 27.08.19.
//

#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"


charm_real_t charm_model_ns_get_mu(p4est_t *p4est, charm_real_t *x, charm_data_t *data);


/*
 *  Volume integrals
 */


static void charm_model_ns_diff_integrals_volume_int_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = charm_get_quad_data(q);
    int                 ibf, igp;
    charm_cons_t        c;
    charm_prim_t        p;
    charm_real_t        phi_x, phi_y, phi_z, phi;
    charm_real_t       *x;
    charm_vec_t         qt;
    charm_tensor_t      tau;

    for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
        for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
            x = data->par.g.quad_gp[igp];
            charm_get_fields(data, x, &c);
            charm_param_cons_to_prim(info->p4est, &p, &c);
            charm_get_visc_tau(data, x, &tau);
            charm_get_heat_q(data, x, qt);

            phi_x = charm_base_func_dx(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_y = charm_base_func_dy(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_z = charm_base_func_dz(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];

            phi = charm_base_func(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];

            data->int_ru[ibf] += (tau.xx*phi_x+tau.xy*phi_y+tau.xz*phi_z);
            data->int_rv[ibf] += (tau.xy*phi_x+tau.yy*phi_y+tau.yz*phi_z);
            data->int_rw[ibf] += (tau.xz*phi_x+tau.yz*phi_y+tau.zz*phi_z);

            data->int_re[ibf] += (tau.xx*p.u+tau.xy*p.v+tau.xz*p.w - qt[0])*phi_x;
            data->int_re[ibf] += (tau.xy*p.u+tau.yy*p.v+tau.yz*p.w - qt[1])*phi_y;
            data->int_re[ibf] += (tau.xz*p.u+tau.yz*p.v+tau.zz*p.w - qt[2])*phi_z;

            data->int_re[ibf] += data->par.model.ns.chem_rhs*phi;
        }
    }
}


/*
 * Surface integrals
 */
static void charm_model_ns_conv_surface_int_iter_bnd (p4est_iter_face_info_t * info, void *user_data) {
    int i, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    charm_real_t n[3];
    charm_real_t qu, qv, qw, qe, qt[3];
    charm_real_t                  fu, fv, fw, ft;
    charm_real_t bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);
    charm_tensor_t          tau[2], ftau;
    int8_t face;
    charm_real_t c[2][3], l[3];
    charm_cons_t cons;
    charm_prim_t prim[2];
    charm_real_t mu, kt;
    charm_vec_t  x;


    CHARM_ASSERT(info->tree_boundary);

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    CHARM_ASSERT(!side[0]->is_hanging);

    charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[0]->treeid);
    if (attr->bnd[side[0]->face]->type == BOUND_WALL_NO_SLIP) {
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

        charm_face_get_center(udata, face, x);
        mu = charm_model_ns_get_mu(p4est, x, udata);
        kt = charm_get_heat_k(info->p4est, x, udata);
        charm_get_fields_avg(udata, &cons);
        charm_param_cons_to_prim(p4est, &(prim[0]), &cons);
        charm_get_heat_q(udata, x, qt);
        charm_bnd_cond(p4est, side[0]->treeid, face, &(prim[0]), &(prim[1]), n);
        charm_real_t vv[3] = {prim[0].u, prim[0].v, prim[0].w};
        charm_real_t un = scalar_prod(vv, n);
        charm_real_t vn[3] = {un*n[0], un*n[1], un*n[2]};
        charm_real_t vt[3] = {vv[0]-vn[0], vv[1]-vn[1], vv[2]-vn[2]};
        charm_real_t ll = vector_length(l);
        qu = -mu*vt[0]/ll;
        qv = -mu*vt[1]/ll;
        qw = -mu*vt[2]/ll;
        qe = -qt[0]*n[0]-qt[1]*n[1]-qt[2]*n[2];
        for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
            if (!side[0]->is.full.is_ghost) {
                bfv = charm_base_func(x, ibf, udata) * udata->par.g.area[face];
                udata->int_ru[ibf] -= qu * bfv;
                udata->int_rv[ibf] -= qv * bfv;
                udata->int_rw[ibf] -= qw * bfv;
                udata->int_re[ibf] -= qe * bfv;
            }
        }
    }
}


static void charm_model_ns_conv_surface_int_iter_inner (p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side, igp, ibf,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    charm_real_t                  n[3];
    charm_real_t                  qu, qv, qw, qe, qt[3];
    charm_real_t                  fu, fv, fw, ft;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
    charm_tensor_t          tau[2], ftau;
    charm_real_t                 *x, gw, gj;
    charm_real_t                  bfv;
    charm_real_t                  c[2][3];
    charm_real_t                  l[3];
    int8_t                  face[2];

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

            for (igp = 0; igp < CHARM_FACE_GP_COUNT; igp++) {
                x  = udata[h_side]->par.g.face_gp[face[h_side]][igp];
                gw = udata[h_side]->par.g.face_gw[face[h_side]][igp];
                gj = udata[h_side]->par.g.face_gj[face[h_side]][igp];
                charm_tensor_zero(&ftau);
                fu = fv = fw = ft = 0.;
                for (i = 0; i < 2; i++) {
                    charm_get_fields(udata[i], x, &(cons[i]));
                    charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                    charm_get_heat_q(udata[i], x, qt);
                    charm_get_visc_tau(udata[i], x, &(tau[i]));
                    charm_tensor_add(&ftau, &(tau[i]));
                    fu += prim[i].u*tau[i].xx + prim[i].v*tau[i].xy + prim[i].w*tau[i].xz;
                    fv += prim[i].u*tau[i].xy + prim[i].v*tau[i].yy + prim[i].w*tau[i].yz;
                    fw += prim[i].u*tau[i].xz + prim[i].v*tau[i].yz + prim[i].w*tau[i].zz;
                    ft += qt[0]*n[0] + qt[1]*n[1] + qt[2]*n[2];
                }
                charm_tensor_mul_scalar(&ftau, 0.5);
                fu *= 0.5;
                fv *= 0.5;
                fw *= 0.5;
                ft *= 0.5;
                qu = ftau.xx*n[0] + ftau.xy*n[1] + ftau.xz*n[2];
                qv = ftau.xy*n[0] + ftau.yy*n[1] + ftau.yz*n[2];
                qw = ftau.xz*n[0] + ftau.yz*n[1] + ftau.zz*n[2];
                qe = fu*n[0] + fv*n[1] + fw*n[2] - ft;
                for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                    for (i = 0; i < 2; i++) {
                        if (i == h_side) {
                            if (!side[i]->is.hanging.is_ghost[j]) {
                                bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                                udata[i]->int_ru[ibf] -= qu * bfv;
                                udata[i]->int_rv[ibf] -= qv * bfv;
                                udata[i]->int_rw[ibf] -= qw * bfv;
                                udata[i]->int_re[ibf] -= qe * bfv;
                            }
                        }
                        else {
                            if (!side[i]->is.full.is_ghost) {
                                bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                                udata[i]->int_ru[ibf] -= qu * bfv;
                                udata[i]->int_rv[ibf] -= qv * bfv;
                                udata[i]->int_rw[ibf] -= qw * bfv;
                                udata[i]->int_re[ibf] -= qe * bfv;
                            }
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

        for (igp = 0; igp < CHARM_FACE_GP_COUNT; igp++) {
            x  = udata[0]->par.g.face_gp[face[0]][igp];
            gw = udata[0]->par.g.face_gw[face[0]][igp];
            gj = udata[0]->par.g.face_gj[face[0]][igp];
            charm_tensor_zero(&ftau);
            fu = fv = fw = ft = 0.;
            for (i = 0; i < 2; i++) {
                charm_get_fields(udata[i], x, &(cons[i]));
                charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                charm_get_heat_q(udata[i], x, qt);
                charm_get_visc_tau(udata[i], x, &(tau[i]));
                charm_tensor_add(&ftau, &(tau[i]));
                fu += prim[i].u*tau[i].xx + prim[i].v*tau[i].xy + prim[i].w*tau[i].xz;
                fv += prim[i].u*tau[i].xy + prim[i].v*tau[i].yy + prim[i].w*tau[i].yz;
                fw += prim[i].u*tau[i].xz + prim[i].v*tau[i].yz + prim[i].w*tau[i].zz;
                ft += qt[0]*n[0] + qt[1]*n[1] + qt[2]*n[2];
            }
            charm_tensor_mul_scalar(&ftau, 0.5);
            fu *= 0.5;
            fv *= 0.5;
            fw *= 0.5;
            ft *= 0.5;
            qu = ftau.xx*n[0] + ftau.xy*n[1] + ftau.xz*n[2];
            qv = ftau.xy*n[0] + ftau.yy*n[1] + ftau.yz*n[2];
            qw = ftau.xz*n[0] + ftau.yz*n[1] + ftau.zz*n[2];
            qe = fu*n[0] + fv*n[1] + fw*n[2] - ft;
            for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                for (i = 0; i < 2; i++) {
                    if (!side[i]->is.full.is_ghost) {
                        bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                        udata[i]->int_ru[ibf] -= qu * bfv;
                        udata[i]->int_rv[ibf] -= qv * bfv;
                        udata[i]->int_rw[ibf] -= qw * bfv;
                        udata[i]->int_re[ibf] -= qe * bfv;
                    }
                }
            }
        }
    }
}


static void charm_model_ns_diff_integrals_surface_int_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        charm_model_ns_conv_surface_int_iter_bnd(info, user_data);
    }
    else {
        charm_model_ns_conv_surface_int_iter_inner(info, user_data);
    }

}


void charm_model_ns_timestep_diff_integrals(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   charm_model_ns_diff_integrals_volume_int_iter_fn,
                   charm_model_ns_diff_integrals_surface_int_iter_fn,
                   NULL, NULL);

}
