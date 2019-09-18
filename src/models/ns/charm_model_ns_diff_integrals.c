//
// Created by zhrv on 27.08.19.
//

#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"




/*
 *  Volume integrals
 */


static void _charm_model_ns_diff_integrals_volume_int_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = charm_get_quad_data(q);
    int                 ibf, igp;
    charm_cons_t        c;
    charm_prim_t        p;
    double              phi_x, phi_y, phi_z, phi;
    double             *x;
    double              lambda  = charm_get_visc_lambda(info->p4est, data);
    double              mu      = charm_get_visc_mu(info->p4est, data);

    double              lp = (lambda+4.*mu/3.);
    double              lm = (lambda-2.*mu/3.);
    charm_tensor_t      tau;

    for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
        for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
            x = data->par.g.quad_gp[igp];

            charm_get_fields(data, x, &c);
            charm_param_cons_to_prim(info->p4est, &p, &c);
            charm_get_visc_tau(data, x, &tau);

            phi_x = charm_base_func_dx(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_y = charm_base_func_dy(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_z = charm_base_func_dz(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];

            data->int_ru[ibf] += (tau.xx*phi_x+tau.xy*phi_y+tau.xz*phi_z);
            data->int_rv[ibf] += (tau.xy*phi_x+tau.yy*phi_y+tau.yz*phi_z);
            data->int_rw[ibf] += (tau.xz*phi_x+tau.yz*phi_y+tau.zz*phi_z);

            data->int_re[ibf] += (tau.xx*p.u+tau.xy*p.v+tau.xz*p.w)*phi_x;
            data->int_re[ibf] += (tau.xy*p.u+tau.yy*p.v+tau.yz*p.w)*phi_y;
            data->int_re[ibf] += (tau.xz*p.u+tau.yz*p.v+tau.zz*p.w)*phi_z;
        }
    }
}


/*
 * Surface integrals
 */


static void _charm_model_ns_conv_surface_int_iter_bnd (p4est_iter_face_info_t * info, void *user_data) {
    int i, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_ctx_t * ctx = charm_get_ctx(p4est);
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    double n[3];
    double qxx, qyy, qzz, qxy, qxz, qyz;
    double fu, fv, fw;
    double bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);
    size_t              c_count = charm_get_comp_count(info->p4est);


    int8_t face;
    double c[2][3], l[3];
    double r_[2], p_[2], u_[2], v_[2], w_[2], e_[2];
    charm_cons_t cons;
    charm_prim_t prim[2];
    double *x, gw, gj;
    double lambda, mu, lp, lm;
    int j;


    CHARM_ASSERT(info->tree_boundary);

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
            CHARM_ASSERT(!side[0]->is_hanging);

    if (side[0]->is.full.is_ghost) {
                CHARM_ASSERT(0);
        udata = &(ghost_data[side[0]->is.full.quadid]);
    } else {
        udata = charm_get_quad_data(side[0]->is.full.quad);
    }
    lambda   = charm_get_visc_lambda(p4est, udata);
    mu       = charm_get_visc_mu(p4est, udata);
    lp = (lambda+4.*mu/3.);
    lm = (lambda-2.*mu/3.);

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

    for (igp = 0; igp < CHARM_FACE_GP_COUNT; igp++) {
        x = udata->par.g.face_gp[face][igp];
        gw = udata->par.g.face_gw[face][igp];
        gj = udata->par.g.face_gj[face][igp];
        charm_get_fields(udata, x, &cons);
        charm_param_cons_to_prim(p4est, &(prim[0]), &cons);
        charm_bnd_cond(p4est, side[0]->treeid, face, &(prim[0]), &(prim[1]), n);
        fu = (prim[0].u+prim[1].u)*0.5;
        fv = (prim[0].v+prim[1].v)*0.5;
        fw = (prim[0].w+prim[1].w)*0.5;
        qxx = lp*fu*n[0] + lm*fv*n[1] + lm*fw*n[2];
        qyy = lm*fu*n[0] + lp*fv*n[1] + lm*fw*n[2];
        qzz = lm*fu*n[0] + lm*fv*n[1] + lp*fw*n[2];
        qxy = mu*(fu*n[1] + fv*n[0]);
        qxz = mu*(fw*n[0] + fu*n[2]);
        qyz = mu*(fw*n[1] + fv*n[2]);
        for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
            if (!side[0]->is.full.is_ghost) {
                bfv = charm_base_func(x, ibf, udata) * gw * gj;

                udata->int_tau_xx[ibf] += qxx * bfv;
                udata->int_tau_yy[ibf] += qyy * bfv;
                udata->int_tau_zz[ibf] += qzz * bfv;
                udata->int_tau_xy[ibf] += qxy * bfv;
                udata->int_tau_xz[ibf] += qxz * bfv;
                udata->int_tau_yz[ibf] += qyz * bfv;
            }
        }
    }
}


static void _charm_model_ns_conv_surface_int_iter_inner (p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side, igp, ibf,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    double                  n[3];
    double                  qu, qv, qw, qe;
    double                  fu, fv, fw;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
    charm_tensor_t          tau[2], ftau;
    double                 *x, gw, gj;
    double                  bfv;
    double                  c[2][3];
    double                  l[3];
    int8_t                  face[2];
    double                  lambda[2], mu[2], flambda, fmu, lp, lm;

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
                lambda[i]   = charm_get_visc_lambda(p4est, udata[i]);
                mu[i]       = charm_get_visc_mu(p4est, udata[i]);
            }

            fmu     = 0.5*(mu[0]+mu[1]);
            flambda = 0.5*(lambda[0]+lambda[1]);
            lp = (flambda+4.*fmu/3.);
            lm = (flambda-2.*fmu/3.);

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
                fu = fv = fw = 0.;
                for (i = 0; i < 2; i++) {
                    charm_get_fields(udata[i], x, &(cons[i]));
                    charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                    charm_get_visc_tau(udata[i], x, &(tau[i]));
                    charm_tensor_add(&ftau, &(tau[i]));
                    fu += prim[i].u*tau[i].xx + prim[i].v*tau[i].xy + prim[i].w*tau[i].xz;
                    fv += prim[i].u*tau[i].xy + prim[i].v*tau[i].yy + prim[i].w*tau[i].yz;
                    fw += prim[i].u*tau[i].xz + prim[i].v*tau[i].yz + prim[i].w*tau[i].zz;
                }
                charm_tensor_mul_scalar(&ftau, 0.5);
                fu *= 0.5;
                fv *= 0.5;
                fw *= 0.5;
                qu = ftau.xx*n[0] + ftau.xy*n[1] + ftau.xz*n[2];
                qv = ftau.xy*n[0] + ftau.yy*n[1] + ftau.yz*n[2];
                qw = ftau.xz*n[0] + ftau.yz*n[1] + ftau.zz*n[2];
                qe = fu*n[0] + fv*n[1] + fw*n[2];
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
    else {

        for (i = 0; i < 2; i++) {
            if (side[i]->is.full.is_ghost) {
                udata[i] = &(ghost_data[side[i]->is.full.quadid]);
            }
            else {
                udata[i] = charm_get_quad_data(side[i]->is.full.quad);//(charm_data_t *) side[i]->is.full.quad->p.user_data;
            }
            lambda[i]   = charm_get_visc_lambda(p4est, udata[i]);
            mu[i]       = charm_get_visc_mu(p4est, udata[i]);
        }

        fmu     = 0.5*(mu[0]+mu[1]);
        flambda = 0.5*(lambda[0]+lambda[1]);
        lp = (flambda+4.*fmu/3.);
        lm = (flambda-2.*fmu/3.);

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
            fu = fv = fw = 0.;
            for (i = 0; i < 2; i++) {
                charm_get_fields(udata[i], x, &(cons[i]));
                charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                charm_get_visc_tau(udata[i], x, &(tau[i]));
                charm_tensor_add(&ftau, &(tau[i]));
                fu += prim[i].u*tau[i].xx + prim[i].v*tau[i].xy + prim[i].w*tau[i].xz;
                fv += prim[i].u*tau[i].xy + prim[i].v*tau[i].yy + prim[i].w*tau[i].yz;
                fw += prim[i].u*tau[i].xz + prim[i].v*tau[i].yz + prim[i].w*tau[i].zz;
            }
            charm_tensor_mul_scalar(&ftau, 0.5);
            fu *= 0.5;
            fv *= 0.5;
            fw *= 0.5;
            qu = ftau.xx*n[0] + ftau.xy*n[1] + ftau.xz*n[2];
            qv = ftau.xy*n[0] + ftau.yy*n[1] + ftau.yz*n[2];
            qw = ftau.xz*n[0] + ftau.yz*n[1] + ftau.zz*n[2];
            qe = fu*n[0] + fv*n[1] + fw*n[2];
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

static void _charm_model_ns_diff_integrals_surface_int_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        //_charm_model_ns_conv_surface_int_iter_bnd(info, user_data);
    }
    else {
        _charm_model_ns_conv_surface_int_iter_inner(info, user_data);
    }

}


static void _charm_model_ns_diff_integrals_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);

    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_xx, data->par.tau.xx);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_yy, data->par.tau.yy);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_zz, data->par.tau.zz);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_xy, data->par.tau.xy);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_xz, data->par.tau.xz);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_yz, data->par.tau.yz);
}


static void _charm_model_ns_diff_integrals_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->int_tau_xx[i] = 0.;
        data->int_tau_yy[i] = 0.;
        data->int_tau_zz[i] = 0.;
        data->int_tau_xy[i] = 0.;
        data->int_tau_xz[i] = 0.;
        data->int_tau_yz[i] = 0.;
    }
}




void charm_model_ns_timestep_diff_integrals(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_model_ns_diff_integrals_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_model_ns_diff_integrals_volume_int_iter_fn,
                   _charm_model_ns_diff_integrals_surface_int_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est, NULL,
                   NULL,
                   _charm_model_ns_diff_integrals_update_quad_iter_fn,
                   NULL, NULL, NULL);


    p4est_ghost_exchange_data (p4est, ghost, ghost_data); /* synchronize the ghost data */
}
