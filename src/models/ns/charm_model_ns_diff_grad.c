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


static void _charm_model_ns_diff_grad_volume_int_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = charm_get_quad_data(q);
    int                 ibf, igp;
    charm_cons_t        c;
    charm_prim_t        p;
    double              phi_x, phi_y, phi_z, phi, tmp_qx, tmp_qy, tmp_qz;
    double             *x;
    double              lambda, mu, kt, lp, lm, h;
    int                 i;
    size_t              c_count = charm_get_comp_count(info->p4est);
    charm_comp_t       *comp;

    for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
        for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
            x = data->par.g.quad_gp[igp];

            kt      = charm_get_heat_k(info->p4est, x, data);
            lambda  = charm_get_visc_lambda(info->p4est, data);
            mu      = charm_get_visc_mu(info->p4est, x, data);
            lp = (lambda+4.*mu/3.);
            lm = (lambda-2.*mu/3.);
            charm_get_fields(data, x, &c);
            charm_param_cons_to_prim(info->p4est, &p, &c);

            phi_x = charm_base_func_dx(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_y = charm_base_func_dy(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_z = charm_base_func_dz(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];

            data->int_tau_xx[ibf] -= (lp*p.u*phi_x+lm*p.v*phi_y+lm*p.w*phi_z);
            data->int_tau_yy[ibf] -= (lm*p.u*phi_x+lp*p.v*phi_y+lm*p.w*phi_z);
            data->int_tau_zz[ibf] -= (lm*p.u*phi_x+lm*p.v*phi_y+lp*p.w*phi_z);

            data->int_tau_xy[ibf] -= mu*(p.u*phi_y+p.v*phi_x);
            data->int_tau_xz[ibf] -= mu*(p.w*phi_x+p.u*phi_z);
            data->int_tau_yz[ibf] -= mu*(p.w*phi_y+p.v*phi_z);

            tmp_qx = kt*p.t*phi_x;
            tmp_qy = kt*p.t*phi_y;
            tmp_qz = kt*p.t*phi_z;
            for (i = 0; i < c_count; i++) {
                comp = charm_get_comp(info->p4est, i);
                h = charm_comp_calc_enthalpy(comp, p.t);
                tmp_qx += h*p.r*data->par.model.ns.d[i]*p.c[i]*phi_x;
                tmp_qy += h*p.r*data->par.model.ns.d[i]*p.c[i]*phi_y;
                tmp_qz += h*p.r*data->par.model.ns.d[i]*p.c[i]*phi_z;
            }
            data->int_q_x[ibf] += tmp_qx;
            data->int_q_y[ibf] += tmp_qy;
            data->int_q_z[ibf] += tmp_qz;
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
    double qxx, qyy, qzz, qxy, qxz, qyz, qtx, qty, qtz;
    double fu, fv, fw, ft;
    double bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);
    size_t              c_count = charm_get_comp_count(info->p4est);


    int8_t face;
    double c[2][3], l[3];
    double r_[2], p_[2], u_[2], v_[2], w_[2], e_[2];
    charm_cons_t cons;
    charm_prim_t prim[2];
    charm_comp_t *comp;
    double *x, gw, gj;
    double lambda, mu, lp, lm, kt, fr, fc, h;
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

    face = side[0]->face;

    charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[0]->treeid);

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

        lambda   = charm_get_visc_lambda(p4est, udata);
        mu       = charm_get_visc_mu(p4est, x, udata);
        lp  = (lambda+4.*mu/3.);
        lm  = (lambda-2.*mu/3.);
        kt  = charm_get_heat_k(p4est, x, udata);

        charm_get_fields(udata, x, &cons);
        charm_param_cons_to_prim(p4est, &(prim[0]), &cons);
        charm_bnd_cond(p4est, side[0]->treeid, face, &(prim[0]), &(prim[1]), n);
        if (attr->bnd[face]->type == BOUND_WALL_NO_SLIP) {
            fu = fv = fw = 0.;                              //скорости на границе с прилипанием равны нулю
            ft = attr->bnd[face]->params[0];
        }
        else {
            fu = (prim[0].u+prim[1].u)*0.5;
            fv = (prim[0].v+prim[1].v)*0.5;
            fw = (prim[0].w+prim[1].w)*0.5;
            ft = (prim[0].t+prim[1].t)*0.5;
        }
        fr = (prim[0].r+prim[1].r)*0.5;
        qxx = lp*fu*n[0] + lm*fv*n[1] + lm*fw*n[2];
        qyy = lm*fu*n[0] + lp*fv*n[1] + lm*fw*n[2];
        qzz = lm*fu*n[0] + lm*fv*n[1] + lp*fw*n[2];
        qxy = mu*(fu*n[1] + fv*n[0]);
        qxz = mu*(fw*n[0] + fu*n[2]);
        qyz = mu*(fw*n[1] + fv*n[2]);
        qtx = kt*ft*n[0];
        qty = kt*ft*n[1];
        qtz = kt*ft*n[2];
        for (i = 0; i < c_count; i++) {
            comp = charm_get_comp(p4est, i);
            h = charm_comp_calc_enthalpy(comp, ft);
            fc = (prim[0].c[i]+prim[1].c[i])*0.5;
            qtx += h*fr*udata->par.model.ns.d[i]*fc*n[0];
            qty += h*fr*udata->par.model.ns.d[i]*fc*n[1];
            qtz += h*fr*udata->par.model.ns.d[i]*fc*n[2];
        }
        for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
            if (!side[0]->is.full.is_ghost) {
                bfv = charm_base_func(x, ibf, udata) * gw * gj;

                udata->int_tau_xx[ibf] += qxx * bfv;
                udata->int_tau_yy[ibf] += qyy * bfv;
                udata->int_tau_zz[ibf] += qzz * bfv;
                udata->int_tau_xy[ibf] += qxy * bfv;
                udata->int_tau_xz[ibf] += qxz * bfv;
                udata->int_tau_yz[ibf] += qyz * bfv;

                udata->int_q_x[ibf] -= qtx * bfv;
                udata->int_q_y[ibf] -= qty * bfv;
                udata->int_q_z[ibf] -= qtz * bfv;
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
    double                  qxx, qyy, qzz, qxy, qxz, qyz, qtx, qty, qtz;
    double                  fu, fv, fw, ft, fr, fc, h;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
    charm_comp_t           *comp;
    double                 *x, gw, gj;
    double                  bfv;
    double                  c[2][3];
    double                  l[3];
    int8_t                  face[2];
    double                  lambda[2], mu[2], flambda, fmu, lp, lm, kt, dm;
    size_t                  c_count = charm_get_comp_count(p4est);

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
                kt = 0.;
                for (i = 0; i < 2; i++) {
                    charm_get_fields(udata[i], x, &(cons[i]));
                    charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                    lambda[i]   = charm_get_visc_lambda(p4est, udata[i]);
                    mu[i]       = charm_get_visc_mu(p4est, x, udata[i]);
                    kt         += charm_get_heat_k(p4est, x, udata[i]);
                }
                fmu     = 0.5*(mu[0]+mu[1]);
                flambda = 0.5*(lambda[0]+lambda[1]);
                lp = (flambda+4.*fmu/3.);
                lm = (flambda-2.*fmu/3.);
                kt *= 0.5;

                fr = (prim[0].r+prim[1].r)*0.5;
                ft = (prim[0].t+prim[1].t)*0.5;
                fu = (prim[0].u+prim[1].u)*0.5;
                fv = (prim[0].v+prim[1].v)*0.5;
                fw = (prim[0].w+prim[1].w)*0.5;
                qxx = lp*fu*n[0] + lm*fv*n[1] + lm*fw*n[2];
                qyy = lm*fu*n[0] + lp*fv*n[1] + lm*fw*n[2];
                qzz = lm*fu*n[0] + lm*fv*n[1] + lp*fw*n[2];
                qxy = fmu*(fu*n[1] + fv*n[0]);
                qxz = fmu*(fw*n[0] + fu*n[2]);
                qyz = fmu*(fw*n[1] + fv*n[2]);
                qtx = kt*ft*n[0];
                qty = kt*ft*n[1];
                qtz = kt*ft*n[2];
                for (i = 0; i < c_count; i++) {
                    comp = charm_get_comp(p4est, i);
                    h = charm_comp_calc_enthalpy(comp, ft);
                    fc = (prim[0].c[i]+prim[1].c[i])*0.5;
                    dm = (udata[0]->par.model.ns.d[i]+udata[1]->par.model.ns.d[i])*0.5;
                    qtx += h*fr*dm*fc*n[0];
                    qty += h*fr*dm*fc*n[1];
                    qtz += h*fr*dm*fc*n[2];
                }
                for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                    for (i = 0; i < 2; i++) {
                        if (i == h_side) {
                            if (!side[i]->is.hanging.is_ghost[j]) {
                                bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                                udata[i]->int_tau_xx[ibf] += qxx * bfv;
                                udata[i]->int_tau_yy[ibf] += qyy * bfv;
                                udata[i]->int_tau_zz[ibf] += qzz * bfv;
                                udata[i]->int_tau_xy[ibf] += qxy * bfv;
                                udata[i]->int_tau_xz[ibf] += qxz * bfv;
                                udata[i]->int_tau_yz[ibf] += qyz * bfv;

                                udata[i]->int_q_x[ibf] -= qtx * bfv;
                                udata[i]->int_q_y[ibf] -= qty * bfv;
                                udata[i]->int_q_z[ibf] -= qtz * bfv;
                            }
                        }
                        else {
                            if (!side[i]->is.full.is_ghost) {
                                bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                                udata[i]->int_tau_xx[ibf] += qxx * bfv;
                                udata[i]->int_tau_yy[ibf] += qyy * bfv;
                                udata[i]->int_tau_zz[ibf] += qzz * bfv;
                                udata[i]->int_tau_xy[ibf] += qxy * bfv;
                                udata[i]->int_tau_xz[ibf] += qxz * bfv;
                                udata[i]->int_tau_yz[ibf] += qyz * bfv;

                                udata[i]->int_q_x[ibf] -= qtx * bfv;
                                udata[i]->int_q_y[ibf] -= qty * bfv;
                                udata[i]->int_q_z[ibf] -= qtz * bfv;
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
            kt = 0.;
            for (i = 0; i < 2; i++) {
                charm_get_fields(udata[i], x, &(cons[i]));
                charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                lambda[i]   = charm_get_visc_lambda(p4est, udata[i]);
                mu[i]       = charm_get_visc_mu(p4est, x, udata[i]);
                kt         += charm_get_heat_k(p4est, x, udata[i]);
            }
            fmu     = 0.5*(mu[0]+mu[1]);
            flambda = 0.5*(lambda[0]+lambda[1]);
            lp = (flambda+4.*fmu/3.);
            lm = (flambda-2.*fmu/3.);
            kt *= 0.5;

            fr = (prim[0].r+prim[1].r)*0.5;
            ft = (prim[0].t+prim[1].t)*0.5;
            fu = (prim[0].u+prim[1].u)*0.5;
            fv = (prim[0].v+prim[1].v)*0.5;
            fw = (prim[0].w+prim[1].w)*0.5;
            qxx = lp*fu*n[0] + lm*fv*n[1] + lm*fw*n[2];
            qyy = lm*fu*n[0] + lp*fv*n[1] + lm*fw*n[2];
            qzz = lm*fu*n[0] + lm*fv*n[1] + lp*fw*n[2];
            qxy = fmu*(fu*n[1] + fv*n[0]);
            qxz = fmu*(fw*n[0] + fu*n[2]);
            qyz = fmu*(fw*n[1] + fv*n[2]);
            qtx = kt*ft*n[0];
            qty = kt*ft*n[1];
            qtz = kt*ft*n[2];
            for (i = 0; i < c_count; i++) {
                comp = charm_get_comp(p4est, i);
                h = charm_comp_calc_enthalpy(comp, ft);
                fc = (prim[0].c[i]+prim[1].c[i])*0.5;
                dm = (udata[0]->par.model.ns.d[i]+udata[1]->par.model.ns.d[i])*0.5;
                qtx += h*fr*dm*fc*n[0];
                qty += h*fr*dm*fc*n[1];
                qtz += h*fr*dm*fc*n[2];
            }
            for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                for (i = 0; i < 2; i++) {
                    if (!side[i]->is.full.is_ghost) {
                        bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                        udata[i]->int_tau_xx[ibf] += qxx * bfv;
                        udata[i]->int_tau_yy[ibf] += qyy * bfv;
                        udata[i]->int_tau_zz[ibf] += qzz * bfv;
                        udata[i]->int_tau_xy[ibf] += qxy * bfv;
                        udata[i]->int_tau_xz[ibf] += qxz * bfv;
                        udata[i]->int_tau_yz[ibf] += qyz * bfv;

                        udata[i]->int_q_x[ibf] -= qtx * bfv;
                        udata[i]->int_q_y[ibf] -= qty * bfv;
                        udata[i]->int_q_z[ibf] -= qtz * bfv;
                    }
                }
            }
        }
    }
}

static void _charm_model_ns_diff_grad_surface_int_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_model_ns_conv_surface_int_iter_bnd(info, user_data);
    }
    else {
        _charm_model_ns_conv_surface_int_iter_inner(info, user_data);
    }

}


static void _charm_model_ns_diff_grad_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);

    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_xx, data->par.model.ns.tau.xx);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_yy, data->par.model.ns.tau.yy);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_zz, data->par.model.ns.tau.zz);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_xy, data->par.model.ns.tau.xy);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_xz, data->par.model.ns.tau.xz);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_tau_yz, data->par.model.ns.tau.yz);

    charm_matr_vect_mult(data->par.g.a_inv, data->int_q_x, data->par.model.ns.q.x);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_q_y, data->par.model.ns.q.y);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_q_z, data->par.model.ns.q.z);
}


static void _charm_model_ns_diff_grad_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
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

        data->int_q_x[i] = 0.;
        data->int_q_y[i] = 0.;
        data->int_q_z[i] = 0.;
    }
}




void charm_model_ns_timestep_diff_grad(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_model_ns_diff_grad_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_model_ns_diff_grad_volume_int_iter_fn,
                   _charm_model_ns_diff_grad_surface_int_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est, NULL,
                   NULL,
                   _charm_model_ns_diff_grad_update_quad_iter_fn,
                   NULL, NULL, NULL);


    p4est_ghost_exchange_data (p4est, ghost, ghost_data); /* synchronize the ghost data */
}
