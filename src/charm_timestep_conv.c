//
// Created by zhrv on 15.09.18.
//

#include <p8est_iterate.h>
#include "charm_timestep_conv.h"
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"

/*
 *  Volume integrals
 */


void charm_convect_volume_int_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    int                 ibf, igp;
    charm_cons_t        c;
    charm_prim_t        p;
    charm_tree_attr_t  *attr;
    double              fr, fu, fv, fw, fe;
    double              gr, gu, gv, gw, ge;
    double              hr, hu, hv, hw, he;
    double              phi_x, phi_y, phi_z;
    double             *x;

    attr = charm_get_tree_attr(info->p4est, info->treeid);

    CHARM_ARR_SET_ZERO(data->int_ro);
    CHARM_ARR_SET_ZERO(data->int_ru);
    CHARM_ARR_SET_ZERO(data->int_rv);
    CHARM_ARR_SET_ZERO(data->int_rw);
    CHARM_ARR_SET_ZERO(data->int_re);


    for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
        for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
            x = data->par.g.quad_gp[igp];

            charm_get_fields(q, x, &c);
            charm_param_cons_to_prim(attr->reg->mat, &p, &c);

            fr = c.ru;
            fu = fr*p.u+p.p;
            fv = fr*p.v;
            fw = fr*p.w;
            fe = fr*p.e_tot+p.p*p.u;

            gr = c.rv;
            gu = gr*p.u;
            gv = gr*p.v+p.p;
            gw = gr*p.w;
            ge = gr*p.e_tot+p.p*p.v;

            hr = c.rw;
            hu = hr*p.u;
            hv = hr*p.v;
            hw = hr*p.w+p.p;
            he = hr*p.e_tot+p.p*p.w;

            phi_x = charm_base_func_dx(x, ibf, q) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_y = charm_base_func_dy(x, ibf, q) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_z = charm_base_func_dz(x, ibf, q) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];

            data->int_ro[ibf] += (fr*phi_x+gr*phi_y+hr*phi_z);
            data->int_ru[ibf] += (fu*phi_x+gu*phi_y+hu*phi_z);
            data->int_rv[ibf] += (fv*phi_x+gv*phi_y+hv*phi_z);
            data->int_rw[ibf] += (fw*phi_x+gw*phi_y+hw*phi_z);
            data->int_re[ibf] += (fe*phi_x+ge*phi_y+he*phi_z);
        }
    }
}


/*
 * Surface integrals
 */


static void _charm_convect_surface_int_iter_bnd (p4est_iter_face_info_t * info, void *user_data)
{
    int                 i, ibf, igp;
    p4est_t            *p4est = info->p4est;
    charm_data_t       *ghost_data = (charm_data_t *) user_data;
    charm_data_t       *udata[2];
    double              n[3];
    double              qr, qu, qv, qw, qe;
    double              bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);
    charm_prim_t        prim1, prim2;

    int bnd = 0;
    int8_t face[2];
    double c[2][3], l[3];
    double r_[2], p_[2], u_[2], v_[2], w_[2], e_[2];
    charm_tree_attr_t *attr;
    charm_cons_t        cons[2];
    charm_prim_t        prim[2];
    double             *x;

    P4EST_ASSERT(info->tree_boundary);

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    P4EST_ASSERT(!side[0]->is_hanging);

    if (side[0]->is.full.is_ghost) {
        udata[0] = &(ghost_data[side[0]->is.full.quadid]);
    }
    else {
        udata[0] = (charm_data_t *) side[0]->is.full.quad->p.user_data;
    }
    face[0] = side[0]->face;
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

    udata[1] = P4EST_ALLOC(charm_data_t, 1);

    for (igp = 0; igp < CHARM_FACE_GP_COUNT; igp++) {

        x = udata[0]->par.g.face_gp[face[0]][igp];

        attr = charm_get_tree_attr(p4est, side[0]->treeid);
        charm_get_fields(side[0]->is.full.quad, x, &(cons[0]));
        charm_param_cons_to_prim(attr->reg->mat, &(prim[0]), &(cons[0]));


        charm_bnd_cond(p4est, side[0]->treeid, face[0], &(prim[0]), &(prim[1]), n);


        for (i = 0; i < 2; i++) {
            r_[i] = prim[i].r;
            u_[i] = prim[i].u;
            v_[i] = prim[i].v;
            w_[i] = prim[i].w;
            p_[i] = prim[i].p;
        }

        /* flux from side 0 to side 1 */
        charm_calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, n);
        for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
            if (!side[0]->is.full.is_ghost) {
                bfv = -1. * charm_base_func(x, ibf, side[0]->is.full.quad) * udata[0]->par.g.face_gw[face[0]][igp] *
                      udata[0]->par.g.face_gj[face[0]][igp];
                udata[0]->int_ro[ibf] += qr * bfv;
                udata[0]->int_ru[ibf] += qu * bfv;
                udata[0]->int_rv[ibf] += qv * bfv;
                udata[0]->int_rw[ibf] += qw * bfv;
                udata[0]->int_re[ibf] += qe * bfv;
            }
        }
    }
    P4EST_FREE(udata[1]);

}

static void _charm_convect_surface_int_iter_inner (p4est_iter_face_info_t * info, void *user_data)
{
    int                 i, j, h_side, igp, ibf;
    p4est_t            *p4est = info->p4est;
    charm_data_t       *ghost_data = (charm_data_t *) user_data;
    charm_data_t       *udata[2];
    double              n[3];
    double              qr, qu, qv, qw, qe;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);
    charm_cons_t        cons[2];
    charm_prim_t        prim[2];
    charm_tree_attr_t  *attr;
    double             *x;
    double              bfv;

    int8_t face[2];
    double c[2][3], l[3];
    double r_[2], p_[2], u_[2], v_[2], w_[2], e_[2];


    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    face[0] = side[0]->face;
    face[1] = side[1]->face;

    h_side = -1;
    if (side[0]->is_hanging || side[1]->is_hanging) { // @todo
        P4EST_ASSERT(0);
//        for (j = 0; j < P4EST_HALF; j++) {
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
//            P4EST_ASSERT(h_side != -1);
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
                udata[i] = (charm_data_t *) side[i]->is.full.quad->p.user_data;
            }
        }
        //facearea = charm_face_get_normal(udata[0], face[0], n);
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

        for (igp = 0; igp < CHARM_BASE_FN_COUNT; igp++) {
            x = udata[0]->par.g.face_gp[face[i]][igp];
            for (i = 0; i < 2; i++) {

                attr = charm_get_tree_attr(p4est, side[i]->treeid);
                charm_get_fields(side[i]->is.full.quad, x, &(cons[i]));
                charm_param_cons_to_prim(attr->reg->mat, &(prim[i]), &(cons[i]));
                r_[i] = prim[i].r;
                u_[i] = prim[i].u;
                v_[i] = prim[i].v;
                w_[i] = prim[i].w;
                p_[i] = prim[i].p;
            }

            /* flux from side 0 to side 1 */
            charm_calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, n);

            for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                for (i = 0; i < 2; i++) {
                    if (!side[i]->is.full.is_ghost) {
                        bfv = (i ? 1. : -1.) * charm_base_func(x, ibf, side[i]->is.full.quad)
                              * udata[i]->par.g.face_gw[face[0]][igp]
                              * udata[i]->par.g.face_gj[face[0]][igp];
                        udata[i]->int_ro[ibf] += qr * bfv;
                        udata[i]->int_ru[ibf] += qu * bfv;
                        udata[i]->int_rv[ibf] += qv * bfv;
                        udata[i]->int_rw[ibf] += qw * bfv;
                        udata[i]->int_re[ibf] += qe * bfv;
                    }
                }
            }
        }

    }

}

void charm_convect_surface_int_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{return;
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_convect_surface_int_iter_bnd(info, user_data);
    }
    else {
        _charm_convect_surface_int_iter_inner(info, user_data);
    }

}
