//
// Created by zhrv on 27.10.17.
//

#include "charm_grad.h"

//static void charm_grad_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
//{
//    int                 i, j;
//    p4est_t            *p4est = info->p4est;
//    charm_ctx_t          *ctx = (charm_ctx_t *) p4est->user_pointer;
//    charm_data_t         *ghost_data = (charm_data_t *) user_data;
//    charm_data_t         *udata;
//    p4est_quadrant_t   *quad;
//    double              vdotn = 0.;
//    double              ro[2];
//    double              ru[2];
//    double              rv[2];
//    double              rw[2];
//    double              re[2];
//    double              qr, qu, qv, qw, qe;
//    double              h, facearea;
//    int                 which_face;
//    p4est_iter_face_side_t *side[2];
//    sc_array_t         *sides = &(info->sides);
//
//    if (sides->elem_count != 2) {
//        side[0] = p4est_iter_fside_array_index_int(sides, 0);
//
//        P4EST_ASSERT(info->tree_boundary);
//
//        /* which of the quadrant's faces the interface touches */
//        which_face = side[0]->face;
//
//        ro[0] = 0;
//        ru[0] = 0;
//        rv[0] = 0;
//        rw[0] = 0;
//        re[0] = 0;
//        if (side[0]->is_hanging) {
//            /* there are 2^(d-1) (P4EST_HALF) subfaces */
//            for (j = 0; j < P4EST_HALF; j++) {
//                if (side[0]->is.hanging.is_ghost[j]) {
//                    udata =
//                            (charm_data_t *) &ghost_data[side[0]->is.hanging.quadid[j]];
//                }
//                else {
//                    udata =
//                            (charm_data_t *) side[0]->is.hanging.quad[j]->p.user_data;
//                }
//                ro[0] += udata->par.c.ro;
//                ru[0] += udata->par.c.ru;
//                rv[0] += udata->par.c.rv;
//                rw[0] += udata->par.c.rw;
//                re[0] += udata->par.c.re;
//            }
//            ro[0] /= P4EST_HALF;
//            ru[0] /= P4EST_HALF;
//            rv[0] /= P4EST_HALF;
//            rw[0] /= P4EST_HALF;
//            re[0] /= P4EST_HALF;
//        }
//        else {
//            if (side[0]->is.full.is_ghost) {
//                udata = (charm_data_t *) &ghost_data[side[0]->is.full.quadid];
//            }
//            else {
//                udata = (charm_data_t *) side[0]->is.full.quad->p.user_data;
//            }
//            ro[0] = udata->par.c.ro;
//            ru[0] = udata->par.c.ru;
//            rv[0] = udata->par.c.rv;
//            rw[0] = udata->par.c.rw;
//            re[0] = udata->par.c.re;
//        }
//
//        /* boundary conditions */
//        double r_, p_, u_, v_;
//        switch (which_face) {
//            case 0:                      /* -x side */
//                vdotn = -1.0;
//
//                ro[1] =  ro[0];
//                ru[1] = -ru[0];
//                rv[1] =  rv[0];
//                rw[1] =  rw[0];
//                re[1] =  re[0];
//                break;
//
//            case 1:                      /* +x side */
//                vdotn = 1.0;
//
//                ro[1] =  ro[0];
//                ru[1] = -ru[0];
//                rv[1] =  rv[0];
//                rw[1] =  rw[0];
//                re[1] =  re[0];
//                break;
//
//            case 2:                      /* -y side */
//                vdotn = -1.0;
//
//                ro[1] =  ro[0];
//                ru[1] =  ru[0];
//                rv[1] = -rv[0];
//                rw[1] =  rw[0];
//                re[1] =  re[0];
//                break;
//
//            case 3:                      /* +y side */
//                vdotn = 1.0;
//
//                ro[1] =  ro[0];
//                ru[1] =  ru[0];
//                rv[1] = -rv[0];
//                rw[1] =  rw[0];
//                re[1] =  re[0];
//                break;
//            case 4:                      /* -z side */
//                vdotn = -1.0;
////
////                r_ = 12.09;
////                u_ = 0.0;
////                v_ = 97.76;
////                p_ = 2.152e+5;
////
////                ro[1] = r_;
////                ru[1] = r_*u_;
////                rv[1] = r_*v_;
////                re[1] = p_/(GAM-1.0)+0.5*r_*(u_*u_+v_*v_);
////
//                ro[1] =  ro[0];
//                ru[1] =  ru[0];
//                rv[1] =  rv[0];
//                rw[1] =  rw[0];
//                re[1] =  re[0];
//                break;
//
//            case 5:                      /* +z side */
//                vdotn = 1.0;
//
//                ro[1] =  ro[0];
//                ru[1] =  ru[0];
//                rv[1] =  rv[0];
//                rw[1] = -rw[0];
//                re[1] =  re[0];
//                break;
//        }
//    }
//    else {
//
//        side[0] = p4est_iter_fside_array_index_int(sides, 0);
//        side[1] = p4est_iter_fside_array_index_int(sides, 1);
//
//        /* which of the quadrant's faces the interface touches */
//        which_face = side[0]->face;
//
//        switch (which_face) {
//            case 0:                      /* -x side */
//                vdotn = -1.0;
//                break;
//            case 1:                      /* +x side */
//                vdotn = 1.0;
//                break;
//            case 2:                      /* -y side */
//                vdotn = -1.0;
//                break;
//            case 3:                      /* +y side */
//                vdotn = 1.0;
//                break;
//            case 4:                      /* -z side */
//                vdotn = -1.0;
//                break;
//            case 5:                      /* +z side */
//                vdotn = 1.0;
//                break;
//        }
//
////        P4EST_ASSERT (vdotn == 1.0);
//
//        for (i = 0; i < 2; i++) {
//            ro[i] = 0;
//            ru[i] = 0;
//            rv[i] = 0;
//            rw[i] = 0;
//            re[i] = 0;
//            if (side[i]->is_hanging) {
//                /* there are 2^(d-1) (P4EST_HALF) subfaces */
//                for (j = 0; j < P4EST_HALF; j++) {
//                    if (side[i]->is.hanging.is_ghost[j]) {
//                        udata =
//                                (charm_data_t *) &ghost_data[side[i]->is.hanging.quadid[j]];
//                    }
//                    else {
//                        udata =
//                                (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
//                    }
//                    ro[i] += udata->par.c.ro;
//                    ru[i] += udata->par.c.ru;
//                    rv[i] += udata->par.c.rv;
//                    rw[i] += udata->par.c.rw;
//                    re[i] += udata->par.c.re;
//                }
//                ro[i] /= P4EST_HALF;
//                ru[i] /= P4EST_HALF;
//                rv[i] /= P4EST_HALF;
//                rw[i] /= P4EST_HALF;
//                re[i] /= P4EST_HALF;
//            }
//            else {
//                if (side[i]->is.full.is_ghost) {
//                    udata = (charm_data_t *) &ghost_data[side[i]->is.full.quadid];
//                }
//                else {
//                    udata = (charm_data_t *) side[i]->is.full.quad->p.user_data;
//                }
//                ro[i] = udata->par.c.ro;
//                ru[i] = udata->par.c.ru;
//                rv[i] = udata->par.c.rv;
//                rw[i] = udata->par.c.rw;
//                re[i] = udata->par.c.re;
//            }
//        }
//
//    }
//
//
//    /* flux from side 0 to side 1 */
//    qr = sqrt(ro[0]*ro[1]);
//
//    double fSB = sqrt( ro[0] );
//    double fSE = sqrt( ro[1] );
//    double fS_ = 1.0 / ( fSB + fSE );
//
//    //qr = fSB * fSE;
//
//    double UB = ru[0]/ro[0];
//    double UE = ru[1]/ro[1];
//    double VB = rv[0]/ro[0];
//    double VE = rv[1]/ro[1];
//    double WB = rw[0]/ro[0];
//    double WE = rw[1]/ro[1];
//
//    double UI = ( fSB * UB + fSE * UE ) * fS_;
//    double VI = ( fSB * UB + fSE * VE ) * fS_;
//    double WI = ( fSB * WB + fSE * WE ) * fS_;
//
//    //double PB = e*ro[0]*
//
//    double UB_MAG = UB*UB+VB*VB+WB*WB;
//    double UE_MAG = UE*UE+VE*VE+WE*WE;
//    double EB = re[0]/ro[0]-0.5*UB_MAG;
//    double EE = re[1]/ro[1]-0.5*UE_MAG;
//    double PB = ro[0]*EB*(GAM-1.0);
//    double PE = ro[1]*EE*(GAM-1.0);
//
//
//    double HB = EB + UB_MAG*0.5 + PB / ro[0];
//    double HE = EE + UE_MAG*0.5 + PE / ro[1];
//
//    double HI = ( fSB * HB + fSE * HE ) * fS_;
//
//    double PI = ( HI - (UI*UI+VI*VI+WI*WI)*0.5 ) * qr * ( GAM - 1.0 ) / GAM;
//    double EI = PI/(qr*(GAM-1.0));
//
//    qu = UI*qr;
//    qv = VI*qr;
//    qw = WI*qr;
//    qe = qr*(EI+0.5*(UI*UI+VI*VI+WI*WI));
//
//    for (i = 0; i < sides->elem_count; i++) {
//        if (side[i]->is_hanging) {
//            /* there are 2^(d-1) (P4EST_HALF) subfaces */
//            for (j = 0; j < P4EST_HALF; j++) {
//                quad = side[i]->is.hanging.quad[j];
//                h = CHARM_GET_H(quad->level);
//                facearea = h * h;
//                if (!side[i]->is.hanging.is_ghost[j]) {
//                    udata = (charm_data_t *) quad->p.user_data;
//                    udata->dro[which_face/2] += qr * facearea * (i ? 1. : -1.) * vdotn;
//                    udata->dru[which_face/2] += qu * facearea * (i ? 1. : -1.) * vdotn;
//                    udata->drv[which_face/2] += qv * facearea * (i ? 1. : -1.) * vdotn;
//                    udata->drw[which_face/2] += qw * facearea * (i ? 1. : -1.) * vdotn;
//                    udata->dre[which_face/2] += qe * facearea * (i ? 1. : -1.) * vdotn;
//                }
//            }
//        }
//        else {
//            quad = side[i]->is.full.quad;
//            h = CHARM_GET_H(quad->level);
//            facearea = h * h;
//            if (!side[i]->is.full.is_ghost) {
//                udata = (charm_data_t *) quad->p.user_data;
//                udata->dro[which_face/2] += qr * facearea * (i ? 1. : -1.) * vdotn;
//                udata->dru[which_face/2] += qu * facearea * (i ? 1. : -1.) * vdotn;
//                udata->drv[which_face/2] += qv * facearea * (i ? 1. : -1.) * vdotn;
//                udata->drw[which_face/2] += qw * facearea * (i ? 1. : -1.) * vdotn;
//                udata->dre[which_face/2] += qe * facearea * (i ? 1. : -1.) * vdotn;
//            }
//        }
//    }
//}
//
//static void
//charm_grad_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
//{
//    p4est_quadrant_t   *q = info->quad;
//    charm_data_t       *data = (charm_data_t *) q->p.user_data;
//    int i;
//
//    double h = CHARM_GET_H(q->level);
//    double vol = h*h*h;
//
//    for (i = 0; i < P4EST_DIM; i++) {
//        data->dro[i] /= vol;
//        data->dru[i] /= vol;
//        data->drv[i] /= vol;
//        data->drw[i] /= vol;
//        data->dre[i] /= vol;
//    }
//}
//
//
//void charm_calc_grad(p4est_t * p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data)
//{
//    int my_ghost = 0;
//    if (!ghost) {
//        ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
//        ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
//        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
//        my_ghost = 1;
//    }
//
//    p4est_iterate (p4est, ghost, (void *) ghost_data,     /* pass in ghost data that we just exchanged */
//                   charm_reset_derivatives,       /* blank the previously calculated derivatives */
//                   charm_grad_face_iter_fn, /* compute the minmod estimate of each cell's derivative */
//                   NULL,
//                   NULL);         /* there is no callback for the corners between quadrants */
//
//
//    p4est_iterate (p4est, ghost, (void *) ghost_data,     /* pass in ghost data that we just exchanged */
//                   charm_grad_quad_iter_fn,       /* blank the previously calculated derivatives */
//                   NULL,
//                   NULL,
//                   NULL);
//
//    if (ghost && my_ghost) {
//        p4est_ghost_destroy (ghost);
//        P4EST_FREE (ghost_data);
//        ghost      = NULL;
//        ghost_data = NULL;
//    }
//}
//
