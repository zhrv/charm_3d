//
// Created by zhrv on 15.09.18.
//

#include "charm_timestep_diff.h"


static void _charm_diffusion_flux_face_iter_bnd (p4est_iter_face_info_t * info, void *user_data)
{
//    int                 i;
//    p4est_t            *p4est = info->p4est;
//    charm_data_t       *ghost_data = (charm_data_t *) user_data;
//    charm_data_t       *udata[2];
//    double              n[3];
//    double              qu, qv, qw, qe;
//    double              facearea;
//    p4est_iter_face_side_t *side[2];
//    sc_array_t         *sides = &(info->sides);
//
//    int8_t face[2];
//    double c[2][3], l[3];
//
//
//    // @todo is diff. flux a zero on boundary ?
//    P4EST_ASSERT(info->tree_boundary);
//
//    side[0] = p4est_iter_fside_array_index_int(sides, 0);
//    P4EST_ASSERT(!side[0]->is_hanging);
//    face[0] = side[0]->face;
//    charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[0]->treeid);
//    if (attr->bnd[face[0]]->type == BOUND_WALL_NO_SLIP) {
//        if (side[0]->is.full.is_ghost) {
//            udata[0] = &(ghost_data[side[0]->is.full.quadid]);
//        }
//        else {
//            udata[0] = (charm_data_t *) side[0]->is.full.quad->p.user_data;
//        }
//        double ml = attr->reg->mat->ml;
//        double lambda = 0.;
//        charm_param_t *p = &(udata[0]->par);
//        charm_param_cons_to_prim(attr->reg->mat, p);
//        facearea = charm_face_get_normal(udata[0], face[0], n);
//        charm_quad_get_center(udata[0], c[0]);
//        charm_face_get_center(udata[0], face[0], c[1]);
//
//        for (i = 0; i < 3; i++) {
//            l[i] = c[1][i] - c[0][i];
//        }
//
//        if (scalar_prod(n, l) < 0) {
//            for (i = 0; i < 3; i++) {
//                n[i] *= -1.0;
//            }
//        }
//        double vv[3] = {p->p.u, p->p.v, p->p.w};
//        double un = scalar_prod(vv, n);
//        double vn[3] = {un*n[0], un*n[1], un*n[2]};
//        double vt[3] = {vv[0]-vn[0], vv[1]-vn[1], vv[2]-vn[2]};
//        double ll = sqrt(scalar_prod(l,l));
//        qu = -ml*vt[0]/ll;
//        qv = -ml*vt[1]/ll;
//        qw = -ml*vt[2]/ll;
//        qe = 0.;
//        if (!side[0]->is.full.is_ghost) {
//            udata[0]->drudt += qu * facearea;
//            udata[0]->drvdt += qv * facearea;
//            udata[0]->drwdt += qw * facearea;
//            udata[0]->dredt += qe * facearea;
//        }
//    }
}

static void _charm_diffusion_flux_face_iter_inner (p4est_iter_face_info_t * info, void *user_data)
{
//    int                 i, j, h_side;
//    p4est_t            *p4est = info->p4est;
//    charm_data_t       *ghost_data = (charm_data_t *) user_data;
//    charm_data_t       *udata[2];
//    double              n[3];
//    double              qu, qv, qw, qe;
//    double              facearea;
//    p4est_iter_face_side_t *side[2];
//    sc_array_t         *sides = &(info->sides);
//
//    int8_t face[2];
//    double c[2][3], l[3];
//
//    double t_xx, t_yy, t_zz, t_xy, t_yz, t_xz;
//
//
//    side[0] = p4est_iter_fside_array_index_int(sides, 0);
//    side[1] = p4est_iter_fside_array_index_int(sides, 1);
//    face[0] = side[0]->face;
//    face[1] = side[1]->face;
//
//    h_side = -1;
//    if (side[0]->is_hanging || side[1]->is_hanging) {
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
//            qu = 0.;
//            qv = 0.;
//            qw = 0.;
//            qe = 0.;
//            for (i = 0; i < 2; i++) {
//                charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[i]->treeid);
//                double ml = attr->reg->mat->ml;
//                double lambda = 0.;
//                charm_param_t *p = &(udata[i]->par);
//                double div_v = p->grad.u[0] + p->grad.v[1] + p->grad.w[2];
//                charm_param_cons_to_prim(attr->reg->mat, p);
//
//                t_xx = (lambda-2.*ml/3.)*div_v+2.*ml*p->grad.u[0];
//                t_yy = (lambda-2.*ml/3.)*div_v+2.*ml*p->grad.v[1];
//                t_zz = (lambda-2.*ml/3.)*div_v+2.*ml*p->grad.w[2];
//                t_xy = ml*(p->grad.u[1] + p->grad.v[0]);
//                t_xz = ml*(p->grad.w[0] + p->grad.u[2]);
//                t_yz = ml*(p->grad.w[1] + p->grad.v[2]);
//
//                qu += t_xx*n[0]+t_xy*n[1]+t_xz*n[2];
//                qv += t_xy*n[0]+t_yy*n[1]+t_yz*n[2];
//                qw += t_xz*n[0]+t_yz*n[1]+t_zz*n[2];
//                qe +=   (t_xx*p->p.u+t_xy*p->p.v+t_xz*p->p.w)*n[0]+
//                        (t_xy*p->p.u+t_yy*p->p.v+t_yz*p->p.w)*n[1]+
//                        (t_xz*p->p.u+t_yz*p->p.v+t_zz*p->p.w)*n[2];
//            }
//
//            qu *= 0.5;
//            qv *= 0.5;
//            qw *= 0.5;
//            qe *= 0.5;
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
//                udata[i]->drudt -= qu * facearea * (i ? 1. : -1.);
//                udata[i]->drvdt -= qv * facearea * (i ? 1. : -1.);
//                udata[i]->drwdt -= qw * facearea * (i ? 1. : -1.);
//                udata[i]->dredt -= qe * facearea * (i ? 1. : -1.);
//
//            }
//
//        }
//    }
//    else {
//
//        for (i = 0; i < 2; i++) {
//            if (side[i]->is.full.is_ghost) {
//                udata[i] = &(ghost_data[side[i]->is.full.quadid]);
//            }
//            else {
//                udata[i] = (charm_data_t *) side[i]->is.full.quad->p.user_data;
//            }
//        }
//        facearea = charm_face_get_normal(udata[0], face[0], n);
//        charm_quad_get_center(udata[0], c[0]);
//        charm_face_get_center(udata[0], face[0], c[1]);
//
//        for (i = 0; i < 3; i++) {
//            l[i] = c[1][i]-c[0][i];
//        }
//
//        if (scalar_prod(n, l) < 0) {
//            for (i = 0; i < 3; i++) {
//                n[i] *= -1.0;
//            }
//        }
//
//        qu = 0.;
//        qv = 0.;
//        qw = 0.;
//        qe = 0.;
//        for (i = 0; i < 2; i++) {
//            charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[i]->treeid);
//            double ml = attr->reg->mat->ml;
//            double lambda = attr->reg->mat->lambda;
//            charm_param_t *p = &(udata[i]->par);
//            double div_v = p->grad.u[0] + p->grad.v[1] + p->grad.w[2];
//            charm_param_cons_to_prim(attr->reg->mat, p);
//            t_xx = (lambda-2.*ml/3.)*div_v+2.*ml*p->grad.u[0];
//            t_yy = (lambda-2.*ml/3.)*div_v+2.*ml*p->grad.v[1];
//            t_zz = (lambda-2.*ml/3.)*div_v+2.*ml*p->grad.w[2];
//            t_xy = ml*(p->grad.u[1] + p->grad.v[0]);
//            t_xz = ml*(p->grad.w[0] + p->grad.u[2]);
//            t_yz = ml*(p->grad.w[1] + p->grad.v[2]);
//
//            qu += t_xx*n[0]+t_xy*n[1]+t_xz*n[2];
//            qv += t_xy*n[0]+t_yy*n[1]+t_yz*n[2];
//            qw += t_xz*n[0]+t_yz*n[1]+t_zz*n[2];
//            qe += (t_xx*p->p.u+t_xy*p->p.v+t_xz*p->p.w)*n[0]+
//                  (t_xy*p->p.u+t_yy*p->p.v+t_yz*p->p.w)*n[1]+
//                  (t_xz*p->p.u+t_yz*p->p.v+t_zz*p->p.w)*n[2];
//        }
//
//        qu *= 0.5;
//        qv *= 0.5;
//        qw *= 0.5;
//        qe *= 0.5;
//
//        for (i = 0; i < 2; i++) {
//            if (side[i]->is.full.is_ghost) {
//                continue;
//            }
//            udata[i]->drudt -= qu * facearea * (i ? 1. : -1.);
//            udata[i]->drvdt -= qv * facearea * (i ? 1. : -1.);
//            udata[i]->drwdt -= qw * facearea * (i ? 1. : -1.);
//            udata[i]->dredt -= qe * facearea * (i ? 1. : -1.);
//        }
//
//    }
//
}

static void charm_diffusion_flux_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_diffusion_flux_face_iter_bnd(info, user_data);
    }
    else {
        _charm_diffusion_flux_face_iter_inner(info, user_data);
    }

}



