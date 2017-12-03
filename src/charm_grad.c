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
    int                 i, j, k, h_side;
    p4est_t            *p4est = info->p4est;
    charm_data_t       *ghost_data = (charm_data_t *) user_data;
    charm_data_t       *udata[2];
    double              n[3];
    double              qr, qu, qv, qw, qp;
    double              facearea;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);

    int8_t face[2];
    double c[2][3], l[3];
    double r_[2], p_[2], u_[2], v_[2], w_[2];



    if (sides->elem_count != 2) {
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
        facearea = charm_face_get_normal(udata[0], face[0], n);
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

        charm_bnd_cond(p4est, side[0]->treeid, face[0], &(udata[0]->par), &(udata[1]->par), n);

        charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[0]->treeid);
        charm_param_cons_to_prim(attr->reg->mat, &(udata[1]->par));


        for (i = 0; i < 2; i++) {
            attr = charm_get_tree_attr(p4est, side[0]->treeid);
            charm_param_cons_to_prim(attr->reg->mat, &(udata[i]->par));
            r_[i] = udata[i]->par.p.r;
            u_[i] = udata[i]->par.p.u;
            v_[i] = udata[i]->par.p.v;
            w_[i] = udata[i]->par.p.w;
            p_[i] = udata[i]->par.p.p;
        }


        qr = r_[1];
        qu = u_[1];
        qv = v_[1];
        qw = w_[1];
        qp = p_[1];

        if (!side[0]->is.full.is_ghost) {
            for (k = 0; k < CHARM_DIM; k++) {
                udata[0]->par.grad.r[k] += qr * n[k] * facearea;
                udata[0]->par.grad.u[k] += qu * n[k] * facearea;
                udata[0]->par.grad.v[k] += qv * n[k] * facearea;
                udata[0]->par.grad.w[k] += qw * n[k] * facearea;
                udata[0]->par.grad.p[k] += qp * n[k] * facearea;
            }
        }
        P4EST_FREE(udata[1]);
    }
    else {
        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        side[1] = p4est_iter_fside_array_index_int(sides, 1);
        face[0] = side[0]->face;
        face[1] = side[1]->face;

        h_side = -1;

        if (side[0]->is_hanging || side[1]->is_hanging) {
            for (j = 0; j < P4EST_HALF; j++) {
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

                P4EST_ASSERT(h_side != -1);

                facearea = charm_face_get_area(udata[h_side], side[h_side]->face);

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

                for (i = 0; i < 2; i++) {
                    charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[i]->treeid);
                    charm_param_cons_to_prim(attr->reg->mat, &(udata[i]->par));
                    r_[i] = udata[i]->par.p.r;
                    u_[i] = udata[i]->par.p.u;
                    v_[i] = udata[i]->par.p.v;
                    w_[i] = udata[i]->par.p.w;
                    p_[i] = udata[i]->par.p.p;
                }

                qr = 0.5*(r_[0]+r_[1]);
                qu = 0.5*(u_[0]+u_[1]);
                qv = 0.5*(v_[0]+v_[1]);
                qw = 0.5*(w_[0]+w_[1]);
                qp = 0.5*(p_[0]+p_[1]);

                for (i = 0; i < 2; i++) {
                    if (side[i]->is_hanging) {
                        if (side[i]->is.hanging.is_ghost[j]) {
                            continue;
                        }

                    }
                    else {
                        if (side[i]->is.full.is_ghost) {
                            continue;
                        }
                    }
                    for (k = 0; k < CHARM_DIM; k++) {
                        udata[i]->par.grad.r[k] += qr * n[k] * facearea * (i ? -1. : 1.);
                        udata[i]->par.grad.u[k] += qu * n[k] * facearea * (i ? -1. : 1.);
                        udata[i]->par.grad.v[k] += qv * n[k] * facearea * (i ? -1. : 1.);
                        udata[i]->par.grad.w[k] += qw * n[k] * facearea * (i ? -1. : 1.);
                        udata[i]->par.grad.p[k] += qp * n[k] * facearea * (i ? -1. : 1.);
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
                    udata[i] = (charm_data_t *) side[i]->is.full.quad->p.user_data;
                }
            }
            facearea = charm_face_get_normal(udata[0], face[0], n);
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

            for (i = 0; i < 2; i++) {
                charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[0]->treeid);
                charm_param_cons_to_prim(attr->reg->mat, &(udata[i]->par));
                r_[i] = udata[i]->par.p.r;
                u_[i] = udata[i]->par.p.u;
                v_[i] = udata[i]->par.p.v;
                w_[i] = udata[i]->par.p.w;
                p_[i] = udata[i]->par.p.p;
            }

            qr = 0.5*(r_[0]+r_[1]);
            qu = 0.5*(u_[0]+u_[1]);
            qv = 0.5*(v_[0]+v_[1]);
            qw = 0.5*(w_[0]+w_[1]);
            qp = 0.5*(p_[0]+p_[1]);

            for (i = 0; i < 2; i++) {
                if (!side[i]->is.full.is_ghost) {
                    //facearea = charm_face_get_area(udata[i], side[i]->face);
                    for (k = 0; k < CHARM_DIM; k++) {
                        udata[i]->par.grad.r[k] += qr * n[k] * facearea * (i ? -1. : 1.);
                        udata[i]->par.grad.u[k] += qu * n[k] * facearea * (i ? -1. : 1.);
                        udata[i]->par.grad.v[k] += qv * n[k] * facearea * (i ? -1. : 1.);
                        udata[i]->par.grad.w[k] += qw * n[k] * facearea * (i ? -1. : 1.);
                        udata[i]->par.grad.p[k] += qp * n[k] * facearea * (i ? -1. : 1.);
                    }
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
    int my_ghost = 0;
    if (!ghost) {
        ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
        ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
        my_ghost = 1;
    }

    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   charm_grad_reset_quad_iter_fn,
                   charm_grad_face_iter_fn,
                   NULL,
                   NULL);


    p4est_iterate (p4est, ghost, (void *) ghost_data,     /* pass in ghost data that we just exchanged */
                   charm_grad_quad_iter_fn,       /* blank the previously calculated derivatives */
                   NULL,
                   NULL,
                   NULL);

    p4est_ghost_exchange_data(p4est, ghost, ghost_data);

    if (ghost && my_ghost) {
        p4est_ghost_destroy (ghost);
        P4EST_FREE (ghost_data);
        ghost      = NULL;
        ghost_data = NULL;
    }

}

