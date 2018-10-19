//
// Created by appmath on 14.10.18.
//

#include "charm_globals.h"
#include "charm_base_func.h"

/**
 *
 * @param info
 * @param user_data
 */
static void _charm_grad_p_volume_int_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 ibf, igp;
    double              press;
    double              phi_x, phi_y, phi_z;
    double             *x;

    for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
        for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
            x = data->par.g.quad_gp[igp];

            press = charm_get_field_p(data, x);

            phi_x = charm_base_func_dx(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_y = charm_base_func_dy(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_z = charm_base_func_dz(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];

            data->int_ru[ibf] -= (press*phi_x);
            data->int_rv[ibf] -= (press*phi_y);
            data->int_rw[ibf] -= (press*phi_z);
        }
    }
}


/*
 * Surface integrals
 */

/**
 *
 * @param info
 * @param user_data
 */
static void _charm_grad_p_surface_int_iter_bnd (p4est_iter_face_info_t * info, void *user_data) {
    int i, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_ctx_t * ctx = charm_get_ctx(p4est);
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    double n[3];
    double p[2];
    double qu, qv, qw;
    double bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);
//    size_t              c_count = charm_get_comp_count(info->p4est);

    int8_t face;
    double c[2][3], l[3];
    double *x, gw, gj;

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

    for (igp = 0; igp < CHARM_FACE_GP_COUNT; igp++) {
        x = udata->par.g.face_gp[face][igp];
        gw = udata->par.g.face_gw[face][igp];
        gj = udata->par.g.face_gj[face][igp];
        p[0] = p[1] = charm_get_field_p(udata, x);
        qu = 0.5*n[0]*(p[0]+p[1]);
        qv = 0.5*n[1]*(p[0]+p[1]);
        qw = 0.5*n[2]*(p[0]+p[1]);
        for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
            if (!side[0]->is.full.is_ghost) {
                bfv = charm_base_func(x, ibf, udata) * gw * gj;
                udata->int_ru[ibf] += qu * bfv;
                udata->int_rv[ibf] += qv * bfv;
                udata->int_rw[ibf] += qw * bfv;
            }
        }
    }
}


/**
 *
 * @param info
 * @param user_data
 */
static void _charm_grad_p_surface_int_iter_inner (p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side, igp, ibf,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    double                  n[3];
    double                  qu, qv, qw;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    double                  p[2];
    double                 *x, gw, gj;
    double                  bfv;
    double                  c[2][3];
    double                  l[3];
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
                for (i = 0; i < 2; i++) {
                    p[i] = charm_get_field_p(udata[i], x);
                }
                qu = 0.5*n[0]*(p[0]+p[1]);
                qv = 0.5*n[1]*(p[0]+p[1]);
                qw = 0.5*n[2]*(p[0]+p[1]);
                for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                    for (i = 0; i < 2; i++) {
                        if (!side[i]->is.full.is_ghost) {
                            bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                            udata[i]->int_ru[ibf] += qu * bfv;
                            udata[i]->int_rv[ibf] += qv * bfv;
                            udata[i]->int_rw[ibf] += qw * bfv;
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
            for (i = 0; i < 2; i++) {
                p[i] = charm_get_field_p(udata[i], x);
            }
            qu = 0.5*n[0]*(p[0]+p[1]);
            qv = 0.5*n[1]*(p[0]+p[1]);
            qw = 0.5*n[2]*(p[0]+p[1]);
            for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                for (i = 0; i < 2; i++) {
                    if (!side[i]->is.full.is_ghost) {
                        bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                        udata[i]->int_ru[ibf] += qu * bfv;
                        udata[i]->int_rv[ibf] += qv * bfv;
                        udata[i]->int_rw[ibf] += qw * bfv;
                    }
                }
            }
        }
    }
}


static void _charm_grad_p_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->int_ru[i] = 0.;
        data->int_rv[i] = 0.;
        data->int_rw[i] = 0.;
    }
}


/**
 *
 * @param info
 * @param user_data
 */
static void _charm_grad_p_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i, j;

    charm_matr_vect_mult(data->par.g.a_inv, data->int_ru, data->par.c.grad_p[0]);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rv, data->par.c.grad_p[1]);
    charm_matr_vect_mult(data->par.g.a_inv, data->int_rw, data->par.c.grad_p[2]);
    for (i = 0; i < CHARM_DIM; i++) {
        for (j = 0; j < CHARM_BASE_FN_COUNT; j++) {
            if (fabs(data->par.c.grad_p[i][j]) < CHARM_EPS) data->par.c.grad_p[i][j] = 0.;
        }
    }
}


static void _charm_grad_p_surface_int_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_grad_p_surface_int_iter_bnd(info, user_data);
    }
    else {
        _charm_grad_p_surface_int_iter_inner(info, user_data);
    }

}


/**
 *
 * @param p4est
 * @param ghost
 * @param ghost_data
 */
void charm_timestep_press_grad(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data)
{
    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_grad_p_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_grad_p_volume_int_iter_fn,
                   _charm_grad_p_surface_int_iter_fn,
                   NULL, NULL);
    p4est_iterate (p4est, NULL, NULL,
                   _charm_grad_p_update_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

}









