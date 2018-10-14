//
// Created by appmath on 14.10.18.
//

#include "charm_globals.h"
#include "charm_base_func.h"
#include "charm_bnd_cond.h"


static void _charm_p_rhs_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->int_rh[i] = 0.;
       // data->par.c.p[i] = 0.;
    }
}

static void _charm_p_rhs_volume_int_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = charm_get_quad_data(q);
    int                 ibf, igp;
    charm_cons_t        c;
    double              fu;
    double              gu;
    double              hu;
    double              phi_x, phi_y, phi_z;
    double             *x;

    for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
        for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
            x = data->par.g.quad_gp[igp];

            charm_get_fields(data, x, &c);

            fu = c.ru;
            gu = c.rv;
            hu = c.rw;

            phi_x = charm_base_func_dx(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_y = charm_base_func_dy(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];
            phi_z = charm_base_func_dz(x, ibf, data) * data->par.g.quad_gj[igp] * data->par.g.quad_gw[igp];

            data->int_rh[ibf] -= (fu*phi_x+gu*phi_y+hu*phi_z);
        }
    }
}

static void _charm_p_rhs_surface_int_iter_bnd (p4est_iter_face_info_t * info, void *user_data)
{
    int i, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_ctx_t * ctx = charm_get_ctx(p4est);
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    double n[3];
    double qh;
    double bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);


    int8_t face;
    double c[2][3], l[3];
    charm_cons_t cons;
    charm_prim_t prim[2];
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
        charm_get_fields(udata, x, &cons);
        charm_param_cons_to_prim(p4est, &(prim[0]), &cons);
        charm_bnd_cond(p4est, side[0]->treeid, face, &(prim[0]), &(prim[1]), n);
        qh = cons.ru*n[0]+cons.rv*n[1]+cons.rw*n[2];
        for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
            if (!side[0]->is.full.is_ghost) {
                bfv = charm_base_func(x, ibf, udata) * gw * gj;
                udata->int_rh[ibf] += qh * bfv;
            }
        }
    }
}


/**
 *
 * @param info
 * @param user_data
 */
static void _charm_p_rhs_surface_int_iter_inner (p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side, igp, ibf,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    double                  n[3];
    double                  qh;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
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
                    charm_get_fields(udata[i], x, &(cons[i]));
                    charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                }
                qh = 0.5*((cons[0].ru+cons[1].ru)*n[0]+(cons[0].rv+cons[1].rv)*n[1]+(cons[0].rw+cons[1].rw)*n[2]);
                for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                    for (i = 0; i < 2; i++) {
                        if (!side[i]->is.full.is_ghost) {
                            bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                            udata[i]->int_rh[ibf] += qh * bfv;
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
                charm_get_fields(udata[i], x, &(cons[i]));
                charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
            }
            qh = 0.5*((cons[0].ru+cons[1].ru)*n[0]+(cons[0].rv+cons[1].rv)*n[1]+(cons[0].rw+cons[1].rw)*n[2]);
            for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                for (i = 0; i < 2; i++) {
                    if (!side[i]->is.full.is_ghost) {
                        bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                        udata[i]->int_rh[ibf] += qh * bfv;
                    }
                }
            }
        }
    }
}

static void _charm_p_rhs_surface_int_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_p_rhs_surface_int_iter_bnd(info, user_data);
    }
    else {
        _charm_p_rhs_surface_int_iter_inner(info, user_data);
    }

}
static void _charm_p_rhs_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    double              dt = *((double *) user_data);
    int                 i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->int_rh[i] /= dt;
    }
}




/**
 *
 * @param p4est
 * @param ghost
 * @param ghost_data
 * @param dt
 * @param stop
 */
void charm_timestep_press_rhs(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data, double *dt)
{
    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_p_rhs_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_p_rhs_volume_int_iter_fn,
                   _charm_p_rhs_surface_int_iter_fn,
                   NULL, NULL);
    p4est_iterate (p4est, NULL,
                   (void *) dt,
                   _charm_p_rhs_update_quad_iter_fn,
                   NULL, NULL, NULL);


    p4est_ghost_exchange_data (p4est, ghost, ghost_data);
}





