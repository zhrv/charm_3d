/**
 *  FVM for du/dt+V1*du/dx+V2*du/dy+V3*du/dz = 0
 *  u = ru[0]
*/

#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"

static const double V[3] = {100., 0., 0.};


void charm_adv_adapt(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data);


/** Compute the timestep.
 *
 */
charm_real_t charm_model_adv_get_dt (p4est_t * p4est)
{
    return charm_get_ctx(p4est)->dt;
}

/*
 * Surface integrals
 */
static void _charm_adv_surface_int_iter_bnd (p4est_iter_face_info_t * info, void *user_data) {
    int i, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_ctx_t * ctx = charm_get_ctx(p4est);
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    charm_real_t n[3], s;
    charm_real_t qu;
    charm_real_t bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);
    size_t              c_count = charm_get_comp_count(info->p4est);


    int8_t face;
    charm_real_t c[2][3], l[3];
    charm_real_t r_[2], p_[2], u_[2], v_[2], w_[2], e_[2];
    charm_cons_t cons;
    charm_prim_t prim[2];
    charm_real_t *x, gw, gj;
    charm_real_t intg[2][5];
    int j;


    CHARM_ASSERT(info->tree_boundary);

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    CHARM_ASSERT(!side[0]->is_hanging);
    CHARM_ASSERT(!side[0]->is.full.is_ghost);

    udata = charm_get_quad_data(side[0]->is.full.quad);

    face = side[0]->face;
    charm_face_get_normal(udata, face, n);
    charm_quad_get_center(udata, c[0]);
    charm_face_get_center(udata, face, c[1]);
    s = charm_face_get_area(udata, face);

    for (i = 0; i < 3; i++) {
        l[i] = c[1][i] - c[0][i];
    }

    if (scalar_prod(n, l) < 0) {
        for (i = 0; i < 3; i++) {
            n[i] *= -1.0;
        }
    }

    charm_real_t u = udata->par.c.rc[0][0];
    qu = u*(V[0]*n[0]+V[1]*n[1]+V[2]*n[2]);
    udata->int_rc[0][0] += qu * s;
}


static void _charm_adv_surface_int_iter_inner (p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side, igp, ibf,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    charm_real_t            n[3];
    charm_real_t            u, qu, bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
    charm_real_t            s;
    charm_real_t            c[2][3];
    charm_real_t            l[3];
    charm_real_t            vn;
    int8_t                  face[2];
    size_t                  c_count = charm_get_comp_count(info->p4est);


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
            s = charm_face_get_area(udata[h_side], face[h_side]);

            for (i = 0; i < 3; i++) {
                l[i] = c[1][i]-c[0][i];
            }

            if (scalar_prod(n, l) < 0) {
                for (i = 0; i < 3; i++) {
                    n[i] *= -1.0;
                }
            }

            vn = V[0]*n[0]+V[1]*n[1]+V[2]*n[2];

            if (vn > 0) {
                u = udata[0]->par.c.rc[0][0];
            }
            else {
                u = udata[1]->par.c.rc[0][0];
            }

            qu = u*vn;

            for (i = 0; i < 2; i++) {
                if (i == h_side) {
                    if (!side[i]->is.hanging.is_ghost[j]) {
                        bfv = (i ? -1. : 1.)  * s;
                        udata[i]->int_rc[0][0] += qu * bfv;
                    }
                }
                else {
                    if (!side[i]->is.full.is_ghost) {
                        bfv = (i ? -1. : 1.) * s;
                        udata[i]->int_rc[0][0] += qu * bfv;
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
        s = charm_face_get_area(udata[0], face[0]);

        for (i = 0; i < 3; i++) {
            l[i] = c[1][i]-c[0][i];
        }

        if (scalar_prod(n, l) < 0) {
            for (i = 0; i < 3; i++) {
                n[i] *= -1.0;
            }
        }

        vn = V[0]*n[0]+V[1]*n[1]+V[2]*n[2];

        if (vn > 0) {
            qu = udata[0]->par.c.rc[0][0] * vn;
        }
        else {
            qu = udata[1]->par.c.rc[0][0] * vn;
        }

        for (i = 0; i < 2; i++) {
            if (!side[i]->is.full.is_ghost) {
                bfv = (i ? -1. : 1.) * s;
                udata[i]->int_rc[0][0] += qu * bfv;
            }
        }
    }
}

static void _charm_adv_surface_int_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_adv_surface_int_iter_bnd(info, user_data);
    }
    else {
        _charm_adv_surface_int_iter_inner(info, user_data);
    }

}



static void _charm_adv_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_real_t        dt = *((charm_real_t *) user_data);
    charm_real_t        v = charm_quad_get_volume(data);

    data->par.c.rc[0][0] -= dt * data->int_rc[0][0] / v;
    data->par.c.rc[1][0] = 1. - data->par.c.rc[0][0];
}


static void _charm_adv_timestep_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i, j;

    data->int_rc[0][0] = 0.;
}



void charm_model_adv_timestep_single(p4est_t * p4est, charm_real_t *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    int                 refine_period = ctx->refine_period;
    int                 repartition_period = ctx->repartition_period;
    int                 write_period = ctx->write_period;
    int                 allowcoarsening = 1;
    int                 mpiret;
    charm_real_t              orig_max_err = ctx->max_err;
    charm_real_t              umax, global_umax;
    int                 ref_flag = 0;
    p4est_ghost_t      *ghost       = *_ghost;
    charm_data_t       *ghost_data  = *_ghost_data;

    /* refine */
    ref_flag = 0;
    if (refine_period) {
        if (!(ctx->timestep % refine_period)) {
            if (ctx->timestep) {
                ctx->amr_fn(p4est, ghost, ghost_data); /* adapt */
                if (ghost) {
                    p4est_ghost_destroy(ghost);
                    CHARM_FREE (ghost_data);
                    ghost = NULL;
                    ghost_data = NULL;
                }

                ref_flag = 1;
            }
            *dt = ctx->get_dt_fn(p4est);

        }
    }
    else {
        *dt = ctx->get_dt_fn(p4est);
    }

    /* repartition */
    if (repartition_period) {
        if (ctx->timestep && !(ctx->timestep % repartition_period)) {

            p4est_partition(p4est, allowcoarsening, NULL);

            if (ghost) {
                p4est_ghost_destroy(ghost);
                CHARM_FREE (ghost_data);
                ghost = NULL;
                ghost_data = NULL;
            }
        }
    }

    /* write out solution */
    if (!(ctx->timestep % write_period)) {
        charm_write_solution (p4est);
        CHARM_GLOBAL_ESSENTIALF (" File for step #%d is saved \n", ctx->timestep);
    }

    /* synchronize the ghost data */
    if (!ghost) {
        ghost = p4est_ghost_new (p4est, CHARM_CONNECT_FULL);
        ghost_data = CHARM_ALLOC (charm_data_t, ghost->ghosts.elem_count);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_adv_timestep_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   NULL,
                   _charm_adv_surface_int_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est, NULL,
                   (void *) dt,
                   _charm_adv_timestep_update_quad_iter_fn,
                   NULL, NULL, NULL);


    p4est_ghost_exchange_data (p4est, ghost, ghost_data); /* synchronize the ghost data */

    *_ghost       = ghost;
    *_ghost_data  = ghost_data;

}


