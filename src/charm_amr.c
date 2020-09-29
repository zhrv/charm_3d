//
// Created by zhrv on 26.10.17.
//

#include "charm_globals.h"
#include "charm_amr.h"
#include "charm_base_func.h"
#include "charm_bnd_cond.h"
#include "charm_fluxes.h"


static charm_real_t charm_error_sqr_estimate (p4est_quadrant_t * q)
{
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    int                 i;
    charm_vector_t      du;
    charm_real_t        vol = charm_quad_get_volume((charm_data_t *)q->p.user_data);
    charm_real_t        diff2 = 0.;

    for (i = 0; i < CHARM_DIM; i++) {
        du[i] = data->par.a.grad_u[i];
    }

    /* use the approximate derivative to estimate the L2 error */
    for (i = 0; i < CHARM_DIM; i++) {
        diff2 += du[i] * du[i];
    }

    return diff2 * (1. / 12.) * exp(log(vol)*5./3.);
}

static int charm_refine_err_estimate (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t              global_err = ctx->max_err;
    charm_real_t              global_err2 = global_err * global_err;
    charm_real_t              vol, err2;

    vol = charm_quad_get_volume((charm_data_t *)q->p.user_data);

    err2 = charm_error_sqr_estimate (q);
    if (err2 > global_err2 * vol) {
        return 1;
    }
    else {
        return 0;
    }
}

static int charm_refine_init_err_estimate (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t              global_err = ctx->max_err;
    charm_real_t              global_err2 = global_err * global_err;
    charm_data_t       *p = charm_get_quad_data(q);
    charm_real_t              vol = charm_quad_get_volume(p);
    charm_real_t              mp[3], err2;

    if (q->level >= ctx->max_level) {
        return 0;
    }

    charm_quad_get_center (q->p.user_data, mp);

#ifdef POGGI
    charm_real_t h = pow(vol, 1./3.);
    if ( ((-h < mp[2]) && (mp[2] < h)) || ((-0.00125 < mp[2]) && (mp[2] < 0.00125)) ) {
        return 1;
    }
    else {
        return 0;
    }

#else
    err2 = charm_error_sqr_estimate (q);
    if (err2 > global_err2 * vol) {
        return 1;
    }
    else {
        return 0;
    }
#endif
}

static int charm_coarsen_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * children[])
{
    return 0;
}

static int charm_coarsen_err_estimate (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * children[])
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_real_t              global_err = ctx->max_err;
    charm_real_t              global_err2 = global_err * global_err;
    charm_data_t       *data;
    charm_cons_t        cons;
    charm_prim_t        prim;
    charm_real_t              vol, err2, childerr2;
    charm_real_t              parentu;
    charm_real_t              diff;
    int                 i;


    if (children[0]->level <= ctx->min_level) {
        return 0;
    }


    /* the quadrant's volume is also its volume fraction */
    vol = charm_quad_get_volume((charm_data_t *)children[0]->p.user_data);

    /* compute the average */
    parentu = 0.;
    for (i = 0; i < CHARM_CHILDREN; i++) {
        data = (charm_data_t *) children[i]->p.user_data;
        charm_get_fields(data, data->par.g.c, &cons);
        charm_param_cons_to_prim(p4est, &prim, &cons);
        parentu += prim.r / CHARM_CHILDREN;
    }

    err2 = 0.;
    for (i = 0; i < CHARM_CHILDREN; i++) {
        data = (charm_data_t *) children[i]->p.user_data;
        childerr2 = charm_error_sqr_estimate (children[i]);
        charm_get_fields(data, data->par.g.c, &cons);
        charm_param_cons_to_prim(p4est, &prim, &cons);

        if (childerr2 > global_err2 * vol) {
            return 0;
        }
        err2 += childerr2;
        diff = (parentu - prim.r) * (parentu - prim.r);
        err2 += diff * vol;
    }
    if (err2 < global_err2 * (vol * CHARM_CHILDREN)) {
        return 1;
    }
    else {
        return 0;
    }

}

static int charm_refine_flag_estimate (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    charm_data_t       *data = (charm_data_t *) q->p.user_data;

    if (data->ref_flag >= 4) {
        return 1;
    }
    else {
        return 0;
    }
}

static void charm_ref_flag_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;

    data->ref_flag = 0;
}


static void charm_ref_flag_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    int                 i, j;
    charm_data_t         *ghost_data = (charm_data_t *) user_data;
    charm_data_t         *udata;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        side[0] = p4est_iter_fside_array_index_int(sides, 0);

        CHARM_ASSERT(info->tree_boundary);
        CHARM_ASSERT(!side[0]->is_hanging);
        CHARM_ASSERT(!side[0]->is.full.is_ghost);

        ((charm_data_t *) side[0]->is.full.quad->p.user_data)->ref_flag++;

    }
    else {

        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        side[1] = p4est_iter_fside_array_index_int(sides, 1);

        for (i = 0; i < 2; i++) {
            if (side[i]->is_hanging) {
                j = (i+1)%2;
                CHARM_ASSERT(!side[j]->is_hanging);
                if (side[j]->is.full.is_ghost) {
                    udata = &ghost_data[side[j]->is.full.quadid];
                }
                else {
                    udata = (charm_data_t *) side[j]->is.full.quad->p.user_data;
                }
                udata->ref_flag++;
            }
        }

    }
}


static void charm_replace_quads (p4est_t * p4est, p4est_topidx_t which_tree,
                                 int num_outgoing, p4est_quadrant_t * outgoing[],
                                 int num_incoming, p4est_quadrant_t * incoming[])
{
    const int           N = CHARM_BASE_FN_COUNT;
    size_t              c_count = charm_get_comp_count(p4est);
    charm_data_t       *parent_data, *child_data;
    int                 i, j, m, n, igp;
    charm_real_t              vol, svol;
    charm_real_t              ar[N][N], cr[N], rhs_ru[N], rhs_rv[N], rhs_rw[N], rhs_re[N], rhs_rc[CHARM_MAX_COMPONETS_COUNT][N];
    charm_real_t              *fld[5];

    if (num_outgoing > 1) {
        charm_geom_quad_calc(p4est, incoming[0], which_tree);
        /* this is coarsening */
        parent_data = (charm_data_t *) incoming[0]->p.user_data;
        child_data  = (charm_data_t *) outgoing[0]->p.user_data;
        parent_data->par.mat_id  = child_data->par.mat_id;
        parent_data->par.grav[0] = child_data->par.grav[0];
        parent_data->par.grav[1] = child_data->par.grav[1];
        parent_data->par.grav[2] = child_data->par.grav[2];
        charm_vect_zero(rhs_ru);
        charm_vect_zero(rhs_rv);
        charm_vect_zero(rhs_rw);
        charm_vect_zero(rhs_re);
        for (j = 0; j < c_count; j++) {
            charm_vect_zero(rhs_rc[j]);
        }
        for (i = 0; i < CHARM_CHILDREN; i++) {
            child_data = (charm_data_t *) outgoing[i]->p.user_data;
            for (m = 0; m < N; m++) {
                for (n = 0; n < N; n++) {
                    ar[m][n] = 0.;
                    for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
                        ar[m][n] += child_data->par.g.quad_gw[igp]*child_data->par.g.quad_gj[igp]
                                            * charm_base_func(child_data->par.g.quad_gp[igp], m, parent_data)
                                            * charm_base_func(child_data->par.g.quad_gp[igp], n, child_data);

                    }
                }
            }
            charm_matr_vect_mult(ar, child_data->par.c.ru, cr);
            charm_vect_add(rhs_ru, cr);
            charm_matr_vect_mult(ar, child_data->par.c.rv, cr);
            charm_vect_add(rhs_rv, cr);
            charm_matr_vect_mult(ar, child_data->par.c.rw, cr);
            charm_vect_add(rhs_rw, cr);
            charm_matr_vect_mult(ar, child_data->par.c.re, cr);
            charm_vect_add(rhs_re, cr);
            for (j = 0; j < c_count; j++) {
                charm_matr_vect_mult(ar, child_data->par.c.rc[j], cr);
                charm_vect_add(rhs_rc[j], cr);
            }
        }
        charm_matr_vect_mult(parent_data->par.g.a_inv, rhs_ru, parent_data->par.c.ru);
        charm_matr_vect_mult(parent_data->par.g.a_inv, rhs_rv, parent_data->par.c.rv);
        charm_matr_vect_mult(parent_data->par.g.a_inv, rhs_rw, parent_data->par.c.rw);
        charm_matr_vect_mult(parent_data->par.g.a_inv, rhs_re, parent_data->par.c.re);
        for (j = 0; j < c_count; j++) {
            charm_matr_vect_mult(parent_data->par.g.a_inv, rhs_rc[j], parent_data->par.c.rc[j]);
        }
    }
    else {
        /* this is refinement */
        parent_data = (charm_data_t *) outgoing[0]->p.user_data;

        for (i = 0; i < CHARM_CHILDREN; i++) {
            charm_geom_quad_calc(p4est, incoming[i], which_tree);
            child_data = (charm_data_t *) incoming[i]->p.user_data;
            child_data->par.mat_id  = parent_data->par.mat_id;
            child_data->par.grav[0] = parent_data->par.grav[0];
            child_data->par.grav[1] = parent_data->par.grav[1];
            child_data->par.grav[2] = parent_data->par.grav[2];
            for (m = 0; m < N; m++) {
                for (n = 0; n < N; n++) {
                    ar[m][n] = 0.;
                    for (igp = 0; igp < CHARM_QUAD_GP_COUNT; igp++) {
                        ar[m][n] += child_data->par.g.quad_gw[igp]*child_data->par.g.quad_gj[igp]
                                    * charm_base_func(child_data->par.g.quad_gp[igp], m, child_data)
                                    * charm_base_func(child_data->par.g.quad_gp[igp], n, parent_data);

                    }
                }
            }
            charm_matr_vect_mult(ar, parent_data->par.c.ru, rhs_ru);
            charm_matr_vect_mult(ar, parent_data->par.c.rv, rhs_rv);
            charm_matr_vect_mult(ar, parent_data->par.c.rw, rhs_rw);
            charm_matr_vect_mult(ar, parent_data->par.c.re, rhs_re);
            for (j = 0; j < c_count; j++) {
                charm_matr_vect_mult(ar, parent_data->par.c.rc[j], rhs_rc[j]);
            }

            charm_matr_vect_mult(child_data->par.g.a_inv, rhs_ru, child_data->par.c.ru);
            charm_matr_vect_mult(child_data->par.g.a_inv, rhs_rv, child_data->par.c.rv);
            charm_matr_vect_mult(child_data->par.g.a_inv, rhs_rw, child_data->par.c.rw);
            charm_matr_vect_mult(child_data->par.g.a_inv, rhs_re, child_data->par.c.re);
            for (j = 0; j < c_count; j++) {
                charm_matr_vect_mult(child_data->par.g.a_inv, rhs_rc[j], child_data->par.c.rc[j]);
            }
        }
    }
}

void charm_adapt_init(p4est_t *p4est)
{
    int                 partforcoarsen;
    int                 recursive;
    charm_data_t         *ghost_data;
    p4est_ghost_t      *ghost;
    charm_ctx_t          *ctx             = (charm_ctx_t *) p4est->user_pointer;

    /* refine and coarsen based on an interpolation error estimate */
    recursive = 1;
    p4est_refine (p4est, recursive, charm_refine_init_err_estimate,
                  charm_init_initial_condition);
    p4est_coarsen (p4est, recursive, charm_coarsen_initial_condition,
                   charm_init_initial_condition);

//    p4est_balance (p4est, CHARM_CONNECT_FACE, charm_init_initial_condition);
//
//    ghost = p4est_ghost_new (p4est, CHARM_CONNECT_FULL);
//    ghost_data = CHARM_ALLOC (charm_data_t, ghost->ghosts.elem_count);
//    p4est_ghost_exchange_data (p4est, ghost, ghost_data);
//
//    p4est_iterate (p4est, ghost, (void *) ghost_data,
//                   charm_ref_flag_quad_iter_fn,
//                   charm_ref_flag_face_iter_fn,
//                   NULL,
//                   NULL);
//    p4est_refine (p4est, recursive, charm_refine_flag_estimate,
//                  charm_init_initial_condition);
//
//    p4est_ghost_destroy (ghost);
//    CHARM_FREE (ghost_data);
//    ghost = NULL;
//    ghost_data = NULL;
//
    partforcoarsen = 1;
    p4est_balance (p4est, CHARM_CONNECT_FACE, charm_init_initial_condition);
    p4est_partition (p4est, partforcoarsen, NULL);

}


static void charm_amr_par_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i;

    for (i = 0; i < CHARM_DIM; i++) {
        data->par.a.grad_u[i] = 0.;
    }
}



static void charm_amr_par_calc_face_iter_bnd (p4est_iter_face_info_t * info, void *user_data)
{
    int i, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    charm_real_t n[3];
    charm_real_t qr, qu, qv, qw, qe;
    charm_real_t bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);


    int8_t face;
    charm_real_t c[2][3], l[3];
    charm_real_t r_[2], p_[2], u_[2], v_[2], w_[2], e_[2];
    charm_cons_t cons;
    charm_prim_t prim[2];
    charm_real_t *x, gw, gj;
    charm_real_t intg[2][5];


    CHARM_ASSERT(info->tree_boundary);

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    CHARM_ASSERT(!side[0]->is_hanging);

    if (side[0]->is.full.is_ghost) {
        CHARM_ASSERT(0);
        udata = &(ghost_data[side[0]->is.full.quadid]);
    } else {
        udata = charm_get_quad_data(side[0]->is.full.quad);//(charm_data_t *) side[0]->is.full.quad->p.user_data;
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
        qr = 0.5*(prim[0].r+prim[1].r);
        bfv = gw * gj;
        for (ibf = 0; ibf < CHARM_DIM; ibf++) {
            udata->par.a.grad_u[ibf] += qr * bfv * n[ibf];
        }
    }
}


static void charm_amr_par_calc_face_iter_inner (p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side, igp, ibf;
    p4est_t                *p4est = info->p4est;
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    charm_real_t                  n[3];
    charm_real_t                  qr, qu, qv, qw, qe;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
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
                for (i = 0; i < 2; i++) {
                    charm_get_fields(udata[i], x, &(cons[i]));
                    charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                }
                for (ibf = 0; ibf < CHARM_DIM; ibf++) {
                    for (i = 0; i < 2; i++) {
                        if (!side[i]->is.full.is_ghost) {
                            udata[i]->par.a.grad_u[ibf] += (i ? -1. : 1.) * gw * gj * 0.5*(prim[0].r+prim[1].r) * n[ibf];
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
                udata[i] = charm_get_quad_data(side[i]->is.full.quad);
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
            qr = 0.5*(prim[0].r+prim[1].r);
            for (ibf = 0; ibf < CHARM_DIM; ibf++) {
                for (i = 0; i < 2; i++) {
                    if (!side[i]->is.full.is_ghost) {
                        bfv = (i ? -1. : 1.) * gw * gj;
                        udata[i]->par.a.grad_u[ibf] += qr * bfv *n[ibf];
                    }
                }
            }
        }
    }
}


static void charm_amr_par_calc_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        charm_amr_par_calc_face_iter_bnd(info, user_data);
    }
    else {
        charm_amr_par_calc_face_iter_inner(info, user_data);
    }
}

static void charm_amr_par_calc_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i;
    charm_real_t              vol = charm_quad_get_volume(data);

    for (i = 0; i < CHARM_DIM; i++) {
        data->par.a.grad_u[i] /= vol;
    }
}

void charm_adapt(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data)
{
    int             recursive       = 0;
    int             callbackorphans = 0;
    charm_ctx_t      *ctx             = (charm_ctx_t *) p4est->user_pointer;

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   charm_amr_par_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   NULL,
                   charm_amr_par_calc_face_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   charm_amr_par_calc_quad_iter_fn,
                   NULL,
                   NULL, NULL);

    p4est_refine_ext (p4est, recursive, ctx->max_level,
                      charm_refine_err_estimate, NULL,
                      charm_replace_quads);
    p4est_coarsen_ext (p4est, recursive, callbackorphans,
                       charm_coarsen_err_estimate, NULL,
                       charm_replace_quads);
    p4est_balance_ext (p4est, CHARM_CONNECT_FACE, NULL,
                       charm_replace_quads);

}
