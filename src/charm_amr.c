//
// Created by zhrv on 26.10.17.
//

#include "charm_amr.h"
#include "charm_geom.h"


static double charm_error_sqr_estimate (p4est_quadrant_t * q)
{
//    charm_data_t       *data = (charm_data_t *) q->p.user_data;
//    int                 i;
//    double              du[P4EST_DIM];
//    double              h = CHARM_GET_H(q->level);
//    double              vol = charm_quad_get_volume((charm_data_t *)q->p.user_data);
//    double              diff2 = 0.;
//
//    for (i = 0; i < CHARM_DIM; i++) {
//        du[i] = data->par.grad.r[i];
//    }
//
//    /* use the approximate derivative to estimate the L2 error */
//    for (i = 0; i < CHARM_DIM; i++) {
//        diff2 += du[i] * du[i];
//    }
//
//    return diff2 * (1. / 12.) * exp(log(vol)*5./3.);
}

static int
charm_refine_err_estimate (p4est_t * p4est, p4est_topidx_t which_tree,
                           p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    double              global_err = ctx->max_err;
    double              global_err2 = global_err * global_err;
    double              h = CHARM_GET_H(q->level);
    double              vol, err2;

    vol = charm_quad_get_volume((charm_data_t *)q->p.user_data);

    err2 = charm_error_sqr_estimate (q);
    if (err2 > global_err2 * vol) {
        return 1;
    }
    else {
        return 0;
    }
}

static int
charm_refine_init_err_estimate (p4est_t * p4est, p4est_topidx_t which_tree,
                                p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
//    double              h = CHARM_GET_H(q->level);
    double              mp[3];

    if (q->level >= ctx->max_level) {
        return 0;
    }

    charm_quad_get_center (q->p.user_data, mp);

//    if (( (-0.001 < mp[2]) && (mp[2] < 0.001) ) || ( (-0.006 < mp[2]) && (mp[2] < -0.004) )) {
    if (( (-0.005 < mp[2]) && (mp[2] < 0.002) )) {
        return 1;
    }
    else {
        return 0;
    }



//        err2 = charm_error_sqr_estimate (q);
//    if (err2 > global_err2 * vol) {
//        return 1;
//    }
//    else {
//        return 0;
//    }
}

static int charm_coarsen_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * children[])
{
    return 0;
}

static int charm_coarsen_err_estimate (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * children[])
{
//    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
//    double              global_err = ctx->max_err;
//    double              global_err2 = global_err * global_err;
//    double              h;
//    charm_data_t       *data;
//    double              vol, err2, childerr2;
//    double              parentu;
//    double              diff;
//    int                 i;
//
//    if (children[0]->level <= ctx->min_level) {
//        return 0;
//    }
//
//
//    h =     CHARM_GET_H(children[0]->level);
//    /* the quadrant's volume is also its volume fraction */
//    vol = h * h * h;
//
//    /* compute the average */
//    parentu = 0.;
//    for (i = 0; i < P4EST_CHILDREN; i++) {
//        data = (charm_data_t *) children[i]->p.user_data;
//        parentu += data->par.c.ro / P4EST_CHILDREN;
//    }
//
//    err2 = 0.;
//    for (i = 0; i < P4EST_CHILDREN; i++) {
//        data = (charm_data_t *) children[i]->p.user_data;
//        childerr2 = charm_error_sqr_estimate (children[i]);
//
//        if (childerr2 > global_err2 * vol) {
//            return 0;
//        }
//        err2 += childerr2;
//        diff = (parentu - data->par.c.ro) * (parentu - data->par.c.ro);
//        err2 += diff * vol;
//    }
//    if (err2 < global_err2 * (vol * P4EST_CHILDREN)) {
//        return 1;
//    }
//    else {
//        return 0;
//    }
//
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

        P4EST_ASSERT(info->tree_boundary);
        P4EST_ASSERT(!side[0]->is_hanging);
        P4EST_ASSERT(!side[0]->is.full.is_ghost);

        ((charm_data_t *) side[0]->is.full.quad->p.user_data)->ref_flag++;

    }
    else {

        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        side[1] = p4est_iter_fside_array_index_int(sides, 1);

        for (i = 0; i < 2; i++) {
            if (side[i]->is_hanging) {
                j = (i+1)%2;
                P4EST_ASSERT(!side[j]->is_hanging);
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
//    charm_data_t       *parent_data, *child_data;
//    int                 i;
//    double              vol, svol;
//
//    if (num_outgoing > 1) {
//        charm_geom_quad_calc(p4est, incoming[0], which_tree);
//        /* this is coarsening */
//        parent_data = (charm_data_t *) incoming[0]->p.user_data;
//        parent_data->par.c.ro = 0.;
//        parent_data->par.c.ru = 0.;
//        parent_data->par.c.rv = 0.;
//        parent_data->par.c.rw = 0.;
//        parent_data->par.c.re = 0.;
//
//        svol = 0.;
//
//        for (i = 0; i < P4EST_CHILDREN; i++) {
//            child_data = (charm_data_t *) outgoing[i]->p.user_data;
//            vol  = charm_quad_get_volume(child_data);
//            svol += vol;
//            parent_data->par.c.ro += child_data->par.c.ro * vol;
//            parent_data->par.c.ru += child_data->par.c.ru * vol;
//            parent_data->par.c.rv += child_data->par.c.rv * vol;
//            parent_data->par.c.rw += child_data->par.c.rw * vol;
//            parent_data->par.c.re += child_data->par.c.re * vol;
//        }
//        parent_data->par.c.ro /= svol;
//        parent_data->par.c.ru /= svol;
//        parent_data->par.c.rv /= svol;
//        parent_data->par.c.rw /= svol;
//        parent_data->par.c.re /= svol;
//    }
//    else {
//        /* this is refinement */
//        parent_data = (charm_data_t *) outgoing[0]->p.user_data;
//
//        for (i = 0; i < P4EST_CHILDREN; i++) {
//            charm_geom_quad_calc(p4est, incoming[i], which_tree);
//            child_data = (charm_data_t *) incoming[i]->p.user_data;
//            child_data->par.c.ro = parent_data->par.c.ro;
//            child_data->par.c.ru = parent_data->par.c.ru;
//            child_data->par.c.rv = parent_data->par.c.rv;
//            child_data->par.c.rw = parent_data->par.c.rw;
//            child_data->par.c.re = parent_data->par.c.re;
//        }
//    }
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

//    p4est_balance (p4est, P4EST_CONNECT_FACE, charm_init_initial_condition);
//
//    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
//    ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
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
//    P4EST_FREE (ghost_data);
//    ghost = NULL;
//    ghost_data = NULL;
//
    partforcoarsen = 1;
    p4est_balance (p4est, P4EST_CONNECT_FACE, charm_init_initial_condition);
    p4est_partition (p4est, partforcoarsen, NULL);

}

void charm_adapt(p4est_t *p4est)
{
    int             recursive       = 0;
    int             callbackorphans = 0;
    charm_ctx_t      *ctx             = (charm_ctx_t *) p4est->user_pointer;

    p4est_refine_ext (p4est, recursive, ctx->max_level,
                      charm_refine_err_estimate, NULL,
                      charm_replace_quads);
    p4est_coarsen_ext (p4est, recursive, callbackorphans,
                       charm_coarsen_err_estimate, NULL,
                       charm_replace_quads);
    p4est_balance_ext (p4est, P4EST_CONNECT_FACE, NULL,
                       charm_replace_quads);

}
