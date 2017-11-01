//
// Created by zhrv on 26.10.17.
//

#include <p8est_iterate.h>
#include "charm_timestep.h"
#include "charm_bnd_cond.h"

static void charm_quad_divergence (p4est_iter_volume_info_t * info, void *user_data);
static void charm_upwind_flux (p4est_iter_face_info_t * info, void *user_data);
static void charm_grad_face_iter_fn (p4est_iter_face_info_t * info, void *user_data);
static void charm_grad_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data);
static double charm_get_timestep (p4est_t * p4est);
static void charm_timestep_update (p4est_iter_volume_info_t * info, void *user_data);
static void charm_reset_derivatives (p4est_iter_volume_info_t * info, void *user_data);


/** Timestep the advection problem.
 *
 * Update the state, refine, repartition, and write the solution to file.
 *
 * \param [in,out] p4est the forest, whose state is updated
 * \param [in] time      the end time
 */
void charm_timestep (p4est_t * p4est, double time)
{
    double              t = 0.;
    double              dt = 0.;
    int                 i;
    charm_data_t       *ghost_data;
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    int                 refine_period = ctx->refine_period;
    int                 repartition_period = ctx->repartition_period;
    int                 write_period = ctx->write_period;
    int                 allowcoarsening = 1;
    int                 mpiret;
    double              orig_max_err = ctx->max_err;
    double              umax, global_umax;
    p4est_ghost_t      *ghost;


    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

//    charm_calc_grad(p4est, ghost, ghost_data);


    for (t = 0., i = 0; t < time; t += dt, i++) {


        /* refine */
        if (!(i % refine_period)) {
            if (i) {
                umax = 1.;
//                p4est_iterate (p4est, NULL, (void *) &umax,     /* pass in ghost data that we just exchanged */
//                               charm_compute_max,       /* blank the previously calculated derivatives */
//                               NULL,    /* there is no callback for the faces between quadrants */
//                               NULL);   /* there is no callback for the corners between quadrants */

                mpiret =
                        sc_MPI_Allreduce (&umax, &global_umax, 1, sc_MPI_DOUBLE, sc_MPI_MAX,
                                          p4est->mpicomm);
                SC_CHECK_MPI (mpiret);
                ctx->max_err = orig_max_err * global_umax;
                P4EST_GLOBAL_PRODUCTIONF ("u_max %f\n", global_umax);

//                charm_adapt(p4est, &ghost, &ghost_data); /* adapt */

            }
            dt = charm_get_timestep (p4est);
        }

        /* repartition */
        if (i && !(i % repartition_period)) {

            p4est_partition (p4est, allowcoarsening, NULL);

            if (ghost) {
                p4est_ghost_destroy (ghost);
                P4EST_FREE (ghost_data);
                ghost = NULL;
                ghost_data = NULL;
            }
        }


        /* write out solution */
        if (!(i % write_period)) {
            charm_write_solution (p4est, i);
            P4EST_GLOBAL_ESSENTIALF ("**************** File for %6d step is saved ***************\n", i);
        }


        /* synchronize the ghost data */
        if (!ghost) {
            ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
            ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
            p4est_ghost_exchange_data (p4est, ghost, ghost_data);
        }

        //charm_calc_grad(p4est, ghost, ghost_data);

        /* compute du/dt */
        p4est_iterate (p4est,                 /* the forest */
                       ghost,                 /* the ghost layer */
                       (void *) ghost_data,   /* the synchronized ghost data */
                       charm_quad_divergence,   /* callback to compute each quad's interior contribution to du/dt */
                       charm_upwind_flux,       /* callback to compute each quads' faces' contributions to du/du */
                       NULL,
                       NULL);                 /* there is no callback for the
                                             corners between quadrants */


        /* update u */
        p4est_iterate (p4est, NULL,               /* ghosts are not needed for this loop */
                       (void *) &dt,              /* pass in dt */
                       charm_timestep_update,       /* update each cell */
                       NULL,                      /* there is no callback for the faces between quadrants */
                       NULL,                      /* there is no callback for the faces between quadrants */
                       NULL);                     /* there is no callback for the corners between quadrants */


        /* synchronize the ghost data */
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);


        /* update du/dx estimate */
//        charm_calc_grad(p4est, ghost, ghost_data);
    }

    P4EST_FREE (ghost_data);
    p4est_ghost_destroy (ghost);
}


static void charm_quad_divergence (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;

    data->drodt = 0.;
    data->drudt = 0.;
    data->drvdt = 0.;
    data->drwdt = 0.;
    data->dredt = 0.;
}




static void charm_upwind_flux (p4est_iter_face_info_t * info, void *user_data)
{
    int                 i, j;
    p4est_t            *p4est = info->p4est;
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *ghost_data = (charm_data_t *) user_data;
    charm_data_t       *udata, edata;
    p4est_quadrant_t   *quad;
    double              vdotn = 0., n[3];
    double              ro_avg[2], ru_avg[2], rv_avg[2], rw_avg[2], re_avg[2];
    double              qr, qu, qv, qw, qe;
    double              h, facearea, dh, vdx, vdy, vdz;
    int                 which_face;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);

    int bnd = 0;
    int8_t face[2];
    double c[2][3], l[3], scp;

    side[0] = p4est_iter_fside_array_index_int(sides, 0);

    face[0] = side[0]->face;
    facearea = charm_face_get_normal(side[0]->is.full.quad, face[0], n);
    charm_quad_get_center(side[0]->is.full.quad, c[0]);
    charm_face_get_center(side[0]->is.full.quad, face[0], c[1]);
    for (i = 0; i < 3; i++) {
        l[i] = c[1][i]-c[0][i];
    }

    if (scalar_prod(n, l) < 0) {
        for (i = 0; i < 3; i++) {
            n[i] *= -1.0;
        }
    }

    if (sides->elem_count != 2) {

        P4EST_ASSERT(!side[0]->is_hanging);

        udata = (charm_data_t *) side[0]->is.full.quad->p.user_data;

        ro_avg[0] = udata->par.c.ro;
        ru_avg[0] = udata->par.c.ru;
        rv_avg[0] = udata->par.c.rv;
        rw_avg[0] = udata->par.c.rw;
        re_avg[0] = udata->par.c.re;

        charm_bnd_cond(p4est, side[0]->treeid, face[0], &udata->par, &edata.par);

        ro_avg[1] = edata.par.c.ro;
        ru_avg[1] = edata.par.c.ru;
        rv_avg[1] = edata.par.c.rv;
        rw_avg[1] = edata.par.c.rw;
        re_avg[1] = edata.par.c.re;

    }
    else {
        side[1] = p4est_iter_fside_array_index_int(sides, 1);
        face[1] = side[1]->face;

        for (i = 0; i < 2; i++) {
            ro_avg[i] = 0;
            ru_avg[i] = 0;
            rv_avg[i] = 0;
            rw_avg[i] = 0;
            re_avg[i] = 0;
            if (side[i]->is_hanging) {
                P4EST_ASSERT(0);
//                /* there are 2^(d-1) (P4EST_HALF) subfaces */
//                for (j = 0; j < P4EST_HALF; j++) {
//                    if (side[i]->is.hanging.is_ghost[j]) {
//                        udata = (charm_data_t *) &ghost_data[side[i]->is.hanging.quadid[j]];
//                    }
//                    else {
//                        udata =
//                                (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
//                    }
//                    ro_avg[i] += udata->par.p.ro;
//                    ru_avg[i] += udata->par.p.ru;
//                    rv_avg[i] += udata->par.p.rv;
//                    rw_avg[i] += udata->par.p.rw;
//                    re_avg[i] += udata->par.p.re;
//                }
//                ro_avg[i] /= P4EST_HALF;
//                ru_avg[i] /= P4EST_HALF;
//                rv_avg[i] /= P4EST_HALF;
//                rw_avg[i] /= P4EST_HALF;
//                re_avg[i] /= P4EST_HALF;
            }
            else {
                if (side[i]->is.full.is_ghost) {
                    udata = (charm_data_t *) &ghost_data[side[i]->is.full.quadid];
                }
                else {
                    udata = (charm_data_t *) side[i]->is.full.quad->p.user_data;
                }
                ro_avg[i] = udata->par.c.ro;
                ru_avg[i] = udata->par.c.ru;
                rv_avg[i] = udata->par.c.rv;
                rw_avg[i] = udata->par.c.rw;
                re_avg[i] = udata->par.c.re;
            }
        }

    }

    double ri, pi, ei, ui, vi, wi;
    double r_[2], p_[2], u_[2], v_[2], w_[2];
    for (i = 0; i < 2; i++) {
        r_[i] = ro_avg[i];
        u_[i] = ru_avg[i]/r_[i];
        v_[i] = rv_avg[i]/r_[i];
        w_[i] = rw_avg[i]/r_[i];
        p_[i] = (re_avg[i]-0.5*r_[i]*(u_[i]*u_[i]+v_[i]*v_[i]+w_[i]*w_[i]))*(GAM-1.0);
        w_[i] = 0.0;
    }

    /* flux from side 0 to side 1 */
    calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, n, bnd);

    for (i = 0; i < sides->elem_count; i++) {
        if (side[i]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
//            for (j = 0; j < P4EST_HALF; j++) {
//                quad = side[i]->is.hanging.quad[j];
//                h = CHARM_GET_H(quad->level);
//                facearea = h * h;
//                if (!side[i]->is.hanging.is_ghost[j]) {
//                    udata = (charm_data_t *) quad->p.user_data;
//                    udata->drodt += qr * facearea * (i ? 1. : -1.);
//                    udata->drudt += qu * facearea * (i ? 1. : -1.);
//                    udata->drvdt += qv * facearea * (i ? 1. : -1.);
//                    udata->drwdt += qw * facearea * (i ? 1. : -1.);
//                    udata->dredt += qe * facearea * (i ? 1. : -1.);
//                }
//            }
        }
        else {
            quad = side[i]->is.full.quad;
            if (!side[i]->is.full.is_ghost) {
                udata = (charm_data_t *) quad->p.user_data;
                udata->drodt += qr * facearea * (i ? 1. : -1.);
                udata->drudt += qu * facearea * (i ? 1. : -1.);
                udata->drvdt += qv * facearea * (i ? 1. : -1.);
                udata->drwdt += qw * facearea * (i ? 1. : -1.);
                udata->dredt += qe * facearea * (i ? 1. : -1.);
            }
        }
    }
}






/** Compute the timestep.
 *
 * Find the smallest quadrant and scale the timestep based on that length and
 * the advection velocity.
 *
 * \param [in] p4est the forest
 * \return the timestep.
 */
static double charm_get_timestep (p4est_t * p4est)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
//    p4est_topidx_t      t, flt, llt;
//    p4est_tree_t       *tree;
//    int                 max_level, global_max_level;
//    int                 mpiret, i;
//    double              min_h, vnorm;
//    double              dt;
//
//    /* compute the timestep by finding the smallest quadrant */
//    flt = p4est->first_local_tree;
//    llt = p4est->last_local_tree;
//
//    max_level = 0;
//    for (t = flt; t <= llt; t++) {
//        tree = p4est_tree_array_index (p4est->trees, t);
//        max_level = SC_MAX (max_level, tree->maxlevel);
//
//    }
//    mpiret =
//            sc_MPI_Allreduce (&max_level, &global_max_level, 1, sc_MPI_INT,
//                              sc_MPI_MAX, p4est->mpicomm);
//    SC_CHECK_MPI (mpiret);
//
//    min_h = CHARM_GET_H(global_max_level);
//
//    vnorm = 0;
//    for (i = 0; i < P4EST_DIM; i++) {
////        vnorm += ctx->v[i] * ctx->v[i];
//    }
//    vnorm = sqrt (vnorm);
//
//    min_h = CHARM_GET_H(ctx->allowed_level);
//
//    dt = ctx->CFL* min_h  / ctx->v_ref;

    return ctx->dt;
}




static void
charm_timestep_update (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    double              dt = *((double *) user_data);
    double              vol;
    double              h = CHARM_GET_H(q->level);

    vol = charm_quad_get_volume(q);

    data->par.c.ro += dt * data->drodt / vol;
    data->par.c.ru += dt * data->drudt / vol;
    data->par.c.rv += dt * data->drvdt / vol;
    data->par.c.rw += dt * data->drwdt / vol;
    data->par.c.re += dt * data->dredt / vol;
}

static void
charm_reset_derivatives (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t         *data = (charm_data_t *) q->p.user_data;
    int                 j;

//    for (j = 0; j < P4EST_DIM; j++) {
//        data->dro[j] = 0.;
//        data->dru[j] = 0.;
//        data->drv[j] = 0.;
//        data->drw[j] = 0.;
//        data->dre[j] = 0.;
//    }
}
