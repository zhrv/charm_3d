//
// Created by zhrv on 26.10.17.
//

#include <p8est_iterate.h>
#include "charm_timestep.h"
#include "charm_bnd_cond.h"

static void charm_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data);
static void charm_convect_flux_face_iter_fn (p4est_iter_face_info_t * info, void *user_data);
static void charm_diffusion_flux_face_iter_fn (p4est_iter_face_info_t * info, void *user_data);
static void charm_grad_face_iter_fn (p4est_iter_face_info_t * info, void *user_data);
static void charm_grad_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data);
static double charm_get_timestep (p4est_t * p4est);
static void charm_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data);


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

    charm_calc_grad(p4est, ghost, ghost_data);


    for (t = 0., i = 0; t < time; t += dt, i++) {

        /* refine */
        if (!(i % refine_period)) {
            if (i) {
                umax = 1.;
//                p4est_iterate (p4est, NULL, (void *) &umax,     /* pass in ghost data that we just exchanged */
//                               charm_compute_max,       /* blank the previously calculated derivatives */
//                               NULL,    /* there is no callback for the faces between quadrants */
//                               NULL);   /* there is no callback for the corners between quadrants */

//                mpiret =
//                        sc_MPI_Allreduce (&umax, &global_umax, 1, sc_MPI_DOUBLE, sc_MPI_MAX,
//                                          p4est->mpicomm);
//                SC_CHECK_MPI (mpiret);
//                ctx->max_err = orig_max_err * global_umax;
//                P4EST_GLOBAL_PRODUCTIONF ("u_max %f\n", global_umax);

                charm_adapt(p4est); /* adapt */
                if (ghost) {
                    p4est_ghost_destroy(ghost);
                    P4EST_FREE (ghost_data);
                    ghost = NULL;
                    ghost_data = NULL;
                }


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
            DBG_CH(p4est->mpirank);

            ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
            ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
            p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        }
        charm_calc_grad(p4est, ghost, ghost_data);
        /* compute du/dt */
        p4est_iterate (p4est,                                /* the forest */
                       ghost,                                /* the ghost layer */
                       (void *) ghost_data,                  /* the synchronized ghost data */
                       charm_zero_quad_iter_fn,              /* callback to compute each quad's interior contribution to du/dt */
                       charm_convect_flux_face_iter_fn,      /* callback to compute each quads' faces' contributions to du/du */
                       NULL,
                       NULL);                                /* there is no callback for the
                                                                corners between quadrants */

        p4est_iterate (p4est,                                /* the forest */
                       ghost,                                /* the ghost layer */
                       (void *) ghost_data,                  /* the synchronized ghost data */
                       NULL,                                 /* callback to compute each quad's interior contribution to du/dt */
                       charm_diffusion_flux_face_iter_fn,    /* callback to compute each quads' faces' contributions to du/du */
                       NULL,
                       NULL);                                /* there is no callback for the
                                                                corners between quadrants */
        /* update u */
        p4est_iterate (p4est, NULL,                              /* ghosts are not needed for this loop */
                       (void *) &dt,                             /* pass in dt */
                       charm_timestep_update_quad_iter_fn,       /* update each cell */
                       NULL,                                     /* there is no callback for the faces between quadrants */
                       NULL,                                     /* there is no callback for the faces between quadrants */
                       NULL);                                    /* there is no callback for the corners between quadrants */

        /* synchronize the ghost data */
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);

        /* update grad */
        charm_calc_grad(p4est, ghost, ghost_data);
    }

    P4EST_FREE (ghost_data);
    p4est_ghost_destroy (ghost);
}


static void charm_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;

    data->drodt = 0.;
    data->drudt = 0.;
    data->drvdt = 0.;
    data->drwdt = 0.;
    data->dredt = 0.;

}




static void charm_convect_flux_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    int                 i, j;
    p4est_t            *p4est = info->p4est;
    charm_data_t       *ghost_data = (charm_data_t *) user_data;
    charm_data_t       *udata[2];
    p4est_quadrant_t   *quad;
    double              n[3];
    double              ro_avg[2], ru_avg[2], rv_avg[2], rw_avg[2], re_avg[2];
    double              qr, qu, qv, qw, qe;
    double              facearea;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);

    int bnd = 0;
    int8_t face[2];
    double c[2][3], l[3];
    int8_t is_hanging = 0;
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

        charm_bnd_cond(p4est, side[0]->treeid, face[0], &(udata[0]->par), &(udata[1]->par));

        charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[0]->treeid);
        charm_param_cons_to_prim(attr->reg->mat, &(udata[1]->par));


        for (i = 0; i < 2; i++) {
            charm_tree_attr_t *attr = charm_get_tree_attr(p4est, side[0]->treeid);
            charm_param_cons_to_prim(attr->reg->mat, &(udata[i]->par));
            r_[i] = udata[i]->par.p.r;
            u_[i] = udata[i]->par.p.u;
            v_[i] = udata[i]->par.p.v;
            w_[i] = udata[i]->par.p.w;
            p_[i] = udata[i]->par.p.p;
        }

        /* flux from side 0 to side 1 */
        calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, n, bnd);

        if (!side[0]->is.full.is_ghost) {
            //facearea = charm_face_get_area(udata[0], side[0]->face);
            udata[0]->drodt -= qr * facearea;
            udata[0]->drudt -= qu * facearea;
            udata[0]->drvdt -= qv * facearea;
            udata[0]->drwdt -= qw * facearea;
            udata[0]->dredt -= qe * facearea;
        }
        P4EST_FREE(udata[1]);
    }
    else {
        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        side[1] = p4est_iter_fside_array_index_int(sides, 1);
        face[0] = side[0]->face;
        face[1] = side[1]->face;

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
                        facearea = charm_face_get_area(udata[i], side[i]->face);
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

                /* flux from side 0 to side 1 */
                calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, n, bnd);

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
                    //facearea = charm_face_get_area(udata[i], side[i]->face);
                    udata[i]->drodt += qr * facearea * (i ? 1. : -1.);
                    udata[i]->drudt += qu * facearea * (i ? 1. : -1.);
                    udata[i]->drvdt += qv * facearea * (i ? 1. : -1.);
                    udata[i]->drwdt += qw * facearea * (i ? 1. : -1.);
                    udata[i]->dredt += qe * facearea * (i ? 1. : -1.);

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

            /* flux from side 0 to side 1 */
            calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, n, bnd);

            for (i = 0; i < 2; i++) {
                if (!side[i]->is.full.is_ghost) {
                    facearea = charm_face_get_area(udata[i], side[i]->face);
                    udata[i]->drodt += qr * facearea * (i ? 1. : -1.);
                    udata[i]->drudt += qu * facearea * (i ? 1. : -1.);
                    udata[i]->drvdt += qv * facearea * (i ? 1. : -1.);
                    udata[i]->drwdt += qw * facearea * (i ? 1. : -1.);
                    udata[i]->dredt += qe * facearea * (i ? 1. : -1.);
                }
            }

        }

    }

}


static void charm_diffusion_flux_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{

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




static void charm_timestep_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    double              dt = *((double *) user_data);
    double              vol;

    vol = charm_quad_get_volume(data);

    data->par.c.ro += dt * data->drodt / vol;
    data->par.c.ru += dt * data->drudt / vol;
    data->par.c.rv += dt * data->drvdt / vol;
    data->par.c.rw += dt * data->drwdt / vol;
    data->par.c.re += dt * data->dredt / vol;
}

