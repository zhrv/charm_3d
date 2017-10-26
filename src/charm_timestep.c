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


static double charm_face_get_normal(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, double* n)
{
    const int ftv[6][4] =
            {{ 0, 2, 4, 6 },
             { 1, 3, 5, 7 },
             { 0, 1, 4, 5 },
             { 2, 3, 6, 7 },
             { 0, 1, 2, 3 },
             { 4, 5, 6, 7 }};

    int i,k;
    double x[4][3], v[2][3], nl;
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    for (i = 0; i < 4; i++) {
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l * (ftv[face][i] % 2), q->y + l * ((ftv[face][i] / 2) % 2), q->z + l * (ftv[face][i] / 4), x[i]);
    }

    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[0][i];
        v[1][i] = x[2][i]-x[0][i];
    }

    vect_prod(v[0], v[1], n);
    nl = vect_length(n);
    for (i = 0; i < 3; i++) {
        n[i] /= nl;
    }

    return nl;
}


static double charm_face_get_area(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face)
{
    const int ftv[6][4] =
            {{ 0, 2, 4, 6 },
             { 1, 3, 5, 7 },
             { 0, 1, 4, 5 },
             { 2, 3, 6, 7 },
             { 0, 1, 2, 3 },
             { 4, 5, 6, 7 }};

    int i,k;
    double x[4][3], v[2][3], n[3], s1, s2;
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    for (i = 0; i < 4; i++) {
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l * (ftv[face][i] % 2), q->y + l * ((ftv[face][i] / 2) % 2), q->z + l * (ftv[face][i] / 4), x[i]);
    }

    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[0][i];
        v[1][i] = x[2][i]-x[0][i];
    }

    vect_prod(v[0], v[1], n);
    s1 = 0.5*vect_length(n);
    for (i = 0; i < 3; i++) {
        v[0][i] = x[1][i]-x[3][i];
        v[1][i] = x[2][i]-x[3][i];
    }

    vect_prod(v[0], v[1], n);
    s2 = 0.5*vect_length(n);

    return s1+s2;
}



static void charm_quad_get_center(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, double* c)
{
    p4est_qcoord_t l2 = P4EST_QUADRANT_LEN(q->level) / 2;
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + l2, q->y + l2, q->z + l2, c);


}


static void charm_face_get_center(p4est_t* p4est, p4est_quadrant_t* q, p4est_topidx_t treeid, int8_t face, double* c)
{
    p4est_qcoord_t l  = P4EST_QUADRANT_LEN(q->level);
    p4est_qcoord_t l2 = l / 2;
    p4est_qcoord_t fc[6][3] = {
            {0, l2, l2},
            {l, l2, l2},
            {l2, 0, l2},
            {l2, l, l2},
            {l2, l2, 0},
            {l2, l2, l},
    };
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, q->x + fc[face][0], q->y + fc[face][1], q->z + fc[face][2], c);
}



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

                mpiret =
                        sc_MPI_Allreduce (&umax, &global_umax, 1, sc_MPI_DOUBLE, sc_MPI_MAX,
                                          p4est->mpicomm);
                SC_CHECK_MPI (mpiret);
                ctx->max_err = orig_max_err * global_umax;
                P4EST_GLOBAL_PRODUCTIONF ("u_max %f\n", global_umax);

                charm_adapt(p4est, &ghost, &ghost_data); /* adapt */

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
        charm_calc_grad(p4est, ghost, ghost_data);
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
    charm_data_t       *udata;
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
    facearea = charm_face_get_normal(p4est, side[0]->is.full.quad, side[0]->treeid, face[0], n);
    charm_quad_get_center(p4est, side[0]->is.full.quad, side[0]->treeid, c[0]);
    charm_face_get_center(p4est, side[0]->is.full.quad, side[0]->treeid, face[0], c[1]);
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

        ro_avg[0] = udata->par.p.ro;
        ru_avg[0] = udata->par.p.ru;
        rv_avg[0] = udata->par.p.rv;
        rw_avg[0] = udata->par.p.rw;
        re_avg[0] = udata->par.p.re;

        charm_bnd_cond(p4est, side[0]->treeid, face[0], ro_avg[0], ru_avg[0], rv_avg[0], rw_avg[0], re_avg[0], &ro_avg[1], &ru_avg[1], &rv_avg[1], &rw_avg[1], &re_avg[1]);

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
                ro_avg[i] = udata->par.p.ro;
                ru_avg[i] = udata->par.p.ru;
                rv_avg[i] = udata->par.p.rv;
                rw_avg[i] = udata->par.p.rw;
                re_avg[i] = udata->par.p.re;
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


static void charm_grad_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    int                 i, j;
    p4est_t            *p4est = info->p4est;
    charm_ctx_t          *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t         *ghost_data = (charm_data_t *) user_data;
    charm_data_t         *udata;
    p4est_quadrant_t   *quad;
    double              vdotn = 0.;
    double              ro[2];
    double              ru[2];
    double              rv[2];
    double              rw[2];
    double              re[2];
    double              qr, qu, qv, qw, qe;
    double              h, facearea;
    int                 which_face;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        side[0] = p4est_iter_fside_array_index_int(sides, 0);

        P4EST_ASSERT(info->tree_boundary);

        /* which of the quadrant's faces the interface touches */
        which_face = side[0]->face;

        ro[0] = 0;
        ru[0] = 0;
        rv[0] = 0;
        rw[0] = 0;
        re[0] = 0;
        if (side[0]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P4EST_HALF; j++) {
                if (side[0]->is.hanging.is_ghost[j]) {
                    udata =
                            (charm_data_t *) &ghost_data[side[0]->is.hanging.quadid[j]];
                }
                else {
                    udata =
                            (charm_data_t *) side[0]->is.hanging.quad[j]->p.user_data;
                }
                ro[0] += udata->par.p.ro;
                ru[0] += udata->par.p.ru;
                rv[0] += udata->par.p.rv;
                rw[0] += udata->par.p.rw;
                re[0] += udata->par.p.re;
            }
            ro[0] /= P4EST_HALF;
            ru[0] /= P4EST_HALF;
            rv[0] /= P4EST_HALF;
            rw[0] /= P4EST_HALF;
            re[0] /= P4EST_HALF;
        }
        else {
            if (side[0]->is.full.is_ghost) {
                udata = (charm_data_t *) &ghost_data[side[0]->is.full.quadid];
            }
            else {
                udata = (charm_data_t *) side[0]->is.full.quad->p.user_data;
            }
            ro[0] = udata->par.p.ro;
            ru[0] = udata->par.p.ru;
            rv[0] = udata->par.p.rv;
            rw[0] = udata->par.p.rw;
            re[0] = udata->par.p.re;
        }

        /* boundary conditions */
        double r_, p_, u_, v_;
        switch (which_face) {
            case 0:                      /* -x side */
                vdotn = -1.0;

                ro[1] =  ro[0];
                ru[1] = -ru[0];
                rv[1] =  rv[0];
                rw[1] =  rw[0];
                re[1] =  re[0];
                break;

            case 1:                      /* +x side */
                vdotn = 1.0;

                ro[1] =  ro[0];
                ru[1] = -ru[0];
                rv[1] =  rv[0];
                rw[1] =  rw[0];
                re[1] =  re[0];
                break;

            case 2:                      /* -y side */
                vdotn = -1.0;

                ro[1] =  ro[0];
                ru[1] =  ru[0];
                rv[1] = -rv[0];
                rw[1] =  rw[0];
                re[1] =  re[0];
                break;

            case 3:                      /* +y side */
                vdotn = 1.0;

                ro[1] =  ro[0];
                ru[1] =  ru[0];
                rv[1] = -rv[0];
                rw[1] =  rw[0];
                re[1] =  re[0];
                break;
            case 4:                      /* -z side */
                vdotn = -1.0;
//
//                r_ = 12.09;
//                u_ = 0.0;
//                v_ = 97.76;
//                p_ = 2.152e+5;
//
//                ro[1] = r_;
//                ru[1] = r_*u_;
//                rv[1] = r_*v_;
//                re[1] = p_/(GAM-1.0)+0.5*r_*(u_*u_+v_*v_);
//
                ro[1] =  ro[0];
                ru[1] =  ru[0];
                rv[1] =  rv[0];
                rw[1] =  rw[0];
                re[1] =  re[0];
                break;

            case 5:                      /* +z side */
                vdotn = 1.0;

                ro[1] =  ro[0];
                ru[1] =  ru[0];
                rv[1] =  rv[0];
                rw[1] = -rw[0];
                re[1] =  re[0];
                break;
        }
    }
    else {

        side[0] = p4est_iter_fside_array_index_int(sides, 0);
        side[1] = p4est_iter_fside_array_index_int(sides, 1);

        /* which of the quadrant's faces the interface touches */
        which_face = side[0]->face;

        switch (which_face) {
            case 0:                      /* -x side */
                vdotn = -1.0;
                break;
            case 1:                      /* +x side */
                vdotn = 1.0;
                break;
            case 2:                      /* -y side */
                vdotn = -1.0;
                break;
            case 3:                      /* +y side */
                vdotn = 1.0;
                break;
            case 4:                      /* -z side */
                vdotn = -1.0;
                break;
            case 5:                      /* +z side */
                vdotn = 1.0;
                break;
        }

//        P4EST_ASSERT (vdotn == 1.0);

        for (i = 0; i < 2; i++) {
            ro[i] = 0;
            ru[i] = 0;
            rv[i] = 0;
            rw[i] = 0;
            re[i] = 0;
            if (side[i]->is_hanging) {
                /* there are 2^(d-1) (P4EST_HALF) subfaces */
                for (j = 0; j < P4EST_HALF; j++) {
                    if (side[i]->is.hanging.is_ghost[j]) {
                        udata =
                                (charm_data_t *) &ghost_data[side[i]->is.hanging.quadid[j]];
                    }
                    else {
                        udata =
                                (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
                    }
                    ro[i] += udata->par.p.ro;
                    ru[i] += udata->par.p.ru;
                    rv[i] += udata->par.p.rv;
                    rw[i] += udata->par.p.rw;
                    re[i] += udata->par.p.re;
                }
                ro[i] /= P4EST_HALF;
                ru[i] /= P4EST_HALF;
                rv[i] /= P4EST_HALF;
                rw[i] /= P4EST_HALF;
                re[i] /= P4EST_HALF;
            }
            else {
                if (side[i]->is.full.is_ghost) {
                    udata = (charm_data_t *) &ghost_data[side[i]->is.full.quadid];
                }
                else {
                    udata = (charm_data_t *) side[i]->is.full.quad->p.user_data;
                }
                ro[i] = udata->par.p.ro;
                ru[i] = udata->par.p.ru;
                rv[i] = udata->par.p.rv;
                rw[i] = udata->par.p.rw;
                re[i] = udata->par.p.re;
            }
        }

    }


    /* flux from side 0 to side 1 */
    qr = sqrt(ro[0]*ro[1]);

    double fSB = sqrt( ro[0] );
    double fSE = sqrt( ro[1] );
    double fS_ = 1.0 / ( fSB + fSE );

    //qr = fSB * fSE;

    double UB = ru[0]/ro[0];
    double UE = ru[1]/ro[1];
    double VB = rv[0]/ro[0];
    double VE = rv[1]/ro[1];
    double WB = rw[0]/ro[0];
    double WE = rw[1]/ro[1];

    double UI = ( fSB * UB + fSE * UE ) * fS_;
    double VI = ( fSB * UB + fSE * VE ) * fS_;
    double WI = ( fSB * WB + fSE * WE ) * fS_;

    //double PB = e*ro[0]*

    double UB_MAG = UB*UB+VB*VB+WB*WB;
    double UE_MAG = UE*UE+VE*VE+WE*WE;
    double EB = re[0]/ro[0]-0.5*UB_MAG;
    double EE = re[1]/ro[1]-0.5*UE_MAG;
    double PB = ro[0]*EB*(GAM-1.0);
    double PE = ro[1]*EE*(GAM-1.0);


    double HB = EB + UB_MAG*0.5 + PB / ro[0];
    double HE = EE + UE_MAG*0.5 + PE / ro[1];

    double HI = ( fSB * HB + fSE * HE ) * fS_;

    double PI = ( HI - (UI*UI+VI*VI+WI*WI)*0.5 ) * qr * ( GAM - 1.0 ) / GAM;
    double EI = PI/(qr*(GAM-1.0));

    qu = UI*qr;
    qv = VI*qr;
    qw = WI*qr;
    qe = qr*(EI+0.5*(UI*UI+VI*VI+WI*WI));

    for (i = 0; i < sides->elem_count; i++) {
        if (side[i]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P4EST_HALF; j++) {
                quad = side[i]->is.hanging.quad[j];
                h = CHARM_GET_H(quad->level);
                facearea = h * h;
                if (!side[i]->is.hanging.is_ghost[j]) {
                    udata = (charm_data_t *) quad->p.user_data;
                    udata->dro[which_face/2] += qr * facearea * (i ? 1. : -1.) * vdotn;
                    udata->dru[which_face/2] += qu * facearea * (i ? 1. : -1.) * vdotn;
                    udata->drv[which_face/2] += qv * facearea * (i ? 1. : -1.) * vdotn;
                    udata->drw[which_face/2] += qw * facearea * (i ? 1. : -1.) * vdotn;
                    udata->dre[which_face/2] += qe * facearea * (i ? 1. : -1.) * vdotn;
                }
            }
        }
        else {
            quad = side[i]->is.full.quad;
            h = CHARM_GET_H(quad->level);
            facearea = h * h;
            if (!side[i]->is.full.is_ghost) {
                udata = (charm_data_t *) quad->p.user_data;
                udata->dro[which_face/2] += qr * facearea * (i ? 1. : -1.) * vdotn;
                udata->dru[which_face/2] += qu * facearea * (i ? 1. : -1.) * vdotn;
                udata->drv[which_face/2] += qv * facearea * (i ? 1. : -1.) * vdotn;
                udata->drw[which_face/2] += qw * facearea * (i ? 1. : -1.) * vdotn;
                udata->dre[which_face/2] += qe * facearea * (i ? 1. : -1.) * vdotn;
            }
        }
    }
}

static void
charm_grad_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    int i;

    double h = CHARM_GET_H(q->level);
    double vol = h*h*h;

    for (i = 0; i < P4EST_DIM; i++) {
        data->dro[i] /= vol;
        data->dru[i] /= vol;
        data->drv[i] /= vol;
        data->drw[i] /= vol;
        data->dre[i] /= vol;
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

    p4est_iterate (p4est, ghost, (void *) ghost_data,     /* pass in ghost data that we just exchanged */
                   charm_reset_derivatives,       /* blank the previously calculated derivatives */
                   charm_grad_face_iter_fn, /* compute the minmod estimate of each cell's derivative */
                   NULL,
                   NULL);         /* there is no callback for the corners between quadrants */


    p4est_iterate (p4est, ghost, (void *) ghost_data,     /* pass in ghost data that we just exchanged */
                   charm_grad_quad_iter_fn,       /* blank the previously calculated derivatives */
                   NULL,
                   NULL,
                   NULL);

    if (ghost && my_ghost) {
        p4est_ghost_destroy (ghost);
        P4EST_FREE (ghost_data);
        ghost      = NULL;
        ghost_data = NULL;
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

    vol = h * h * h;

    data->par.p.ro += dt * data->drodt / vol;
    data->par.p.ru += dt * data->drudt / vol;
    data->par.p.rv += dt * data->drvdt / vol;
    data->par.p.rw += dt * data->drwdt / vol;
    data->par.p.re += dt * data->dredt / vol;
}

static void
charm_reset_derivatives (p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t         *data = (charm_data_t *) q->p.user_data;
    int                 j;

    for (j = 0; j < P4EST_DIM; j++) {
        data->dro[j] = 0.;
        data->dru[j] = 0.;
        data->drv[j] = 0.;
        data->drw[j] = 0.;
        data->dre[j] = 0.;
    }
}
