#define P4EST_ENABLE_DEBUG

#include <p4est_to_p8est.h>

#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>


#include "charm_globals.h"
#include "charm_connectivity.h"



static void charm_initial_condition (double x[], double u[FLD_COUNT], double du[FLD_COUNT][P4EST_DIM], charm_ctx_t * ctx)
{
    int                 i;
    double r_,p_,u_,v_,w_,e_;

    double pi = 4.0*atan(1.0);
    double pi2 = pi*2.0;

    if (x[2] < -0.005) {
        r_ = 12.09;
        u_ = 0.0;
        v_ = 0.0;
        w_ = 97.76;
        p_ = 2.152e+5;
    }
//    else if (x[2] > -0.00005*(1.0-cos(pi2*x[0]/0.001))) {
    else if (x[2] > -0.0005*(1.0-cos(pi2*x[0]/0.001))*(1.0-cos(pi2*x[1]/0.001))) {
        r_ = 1.198;
        u_ = 0.0;
        v_ = 0.0;
        w_ = 0.0;
        p_ = 1.0e+5;
    }
    else {
        r_ = 6.037;
        u_ = 0.0;
        v_ = 0.0;
        w_ = 0.0;
        p_ = 1.0e+5;
    }
    u[0] = r_;
    u[1] = r_*u_;
    u[2] = r_*v_;
    u[3] = r_*w_;
    u[4] = p_/0.4+0.5*r_*(u_*u_+v_*v_+w_*w_);

    if (du) {
        for (i = 0; i < P4EST_DIM; i++) {
            du[0][i] = 0.0;
            du[1][i] = 0.0;
            du[2][i] = 0.0;
            du[3][i] = 0.0;
            du[4][i] = 0.0;
        }
    }
}

static void charm_get_midpoint (p4est_t * p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t * q, double xyz[3])
{
    p4est_qcoord_t      half_length = P4EST_QUADRANT_LEN (q->level) / 2;

    p4est_qcoord_to_vertex (p4est->connectivity, which_tree,
                            q->x + half_length, q->y + half_length, q->z + half_length,
                            xyz);
}

static void charm_init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree,
                              p4est_quadrant_t * q)
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    double              midpoint[3];
    int i;

    double du[FLD_COUNT][P4EST_DIM], u[FLD_COUNT];

    charm_get_midpoint (p4est, which_tree, q, midpoint);
    charm_initial_condition (midpoint, u, du, ctx);

    data->par.p.ro = u[0];
    data->par.p.ru = u[1];
    data->par.p.rv = u[2];
    data->par.p.rw = u[3];
    data->par.p.re = u[4];

    for (i = 0; i < P4EST_DIM; i++) {
        data->dro[i] = du[0][i];
        data->dru[i] = du[1][i];
        data->drv[i] = du[2][i];
        data->drw[i] = du[3][i];
        data->dre[i] = du[4][i];
    }
}


static double
charm_error_sqr_estimate (p4est_quadrant_t * q)
{
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    int                 i;
    double              diff2;
    double              du[P4EST_DIM];
    double              h = CHARM_GET_H(q->level);
    double              vol;

    for (i = 0; i < P4EST_DIM; i++) {
        du[i] = data->dro[i];
    }

    vol = h * h * h;

    diff2 = 0.;
    /* use the approximate derivative to estimate the L2 error */
    for (i = 0; i < P4EST_DIM; i++) {
        diff2 += du[i] * du[i] * (1. / 12.) * h * h * vol;
    }

    return diff2;
}

static int
charm_refine_flag_estimate (p4est_t * p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t * q)
{
    charm_data_t       *data = (charm_data_t *) q->p.user_data;

    if (data->ref_flag >= 4) {
        return 1;
    }
    else {
        return 0;
    }
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

    vol = h * h * h;

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

    if (q->level >= ctx->allowed_level) {
        return 0;
    }

    charm_get_midpoint (p4est, which_tree, q, mp);

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

static int
charm_coarsen_initial_condition (p4est_t * p4est,
                                 p4est_topidx_t which_tree,
                                 p4est_quadrant_t * children[])
{
    return 0;
}

static int charm_coarsen_err_estimate (p4est_t * p4est,
                                     p4est_topidx_t which_tree,
                                     p4est_quadrant_t * children[])
{
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    double              global_err = ctx->max_err;
    double              global_err2 = global_err * global_err;
    double              h;
    charm_data_t       *data;
    double              vol, err2, childerr2;
    double              parentu;
    double              diff;
    int                 i;

    if (children[0]->level <= ctx->min_level) {
        return 0;
    }


    h =     CHARM_GET_H(children[0]->level);
    /* the quadrant's volume is also its volume fraction */
    vol = h * h * h;

    /* compute the average */
    parentu = 0.;
    for (i = 0; i < P4EST_CHILDREN; i++) {
        data = (charm_data_t *) children[i]->p.user_data;
        parentu += data->par.p.ro / P4EST_CHILDREN;
    }

    err2 = 0.;
    for (i = 0; i < P4EST_CHILDREN; i++) {
        childerr2 = charm_error_sqr_estimate (children[i]);

        if (childerr2 > global_err2 * vol) {
            return 0;
        }
        err2 += childerr2;
        diff = (parentu - data->par.p.ro) * (parentu - data->par.p.ro);
        err2 += diff * vol;
    }
    if (err2 < global_err2 * (vol * P4EST_CHILDREN)) {
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
                    udata = (charm_data_t *) &ghost_data[side[j]->is.full.quadid];
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
    charm_data_t       *parent_data, *child_data;
    int                 i, j;
    double              h;
    double              du_old, du_est;

    if (num_outgoing > 1) {
        /* this is coarsening */
        parent_data = (charm_data_t *) incoming[0]->p.user_data;
        parent_data->par.p.ro = 0.;
        parent_data->par.p.ru = 0.;
        parent_data->par.p.rv = 0.;
        parent_data->par.p.rw = 0.;
        parent_data->par.p.re = 0.;

        for (j = 0; j < P4EST_DIM; j++) {
            parent_data->dro[j] = 0.0;
            parent_data->dru[j] = 0.0;
            parent_data->drv[j] = 0.0;
            parent_data->drw[j] = 0.0;
            parent_data->dre[j] = 0.0;
        }
        for (i = 0; i < P4EST_CHILDREN; i++) {
            child_data = (charm_data_t *) outgoing[i]->p.user_data;
            parent_data->par.p.ro += child_data->par.p.ro / P4EST_CHILDREN;
            parent_data->par.p.ru += child_data->par.p.ru / P4EST_CHILDREN;
            parent_data->par.p.rv += child_data->par.p.rv / P4EST_CHILDREN;
            parent_data->par.p.rw += child_data->par.p.rw / P4EST_CHILDREN;
            parent_data->par.p.re += child_data->par.p.re / P4EST_CHILDREN;
            for (j = 0; j < P4EST_DIM; j++) {
                parent_data->dro[j] += child_data->dro[j] / P4EST_CHILDREN;
                parent_data->dru[j] += child_data->dru[j] / P4EST_CHILDREN;
                parent_data->drv[j] += child_data->drv[j] / P4EST_CHILDREN;
                parent_data->drw[j] += child_data->drw[j] / P4EST_CHILDREN;
                parent_data->dre[j] += child_data->dre[j] / P4EST_CHILDREN;
            }
        }

    }
    else {
        /* this is refinement */
        parent_data = (charm_data_t *) outgoing[0]->p.user_data;
        h = CHARM_GET_H(outgoing[0]->level);

        for (i = 0; i < P4EST_CHILDREN; i++) {
            child_data = (charm_data_t *) incoming[i]->p.user_data;
            child_data->par.p.ro = parent_data->par.p.ro;
            child_data->par.p.ru = parent_data->par.p.ru;
            child_data->par.p.rv = parent_data->par.p.rv;
            child_data->par.p.rw = parent_data->par.p.rw;
            child_data->par.p.re = parent_data->par.p.re;
            for (j = 0; j < P4EST_DIM; j++) {
                child_data->dro[j] = parent_data->dro[j];
                child_data->dru[j] = parent_data->dru[j];
                child_data->drv[j] = parent_data->drv[j];
                child_data->drw[j] = parent_data->drw[j];
                child_data->dre[j] = parent_data->dre[j];
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

    p4est_balance (p4est, P4EST_CONNECT_FACE, charm_init_initial_condition);

    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    ghost_data = P4EST_ALLOC (charm_data_t, ghost->ghosts.elem_count);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   charm_ref_flag_quad_iter_fn,
                   charm_ref_flag_face_iter_fn,
                   NULL,
                   NULL);
    p4est_refine (p4est, recursive, charm_refine_flag_estimate,
                  charm_init_initial_condition);

    p4est_ghost_destroy (ghost);
    P4EST_FREE (ghost_data);
    ghost = NULL;
    ghost_data = NULL;

    partforcoarsen = 1;
    p4est_balance (p4est, P4EST_CONNECT_FACE, charm_init_initial_condition);
    p4est_partition (p4est, partforcoarsen, NULL);

}

void charm_adapt(p4est_t *p4est, p4est_ghost_t **ghost, charm_data_t **ghost_data)
{
    int             recursive       = 0;
    int             callbackorphans = 0;
    charm_ctx_t      *ctx             = (charm_ctx_t *) p4est->user_pointer;

    p4est_refine_ext (p4est, recursive, ctx->allowed_level,
                      charm_refine_err_estimate, NULL,
                      charm_replace_quads);
    p4est_coarsen_ext (p4est, recursive, callbackorphans,
                       charm_coarsen_err_estimate, NULL,
                       charm_replace_quads);
    p4est_balance_ext (p4est, P4EST_CONNECT_FACE, NULL,
                       charm_replace_quads);

    if (*ghost) {
        p4est_ghost_destroy(*ghost);
        P4EST_FREE (*ghost_data);
    }
    *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
    *ghost_data = P4EST_ALLOC (charm_data_t, (*ghost)->ghosts.elem_count);
    p4est_ghost_exchange_data (p4est, *ghost, *ghost_data);



    p4est_iterate (p4est, *ghost, (void *) *ghost_data,
                   charm_ref_flag_quad_iter_fn,
                   charm_ref_flag_face_iter_fn,
                   NULL,
                   NULL);
    p4est_refine_ext (p4est, recursive, ctx->allowed_level,
                      charm_refine_flag_estimate, NULL,
                      charm_replace_quads);

    p4est_balance_ext (p4est, P4EST_CONNECT_FACE, NULL,
                       charm_replace_quads);


    p4est_ghost_destroy (*ghost);
    P4EST_FREE (*ghost_data);
    *ghost = NULL;
    *ghost_data = NULL;


}

static void charm_interpolate_solution (p4est_iter_volume_info_t * info, void *user_data)
{
    sc_array_t         **u_interp = (sc_array_t **) user_data;      /* we passed the array of values to fill as the user_data in the call to p4est_iterate */
    p4est_t            *p4est = info->p4est;
    p4est_quadrant_t   *q = info->quad;
    p4est_topidx_t      which_tree = info->treeid;
    p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
    p4est_tree_t       *tree;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    double              h;
    p4est_locidx_t      arrayoffset;
    double              this_u[FLD_COUNT];
    double             *this_u_ptr;
    int                 i, j;
    p4est_locidx_t      fld_shift;

    fld_shift = p4est->local_num_quadrants * P4EST_CHILDREN;

    tree = p4est_tree_array_index (p4est->trees, which_tree);
    local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
    arrayoffset = P4EST_CHILDREN * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in u_interp */
    h =  CHARM_GET_H(q->level);

    for (i = 0; i < P4EST_CHILDREN; i++) {
        this_u[0] = data->par.p.ro;
        this_u[1] = data->par.p.ru;
        this_u[2] = data->par.p.rv;
        this_u[3] = data->par.p.rw;
        this_u[4] = data->par.p.re;
        /* loop over the derivative components and linearly interpolate from the
         * midpoint to the corners */
        for (j = 0; j < P4EST_DIM; j++) {
            /* In order to know whether the direction from the midpoint to the corner is
             * negative or positive, we take advantage of the fact that the corners
             * are in z-order.  If i is an odd number, it is on the +x side; if it
             * is even, it is on the -x side.  If (i / 2) is an odd number, it is on
             * the +y side, etc. */
//            this_u[0] += (h / 2) * data->dro[j] * ((i & (1 << j)) ? 1. : -1.);
//            this_u[1] += (h / 2) * data->dru[j] * ((i & (1 << j)) ? 1. : -1.);
//            this_u[2] += (h / 2) * data->drv[j] * ((i & (1 << j)) ? 1. : -1.);
//            this_u[3] += (h / 2) * data->dre[j] * ((i & (1 << j)) ? 1. : -1.);
        }
        for (j = 0; j < 5; j++) {
            this_u_ptr = (double *) sc_array_index (u_interp[j], arrayoffset + i);
            this_u_ptr[0] = this_u[j];
        }
//        u_interp[arrayoffset + i + fld_shift*0] = data->ro;
//        u_interp[arrayoffset + i + fld_shift*1] = data->ru;
//        u_interp[arrayoffset + i + fld_shift*2] = data->rv;
//        u_interp[arrayoffset + i + fld_shift*3] = data->re;
    }

}




static void charm_write_solution (p4est_t * p4est, int timestep)
{
    char                filename[BUFSIZ] = { '\0' };
    sc_array_t         **u_interp;
    p4est_locidx_t      numquads;
    int                 fld_count = 5, i;

    snprintf (filename, 33, P4EST_STRING "_charm_%08d", timestep);

    numquads = p4est->local_num_quadrants;

    /* create a vector with one value for the corner of every local quadrant
     * (the number of children is always the same as the number of corners) */
    u_interp = P4EST_ALLOC(sc_array_t*, 5);
    for (i = 0; i < 5; i++) {
        u_interp[i] = sc_array_new_size (sizeof(double), numquads * P4EST_CHILDREN);
    }

    /* Use the iterator to visit every cell and fill in the solution values.
     * Using the iterator is not absolutely necessary in this case: we could
     * also loop over every tree (there is only one tree in this case) and loop
     * over every quadrant within every tree, but we are trying to demonstrate
     * the usage of p4est_iterate in this example */
    p4est_iterate (p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                   (void *) u_interp,     /* pass in u_interp so that we can fill it */
                   charm_interpolate_solution,    /* callback function that interpolates from the cell center to the cell corners, defined above */
                   NULL,          /* there is no callback for the faces between quadrants */
#ifdef P4_TO_P8
                   NULL,          /* there is no callback for the edges between quadrants */
#endif
                   NULL);         /* there is no callback for the corners between quadrants */

    /* create VTK output context and set its parameters */
    p4est_vtk_context_t *context = p4est_vtk_context_new (p4est, filename);
    p4est_vtk_context_set_scale (context, 0.9999);  /* quadrant at almost full scale */

    /* begin writing the output files */
    context = p4est_vtk_write_header (context);
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing vtk header");

    /* do not write the tree id's of each quadrant
     * (there is only one tree in this example) */
    context = p4est_vtk_write_cell_dataf (context, 0, 1,  /* do write the refinement level of each quadrant */
                                          1,      /* do write the mpi process id of each quadrant */
                                          0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                          0,      /* there is no custom cell scalar data. */
                                          0,      /* there is no custom cell vector data. */
                                          context);       /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");

    /* write one scalar field: the solution value */
    context = p4est_vtk_write_point_dataf (context, fld_count, 0, /* write no vector fields */
                         "RO", u_interp[0],
                         "RU", u_interp[1],
                         "RV", u_interp[2],
                         "RW", u_interp[3],
                         "RE", u_interp[4],
                         context);        /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");

    const int           retval = p4est_vtk_write_footer (context);
    SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
    for (i = 0; i < 5; i++) {
        sc_array_destroy(u_interp[i]);
    }
    P4EST_FREE(u_interp);

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




void rim_orig(  double* RI, double* EI, double* PI, double* UI, double* VI, double* WI,
                double RB, double PB, double UB, double VB, double WB,
                double RE, double PE, double UE, double VE, double WE, double gam) {

    double AGAM = (gam - 1.0);
    double BGAM = (2.0 * sqrt(gam / AGAM));
    double CGAM = (1.0 / gam);
    double DGAM = (2.0 / AGAM);
    double EGAM = (AGAM / (gam + 1.0));
    double GGAM = (sqrt(gam * AGAM));
    double HGAM = (AGAM / 2.0);
    double FGAM = (3.0 * gam - 1.0);
    double OGAM = (AGAM / (2.0 * gam));
    double QGAM = (gam + 1.0);
    double PGAM = (QGAM / (2.0 * gam));
    double RGAM = (4.0 * gam);
    double SGAM = (gam * AGAM);
    double TGAM = (QGAM / 2.0);
    double UGAM = (sqrt(AGAM / gam));

    double RF, RS, EF, ES, SBL, SFL, SSL, SEL, D, FS1, F1, ZNB, PKB, ZFB, PKE, ZNE, F2, FS2, ZFE, DP, UBD, RUBD, UF, UED, RUED, US, PPE, PPB, P;
    double eps = RIM_EPS;
    double CB = sqrt(gam * PB / RB);
    double CE = sqrt(gam * PE / RE);
    double EB = CB * CB / SGAM;
    double EE = CE * CE / SGAM;
    double RCB = RB * CB;
    double RCE = RE * CE;
    double DU = UB - UE;
    if (DU < -2.0 * (CB + CE) / AGAM) {
        printf(" ATTENTION!!!  VACUUM \n");
        RF = 0.0;
        RS = 0.0;
        EF = 0.0;
        ES = 0.0;
        SBL = UB - CB;
        SFL = UB + 2.0 * CB / AGAM;
        SSL = UE - 2.0 * CE / AGAM;
        SEL = UE + CE;
    } else {
        P = (PB * RCE + PE * RCB + DU * RCB * RCE) / (RCB + RCE);
        do {

            if (P < eps) P = eps;

            PPB = P / PB;
            if (PB <= P) {
                PKB = PGAM * PPB + OGAM;
                ZNB = RCB * sqrt(PKB);
                F1 = (P - PB) / ZNB;
                FS1 = (QGAM * PPB + FGAM) / (RGAM * ZNB * PKB);
            }
            else {
                ZFB = CB * exp(log(PPB) * OGAM);
                F1 = DGAM * (ZFB - CB);
                FS1 = ZFB / (gam * P);
            }
            PPE = P / PE;
            if (PE <= P) {
                PKE = PGAM * PPE + OGAM;
                ZNE = RCE * sqrt(PKE);
                F2 = (P - PE) / ZNE;
                FS2 = (QGAM * PPE + FGAM) / (RGAM * ZNE * PKE);
            }
            else {
                ZFE = CE * exp(log(PPE) * OGAM);
                F2 = DGAM * (ZFE - CE);
                FS2 = ZFE / (gam * P);
            }
            DP = (DU - F1 - F2) / (FS1 + FS2);
            P = P + DP;
        } while (fabs(DU - F1 - F2) > eps);


        PPB = P / PB;
        PPE = P / PE;

//       ZFB=CB*PPB**OGAM;
//       ZFE=CE*PPE**OGAM;
        ZFB = CB * exp(log(PPB) * OGAM);
        ZFE = CE * exp(log(PPE) * OGAM);
        if (PB <= P) {
            D = UB - sqrt((TGAM * P + HGAM * PB) / RB);
            UBD = UB - D;
            RUBD = RB * UBD;
            RF = RUBD * RUBD / (PB - P + RUBD * UBD);
            UF = D + RUBD / RF;
            EF = P / (AGAM * RF);
            SBL = D;
            SFL = D;
        }
        else {
            EF = ZFB * ZFB / SGAM;
            UF = UB + DGAM * (CB - ZFB);
            RF = P / (AGAM * EF);
            SBL = UB - CB;
            SFL = UF - ZFB;
        }
        if (PE <= P) {
            D = UE + sqrt((TGAM * P + HGAM * PE) / RE);
            UED = UE - D;
            RUED = RE * UED;
            RS = RUED * RUED / (PE - P + RUED * UED);
            US = D + RUED / RS;
            ES = P / (AGAM * RS);
            SEL = D;
            SSL = D;
        }
        else {
            ES = ZFE * ZFE / SGAM;
            US = UE - DGAM * (CE - ZFE);
            RS = P / (AGAM * ES);
            SSL = US + ZFE;
            SEL = UE + CE;
        }
    }
//
// C     compute the interpolation value
    if (SEL<=0.0) {
        *RI= RE;
        *EI= EE;
        *UI= UE;
        *VI= VE;
        *WI= WE;
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if (SBL>=0.0) {
        *RI= RB;
        *EI= EB;
        *UI= UB;
        *VI= VB;
        *WI= WB;
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if ((SSL>=0.0)&&(SFL<=0.0)) {
        if (US>=0.0) {
            *RI= RF;
            *EI= EF;
            *UI= UF;
            *VI= VB;
            *WI= WB;
        } else {
            *RI= RS;
            *EI= ES;
            *UI= US;
            *VI= VE;
            *WI= WE;
        }
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if (SFL>0.0) {
        *UI= (UB+DGAM*GGAM*sqrt(EB))/(1+DGAM);
        *VI= VB;
        *WI= WB;
        *EI= ((*UI)*(*UI))/SGAM;
        *RI= RB*exp(log((*EI)/EB)*(1/AGAM));
    } else {
        *UI= (UE-DGAM*GGAM*sqrt(EE))/(1+DGAM);
        *VI= VE;
        *WI= WE;
        *EI= ((*UI)*(*UI))/SGAM;
        *RI= RE*exp(log((*EI)/EE)*(1/AGAM)) ;
    }

    *PI= AGAM*(*EI)*(*RI);

    return;
}


void calc_flux(double r_[2], double u_[2], double v_[2], double w_[2], double p_[2], double* qr, double* qu, double* qv, double* qw, double* qe, int which_face, int bnd)
{
#ifdef FLUX_RIM
    double ri, ei, pi, ui, vi, wi;
    switch (which_face) {
        case 0:                      /* -x side */
            if (bnd) {
                rim_orig(  &ri, &ei, &pi, &ui, &vi, &wi,
                           r_[1], p_[1], u_[1], v_[1], w_[1],
                           r_[0], p_[0], u_[0], v_[0], w_[0], GAM);
                *qr = ri*ui;
                *qu = (*qr)*ui+pi;
                *qv = (*qr)*vi;
                *qw = (*qr)*wi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*ui;

            }
            else {
                rim_orig(  &ri, &ei, &pi, &ui, &vi, &wi,
                           r_[0], p_[0], u_[0], v_[0], w_[0],
                           r_[1], p_[1], u_[1], v_[1], w_[1], GAM);
                *qr = ri*ui;
                *qu = (*qr)*ui+pi;
                *qv = (*qr)*vi;
                *qw = (*qr)*wi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*ui;
            }
            break;
        case 1:                      /* +x side */

            rim_orig(  &ri, &ei, &pi, &ui, &vi, &wi,
                       r_[0], p_[0], u_[0], v_[0], w_[0],
                       r_[1], p_[1], u_[1], v_[1], w_[1], GAM);
            *qr = ri*ui;
            *qu = (*qr)*ui+pi;
            *qv = (*qr)*vi;
            *qw = (*qr)*wi;
            *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*ui;
            break;
        case 2:                      /* -y side */
            if (bnd) {
                rim_orig(  &ri, &ei, &pi, &vi, &wi, &ui,
                           r_[1], p_[1], v_[1], w_[1], u_[1],
                           r_[0], p_[0], v_[0], w_[0], u_[0], GAM);
                *qr = ri*vi;
                *qu = (*qr)*ui;
                *qv = (*qr)*vi+pi;
                *qw = (*qr)*wi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*vi;
            }
            else {
                rim_orig(  &ri, &ei, &pi, &vi, &wi, &ui,
                           r_[0], p_[0], v_[0], w_[0], u_[0],
                           r_[1], p_[1], v_[1], w_[1], u_[1], GAM);
                *qr = ri*vi;
                *qu = (*qr)*ui;
                *qv = (*qr)*vi+pi;
                *qw = (*qr)*wi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*vi;
            }
            break;
        case 3:                      /* +y side */

            rim_orig(  &ri, &ei, &pi, &vi, &wi, &ui,
                       r_[0], p_[0], v_[0], w_[0], u_[0],
                       r_[1], p_[1], v_[1], w_[1], u_[1], GAM);
            *qr = ri*vi;
            *qu = (*qr)*ui;
            *qv = (*qr)*vi+pi;
            *qw = (*qr)*wi;
            *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*vi;
            break;
        case 4:                      /* -z side */
            if (bnd) {
                rim_orig(  &ri, &ei, &pi, &wi, &ui, &vi,
                           r_[1], p_[1], w_[1], u_[1], v_[1],
                           r_[0], p_[0], w_[0], u_[0], v_[0], GAM);
                *qr = ri*wi;
                *qu = (*qr)*ui;
                *qv = (*qr)*vi;
                *qw = (*qr)*wi+pi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*wi;
            }
            else {
                rim_orig(  &ri, &ei, &pi, &wi, &ui, &vi,
                           r_[0], p_[0], w_[0], u_[0], v_[0],
                           r_[1], p_[1], w_[1], u_[1], v_[1], GAM);
                *qr = ri*wi;
                *qu = (*qr)*ui;
                *qv = (*qr)*vi;
                *qw = (*qr)*wi+pi;
                *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*wi;
            }
            break;
        case 5:                      /* +z side */

            rim_orig(  &ri, &ei, &pi, &wi, &ui, &vi,
                       r_[0], p_[0], w_[0], u_[0], v_[0],
                       r_[1], p_[1], w_[1], u_[1], v_[1], GAM);
            *qr = ri*wi;
            *qu = (*qr)*ui;
            *qv = (*qr)*vi;
            *qw = (*qr)*wi+pi;
            *qe = (ri*(ei+0.5*(ui*ui+vi*vi+wi*wi))+pi)*wi;
            break;
    }
#else
#ifdef FLUX_LF
    double fr[2], fu[2], fv[2], fe[2];
    double ro[2], ru[2], rv[2], re[2];
    double alpha;
    double vn = 1.0;
    int i;

    for (i = 0; i < 2; i++) {
        ro[i] = r_[i];
        ru[i] = r_[i]*u_[i];
        rv[i] = r_[i]*v_[i];
        re[i] = p_[i]/(GAM-1.0)+0.5*r_[i]*(u_[i]*u_[i]+v_[i]*v_[i]);
        switch (which_face) {
            case 0:                      /* -x side */
            case 1:                      /* +x side */
                fr[i] = ru[i];
                fu[i] = fr[i]*u_[i]+p_[i];
                fv[i] = fr[i]*v_[i];
                fe[i] = (re[i]+p_[i])*u_[i];
                break;
            case 2:                      /* -y side */
            case 3:                      /* +y side */
                fr[i] = rv[i];
                fu[i] = fr[i]*u_[i];
                fv[i] = fr[i]*v_[i]+p_[i];
                fe[i] = (re[i]+p_[i])*v_[i];
                break;
        }
    }
    alpha = _MAX_(fabs(u_[0])+sqrt(GAM*p_[0]/r_[0]), fabs(u_[1])+sqrt(GAM*p_[1]/r_[1]));

    if (bnd) {
        switch (which_face) {
            case 0:                      /* -x side */
                vn = -1.0;
                break;
            case 1:                      /* +x side */
                vn = 1.0;
                break;
            case 2:                      /* -y side */
                vn = -1.0;
                break;
            case 3:                      /* +y side */
                vn = 1.0;
                break;
        }

        alpha *= vn;
    }

    *qr = 0.5*(fr[0]+fr[1]-alpha*(ro[1]-ro[0]));
    *qu = 0.5*(fu[0]+fu[1]-alpha*(ru[1]-ru[0]));
    *qv = 0.5*(fv[0]+fv[1]-alpha*(rv[1]-rv[0]));
    *qe = 0.5*(fe[0]+fe[1]-alpha*(re[1]-re[0]));
#endif // FLUX_LF
#endif // FLUX_RIM
}

static void charm_upwind_flux (p4est_iter_face_info_t * info, void *user_data)
{
    int                 i, j;
    p4est_t            *p4est = info->p4est;
    charm_ctx_t        *ctx = (charm_ctx_t *) p4est->user_pointer;
    charm_data_t       *ghost_data = (charm_data_t *) user_data;
    charm_data_t       *udata;
    p4est_quadrant_t   *quad;
    double              vdotn = 0.;
    double              ro_avg[2], ru_avg[2], rv_avg[2], rw_avg[2], re_avg[2];
    double              qr, qu, qv, qw, qe;
    double              h, facearea, dh, vdx, vdy, vdz;
    int                 which_face;
    p4est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);

    int bnd = 0;

    if (sides->elem_count != 2) {
        side[0] = p4est_iter_fside_array_index_int(sides, 0);

        bnd = 1;

        P4EST_ASSERT(info->tree_boundary);

        /* which of the quadrant's faces the interface touches */
        which_face = side[0]->face;

        switch (which_face) {
            case 0:                      /* -x side */
                vdx = -1.0;
                vdy =  0.0;
                vdz =  0.0;
                break;
            case 1:                      /* +x side */
                vdx = 1.0;
                vdy = 0.0;
                vdz = 0.0;
                break;
            case 2:                      /* -y side */
                vdx =  0.0;
                vdy = -1.0;
                vdz =  0.0;
                break;
            case 3:                      /* +y side */
                vdx = 0.0;
                vdy = 1.0;
                vdz = 0.0;
                break;
            case 4:                      /* -z side */
                vdx =  0.0;
                vdy =  0.0;
                vdz = -1.0;
                break;
            case 5:                      /* +z side */
                vdx = 0.0;
                vdy = 0.0;
                vdz = 1.0;
                break;
        }

        ro_avg[0] = 0;
        ru_avg[0] = 0;
        rv_avg[0] = 0;
        rw_avg[0] = 0;
        re_avg[0] = 0;
        if (side[0]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P4EST_HALF; j++) {
                if (side[0]->is.hanging.is_ghost[j]) {
                    udata =
                            (charm_data_t *) &ghost_data[side[0]->is.hanging.quadid[j]];
                    dh = 0.0;
                }
                else {
                    udata =
                            (charm_data_t *) side[0]->is.hanging.quad[j]->p.user_data;
#ifdef SECOND_ORDER
                    dh = 0.5 * CHARM_GET_H(side[0]->is.hanging.quad[j]->level);
#else
                    dh = 0.0;
#endif
                }
                ro_avg[0] += udata->par.p.ro+dh*(udata->dro[0]*vdx+udata->dro[1]*vdy+udata->dro[2]*vdz);
                ru_avg[0] += udata->par.p.ru+dh*(udata->dru[0]*vdx+udata->dru[1]*vdy+udata->dru[2]*vdz);
                rv_avg[0] += udata->par.p.rv+dh*(udata->drv[0]*vdx+udata->drv[1]*vdy+udata->drv[2]*vdz);
                rw_avg[0] += udata->par.p.rw+dh*(udata->drw[0]*vdx+udata->drw[1]*vdy+udata->drw[2]*vdz);
                re_avg[0] += udata->par.p.re+dh*(udata->dre[0]*vdx+udata->dre[1]*vdy+udata->dre[2]*vdz);
            }
            ro_avg[0] /= P4EST_HALF;
            ru_avg[0] /= P4EST_HALF;
            rv_avg[0] /= P4EST_HALF;
            rw_avg[0] /= P4EST_HALF;
            re_avg[0] /= P4EST_HALF;
        }
        else {
            if (side[0]->is.full.is_ghost) {
                udata = (charm_data_t *) &ghost_data[side[0]->is.full.quadid];
                dh = 0.0;
            }
            else {
                udata = (charm_data_t *) side[0]->is.full.quad->p.user_data;
#ifdef SECOND_ORDER
                dh = 0.5 * CHARM_GET_H(side[0]->is.full.quad->level);
#else
                dh = 0.0;
#endif
            }
            ro_avg[0] = udata->par.p.ro+dh*(udata->dro[0]*vdx+udata->dro[1]*vdy+udata->dro[2]*vdz);;
            ru_avg[0] = udata->par.p.ru+dh*(udata->dru[0]*vdx+udata->dru[1]*vdy+udata->dru[2]*vdz);
            rv_avg[0] = udata->par.p.rv+dh*(udata->drv[0]*vdx+udata->drv[1]*vdy+udata->drv[2]*vdz);
            rw_avg[0] = udata->par.p.rw+dh*(udata->drw[0]*vdx+udata->drw[1]*vdy+udata->drw[2]*vdz);
            re_avg[0] = udata->par.p.re+dh*(udata->dre[0]*vdx+udata->dre[1]*vdy+udata->dre[2]*vdz);
        }

        /* boundary conditions */
        double r_, p_, u_, v_, w_;
        switch (which_face) {
            case 0:                      /* -x side */
                vdotn = -1.0;

                ro_avg[1] =  ro_avg[0];
                ru_avg[1] = -ru_avg[0];
                rv_avg[1] =  rv_avg[0];
                rw_avg[1] =  rw_avg[0];
                re_avg[1] =  re_avg[0];
                break;

            case 1:                      /* +x side */
                vdotn = 1.0;

                ro_avg[1] =  ro_avg[0];
                ru_avg[1] = -ru_avg[0];
                rv_avg[1] =  rv_avg[0];
                rw_avg[1] =  rw_avg[0];
                re_avg[1] =  re_avg[0];
                break;

            case 2:                      /* -y side */
                vdotn = -1.0;

                ro_avg[1] =  ro_avg[0];
                ru_avg[1] =  ru_avg[0];
                rv_avg[1] = -rv_avg[0];
                rw_avg[1] =  rw_avg[0];
                re_avg[1] =  re_avg[0];
                break;

            case 3:                      /* +y side */
                vdotn = 1.0;

                ro_avg[1] =  ro_avg[0];
                ru_avg[1] =  ru_avg[0];
                rv_avg[1] = -rv_avg[0];
                re_avg[1] =  re_avg[0];
                break;
            case 4:                      /* -z side */
                vdotn = -1.0;

//                r_ = 12.09;
//                u_ = 0.0;
//                v_ = 97.76;
//                w_ = 0.0;
//                p_ = 2.152e+5;
//
//                ro_avg[1] = r_;
//                ru_avg[1] = r_*u_;
//                rv_avg[1] = r_*v_;
//                rw_avg[1] = r_*w_;
//                re_avg[1] = p_/(GAM-1.0)+0.5*r_*(u_*u_+v_*v_+w_*w_);


                ro_avg[1] =  ro_avg[0];
                ru_avg[1] =  ru_avg[0];
                rv_avg[1] =  rv_avg[0];
                rw_avg[1] =  rw_avg[0];
                re_avg[1] =  re_avg[0];
                break;

            case 5:                      /* +z side */
                vdotn = 1.0;

                ro_avg[1] =  ro_avg[0];
                ru_avg[1] =  ru_avg[0];
                rv_avg[1] =  rv_avg[0];
                rw_avg[1] = -rw_avg[0];
                re_avg[1] =  re_avg[0];
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
                vdx = -1.0;
                vdy =  0.0;
                vdz =  0.0;
                break;
            case 1:                      /* +x side */
                vdotn = 1.0;
                vdx =  1.0;
                vdy =  0.0;
                vdz =  0.0;
                break;
            case 2:                      /* -y side */
                vdotn = -1.0;
                vdx =  0.0;
                vdy = -1.0;
                vdz =  0.0;
                break;
            case 3:                      /* +y side */
                vdotn = 1.0;
                vdx =  0.0;
                vdy =  1.0;
                vdz =  0.0;
                break;
            case 4:                      /* -z side */
                vdotn = -1.0;
                vdx =  0.0;
                vdy =  0.0;
                vdz = -1.0;
                break;
            case 5:                      /* +z side */
                vdotn = 1.0;
                vdx =  0.0;
                vdy =  0.0;
                vdz =  1.0;
                break;
        }

        P4EST_ASSERT (vdotn == 1.0);

        for (i = 0; i < 2; i++) {
            ro_avg[i] = 0;
            ru_avg[i] = 0;
            rv_avg[i] = 0;
            rw_avg[i] = 0;
            re_avg[i] = 0;
            if (side[i]->is_hanging) {
                /* there are 2^(d-1) (P4EST_HALF) subfaces */
                for (j = 0; j < P4EST_HALF; j++) {
                    if (side[i]->is.hanging.is_ghost[j]) {
                        udata = (charm_data_t *) &ghost_data[side[i]->is.hanging.quadid[j]];
                        dh = 0.0;
                    }
                    else {
                        udata =
                                (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
#ifdef SECOND_ORDER
                        dh = 0.5 * CHARM_GET_H(side[i]->is.hanging.quad[j]->level);
#else
                        dh = 0.0;
#endif
                    }
                    ro_avg[i] += udata->par.p.ro+dh*(udata->dro[0]*vdx+udata->dro[1]*vdy+udata->dro[2]*vdz);
                    ru_avg[i] += udata->par.p.ru+dh*(udata->dru[0]*vdx+udata->dru[1]*vdy+udata->dru[2]*vdz);
                    rv_avg[i] += udata->par.p.rv+dh*(udata->drv[0]*vdx+udata->drv[1]*vdy+udata->drv[2]*vdz);
                    rw_avg[i] += udata->par.p.rw+dh*(udata->drw[0]*vdx+udata->drw[1]*vdy+udata->drw[2]*vdz);
                    re_avg[i] += udata->par.p.re+dh*(udata->dre[0]*vdx+udata->dre[1]*vdy+udata->dre[2]*vdz);
                }
                ro_avg[i] /= P4EST_HALF;
                ru_avg[i] /= P4EST_HALF;
                rv_avg[i] /= P4EST_HALF;
                rw_avg[i] /= P4EST_HALF;
                re_avg[i] /= P4EST_HALF;
            }
            else {
                if (side[i]->is.full.is_ghost) {
                    udata = (charm_data_t *) &ghost_data[side[i]->is.full.quadid];
                    dh = 0.0;
                }
                else {
                    udata = (charm_data_t *) side[i]->is.full.quad->p.user_data;
#ifdef SECOND_ORDER
                    dh = 0.5 * CHARM_GET_H(side[0]->is.full.quad->level);
#else
                    dh = 0.0;
#endif
                }
                ro_avg[i] = udata->par.p.ro+dh*(udata->dro[0]*vdx+udata->dro[1]*vdy+udata->dro[2]*vdz);
                ru_avg[i] = udata->par.p.ru+dh*(udata->dru[0]*vdx+udata->dru[1]*vdy+udata->dru[2]*vdz);
                rv_avg[i] = udata->par.p.rv+dh*(udata->drv[0]*vdx+udata->drv[1]*vdy+udata->drv[2]*vdz);
                rw_avg[i] = udata->par.p.rw+dh*(udata->drw[0]*vdx+udata->drw[1]*vdy+udata->drw[2]*vdz);
                re_avg[i] = udata->par.p.re+dh*(udata->dre[0]*vdx+udata->dre[1]*vdy+udata->dre[2]*vdz);
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
    calc_flux(r_, u_, v_, w_, p_, &qr, &qu, &qv, &qw, &qe, which_face, bnd);

    for (i = 0; i < sides->elem_count; i++) {
        if (side[i]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P4EST_HALF; j++) {
                quad = side[i]->is.hanging.quad[j];
                h = CHARM_GET_H(quad->level);
                facearea = h * h;
                if (!side[i]->is.hanging.is_ghost[j]) {
                    udata = (charm_data_t *) quad->p.user_data;
                    udata->drodt += qr * facearea * (i ? 1. : -1.) * vdotn;
                    udata->drudt += qu * facearea * (i ? 1. : -1.) * vdotn;
                    udata->drvdt += qv * facearea * (i ? 1. : -1.) * vdotn;
                    udata->drwdt += qw * facearea * (i ? 1. : -1.) * vdotn;
                    udata->dredt += qe * facearea * (i ? 1. : -1.) * vdotn;
                }
            }
        }
        else {
            quad = side[i]->is.full.quad;
            h = CHARM_GET_H(quad->level);
            facearea = h * h;
            if (!side[i]->is.full.is_ghost) {
                udata = (charm_data_t *) quad->p.user_data;
                udata->drodt += qr * facearea * (i ? 1. : -1.) * vdotn;
                udata->drudt += qu * facearea * (i ? 1. : -1.) * vdotn;
                udata->drvdt += qv * facearea * (i ? 1. : -1.) * vdotn;
                udata->drwdt += qw * facearea * (i ? 1. : -1.) * vdotn;
                udata->dredt += qe * facearea * (i ? 1. : -1.) * vdotn;
            }
        }
    }
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

static void
charm_grad_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
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

        P4EST_ASSERT (vdotn == 1.0);

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



/** Timestep the advection problem.
 *
 * Update the state, refine, repartition, and write the solution to file.
 *
 * \param [in,out] p4est the forest, whose state is updated
 * \param [in] time      the end time
 */
static void charm_timestep (p4est_t * p4est, double time)
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







void charm_init_context(charm_ctx_t *ctx)
{
    ctx->max_err                = 1.e-2;

    ctx->refine_period          = 10;
    ctx->repartition_period     = 100;

    ctx->min_level              = 1;
    ctx->allowed_level          = 1;

    ctx->write_period           = 100;

    ctx->v_ref                  = 20.;
    ctx->CFL                    = 0.1;

    ctx->dt                     = 1.e-8;

    ctx->time                   = 3.5e-3;
}


int main (int argc, char **argv)
{

    int                 mpiret;
    sc_MPI_Comm         mpicomm;
    p4est_t            *p4est;
    p4est_connectivity_t *conn;
    charm_ctx_t         ctx;

    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;

    sc_init (mpicomm, 0, 0, NULL, SC_LP_ESSENTIAL);
    p4est_init (NULL, SC_LP_ESSENTIAL);


    charm_init_context(&ctx);


    /* Create a forest that consists of just one periodic quadtree/octree. */
//    conn = charm_conn_create_n(40, 400);
    conn = charm_conn_create(&ctx);


    P4EST_ASSERT(p4est_connectivity_is_valid(conn));

    p4est = p4est_new_ext (mpicomm,              /* communicator */
                           conn,                 /* connectivity */
                           0,                    /* minimum quadrants per MPI process */
                           ctx.min_level,        /* minimum level of refinement */
                           1,                    /* fill uniform */
                           sizeof (charm_data_t),         /* data size */
                           charm_init_initial_condition,  /* initializes data */
                           (void *) (&ctx));            /* context */

    charm_write_solution (p4est, 0);
    charm_calc_grad(p4est, NULL, NULL);

    charm_adapt_init(p4est);



    charm_write_solution (p4est, 0);



    /* time step */
    charm_timestep (p4est, ctx.time);

    /* Destroy the p4est and the connectivity structure. */
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);

    /* Verify that allocations internal to p4est and sc do not leak memory.
     * This should be called if sc_init () has been called earlier. */
    sc_finalize ();

    /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
    return 0;
}
