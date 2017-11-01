//
// Created by zhrv on 26.10.17.
//

#include "charm_output.h"


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
        this_u[0] = data->par.c.ro;
        this_u[1] = data->par.c.ru;
        this_u[2] = data->par.c.rv;
        this_u[3] = data->par.c.rw;
        this_u[4] = data->par.c.re;
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



static void charm_interpolate_cell_solution (p4est_iter_volume_info_t * info, void *user_data)
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
    double              this_u[7];
    double             *this_u_ptr;
    int                 i, j;
    p4est_locidx_t      fld_shift;

    fld_shift = p4est->local_num_quadrants * P4EST_CHILDREN;

    tree = p4est_tree_array_index (p4est->trees, which_tree);
    local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */

    charm_tree_attr_t * attr = (charm_tree_attr_t *)&(p4est->connectivity->tree_to_attr[which_tree*sizeof(charm_tree_attr_t)]);
    charm_param_cons_to_prim(attr->reg->mat, &(data->par));

    this_u[0] = data->par.p.r;
    this_u[1] = data->par.p.p;
    this_u[2] = data->par.p.e;
    this_u[3] = data->par.p.e_tot;
    this_u[4] = data->par.p.u;
    this_u[5] = data->par.p.v;
    this_u[6] = data->par.p.w;
    /* loop over the derivative components and linearly interpolate from the
     * midpoint to the corners */
    for (j = 0; j < 7; j++) {
        this_u_ptr = (double *) sc_array_index (u_interp[j], local_id);
        this_u_ptr[0] = this_u[j];
    }
}




void charm_write_solution (p4est_t * p4est, int timestep)
{
    char                filename[BUFSIZ] = { '\0' };
    sc_array_t        **u_interp;
    size_t              numquads;
    int                 i;

    snprintf (filename, 33, P4EST_STRING "_charm_%08d", timestep);

    numquads = p4est->local_num_quadrants;

    /* create a vector with one value for the corner of every local quadrant
     * (the number of children is always the same as the number of corners) */
    u_interp = P4EST_ALLOC(sc_array_t*, 7);
    for (i = 0; i < 7; i++) {
        u_interp[i] = sc_array_new_size (sizeof(double), numquads);
    }

    /* Use the iterator to visit every cell and fill in the solution values.
     * Using the iterator is not absolutely necessary in this case: we could
     * also loop over every tree (there is only one tree in this case) and loop
     * over every quadrant within every tree, but we are trying to demonstrate
     * the usage of p4est_iterate in this example */
    p4est_iterate (p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                   (void *) u_interp,     /* pass in u_interp so that we can fill it */
                   charm_interpolate_cell_solution,    /* callback function that interpolates from the cell center to the cell corners, defined above */
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
    context = p4est_vtk_write_cell_dataf (context, 1, 1,  /* do write the refinement level of each quadrant */
                                          1,      /* do write the mpi process id of each quadrant */
                                          0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                          7,      /* there is no custom cell scalar data. */
                                          0,      /* there is no custom cell vector data. */
                                          "R",     u_interp[0],
                                          "P",     u_interp[1],
                                          "E",     u_interp[2],
                                          "E_TOT", u_interp[3],
                                          "U",     u_interp[4],
                                          "V",     u_interp[5],
                                          "W",     u_interp[6],
                                          context);       /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");

    /* write one scalar field: the solution value */
//    context = p4est_vtk_write_point_dataf (context, fld_count, 0, /* write no vector fields */
//                                           "RO", u_interp[0],
//                                           "RU", u_interp[1],
//                                           "RV", u_interp[2],
//                                           "RW", u_interp[3],
//                                           "RE", u_interp[4],
//                                           context);        /* mark the end of the variable cell data. */
//    SC_CHECK_ABORT (context != NULL,
//                    P4EST_STRING "_vtk: Error writing cell data");

    const int           retval = p4est_vtk_write_footer (context);
    SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
    for (i = 0; i < 5; i++) {
        sc_array_destroy(u_interp[i]);
    }
    P4EST_FREE(u_interp);

}



