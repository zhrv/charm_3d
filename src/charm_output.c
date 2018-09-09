//
// Created by zhrv on 26.10.17.
//

#include "charm_output.h"
#include "charm_base_func.h"


static void charm_interpolate_cell_solution (p4est_iter_volume_info_t * info, void *user_data)
{
    sc_array_t         **u_interp = (sc_array_t **) user_data;      /* we passed the array of values to fill as the user_data in the call to p4est_iterate */
    p4est_t            *p4est = info->p4est;
    p4est_quadrant_t   *q = info->quad;
    p4est_topidx_t      which_tree = info->treeid;
    p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
    p4est_tree_t       *tree;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    double              this_u[7];
    double             *this_u_ptr;
    int                 j;
    charm_cons_t        cons;
    charm_primitive_t   prim;

    tree = p4est_tree_array_index (p4est->trees, which_tree);
    local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */

    charm_tree_attr_t * attr = (charm_tree_attr_t *)&(p4est->connectivity->tree_to_attr[which_tree*sizeof(charm_tree_attr_t)]);
    charm_get_fields(q, data->par.g.c, &cons);
    /*
     * @todo TODO TODO TODO
     */
    memset(&cons, 0, sizeof(cons));
    cons.ro = data->par.c.ro[0];
    cons.ru = data->par.c.ru[0];
    cons.rv = data->par.c.rv[0];
    cons.rw = data->par.c.rw[0];
    cons.re = data->par.c.re[0];
    charm_param_cons_to_prim(attr->reg->mat, &prim, &cons);

    if (prim.r != prim.r) {
        int iiii=0;
    }
    this_u[0] = prim.r;
    this_u[1] = prim.p;
    this_u[2] = prim.e;
    this_u[3] = prim.e_tot;
    this_u[4] = prim.u;
    this_u[5] = prim.v;
    this_u[6] = prim.w;
    for (j = 0; j < 7; j++) {
        this_u_ptr = (double *) sc_array_index (u_interp[j], (size_t)local_id);
        this_u_ptr[0] = this_u[j];
    }
}




void charm_write_solution (p4est_t * p4est, int timestep)
{
    char                filename[BUFSIZ] = { '\0' };
    sc_array_t        **u_interp;
    size_t              numquads;
    int                 i;

    snprintf (filename, 33, CHARM_STRING "_%08d", timestep);

    numquads = (size_t)p4est->local_num_quadrants;

    u_interp = P4EST_ALLOC(sc_array_t*, 7);
    for (i = 0; i < 7; i++) {
        u_interp[i] = sc_array_new_size (sizeof(double), numquads);
    }

    p4est_iterate (p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                   (void *) u_interp,     /* pass in u_interp so that we can fill it */
                   charm_interpolate_cell_solution,    /* callback function that interpolates from the cell center to the cell corners, defined above */
                   NULL,          /* there is no callback for the faces between quadrants */
                    NULL,          /* there is no callback for the edges between quadrants */
                   NULL);         /* there is no callback for the corners between quadrants */

    p4est_vtk_context_t *context = p4est_vtk_context_new (p4est, filename);
    p4est_vtk_context_set_scale (context, 1.);  /* quadrant at almost full scale */
    p4est_vtk_context_set_continuous (context, 1);

    context = p4est_vtk_write_header (context);
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing vtk header");

    context = p4est_vtk_write_cell_dataf (context, 1, 1,      /* do write the refinement level of each quadrant */
                                          1,                  /* do write the mpi process id of each quadrant */
                                          0,                  /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                          7,                  /* there is no custom cell scalar data. */
                                          0,                  /* there is no custom cell vector data. */
                                          "R",     u_interp[0],
                                          "P",     u_interp[1],
                                          "E",     u_interp[2],
                                          "E_TOT", u_interp[3],
                                          "U",     u_interp[4],
                                          "V",     u_interp[5],
                                          "W",     u_interp[6],
                                          context);           /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");

    const int           retval = p4est_vtk_write_footer (context);
    SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
    for (i = 0; i < 5; i++) {
        sc_array_destroy(u_interp[i]);
    }
    P4EST_FREE(u_interp);

}


void charm_log_statistics(p4est_t * p4est, int timestep, double time)
{
    CHARM_GLOBAL_ESSENTIAL ("****************************************************\n");
    CHARM_GLOBAL_ESSENTIALF("* STEP = %8d, TIME = %24.16e *\n", timestep, time);
    CHARM_GLOBAL_ESSENTIAL ("****************************************************\n\n");
}
