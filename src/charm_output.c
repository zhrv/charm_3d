//
// Created by zhrv on 26.10.17.
//

#include "charm_output.h"
#include "charm_base_func.h"


static void charm_interpolate_cell_solution (p4est_iter_volume_info_t * info, void *user_data)
{
    sc_array_t         **u_interp = (sc_array_t **) user_data;      /* we passed the array of values to fill as the user_data in the call to p4est_iterate */
    p4est_t            *p4est = info->p4est;
    charm_ctx_t        *ctx = charm_get_ctx(p4est);
    p4est_quadrant_t   *q = info->quad;
    p4est_topidx_t      which_tree = info->treeid;
    p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
    p4est_tree_t       *tree;
    charm_data_t       *data = (charm_data_t *) q->p.user_data;
    size_t              c_count = charm_get_comp_count(p4est);
    double              this_u[9+c_count];
    double             *this_u_ptr;
    int                 j;
    charm_cons_t        cons;
    charm_prim_t        prim;

    tree = p4est_tree_array_index (p4est->trees, which_tree);
    local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */

    charm_tree_attr_t * attr = (charm_tree_attr_t *)&(p4est->connectivity->tree_to_attr[which_tree*sizeof(charm_tree_attr_t)]);
    charm_get_fields(data, data->par.g.c, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);

    this_u[0] = prim.r;
    this_u[1] = prim.p0;
    this_u[2] = prim.p;
    this_u[3] = prim.p0 + prim.p;
    this_u[4] = prim.e;
    this_u[5] = prim.e+_MAG_(prim.u, prim.v, prim.w);
    this_u[6] = prim.u;
    this_u[7] = prim.v;
    this_u[8] = prim.w;
    for (j = 0; j < c_count; j++) {
        this_u[9+j] = prim.c[j];
    }
    for (j = 0; j < 9 + c_count; j++) {
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
    int                 num_cell_scalars;
    int                 num_cell_vectors;
    charm_ctx_t        *ctx;
    char*         names9[] = {"R", "P_0", "P_1", "P", "E", "E_TOT", "U", "V", "W"};
    charm_comp_t       *comp;

    ctx = charm_get_ctx(p4est);
    snprintf (filename, 33, CHARM_STRING "_%08d", timestep);

    numquads = (size_t)p4est->local_num_quadrants;

    num_cell_vectors = 0;
    num_cell_scalars = 9 + (int)ctx->comp->elem_count;
    u_interp = CHARM_ALLOC(sc_array_t*, num_cell_scalars + num_cell_vectors);
    for (i = 0; i < num_cell_scalars + num_cell_vectors; i++) {
        u_interp[i] = sc_array_new_size (sizeof(double), numquads);
    }

    char** names = CHARM_ALLOC (char *, num_cell_scalars + num_cell_vectors);
    for (i = 0; i < 9; i++) {
        names[i] = names9[i];
    }
    for (i = 0; i < ctx->comp->elem_count; i++) {
        comp = charm_get_comp(p4est, i);
        names[i+9] = CHARM_ALLOC (char, 128);
        strcpy(names[i+9], "C_");
        strcat(names[i+9], comp->name);
    }

    p4est_iterate (p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                   (void *) u_interp,     /* pass in u_interp so that we can fill it */
                   charm_interpolate_cell_solution,    /* callback function that interpolates from the cell center to the cell corners, defined above */
                   NULL,          /* there is no callback for the faces between quadrants */
                    NULL,          /* there is no callback for the edges between quadrants */
                   NULL);         /* there is no callback for the corners between quadrants */

    charm_vtk_context_t *context = charm_vtk_context_new (p4est, filename);
    charm_vtk_context_set_scale (context, 1.);  /* quadrant at almost full scale */
    charm_vtk_context_set_continuous (context, 1);

    context = charm_vtk_write_header (context);
    SC_CHECK_ABORT (context != NULL,
                    CHARM_STRING "_vtk: Error writing vtk header");
    context = charm_vtk_write_cell_data (context, 1, 1,      /* do write the refinement level of each quadrant */
                                          1,                  /* do write the mpi process id of each quadrant */
                                          0,                  /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                          num_cell_scalars,                  /* there is no custom cell scalar data. */
                                          num_cell_vectors,
                                         (const char**)names,
                                          u_interp);           /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    CHARM_STRING "_vtk: Error writing cell data");

    const int           retval = charm_vtk_write_footer (context);
    SC_CHECK_ABORT (!retval, CHARM_STRING "_vtk: Error writing footer");
    for (i = 0; i < 5; i++) {
        sc_array_destroy(u_interp[i]);
    }
    CHARM_FREE(u_interp);
    for (i = 0; i < ctx->comp->elem_count; i++) {
        CHARM_FREE (names[i+9]);
    }

}


void charm_log_statistics(p4est_t * p4est, int timestep, double time, double dt, double calc_time)
{
    CHARM_GLOBAL_ESSENTIALF(" STEP = %8d, TIME = %16.8e , DT = %16.8e, ELAPSED_TIME = %16.6e \n", timestep, time, dt, calc_time);
}
