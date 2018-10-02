#include "charm_globals.h"
#include "charm_connectivity.h"
#include "charm_timestep.h"
#include "charm_init.h"
#include "charm_amr.h"


int main (int argc, char **argv)
{
    int                   mpiret;
    sc_MPI_Comm           mpicomm;
    p4est_t              *p4est;
    p4est_connectivity_t *conn;
    charm_ctx_t           ctx;

    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;
    sc_init (mpicomm, 0, 0, NULL, CHARM_LOG_LEVEL);
    p4est_init (NULL, CHARM_LOG_LEVEL);

    charm_package_id = sc_package_register(NULL, CHARM_LOG_LEVEL, CHARM_STRING, "Chemistry on AMR");

    CHARM_GLOBAL_ESSENTIAL("charm_3d started...\n");

    charm_init_context(&ctx);

    conn = charm_conn_create(&ctx);
    if (!conn) {
        sc_finalize ();

        /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
        mpiret = sc_MPI_Finalize ();
        return 1;
    }
    CHARM_ASSERT(p4est_connectivity_is_valid(conn));

    p4est = p4est_new_ext (mpicomm,              /* communicator */
                           conn,                 /* connectivity */
                           0,                    /* minimum quadrants per MPI process */
                           ctx.min_level,        /* minimum level of refinement */
                           1,                    /* fill uniform */
                           sizeof (charm_data_t),         /* data size */
                           charm_init_initial_condition,  /* initializes data */
                           (void *) (&ctx));            /* context */

    charm_write_solution (p4est, 0);
//    charm_calc_grad(p4est, NULL, NULL);

    charm_adapt_init(p4est);



   // charm_write_solution (p4est, 0);



    /* time step */
    charm_timesteps (p4est, ctx.time);

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
