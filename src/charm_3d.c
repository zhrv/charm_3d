#include "charm_globals.h"
#include "charm_amr.h"

void charm_init_context_yaml(charm_ctx_t *ctx);

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

    charm_init_context_yaml(&ctx);

    conn = charm_conn_create(&ctx);
    if (!conn) {
        charm_abort(NULL, 1);
    }
    CHARM_ASSERT(p4est_connectivity_is_valid(conn));

    p4est = p4est_new_ext (mpicomm,                         /* communicator */
                           conn,                            /* connectivity */
                           0,                               /* minimum quadrants per MPI process */
                           ctx.min_level,                   /* minimum level of refinement */
                           1,                               /* fill uniform */
                           sizeof (charm_data_t),           /* data size */
                           charm_init_initial_condition,    /* initializes data */
                           (void *) (&ctx));                /* context */
    charm_set_p4est(p4est);

    charm_write_solution (p4est);
    charm_adapt_init(p4est);

    charm_timesteps (p4est);

    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
    return 0;
}
