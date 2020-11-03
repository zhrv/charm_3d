//
// Created by zhrv on 03.11.2020.
//
#include <fstream>
#include "charm_globals.h"
#include "yaml-cpp/yaml.h"

void charm_write_data(p4est_t *p4est)
{
    char                fname[BUFSIZ] = { '\0' };
    charm_ctx_t        *ctx  = charm_get_ctx(p4est);
    std::string         str;
    try {
        snprintf (fname, 64, CHARM_STRING "_%08d.chrm", ctx->timestep);
        p4est_save(fname, p4est, 1);

        YAML::Node config = YAML::LoadFile("task.yaml");
        config["control"]["STEP_START"] = ctx->timestep;
        config["control"]["TIME_START"] = ctx->t;
        std::ofstream fout("task.yaml");
        fout << config;
    }
    catch(YAML::Exception &e) {
        std::cout << "Error (YAML): " << e.msg << std::endl;
        exit(1);
    }
}


p4est_t* charm_load_data(charm_ctx_t *ctx, p4est_connectivity_t **conn)
{
    p4est_t            *p4est;
    char                fname[BUFSIZ] = { '\0' };
    snprintf (fname, 64, CHARM_STRING "_%08d.chrm", ctx->timestep);
    p4est = p4est_load(fname, sc_MPI_COMM_WORLD, sizeof(charm_data_t), 1, ctx, conn);
    CHARM_ESSENTIALF("Restarting from file '%s':\n\t\ttime = %12.4e\n\t\tstep = %12d\n",
                     fname, ctx->t, ctx->timestep);
    return p4est;
}
