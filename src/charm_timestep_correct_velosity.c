//
// Created by appmath on 14.10.18.
//

#include "charm_timestep_correct_velosity.h"
#include "charm_globals.h"


static void charm_timestep_correct_velosity_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i;
    double              dt = *((double *) user_data);

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ru[i] -= dt*data->par.c.grad_p[0][i];
        data->par.c.rv[i] -= dt*data->par.c.grad_p[1][i];
        data->par.c.rw[i] -= dt*data->par.c.grad_p[2][i];
    }
}


void charm_timestep_correct_velosity(p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data, double *dt)
{
    p4est_iterate (p4est,
                  NULL,
                   (void *) dt,
                   charm_timestep_correct_velosity_iter_fn,
                   NULL, NULL, NULL);

    p4est_ghost_exchange_data (p4est, ghost, ghost_data);
}