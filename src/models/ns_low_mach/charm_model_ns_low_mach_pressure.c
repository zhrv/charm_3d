//
// Created by zhrv on 09.10.19.
//
#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"

void charm_model_ns_low_mach_pressure_div_vel(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_low_mach_pressure_grad(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_low_mach_pressure_integrals(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);

static void _charm_model_ns_low_mach_pressure_update_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    double              dt = *((double *) user_data);
    int                 i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.p[i] += _NORM_(dt * data->int_pi[i]);
    }
}


static void _charm_model_ns_low_mach_pressure_zero_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->int_pi[i] = data->div_vel[i];
    }
}


static void _charm_model_ns_low_mach_pressure_correct_vel_quad_iter_fn (p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    double              dt = *((double *) user_data);
    int                 i, j;

    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        data->par.c.ru[i] -= _NORM_(dt * data->par.grad_p.x[i]);
        data->par.c.rv[i] -= _NORM_(dt * data->par.grad_p.y[i]);
        data->par.c.rw[i] -= _NORM_(dt * data->par.grad_p.z[i]);
    }
}


void charm_model_ns_low_mach_pressure(p4est_t * p4est, double *dt, p4est_ghost_t * ghost, charm_data_t * ghost_data) {
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    int i;
    charm_model_ns_low_mach_pressure_div_vel(p4est, ghost, ghost_data);

    for (i = 0; i < ctx->pressure_iterations; i++) {
        p4est_iterate (p4est,
                       ghost,
                       (void *) ghost_data,
                       _charm_model_ns_low_mach_pressure_zero_quad_iter_fn,
                       NULL, NULL, NULL);
        charm_model_ns_low_mach_pressure_grad(p4est, ghost, ghost_data);
        charm_model_ns_low_mach_pressure_integrals(p4est, ghost, ghost_data);
        p4est_iterate(p4est, NULL,
                      (void *) &(ctx->pressure_tau),
                      _charm_model_ns_low_mach_pressure_update_quad_iter_fn,
                      NULL, NULL, NULL);
        p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }

    p4est_iterate(p4est, NULL,
                  (void *) dt,
                  _charm_model_ns_low_mach_pressure_correct_vel_quad_iter_fn,
                  NULL, NULL, NULL);
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);
}