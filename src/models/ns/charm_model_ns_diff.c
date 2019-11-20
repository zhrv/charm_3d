//
// Created by zhrv on 27.08.19.
//

#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"

void charm_model_ns_timestep_diff_grad(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
void charm_model_ns_timestep_diff_integrals(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void _charm_model_ns_timestep_diffusion_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_cons_t        cons;
    charm_prim_t        prim;
    charm_comp_t       *ci, *cj;
    int                 c_count = charm_get_comp_count(info->p4est);
    int                 i, j;
    double              s_xd, dij, pabs, td, wd, *x, s;

    x = CHARM_ALLOC(double, c_count);
    charm_get_fields(data, data->par.g.c, &cons);
    charm_param_cons_to_prim(info->p4est, &prim, &cons);
    pabs = prim.p/101325.0;
    memset(data->par.model.ns.d, 0, sizeof(double)*CHARM_MAX_COMPONETS_COUNT);

    s = 0.;
    for (i = 0; i < c_count; i++) {
        ci = charm_get_comp(info->p4est, i);
        x[i] = prim.c[i]/ci->m;
        s += x[i];
    }
    for (i = 0; i < c_count; i++) {
        x[i] /= s;
    }
    for (i = 0; i < c_count; i++) {
        if (fabs(prim.c[i]-1.) <= CHARM_EPS) {
            data->par.model.ns.d[i] = 0.;
        }
        else {
            ci = charm_get_comp(info->p4est, i);
            s_xd = 0.;
            for (j = 0; j < c_count; j++) {
                cj = charm_get_comp(info->p4est, j);
                td = prim.t / sqrt(ci->ek * cj->ek);
                wd = 1.06036 / pow(td, 0.1561) + 0.1930 / exp(0.47635 * td) + 1.03587 / exp(1.52996 * td) +
                     1.76474 / exp(3.89411 * td);
                dij = (2.628e-7) * sqrt(0.5 * (1.0 / (1000.0 * ci->m) + 1.0 / (1000.0 * cj->m))) *
                      sqrt(prim.t * prim.t * prim.t) /
                      (pabs * pow((0.5 * (ci->sig + cj->sig)), 2.0) * wd);
                s_xd += i == j ? 0. : x[j] / dij;
            }
            data->par.model.ns.d[i] = (1. - x[i]) / s_xd;
        }
    }
    CHARM_FREE(x);
}


void charm_model_ns_timestep_diffusion(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    size_t c_count = charm_get_comp_count(p4est);
    if (c_count < 2) return;


    p4est_iterate (p4est, NULL, NULL,
                   _charm_model_ns_timestep_diffusion_quad_iter_fn,
                   NULL, NULL, NULL);

}


void charm_model_ns_timestep_diff(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    if (!ctx->model.ns.use_visc) return;
    charm_model_ns_timestep_diffusion(p4est, ghost, ghost_data);
    charm_model_ns_timestep_diff_grad(p4est, ghost, ghost_data);
    charm_model_ns_timestep_diff_integrals(p4est, ghost, ghost_data);
}
