//
// Created by zhrv on 08.11.19.
//
#include <p8est_iterate.h>
#include <charm_globals.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"
#include "charm_eos.h"


static charm_real_t *chem_k, *chem_w, *chem_phi, *chem_psi, *chem_c, *chem_c_, *chem_c_hat;


static void _charm_model_ns_chem_rhs(p4est_t * p4est, charm_data_t *data)
{
    charm_real_t       R   = charm_eos_get_r();
    charm_cons_t cons;
    charm_prim_t prim;
    size_t r_count = charm_get_reactions_count(p4est);
    size_t c_count = charm_get_comp_count(p4est);
    charm_reaction_t * r;
    charm_comp_t *comp;
    int i, j, m, n, ic;
    charm_real_t w;

    charm_get_fields_avg(data, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);
    for (n = 0; n < r_count; n++) {
        r = charm_get_reaction(p4est, n);
        chem_k[n] = r->a*pow(prim.t, r->n)*exp(-r->e/prim.t/R);
        chem_w[n] = chem_k[n];
        for (m = 0; m < 3; m++) {
            i = r->left_comps[m];
            chem_w[n] *= (i < 0 ? 1. : chem_c_[i]);
        }
    }

    for (i = 0; i < c_count; i++) {
        chem_phi[i] = 0.;
        chem_psi[i] = 0.;
        for (n = 0; n < r_count; n++) {
            r = charm_get_reaction(p4est, n);
            for (m = 0; m < 3; m++) {
                if (i == r->left_comps[m]) {
                    w = chem_k[n];
                    for (j = 1; j < 3; j++) {
                        ic = r->left_comps[(m+j) % 3];
                        w *= (ic < 0 ? 1. : chem_c_[ic]);
                    }
                    chem_phi[i] += w;
                }
                if (i == r->right_comps[m]) {
                    chem_psi[i] += chem_w[n];
                }
            }
        }
    }

}


static void _charm_model_ns_chem_rhs_save(p4est_t * p4est, charm_data_t *data)
{
    charm_real_t R = charm_eos_get_r();
    charm_cons_t cons;
    charm_prim_t prim;
    size_t r_count = charm_get_reactions_count(p4est);
    size_t c_count = charm_get_comp_count(p4est);
    charm_reaction_t * r;
    charm_comp_t *comp;
    int i, j, m, n, ic;
    charm_real_t w, h;

    charm_get_fields_avg(data, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);
    for (n = 0; n < r_count; n++) {
        r = charm_get_reaction(p4est, n);
        chem_k[n] = r->a*pow(prim.t, r->n)*exp(-r->e/prim.t/R);
        chem_w[n] = chem_k[n];
        for (m = 0; m < 3; m++) {
            i = r->left_comps[m];
            chem_w[n] *= (i < 0 ? 1. : chem_c_[i]);
        }
    }

    data->par.model.ns.chem_rhs = 0.;
    for (i = 0; i < c_count; i++) {
        w = 0.;
        for (n = 0; n < r_count; n++) {
            r = charm_get_reaction(p4est, n);
            for (m = 0; m < 3; m++) {
                if (i == r->left_comps[m]) {
                    w -= chem_w[i];
                }
                if (i == r->right_comps[m]) {
                    w += chem_w[n];
                }
            }
        }
        comp = charm_get_comp(p4est, i);
        data->par.model.ns.chem_rhs += w*charm_comp_calc_enthalpy(comp, prim.t);
    }

}


static void _charm_model_ns_chem_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_t *p4est = info->p4est;
    charm_data_t *data = charm_get_quad_data(info->quad);
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    size_t c_count = charm_get_comp_count(p4est);
    charm_cons_t cons;
    charm_prim_t prim;
    charm_comp_t *comp;
    charm_reaction_t *r;
    int i, j;

    charm_get_fields_avg(data, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);
    for (i = 0; i < c_count; i++) {
        comp = charm_get_comp(p4est, i);
        chem_c[i] = (cons.rc[i] / comp->m);
        chem_c_hat[i] = chem_c[i];
    }

    for (j = 0; j < 2; j++) {
        for (i = 0; i < c_count; i++) {
            chem_c_[i] = (chem_c_hat[i] + chem_c[i]) * 0.5;
        }
        _charm_model_ns_chem_rhs(p4est, data);

        for (i = 0; i < c_count; i++) {
            chem_c_hat[i]  = (chem_c[i] + ctx->dt * chem_psi[i] * (1. + ctx->dt * chem_phi[i] * 0.5));
            chem_c_hat[i] /= (1. + ctx->dt * chem_phi[i] + _SQR_(ctx->dt * chem_phi[i]) * 0.5);
        }
    }

    for (i = 0; i < c_count; i++) {
        chem_c_[i] = chem_c_hat[i];
    }
    _charm_model_ns_chem_rhs_save(p4est, data);

    // TODO
    for (i = 0; i < c_count; i++) {
        comp = charm_get_comp(p4est, i);
        memset(data->par.c.rc[i], 0, sizeof(charm_real_t)*CHARM_BASE_FN_COUNT);
        data->par.c.rc[i][0] = chem_c_hat[i] * comp->m;
    }


}


static void _charm_model_ns_chem_zero_rhs_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_get_quad_data(info->quad)->par.model.ns.chem_rhs = 0.;
}


void charm_model_ns_timestep_chem(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    if (!ctx->reactions) {
        p4est_iterate (p4est,
                       ghost,
                       (void *) ghost_data,
                       _charm_model_ns_chem_zero_rhs_quad_iter_fn,
                       NULL, NULL, NULL);

        return;
    }

    size_t c_count = charm_get_comp_count(p4est);
    size_t r_count = charm_get_reactions_count(p4est);

    chem_k     = CHARM_ALLOC(charm_real_t, r_count);
    chem_w     = CHARM_ALLOC(charm_real_t, r_count);
    chem_phi   = CHARM_ALLOC(charm_real_t, c_count);
    chem_psi   = CHARM_ALLOC(charm_real_t, c_count);
    chem_c     = CHARM_ALLOC(charm_real_t, c_count);
    chem_c_    = CHARM_ALLOC(charm_real_t, c_count);
    chem_c_hat = CHARM_ALLOC(charm_real_t, c_count);

    p4est_iterate (p4est,
                   NULL, NULL,
                   _charm_model_ns_chem_iter_fn,
                   NULL, NULL, NULL);


    CHARM_FREE(chem_c);
    CHARM_FREE(chem_c_);
    CHARM_FREE(chem_c_hat);
    CHARM_FREE(chem_k);
    CHARM_FREE(chem_w);
    CHARM_FREE(chem_phi);
    CHARM_FREE(chem_psi);
}


