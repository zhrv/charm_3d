//
// Created by zhrv on 08.11.19.
//
#include <p8est_iterate.h>
#include "charm_base_func.h"
#include "charm_fluxes.h"
#include "charm_bnd_cond.h"
#include "charm_globals.h"
#include "charm_limiter.h"
#include "charm_amr.h"
#include "charm_xml.h"
#include "charm_eos.h"


static double *chem_k, *chem_w, *chem_phi, *chem_psi, *chem_c;


static void _charm_model_ns_chem_init_constants(p4est_t * p4est, charm_data_t *data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    size_t r_count = charm_get_reactions_count(p4est);
}


static void _charm_model_ns_chem_rhs(p4est_t * p4est, charm_data_t *data)
{
    double       R   = charm_eos_get_r();
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    size_t r_count = charm_get_reactions_count(p4est);
    size_t c_count = charm_get_comp_count(p4est);
    charm_cons_t cons;
    charm_prim_t prim;
    charm_comp_t *comp;
    charm_reaction_t * r;
    int i, j1, j2, m, n;
    double c1, c2, c3;

    memset(chem_phi, 0, sizeof(double)*c_count);
    memset(chem_psi, 0, sizeof(double)*c_count);
    for (i = 0; i < c_count; i++) {
        chem_phi[i] = 0.;
        chem_psi[i] = 0.;
        for (n = 0; n < r_count; n++) {
            r = charm_get_reaction(p4est, n);
            for (m = 0; m < 3; m++) {
                if (i == r->left_comps[m]) {
                    j1 = (m + 1) % 3;
                    j2 = (m + 2) % 3;
                    c1 = r->left_comps[j1] == 0 ? 1 : chem_c[r->left_comps[j1]];
                    c2 = r->left_comps[j2] == 0 ? 1 : chem_c[r->left_comps[j2]];
                    chem_phi[i] += r->a*exp(-r->e/prim.t/R)*c1*c2;
                }
            }
        }
    }
}


static void _charm_model_ns_chem_solve(p4est_t * p4est, charm_data_t *data)
{
    double       R   = charm_eos_get_r();
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    size_t r_count = charm_get_reactions_count(p4est);
    size_t c_count = charm_get_comp_count(p4est);
    charm_cons_t cons;
    charm_prim_t prim;
    charm_comp_t *comp;
    charm_reaction_t * r;
    int i, j, m, n;
    double ci;

    charm_get_fields_avg(data, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);
    for (i = 0; i < c_count; i++) {
        comp = charm_get_comp(p4est, i);
        chem_c[i] = (cons.rc[i]/comp->m);
    }

//    for (n = 0; n < r_count; n++) {
//        r = charm_get_reaction(p4est, n);
//        chem_w[n] = r->a*exp(-r->e/prim.t/R);
//        for (m = 0; m < 3; m++) {
//            i = r->left_comps[m];
//            comp = charm_get_comp(p4est, i);
//            chem_w[n] *= (cons.rc[i]/comp->m);
//        }
//    }
}


static void _charm_model_ns_chem_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_quadrant_t   *q = info->quad;
    charm_data_t       *data = charm_get_quad_data(q);
    int                 ibf, igp;
    charm_reaction_t    r;

    _charm_model_ns_chem_init_constants(info->p4est, data);
    _charm_model_ns_chem_solve(info->p4est, data);
}

void charm_model_ns_timestep_chem(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    if (!ctx->reactions) return;




    size_t c_count = charm_get_comp_count(p4est);
    size_t r_count = charm_get_reactions_count(p4est);

    chem_k = CHARM_ALLOC(double, r_count);
    chem_w = CHARM_ALLOC(double, r_count);
    chem_phi = CHARM_ALLOC(double, c_count);
    chem_psi = CHARM_ALLOC(double, c_count);
    chem_c   = CHARM_ALLOC(double, c_count);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_model_ns_chem_iter_fn,
                   NULL, NULL, NULL);


    CHARM_FREE(chem_c);
    CHARM_FREE(chem_k);
    CHARM_FREE(chem_w);
    CHARM_FREE(chem_phi);
    CHARM_FREE(chem_psi);
}


static void _charm_model_ns_chem_init_fetch_reaction(charm_ctx_t *ctx, mxml_node_t *node, charm_reaction_t *r)
{
    int i, comp;
    int c_count = ctx->comp->elem_count;
    double lg_a;
    mxml_node_t *left, *right, *n;
    charm_xml_node_child_param_dbl(node, "a", &lg_a);
    charm_xml_node_child_param_dbl(node, "e", &(r->e));
    r->a = pow(10., lg_a);
    left = charm_xml_node_get_child(node, "left");
    for (n = charm_xml_node_get_child(left, "comp"), i = 0;
         n != NULL;
         n = charm_xml_node_get_next_child(n, left, "comp"), i++) {
        if (i >= 3) {
            CHARM_LERROR("More than 3 components in reaction!\n");
            charm_abort(NULL, 1);
        }
        charm_xml_node_value_int(n, &comp);
        if (!( -1 <= comp && comp < c_count )) {
            CHARM_LERRORF("Wrong component number in reaction: %d, but total components count: %d\n", comp, c_count);
            charm_abort(NULL, 1);
        }
        r->left_comps[i] = comp;
    }

    right = charm_xml_node_get_child(node, "right");
    for (n = charm_xml_node_get_child(right, "comp"), i = 0;
         n != NULL;
         n = charm_xml_node_get_next_child(n, right, "comp"), i++) {
        if (i >= 3) {
            CHARM_LERROR("More than 3 components in reaction!\n");
            charm_abort(NULL, 1);
        }
        charm_xml_node_value_int(n, &comp);
        if (!( -1 <= comp && comp < c_count )) {
            CHARM_LERRORF("Wrong component number in reaction: %d, but total components count: %d\n", comp, c_count);
            charm_abort(NULL, 1);
        }
        r->right_comps[i] = comp;
    }

}


void charm_model_ns_chem_init(charm_ctx_t *ctx, mxml_node_t *node_task)
{
    mxml_node_t *r_node, *node;
    r_node = charm_xml_node_get_child_or_null(node_task, "reactions");
    if (r_node == NULL) {
        ctx->reactions = NULL;
        return;
    }
    charm_reaction_t *r;

    ctx->reactions = sc_array_new(sizeof(charm_reaction_t));

    for (node = charm_xml_node_get_child(r_node, "reaction");
         node != NULL;
         node = charm_xml_node_get_next_child(node, r_node, "reaction")) {
        r = (charm_reaction_t *) sc_array_push(ctx->reactions);
        _charm_model_ns_chem_init_fetch_reaction(ctx, node, r);
    }
}