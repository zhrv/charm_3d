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



void charm_model_ns_timestep_chem(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    if (!ctx->reactions) return;


}


static void _charm_model_ns_chem_init_fetch_reaction(charm_ctx_t *ctx, mxml_node_t *node, charm_reaction_t *r)
{
    int i, comp;
    int c_count = ctx->comp->elem_count;
    mxml_node_t *left, *right, *n;
    charm_xml_node_child_param_dbl(node, "a", &(r->lg_a));
    charm_xml_node_child_param_dbl(node, "e", &(r->e));
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