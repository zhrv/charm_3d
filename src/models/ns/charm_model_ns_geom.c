//
// Created by zhrv on 13.06.2020.
//
#include <charm_globals.h>

static void _charm_model_ns_geom_calc_face_iter_fn (p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t             *fc = user_data;
    sc_array_t             *sides = &(info->sides);
    p4est_iter_face_side_t *side;
    charm_data_t           *udata;
    int8_t                  face;
    charm_real_t           *c;

    if (sides->elem_count != 2) {
        side = p4est_iter_fside_array_index_int(sides, 0);
        face = side->face;
        charm_tree_attr_t *attr = charm_get_tree_attr(info->p4est, side->treeid);
        charm_bnd_t *bnd = attr->bnd[face];

        if (bnd->type == BOUND_WALL_NO_SLIP) {
            if (side->is.full.is_ghost) {
                CHARM_ASSERT(0);
            } else {
                udata = charm_get_quad_data(side->is.full.quad);//(charm_data_t *) side[0]->is.full.quad->p.user_data;
            }
            c = sc_array_push(fc);
            c[0] = udata->par.g.fc[face][0];
            c[1] = udata->par.g.fc[face][1];
            c[2] = udata->par.g.fc[face][2];
        }

    }
}


static void _charm_model_ns_geom_calc_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    sc_array_t         *fc = user_data;
    charm_data_t       *data = charm_get_quad_data(info->quad);
    int                 i, j;
    charm_real_t       *c;
    charm_real_t        l;

    data->par.g.y = DBL_MAX;
    for (i = 0; i < fc->elem_count; i++) {
        c = sc_array_index(fc, i);
        l =vect_dist(c, data->par.g.c);
        if (data->par.g.y > l) {
            data->par.g.y = l;
        }
    }
}


void charm_model_ns_geom_calc(p4est_t *p4est)
{
    sc_array_t *fc;
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    if (ctx->model.ns.turb.model_type == TURB_MODEL_UNKNOWN) return;

    fc = sc_array_new(3*sizeof(charm_real_t));
    p4est_iterate (p4est,
                   NULL,
                   (void *) fc,
                   NULL,
                   _charm_model_ns_geom_calc_face_iter_fn,
                   NULL, NULL);
    p4est_iterate (p4est,
                   NULL,
                   (void *) fc,
                   _charm_model_ns_geom_calc_quad_iter_fn,
                   NULL, NULL, NULL);
    sc_array_destroy(fc);
}