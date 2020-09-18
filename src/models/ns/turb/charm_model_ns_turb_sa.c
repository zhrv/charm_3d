#include <charm_globals.h>
#include "charm_globals.h"

#include "charm_globals.h"
#include "charm_base_func.h"



static void _charm_model_ns_turb_sa_params(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{

}



static void _charm_model_ns_turb_sa_grad_zero_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t *data = (charm_data_t *) info->quad->p.user_data;
    data->par.model.ns.turb.model.sa.grad_nu[0] = 0.;
    data->par.model.ns.turb.model.sa.grad_nu[1] = 0.;
    data->par.model.ns.turb.model.sa.grad_nu[2] = 0.;
}


static void _charm_model_ns_turb_sa_grad_update_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    charm_data_t *data = (charm_data_t *) info->quad->p.user_data;
    charm_real_t volume = data->par.g.volume;
    data->par.model.ns.turb.model.sa.grad_nu[0] /= volume;
    data->par.model.ns.turb.model.sa.grad_nu[1] /= volume;
    data->par.model.ns.turb.model.sa.grad_nu[2] /= volume;
}


static void _charm_model_ns_turb_sa_grad_surface_int_iter_bnd(p4est_iter_face_info_t * info, void *user_data)
{
    int i;
    p4est_t                    *p4est = info->p4est;
    charm_ctx_t                *ctx = charm_get_ctx(p4est);
    charm_data_t               *ghost_data = (charm_data_t *) user_data;
    charm_data_t               *udata;
    charm_real_t                n[3];
    charm_real_t                qu;
    p4est_iter_face_side_t     *side[2];
    sc_array_t                 *sides = &(info->sides);
    size_t                      c_count = charm_get_comp_count(info->p4est);


    int8_t face;
    charm_real_t c[2][3], l[3];
    charm_real_t nu_[2], *int_nu;
    charm_cons_t cons;
    charm_prim_t prim[2];
    charm_real_t *x, s;
    charm_real_t intg[2][5];
    int j;


    CHARM_ASSERT(info->tree_boundary);


    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    CHARM_ASSERT(!side[0]->is_hanging);

    if (side[0]->is.full.is_ghost) {
        CHARM_ASSERT(0);
        udata = &(ghost_data[side[0]->is.full.quadid]);
    } else {
        udata = charm_get_quad_data(side[0]->is.full.quad);
    }

    face = side[0]->face;
    charm_face_get_normal(udata, face, n);
    charm_quad_get_center(udata, c[0]);
    charm_face_get_center(udata, face, c[1]);

    int_nu = &(udata->par.model.ns.turb.model.sa.int_nu);

    for (i = 0; i < 3; i++) {
        l[i] = c[1][i] - c[0][i];
    }

    if (scalar_prod(n, l) < 0) {
        for (i = 0; i < 3; i++) {
            n[i] *= -1.0;
        }
    }

    {
        x = udata->par.g.fc[face];
        s = udata->par.g.area[face];
        nu_[0] = udata->par.model.ns.turb.model.sa.nu;
        if (!side[0]->is.full.is_ghost) {

            *int_nu += qu * s;
        }

    }
}


static void _charm_model_ns_turb_sa_grad_surface_int_iter_inner(p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    charm_real_t            n[3];
    charm_real_t            qu;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
    charm_real_t           *x, s;
    charm_real_t            c[2][3];
    charm_real_t            l[3];
    int8_t                  face[2];
    charm_real_t            nu_[2], *int_nu[2], bfv;



    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    face[0] = side[0]->face;
    face[1] = side[1]->face;

    h_side = -1;
    if (side[0]->is_hanging || side[1]->is_hanging) { // @todo
        for (j = 0; j < CHARM_HALF; j++) {
            for (i = 0; i < 2; i++) {
                if (side[i]->is_hanging) {
                    if (side[i]->is.hanging.is_ghost[j]) {
                        udata[i] = &(ghost_data[side[i]->is.hanging.quadid[j]]);
                    }
                    else {
                        udata[i] = (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
                    }
                    h_side = i;
                }
                else {
                    if (side[i]->is.full.is_ghost) {
                        udata[i] = &ghost_data[side[i]->is.full.quadid];
                    }
                    else {
                        udata[i] = (charm_data_t *) side[i]->is.full.quad->p.user_data;
                    }
                }
            }

            CHARM_ASSERT(h_side != -1);

            charm_face_get_normal(udata[0], face[0], n);
            charm_quad_get_center(udata[0], c[0]);
            charm_face_get_center(udata[0], face[0], c[1]);

            for (i = 0; i < 3; i++) {
                l[i] = c[1][i]-c[0][i];
            }

            if (scalar_prod(n, l) < 0) {
                for (i = 0; i < 3; i++) {
                    n[i] *= -1.0;
                }
            }

            charm_face_get_center(udata[h_side], face[h_side], x);
            s = charm_face_get_area(udata[h_side], face[h_side]);
            qu = 0.;
            for (i = 0; i < 2; i++) {
                qu += udata[i]->par.model.ns.turb.model.sa.nu;
            }
            qu *= 0.5;
            for (i = 0; i < 2; i++) {
                if (i == h_side) {
                    if (!side[i]->is.hanging.is_ghost[j]) {
                        udata[i]->par.model.ns.turb.model.sa.grad_nu[0] += qu * (i ? -1. : 1.) * s * n[0];
                        udata[i]->par.model.ns.turb.model.sa.grad_nu[1] += qu * (i ? -1. : 1.) * s * n[1];
                        udata[i]->par.model.ns.turb.model.sa.grad_nu[2] += qu * (i ? -1. : 1.) * s * n[2];
                    }
                }
                else {
                    if (!side[i]->is.full.is_ghost) {
                        udata[i]->par.model.ns.turb.model.sa.grad_nu[0] += qu * (i ? -1. : 1.) * s * n[0];
                        udata[i]->par.model.ns.turb.model.sa.grad_nu[1] += qu * (i ? -1. : 1.) * s * n[1];
                        udata[i]->par.model.ns.turb.model.sa.grad_nu[2] += qu * (i ? -1. : 1.) * s * n[2];
                    }
                }
            }
        }
    }
    else {

        for (i = 0; i < 2; i++) {
            if (side[i]->is.full.is_ghost) {
                udata[i] = &(ghost_data[side[i]->is.full.quadid]);
            }
            else {
                udata[i] = charm_get_quad_data(side[i]->is.full.quad);//(charm_data_t *) side[i]->is.full.quad->p.user_data;
            }
        }
        charm_face_get_normal(udata[0], face[0], n);
        charm_quad_get_center(udata[0], c[0]);
        charm_face_get_center(udata[0], face[0], c[1]);

        for (i = 0; i < 3; i++) {
            l[i] = c[1][i]-c[0][i];
        }

        if (scalar_prod(n, l) < 0) {
            for (i = 0; i < 3; i++) {
                n[i] *= -1.0;
            }
        }

        charm_face_get_center(udata[0], face[0], x);
        s = charm_face_get_area(udata[0], face[0]);
        qu = 0.;
        for (i = 0; i < 2; i++) {
            qu += udata[i]->par.model.ns.turb.model.sa.nu;
        }
        qu *= 0.5;

        for (i = 0; i < 2; i++) {
            if (!side[i]->is.full.is_ghost) {
                udata[i]->par.model.ns.turb.model.sa.grad_nu[0] += qu * (i ? -1. : 1.) * s * n[0];
                udata[i]->par.model.ns.turb.model.sa.grad_nu[1] += qu * (i ? -1. : 1.) * s * n[1];
                udata[i]->par.model.ns.turb.model.sa.grad_nu[2] += qu * (i ? -1. : 1.) * s * n[2];
            }
        }
    }
}


static void _charm_model_ns_turb_sa_grad_surface_int_iter_fn(p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_model_ns_turb_sa_grad_surface_int_iter_bnd(info, user_data);
    }
    else {
        _charm_model_ns_turb_sa_grad_surface_int_iter_inner(info, user_data);
    }

}







static void _charm_model_ns_turb_sa_zero_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
}


static void _charm_model_ns_turb_sa_update_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
}


static void _charm_model_ns_turb_sa_surface_int_iter_bnd(p4est_iter_face_info_t * info, void *user_data)
{
    int i;
    p4est_t *p4est = info->p4est;
    charm_ctx_t * ctx = charm_get_ctx(p4est);
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    charm_real_t n[3];
    charm_real_t qu;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);
    size_t              c_count = charm_get_comp_count(info->p4est);


    int8_t face;
    charm_real_t c[2][3], l[3];
    charm_real_t nu_[2], *int_nu;
    charm_cons_t cons;
    charm_prim_t prim[2];
    charm_real_t *x, s;
    charm_real_t intg[2][5];
    int j;


    CHARM_ASSERT(info->tree_boundary);


    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    CHARM_ASSERT(!side[0]->is_hanging);

    if (side[0]->is.full.is_ghost) {
        CHARM_ASSERT(0);
        udata = &(ghost_data[side[0]->is.full.quadid]);
    } else {
        udata = charm_get_quad_data(side[0]->is.full.quad);
    }

    face = side[0]->face;
    charm_face_get_normal(udata, face, n);
    charm_quad_get_center(udata, c[0]);
    charm_face_get_center(udata, face, c[1]);

    int_nu = &(udata->par.model.ns.turb.model.sa.int_nu);

    for (i = 0; i < 3; i++) {
        l[i] = c[1][i] - c[0][i];
    }

    if (scalar_prod(n, l) < 0) {
        for (i = 0; i < 3; i++) {
            n[i] *= -1.0;
        }
    }

    {
        x = udata->par.g.fc[face];
        s = udata->par.g.area[face];
        nu_[0] = udata->par.model.ns.turb.model.sa.nu;
        if (!side[0]->is.full.is_ghost) {

            *int_nu += qu * s;
        }

    }
}


static void _charm_model_ns_turb_sa_surface_int_iter_inner(p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    charm_real_t            n[3];
    charm_real_t            qu;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
    charm_real_t           *x, s;
    charm_real_t            c[2][3];
    charm_real_t            l[3];
    int8_t                  face[2];
    charm_real_t            nu_[2], *int_nu[2], bfv;



    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);
    face[0] = side[0]->face;
    face[1] = side[1]->face;

    h_side = -1;
    if (side[0]->is_hanging || side[1]->is_hanging) { // @todo
        for (j = 0; j < CHARM_HALF; j++) {
            for (i = 0; i < 2; i++) {
                if (side[i]->is_hanging) {
                    if (side[i]->is.hanging.is_ghost[j]) {
                        udata[i] = &(ghost_data[side[i]->is.hanging.quadid[j]]);
                    }
                    else {
                        udata[i] = (charm_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
                    }
                    h_side = i;
                }
                else {
                    if (side[i]->is.full.is_ghost) {
                        udata[i] = &ghost_data[side[i]->is.full.quadid];
                    }
                    else {
                        udata[i] = (charm_data_t *) side[i]->is.full.quad->p.user_data;
                    }
                }
            }

            CHARM_ASSERT(h_side != -1);

            charm_face_get_normal(udata[0], face[0], n);
            charm_quad_get_center(udata[0], c[0]);
            charm_face_get_center(udata[0], face[0], c[1]);

            for (i = 0; i < 3; i++) {
                l[i] = c[1][i]-c[0][i];
            }

            if (scalar_prod(n, l) < 0) {
                for (i = 0; i < 3; i++) {
                    n[i] *= -1.0;
                }
            }

            charm_face_get_center(udata[h_side], face[h_side], x);
            s = charm_face_get_area(udata[h_side], face[h_side]);
            for (i = 0; i < 2; i++) {
                nu_[i] = udata[i]->par.model.ns.turb.model.sa.nu;
                int_nu[i] = &(udata[i]->par.model.ns.turb.model.sa.int_nu);
            }
            qu = 0.; // TODO
            for (i = 0; i < 2; i++) {
                if (i == h_side) {
                    if (!side[i]->is.hanging.is_ghost[j]) {
                        *(int_nu[i]) += qu * (i ? -1. : 1.) * s;
                    }
                }
                else {
                    if (!side[i]->is.full.is_ghost) {
                        *(int_nu[i]) += qu * (i ? -1. : 1.) * s;
                    }
                }
            }
        }
    }
    else {

        for (i = 0; i < 2; i++) {
            if (side[i]->is.full.is_ghost) {
                udata[i] = &(ghost_data[side[i]->is.full.quadid]);
            }
            else {
                udata[i] = charm_get_quad_data(side[i]->is.full.quad);//(charm_data_t *) side[i]->is.full.quad->p.user_data;
            }
        }
        charm_face_get_normal(udata[0], face[0], n);
        charm_quad_get_center(udata[0], c[0]);
        charm_face_get_center(udata[0], face[0], c[1]);

        for (i = 0; i < 3; i++) {
            l[i] = c[1][i]-c[0][i];
        }

        if (scalar_prod(n, l) < 0) {
            for (i = 0; i < 3; i++) {
                n[i] *= -1.0;
            }
        }

        charm_face_get_center(udata[0], face[0], x);
        s = charm_face_get_area(udata[0], face[0]);
        for (i = 0; i < 2; i++) {
            nu_[i] = udata[i]->par.model.ns.turb.model.sa.nu;
            int_nu[i] = &(udata[i]->par.model.ns.turb.model.sa.int_nu);
        }
        qu = 0.; // TODO

        for (i = 0; i < 2; i++) {
            if (!side[i]->is.full.is_ghost) {
                *(int_nu[i]) += qu * (i ? -1. : 1.) * s;
            }
        }
    }
}


static void _charm_model_ns_turb_sa_surface_int_iter_fn(p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        _charm_model_ns_turb_sa_surface_int_iter_bnd(info, user_data);
    }
    else {
        _charm_model_ns_turb_sa_surface_int_iter_inner(info, user_data);
    }

}


void charm_model_ns_turb_sa(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    _charm_model_ns_turb_sa_params(p4est, ghost, ghost_data);

    // calc gradients

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_model_ns_turb_sa_grad_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   NULL,
                   _charm_model_ns_turb_sa_grad_surface_int_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est, NULL,
                   NULL,
                   _charm_model_ns_turb_sa_grad_update_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_ghost_exchange_data (p4est, ghost, ghost_data); /* synchronize the ghost data */


    // calc main equation

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   _charm_model_ns_turb_sa_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   NULL,
                   _charm_model_ns_turb_sa_surface_int_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est, NULL,
                   NULL,
                   _charm_model_ns_turb_sa_update_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_ghost_exchange_data (p4est, ghost, ghost_data); /* synchronize the ghost data */

}

