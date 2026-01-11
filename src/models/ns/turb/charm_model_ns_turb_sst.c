#include "charm_globals.h"
#include "charm_base_func.h"



static void charm_model_ns_turb_sst_params(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{

}

static void charm_model_ns_turb_sst_zero_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
}


static void charm_model_ns_turb_sst_update_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
}


static void charm_model_ns_turb_sst_surface_int_iter_bnd(p4est_iter_face_info_t * info, void *user_data)
{
    int i, ibf, igp;
    p4est_t *p4est = info->p4est;
    charm_ctx_t * ctx = charm_get_ctx(p4est);
    charm_data_t *ghost_data = (charm_data_t *) user_data;
    charm_data_t *udata;
    charm_real_t n[3];
    charm_real_t qu, qv, qw, qe, *qc;
    charm_real_t bfv;
    p4est_iter_face_side_t *side[2];
    sc_array_t *sides = &(info->sides);
    size_t              c_count = charm_get_comp_count(info->p4est);


    int8_t face;
    charm_real_t c[2][3], l[3];
    charm_real_t r_[2], p_[2], u_[2], v_[2], w_[2], e_[2];
    charm_cons_t cons;
    charm_prim_t prim[2];
    charm_real_t *x, gw, gj;
    charm_real_t intg[2][5];
    int j;


    CHARM_ASSERT(info->tree_boundary);

    qc = CHARM_ALLOC(charm_real_t, c_count);

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

    for (i = 0; i < 3; i++) {
        l[i] = c[1][i] - c[0][i];
    }

    if (scalar_prod(n, l) < 0) {
        for (i = 0; i < 3; i++) {
            n[i] *= -1.0;
        }
    }

    for (igp = 0; igp < CHARM_FACE_GP_COUNT; igp++) {
        x = udata->par.g.face_gp[face][igp];
        gw = udata->par.g.face_gw[face][igp];
        gj = udata->par.g.face_gj[face][igp];
        charm_get_fields(udata, x, &cons);
        charm_param_cons_to_prim(p4est, &(prim[0]), &cons);
        // TODO charm_bnd_cond(p4est, side[0]->treeid, face, &(prim[0]), &(prim[1]), n);
        ctx->flux_fn(p4est, prim, &qu, &qv, &qw, &qe, qc, n); /* flux from side 0 to side 1 */
        for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
            if (!side[0]->is.full.is_ghost) {
                bfv = charm_base_func(x, ibf, udata) * gw * gj;
                for (j = 0; j < c_count; j++) {
                    udata->int_rc[j][ibf] += qc[j] * bfv;
                }

                udata->int_ru[ibf] += qu * bfv;
                udata->int_rv[ibf] += qv * bfv;
                udata->int_rw[ibf] += qw * bfv;
                udata->int_re[ibf] += qe * bfv;
            }
        }
    }
    CHARM_FREE(qc);

}


static void charm_model_ns_turb_sst_surface_int_iter_inner(p4est_iter_face_info_t * info, void *user_data)
{
    int                     i, j, h_side, igp, ibf,cj;
    p4est_t                *p4est = info->p4est;
    charm_ctx_t            *ctx = charm_get_ctx(p4est);
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    charm_real_t                  n[3];
    charm_real_t                  qu, qv, qw, qe, *qc;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
    charm_real_t                 *x, gw, gj;
    charm_real_t                  bfv;
    charm_real_t                  c[2][3];
    charm_real_t                  l[3];
    int8_t                  face[2];
    size_t                  c_count = charm_get_comp_count(info->p4est);


    qc = CHARM_ALLOC(charm_real_t, c_count);

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

            for (igp = 0; igp < CHARM_FACE_GP_COUNT; igp++) {
                x  = udata[h_side]->par.g.face_gp[face[h_side]][igp];
                gw = udata[h_side]->par.g.face_gw[face[h_side]][igp];
                gj = udata[h_side]->par.g.face_gj[face[h_side]][igp];
                for (i = 0; i < 2; i++) {
                    charm_get_fields(udata[i], x, &(cons[i]));
                    charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
                }
                ctx->flux_fn(p4est, prim, &qu, &qv, &qw, &qe, qc, n);  // flux from side 0 to side 1
                for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                    for (i = 0; i < 2; i++) {
                        if (i == h_side) {
                            if (!side[i]->is.hanging.is_ghost[j]) {
                                bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                                for (cj = 0; cj < c_count; cj++) {
                                    udata[i]->int_rc[cj][ibf] += qc[cj] * bfv;
                                }
                                udata[i]->int_ru[ibf] += qu * bfv;
                                udata[i]->int_rv[ibf] += qv * bfv;
                                udata[i]->int_rw[ibf] += qw * bfv;
                                udata[i]->int_re[ibf] += qe * bfv;
                            }
                        }
                        else {
                            if (!side[i]->is.full.is_ghost) {
                                bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                                for (cj = 0; cj < c_count; cj++) {
                                    udata[i]->int_rc[cj][ibf] += qc[cj] * bfv;
                                }
                                udata[i]->int_ru[ibf] += qu * bfv;
                                udata[i]->int_rv[ibf] += qv * bfv;
                                udata[i]->int_rw[ibf] += qw * bfv;
                                udata[i]->int_re[ibf] += qe * bfv;
                            }
                        }
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

        for (igp = 0; igp < CHARM_FACE_GP_COUNT; igp++) {
            x  = udata[0]->par.g.face_gp[face[0]][igp];
            gw = udata[0]->par.g.face_gw[face[0]][igp];
            gj = udata[0]->par.g.face_gj[face[0]][igp];
            for (i = 0; i < 2; i++) {
                charm_get_fields(udata[i], x, &(cons[i]));
                charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));
            }
            ctx->flux_fn(p4est, prim, &qu, &qv, &qw, &qe, qc, n);  // flux from side 0 to side 1
            for (ibf = 0; ibf < CHARM_BASE_FN_COUNT; ibf++) {
                for (i = 0; i < 2; i++) {
                    if (!side[i]->is.full.is_ghost) {
                        bfv = (i ? -1. : 1.) * charm_base_func(x, ibf, udata[i]) * gw * gj;
                        for (j = 0; j < c_count; j++) {
                            udata[i]->int_rc[j][ibf] += qc[j] * bfv;
                        }
                        udata[i]->int_ru[ibf] += qu * bfv;
                        udata[i]->int_rv[ibf] += qv * bfv;
                        udata[i]->int_rw[ibf] += qw * bfv;
                        udata[i]->int_re[ibf] += qe * bfv;
                    }
                }
            }
        }
    }
    CHARM_FREE(qc);
}


static void charm_model_ns_turb_sst_surface_int_iter_fn(p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        charm_model_ns_turb_sst_surface_int_iter_bnd(info, user_data);
    }
    else {
        charm_model_ns_turb_sst_surface_int_iter_inner(info, user_data);
    }

}


void charm_model_ns_turb_sst(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_model_ns_turb_sst_params(p4est, ghost, ghost_data);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   charm_model_ns_turb_sst_zero_quad_iter_fn,
                   NULL, NULL, NULL);

    p4est_iterate (p4est,
                   ghost,
                   (void *) ghost_data,
                   NULL,
                   charm_model_ns_turb_sst_surface_int_iter_fn,
                   NULL, NULL);

    p4est_iterate (p4est, NULL,
                   NULL,
                   charm_model_ns_turb_sst_update_quad_iter_fn,
                   NULL, NULL, NULL);


    p4est_ghost_exchange_data (p4est, ghost, ghost_data); /* synchronize the ghost data */

}
