#include "charm_globals.h"
#include "charm_base_func.h"
#include "charm_bnd_cond.h"


static charm_real_t sigma;
static charm_real_t kappa;
static charm_real_t cb1;
static charm_real_t cb2;
static charm_real_t cw1;
static charm_real_t cw2;
static charm_real_t cw3;
static charm_real_t cv1;
static charm_real_t ct1;
static charm_real_t ct2;
static charm_real_t ct3;
static charm_real_t ct4;

charm_real_t charm_model_ns_get_visc_mu(p4est_t* p4est, charm_real_t *x, charm_data_t* data);
void charm_model_ns_turb_sa_grad(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);


static void charm_model_ns_turb_sa_params(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t                *ctx = charm_get_ctx(p4est);
    sigma = ctx->model.ns.turb.param.sa.sigma;
    kappa = ctx->model.ns.turb.param.sa.kappa;
    cb1 = ctx->model.ns.turb.param.sa.cb1;
    cb2 = ctx->model.ns.turb.param.sa.cb2;
    cw1 = ctx->model.ns.turb.param.sa.cw1;
    cw2 = ctx->model.ns.turb.param.sa.cw2;
    cw3 = ctx->model.ns.turb.param.sa.cw3;
    cv1 = ctx->model.ns.turb.param.sa.cv1;
    ct1 = ctx->model.ns.turb.param.sa.ct1;
    ct2 = ctx->model.ns.turb.param.sa.ct2;
    ct3 = ctx->model.ns.turb.param.sa.ct3;
    ct4 = ctx->model.ns.turb.param.sa.ct4;

}








static void charm_model_ns_turb_sa_zero_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_t *p4est = info->p4est;
    charm_data_t *data = (charm_data_t *) info->quad->p.user_data;
    charm_real_t nu_ = data->par.model.ns.turb.model.sa.nu_;
    charm_real_t volume;
    charm_real_t d = data->par.g.y;
    charm_real_t *grad_nu_ = data->par.model.ns.turb.model.sa.grad_nu_;
    charm_real_t *c;
    charm_cons_t cons;
    charm_prim_t prim;
    charm_real_t mu, nu;
    charm_real_t xi, xi3, fv1, fv2, Omega, S_, Pnu, Dnu, r, g, kd2, ft1, ft2, cw36, fw, Oij;
    charm_int_t i,j;

    for (i = 0; i < CHARM_DIM; i++) {
        for (j = 0; j < CHARM_DIM; j++) {
            Oij = 2.*(data->par.model.ns.turb.model.sa.grad_u[i][j]
                     -data->par.model.ns.turb.model.sa.grad_u[j][i]);
            Omega += Oij*Oij;
        }
    }
    Omega = sqrt(Omega);
    charm_quad_get_center(data, c);
    volume = charm_quad_get_volume(data);
    charm_get_fields(data, c, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);
    mu = charm_model_ns_get_visc_mu(p4est, c, data);
    nu = mu/prim.r;
    xi = nu_/nu;
    xi3 = pow(xi, 3.);
    fv1 = xi3/(xi3+pow(cv1, 3.));
    fv2 = 1-xi/(1+xi*fv1);
    //ft2 = ct3*exp(-ct4*xi*xi);
    kd2 = _SQR_(kappa*d);
    S_= Omega+fv2*nu_/kd2;
    r = nu_/(S_*kd2);
    g = r+cw2*(pow(r,6.)-r);
    cw36 = pow(cw3, 6.);
    fw = g*pow((1+cw36)/(pow(g, 6.)+cw36), 1./6.);
    Pnu = cb1*(1.-ft2)*S_*nu_;
    Dnu = (cw1*fw-cb1*ft2/_SQR_(kappa))*_SQR_(nu_/d);

    data->par.model.ns.turb.model.sa.int_nu_ = -volume*(Pnu-Dnu+(scalar_prod(grad_nu_, grad_nu_))*cb2/sigma);
}


static void charm_model_ns_turb_sa_update_quad_iter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
    p4est_t            *p4est = info->p4est;
    charm_data_t       *data = charm_get_quad_data(info->quad);
    charm_ctx_t        *ctx = (charm_ctx_t*)info->p4est->user_pointer;
    charm_real_t        dt = *((charm_real_t *) user_data);
    charm_real_t        nut, fv1, xi, nu, xi3;
    charm_point_t       c;
    charm_cons_t        cons;
    charm_prim_t        prim;
    int                 i, j;

    charm_quad_get_center(data, c);
    charm_get_fields(data, c, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);
    data->par.model.ns.turb.model.sa.nu_ -= _NORM_(dt * data->par.model.ns.turb.model.sa.int_nu_);
    nu = charm_model_ns_get_visc_mu(p4est, c, data)/prim.r;
    xi = data->par.model.ns.turb.model.sa.nu_/nu;
    xi3 = pow(xi, 3.);
    fv1 = xi3/(xi3+pow(cv1, 3.));
    nut = fv1*data->par.model.ns.turb.model.sa.nu_;
    data->par.model.ns.turb.mu_t = nut*prim.r;
}


static void charm_model_ns_turb_sa_surface_int_iter_bnd(p4est_iter_face_info_t * info, void *user_data)
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

    charm_bnd_types_t bnd_type;
    int8_t face;
    charm_real_t c[2][3], l[3];
    charm_real_t nu_[2], nu[2], nut[2], *int_nu, un[2], dnu_dn[2], mu[2];
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

    bnd_type = charm_bnd_get_type(p4est, side[0]->treeid, face);

    int_nu = &(udata->par.model.ns.turb.model.sa.int_nu_);

    for (i = 0; i < 3; i++) {
        l[i] = c[1][i] - c[0][i];
    }

    if (scalar_prod(n, l) < 0) {
        for (i = 0; i < 3; i++) {
            n[i] *= -1.0;
        }
    }

    x = c[1];
    s = charm_face_get_area(udata, face);
    charm_get_fields(udata, x, &cons);
    charm_param_cons_to_prim(p4est, &(prim[0]), &cons);

    un[0]       = prim[i].u*n[0]+prim[i].v*n[1]+prim[i].w*n[2];
    mu[0]       = charm_model_ns_get_visc_mu(p4est, x, udata);
    nu_[0]      = udata->par.model.ns.turb.model.sa.nu_;
    nu[0]       = mu[i]/prim[i].r;
    int_nu      = &(udata->par.model.ns.turb.model.sa.int_nu_);
    dnu_dn[0]   = udata->par.model.ns.turb.model.sa.grad_nu_[0]*n[0]+
                    udata->par.model.ns.turb.model.sa.grad_nu_[1]*n[1]+
                    udata->par.model.ns.turb.model.sa.grad_nu_[2]*n[2];

    if (bnd_type == BOUND_WALL_NO_SLIP) {
        qu = 0.;
    }
    else {
        qu = 0.5*(
                nu_[0]*un[0]-(nu[0]+nu_[0])*dnu_dn[0]/sigma+
                nu_[1]*un[1]-(nu[1]+nu_[1])*dnu_dn[1]/sigma
        );
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


static void charm_model_ns_turb_sa_surface_int_iter_inner(p4est_iter_face_info_t * info, void *user_data)
{
    p4est_t                *p4est = info->p4est;
    int                     i, j, h_side;
    charm_data_t           *ghost_data = (charm_data_t *) user_data;
    charm_data_t           *udata[2];
    charm_real_t            n[3];
    charm_real_t            qu;
    p4est_iter_face_side_t *side[2];
    sc_array_t             *sides = &(info->sides);
    charm_real_t           *x, s;
    charm_real_t            c[2][3];
    charm_real_t            l[3];
    charm_cons_t            cons[2];
    charm_prim_t            prim[2];
    int8_t                  face[2];
    charm_real_t            nu_[2], nu[2], nut[2], *int_nu[2], un[2], dnu_dn[2], mu[2];



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
                charm_get_fields(udata[i], x, &(cons[i]));
                charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));

                un[i] = prim[i].u*n[0]+prim[i].v*n[1]+prim[i].w*n[2];
                mu[i]       = charm_model_ns_get_visc_mu(p4est, x, udata[i]);
                nu_[i] = udata[i]->par.model.ns.turb.model.sa.nu_;
                nu[i] = mu[i]/prim[i].r;
                int_nu[i] = &(udata[i]->par.model.ns.turb.model.sa.int_nu_);
                dnu_dn[i] = udata[i]->par.model.ns.turb.model.sa.grad_nu_[0]*n[0]+
                            udata[i]->par.model.ns.turb.model.sa.grad_nu_[1]*n[1]+
                            udata[i]->par.model.ns.turb.model.sa.grad_nu_[2]*n[2];
            }
            qu = 0.5*(
                        nu_[0]*un[0]-(nu[0]+nu_[0])*dnu_dn[0]/sigma+
                        nu_[1]*un[1]-(nu[1]+nu_[1])*dnu_dn[1]/sigma
                    );
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
            charm_get_fields(udata[i], x, &(cons[i]));
            charm_param_cons_to_prim(p4est, &(prim[i]), &(cons[i]));

            un[i] = prim[i].u*n[0]+prim[i].v*n[1]+prim[i].w*n[2];
            mu[i]       = charm_model_ns_get_visc_mu(p4est, x, udata[i]);
            nu_[i] = udata[i]->par.model.ns.turb.model.sa.nu_;
            nu[i] = mu[i]/prim[i].r;
            int_nu[i] = &(udata[i]->par.model.ns.turb.model.sa.int_nu_);
            dnu_dn[i] = udata[i]->par.model.ns.turb.model.sa.grad_nu_[0]*n[0]+
                        udata[i]->par.model.ns.turb.model.sa.grad_nu_[1]*n[1]+
                        udata[i]->par.model.ns.turb.model.sa.grad_nu_[2]*n[2];
        }
        qu = 0.5*(
                nu_[0]*un[0]-(nu[0]+nu_[0])*dnu_dn[0]/sigma+
                nu_[1]*un[1]-(nu[1]+nu_[1])*dnu_dn[1]/sigma
        );

        for (i = 0; i < 2; i++) {
            if (!side[i]->is.full.is_ghost) {
                *(int_nu[i]) += qu * (i ? -1. : 1.) * s;
            }
        }
    }
}


static void charm_model_ns_turb_sa_surface_int_iter_fn(p4est_iter_face_info_t * info, void *user_data)
{
    sc_array_t         *sides = &(info->sides);

    if (sides->elem_count != 2) {
        charm_model_ns_turb_sa_surface_int_iter_bnd(info, user_data);
    }
    else {
        charm_model_ns_turb_sa_surface_int_iter_inner(info, user_data);
    }

}

static void charm_model_ns_turb_sa_main(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_ctx_t    *ctx = charm_get_ctx(p4est);
    charm_real_t    dt = ctx->get_dt_fn(p4est);

    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   charm_model_ns_turb_sa_zero_quad_iter_fn, NULL, NULL, NULL);

    p4est_iterate (p4est, ghost, (void *) ghost_data,
                   NULL, charm_model_ns_turb_sa_surface_int_iter_fn, NULL, NULL);

    p4est_iterate (p4est, NULL, (void *) &dt,
                   charm_model_ns_turb_sa_update_quad_iter_fn, NULL, NULL, NULL);

    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

}



void charm_model_ns_turb_sa(p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data)
{
    charm_model_ns_turb_sa_params(p4est, ghost, ghost_data);

    charm_model_ns_turb_sa_grad(p4est, ghost, ghost_data);
    charm_model_ns_turb_sa_main(p4est, ghost, ghost_data);
}

