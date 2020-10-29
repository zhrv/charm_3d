//
// Created by zhrv on 26.10.17.
//
#define GLOBALS_H_FILE

#include <charm_base_func.h>
#include "charm_globals.h"


p4est_t * g_p4est = NULL;
int charm_package_id = -1;

const char *charm_bnd_types[] ={
        "BOUND_INLET",
        "BOUND_OUTLET",
        "BOUND_WALL_SLIP",
        "BOUND_WALL_NO_SLIP",
        "BOUND_MASS_FLOW",
        "BOUND_SYMMETRY",
        "BOUND_FREE_STREAM",
        "BOUND_PRESSURE",
        NULL
};

const char *charm_turb_models[] ={
        "SA",
        "SST",
        NULL
};


p4est_t * charm_get_p4est() { return g_p4est; }
void charm_set_p4est(p4est_t *p4est) { g_p4est = p4est; }



charm_real_t scalar_prod(charm_vec_t v1, charm_vec_t v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}


charm_real_t vector_length(charm_vec_t v)
{
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

void vector_sub(charm_vec_t v1, charm_vec_t v2, charm_vec_t res)
{
    res[0] = v1[0]-v2[0];
    res[1] = v1[1]-v2[1];
    res[2] = v1[2]-v2[2];
}

charm_real_t vector_dist(charm_vec_t v1, charm_vec_t v2)
{
    charm_vec_t res;
    vector_sub(v1, v2, res);
    return vector_length(res);
}

void vector_prod(charm_vec_t v1, charm_vec_t v2, charm_vec_t res)
{
    res[0] =  v1[1]*v2[2]-v1[2]*v2[1];
    res[1] = -v1[0]*v2[2]+v1[2]*v2[0];
    res[2] =  v1[0]*v2[1]-v1[1]*v2[0];
}


charm_real_t charm_face_get_area(charm_data_t *d, int8_t face)
{
    return d->par.g.area[face];
}

charm_real_t charm_face_get_normal(charm_data_t *d, int8_t face, charm_vec_t n)
{
    memcpy(n, d->par.g.n[face], sizeof(charm_vec_t));
    return d->par.g.area[face];
}

void charm_quad_get_center(charm_data_t *d, charm_vec_t c)
{
    memcpy(c, d->par.g.c, sizeof(charm_vec_t));
}

void charm_face_get_center(charm_data_t *d, int8_t face, charm_vec_t c)
{
    memcpy(c, d->par.g.fc[face], sizeof(charm_vec_t));
}

charm_real_t charm_quad_get_volume(charm_data_t *d)
{
    return d->par.g.volume;
}


charm_mat_t * charm_mat_find_by_id(charm_ctx_t *ctx, int id)
{
    size_t i;
    sc_array_t *arr = ctx->mat;
    charm_mat_t * mat;


    for (i = 0; i < arr->elem_count; i++) {
        mat = sc_array_index(arr, i);
        if (mat->id == id) {
            return mat;
        }
    }

    return NULL;
}

int charm_mat_index_find_by_id(charm_ctx_t *ctx, int id, size_t *index)
{
    size_t i;
    sc_array_t *arr = ctx->mat;
    charm_mat_t * mat;


    for (i = 0; i < arr->elem_count; i++) {
        mat = sc_array_index(arr, i);
        if (mat->id == id) {
            *index = i;
            return 1;
        }
    }

    return 0;
}

charm_comp_t * charm_comp_find_by_id(charm_ctx_t *ctx, int id)
{
    size_t i;
    sc_array_t *arr = ctx->comp;
    charm_comp_t * comp;


    for (i = 0; i < arr->elem_count; i++) {
        comp = sc_array_index(arr, i);
        if (comp->id == id) {
            return comp;
        }
    }

    return NULL;
}

int charm_comp_index_find_by_id(charm_ctx_t *ctx, int id, size_t *index)
{
    size_t i;
    sc_array_t *arr = ctx->comp;
    charm_comp_t * comp;


    for (i = 0; i < arr->elem_count; i++) {
        comp = sc_array_index(arr, i);
        if (comp->id == id) {
            *index = i;
            return 1;
        }
    }

    return 0;
}

charm_reg_t * charm_reg_find_by_id(charm_ctx_t *ctx, int id)
{
  size_t i;
  sc_array_t *arr = ctx->reg;
  charm_reg_t * reg;


  for (i = 0; i < arr->elem_count; i++) {
    reg = sc_array_index(arr, i);
    if (reg->id == id) {
      return reg;
    }
  }

  return NULL;
}

charm_mesh_type_t charm_mesh_get_type_by_str(char *str)
{
    if (strcmp(str, "gmsh_msh") == 0) {
        return CHARM_MESH_GMSH_MSH;
    }
    else if (strcmp(str, "gmsh_unv") == 0) {
        return CHARM_MESH_GMSH_UNV;
    }
    else if (strcmp(str, "gmsh_inp") == 0) {
        return CHARM_MESH_GMSH_INP;
    }
    else {
        return CHARM_MESH_UNKNOWN;
    }
}

charm_tree_attr_t * charm_get_tree_attr(p4est_t * p4est, p4est_topidx_t which_tree)
{
    p4est_connectivity_t * conn = p4est->connectivity;
    return (charm_tree_attr_t *)&(conn->tree_to_attr[sizeof(charm_tree_attr_t)*which_tree]);
}

void charm_param_cons_to_prim(p4est_t * p4est, charm_prim_t * p, charm_cons_t * c)
{
    charm_ctx_t * ctx       = charm_get_ctx(p4est);
    size_t        c_count   = ctx->comp->elem_count;
    charm_mat_t * mat       = charm_mat_find_by_id(ctx, c->mat_id);
    size_t i;

    p->mat_id = c->mat_id;
    p->r      = 0.;
    for (i = 0; i < c_count; i++) {
        p->r += c->rc[i];
    }
    p->u      = c->ru/p->r;
    p->v      = c->rv/p->r;
    p->w      = c->rw/p->r;
    p->e_tot  = c->re/p->r;
    p->e      = p->e_tot-0.5*_MAG_(p->u, p->v, p->w);

    for (i = 0; i < c_count; i++) {
        p->c[i] = c->rc[i] / p->r;
        if (p->c[i] < 0.) p->c[i] = 0.; // @todo
        if (p->c[i] > 1.) p->c[i] = 1.;
    }

    mat->eos_fn(p4est, p, 4);  // {p,cz, t}=EOS(r,e)
}


void charm_param_prim_to_cons(p4est_t * p4est, charm_cons_t * c, charm_prim_t * p)
{
    size_t c_count = charm_get_comp_count(p4est);
    size_t i;
    c->mat_id = p->mat_id;
    c->ru = p->r * p->u;
    c->rv = p->r * p->v;
    c->rw = p->r * p->w;
    c->re = p->r * (p->e + 0.5 * _MAG_(p->u, p->v, p->w));
    for (i = 0; i < c_count; i++) {
        c->rc[i] = p->r * p->c[i];
    }
}

void charm_prim_cpy(charm_prim_t * dest, charm_prim_t * src)
{
    dest->mat_id = src->mat_id;
    dest->r      = src->r;
    dest->p      = src->p;
    dest->u      = src->u;
    dest->v      = src->v;
    dest->w      = src->w;
    dest->t      = src->t;
    dest->cz     = src->cz;
    dest->e      = src->e;
    dest->e_tot  = src->e_tot;
    dest->cp     = src->cp;
    dest->cv     = src->cv;
    dest->gam    = src->gam;
    dest->m      = src->m;
    memcpy(dest->c, src->c, CHARM_MAX_COMPONETS_COUNT*sizeof(charm_real_t));
}


charm_real_t charm_prim_vel_mag(charm_prim_t * prim)
{
    return sqrt(prim->u*prim->u+prim->v*prim->v+prim->w*prim->w);
}


charm_real_t charm_matr3_det(charm_real_t a[3][3])
{
    charm_real_t det_a = 0.;

    det_a += a[0][0]*a[1][1]*a[2][2];
    det_a += a[0][2]*a[1][0]*a[2][1];
    det_a += a[2][0]*a[0][1]*a[1][2];
    det_a -= a[0][2]*a[1][1]*a[2][0];
    det_a -= a[0][0]*a[1][2]*a[2][1];
    det_a -= a[0][1]*a[1][0]*a[2][2];

    return det_a;
}

void charm_matr3_inv(charm_real_t a[3][3], charm_real_t a_inv[3][3])
{
    charm_real_t a_[3][3];
    int i, j;
    charm_real_t det_a = charm_matr3_det(a);

    CHARM_ASSERT(det_a != 0.);

    a_[0][0] =  a[1][1]*a[2][2]-a[1][2]*a[2][1];
    a_[0][1] = -a[0][1]*a[2][2]+a[0][2]*a[2][1];
    a_[0][2] =  a[0][1]*a[1][2]-a[0][2]*a[1][1];

    a_[1][0] = -a[1][0]*a[2][2]+a[1][2]*a[2][0];
    a_[1][1] =  a[0][0]*a[2][2]-a[0][2]*a[2][0];
    a_[1][2] = -a[0][0]*a[1][2]+a[0][2]*a[1][0];

    a_[2][0] =  a[1][0]*a[2][1]-a[1][1]*a[2][0];
    a_[2][1] = -a[0][0]*a[2][1]+a[0][1]*a[2][0];
    a_[2][2] =  a[0][0]*a[1][1]-a[0][1]*a[1][0];

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            a_inv[i][j] = a_[i][j]/det_a;
        }
    }
}

void charm_matr_inv(charm_matr_t a_src, charm_matr_t am)
{
    int	       *mask;
    charm_real_t	    fmaxval;
    int		    maxind;
    int		    tmpi;
    charm_real_t	    tmp;
    charm_real_t	  **a;
    int         N = CHARM_BASE_FN_COUNT;
    int         i, j, ni, nj;

    mask = (int*)malloc(N*sizeof(int));//   new int[N];
    a    = (charm_real_t**)malloc(N*sizeof(charm_real_t*)); //new charm_real_t*[N];
    for (i = 0; i < N; i++)
    {
        a[i] = (charm_real_t*)malloc(N*sizeof(charm_real_t)); //new charm_real_t[N];
        for (j = 0; j < N; j++)
        {
            a[i][j] = a_src[i][j];
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i == j)
            {
                am[i][j] = 1.0;
            }
            else {
                am[i][j] = 0.0;
            }
        }
    }
    for (i = 0; i < N; i++)
    {
        mask[i] = i;
    }
    for (i = 0; i < N; i++)
    {
        maxind = i;
        fmaxval = fabs(a[i][i]);
        for (ni = i + 1; ni < N; ni++)
        {
            if (fabs(fmaxval) <= fabs(a[ni][i]))
            {
                fmaxval = fabs(a[ni][i]);
                maxind = ni;
            }
        }
        fmaxval = a[maxind][i];
        CHARM_ASSERT(fmaxval != 0);
        if (i != maxind)
        {
            for (nj = 0; nj < N; nj++)
            {
                tmp = a[i][nj];
                a[i][nj] = a[maxind][nj];
                a[maxind][nj] = tmp;

                tmp = am[i][nj];
                am[i][nj] = am[maxind][nj];
                am[maxind][nj] = tmp;
            }
            tmpi = mask[i];
            mask[i] = mask[maxind];
            mask[maxind] = tmpi;
        }
        charm_real_t aii = a[i][i];
        for (j = 0; j < N; j++)
        {
            a[i][j] = a[i][j] / aii;
            am[i][j] = am[i][j] / aii;
        }
        for (ni = 0; ni < N; ni++)
        {
            if (ni != i)
            {
                charm_real_t fconst = a[ni][i];
                for (nj = 0; nj < N; nj++)
                {
                    a[ni][nj] = a[ni][nj] - fconst *  a[i][nj];
                    am[ni][nj] = am[ni][nj] - fconst * am[i][nj];
                }
            }
        }
    }
    /**/
    for (i = 0; i < N; i++)
    {
        if (mask[i] != i)
        {
            for (j = 0; j < N; j++)
            {
                tmp				= a[i][j];
                a[i][j]			= a[mask[i]][j];
                a[mask[i]][j]	= tmp;
            }
        }
    }
    /**/
    for (i = 0; i < N; i++)
    {
        free(a[i]);
    }
    free(a);
    free(mask);

}


void charm_matr_vect_mult(charm_matr_t a, charm_vect_t b, charm_vect_t res)
{
    int i, j;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        res[i] = 0.;
        for (j = 0; j < CHARM_BASE_FN_COUNT; j++) {
            res[i] += a[i][j]*b[j];
        }
    }
}


void charm_matr_add(charm_matr_t a, charm_matr_t b)
{
    int i, j;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        for (j = 0; j < CHARM_BASE_FN_COUNT; j++) {
            a[i][j] += b[i][j];
        }
    }

}


void charm_matr_zero(charm_matr_t a)
{
    int i, j;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        for (j = 0; j < CHARM_BASE_FN_COUNT; j++) {
            a[i][j] = 0.;
        }
    }

}


void charm_vect_add(charm_vect_t a, charm_vect_t b)
{
    int i, j;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        a[i] += b[i];
    }

}


void charm_vect_zero(charm_vect_t a)
{
    int i;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        a[i] = 0.;
    }

}


charm_data_t * charm_get_quad_data(p4est_quadrant_t *q)
{
    return (charm_data_t *) q->p.user_data;
}


void charm_abort(p4est_t *p4est, int err_code)
{
    int mpiret;

    if (p4est != NULL) {
        charm_write_solution(p4est);
    }

    sc_finalize ();
    MPI_Finalize ();
    exit(1);
}


charm_ctx_t* charm_get_ctx(p4est_t* p4est)
{
    return (charm_ctx_t*)(p4est->user_pointer);
}


charm_comp_t * charm_get_comp(p4est_t * p4est, int i)
{
    charm_ctx_t * ctx = charm_get_ctx(p4est);
    charm_comp_t * comp = sc_array_index(ctx->comp, i);
    return comp;
}


size_t charm_get_comp_count(p4est_t* p4est)
{
    charm_ctx_t * ctx       = charm_get_ctx(p4est);
    return ctx->comp->elem_count;
}


charm_reaction_t * charm_get_reaction(p4est_t * p4est, int i)
{
    charm_ctx_t * ctx = charm_get_ctx(p4est);
    charm_reaction_t * r = sc_array_index(ctx->reactions, i);
    return r;
}


size_t charm_get_reactions_count(p4est_t* p4est)
{
    charm_ctx_t * ctx       = charm_get_ctx(p4est);
    return ctx->reactions == NULL ? 0 : ctx->reactions->elem_count;
}


charm_real_t charm_get_heat_k(p4est_t* p4est, charm_real_t *x, charm_data_t* data)
{
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    size_t c_count = charm_get_comp_count(p4est);
    charm_comp_t *comp;
    charm_cons_t cons;
    charm_prim_t prim;
    charm_real_t kt;
    int i;

    charm_get_fields(data, x, &cons);
    charm_param_cons_to_prim(p4est, &prim, &cons);
    kt = 0.;
    for (i = 0; i < c_count; i++) {
        comp = charm_get_comp(p4est, i);
        kt += prim.c[i]*charm_comp_calc_kp(comp, prim.t);
    }
    return kt;


}


void charm_tensor_add(charm_tensor_t * dest, charm_tensor_t *src)
{
    dest->xx += src->xx;
    dest->yy += src->yy;
    dest->zz += src->zz;
    dest->xy += src->xy;
    dest->yz += src->yz;
    dest->xz += src->xz;
}


void charm_tensor_sum(charm_tensor_t * t1, charm_tensor_t *t2, charm_tensor_t * result)
{
    result->xx = t1->xx + t2->xx;
    result->yy = t1->yy + t2->yy;
    result->zz = t1->zz + t2->zz;
    result->xy = t1->xy + t2->xy;
    result->yz = t1->yz + t2->yz;
    result->xz = t1->xz + t2->xz;
}


void charm_tensor_mul_scalar(charm_tensor_t * dest, charm_real_t x)
{
    dest->xx *= x;
    dest->yy *= x;
    dest->zz *= x;
    dest->xy *= x;
    dest->yz *= x;
    dest->xz *= x;
}


void charm_tensor_zero(charm_tensor_t * t)
{
    memset(t, 0, sizeof(charm_tensor_t));
}


charm_real_t charm_comp_calc_cp(charm_comp_t * comp, charm_real_t t)
{
    int i;
    charm_real_t res = 0.;
    charm_real_t tt  = 1.;
    charm_real_t *cp;
    if (comp->cp_type == COMP_ML_CONST) {
        cp = sc_array_index(comp->cp, 0);
        return *cp;
    }
    else if (comp->ml_type == COMP_CP_POLYNOM) {
        for (i= 0; i < comp->cp->elem_count; i++) {
            cp = sc_array_index(comp->cp, i);
            res += tt*(*cp);
            tt *= t;
        }
        return res;
    }
    else {
        cp = sc_array_index(comp->cp, 0);
        return *cp;
    }
}


charm_real_t charm_comp_calc_cp_dt(charm_comp_t * comp, charm_real_t t)
{
    CHARM_ASSERT(comp->cp_type == COMP_CP_POLYNOM);

    int i;
    charm_real_t res = 0.;
    charm_real_t tt  = 1.;
    charm_real_t *cp;
    if (comp->cp_type == COMP_CP_POLYNOM) {
        for (i= 1; i < comp->cp->elem_count; i++) {
            cp = sc_array_index(comp->cp, i);
            res += tt*(*cp)*i;
            tt *= t;
        }
        return res;
    }
    else {
        CHARM_GLOBAL_LERROR("Wrong call of function 'CALC_CP_DT'\n");
        charm_abort(NULL, 1);
    }
}


charm_real_t charm_comp_calc_ml(charm_comp_t * comp, charm_real_t t)
{
    if (comp->ml_type == COMP_ML_CONST) {
        return comp->ml0;
    }
    else if (comp->ml_type == COMP_ML_SATHERLAND) {
        return comp->ml0 * sqrt( pow( t / comp->t0, 3 ) ) * ( comp->t0 + comp->ts ) / ( t + comp->ts );
    }
    else {
        return comp->ml0;
    }
}


charm_real_t charm_comp_calc_kp(charm_comp_t * comp, charm_real_t t)
{
    if (comp->kp_type == COMP_KP_CONST) {
        return comp->kp0;
    }
    else if (comp->kp_type == COMP_KP_SATHERLAND) {
        return comp->kp0 * sqrt( pow( t / comp->t0, 3 ) ) * ( comp->t0 + comp->ts ) / ( t + comp->ts );
    }
    else {
        return comp->kp0;
    }
}


charm_real_t charm_comp_calc_enthalpy(charm_comp_t * comp, charm_real_t t)
{
    p4est_t *p4est = charm_get_p4est();
    charm_ctx_t *ctx = charm_get_ctx(p4est);
    charm_real_t t1,t2, *cp;
    charm_real_t t_ref = ctx->model.ns.t_ref;
    int i;
    charm_real_t h = comp->h0;
    for (i = 0; i < comp->cp->elem_count; i++){
        cp = (charm_real_t*) sc_array_index(comp->cp, i);
        t2 = pow(t, i+1)/(i+1);
        t1 = pow(t_ref,(i+1))/(i+1);
        h += (*cp)*(t2-t1);
    }
    return h;
}


