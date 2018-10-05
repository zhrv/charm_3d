//
// Created by zhrv on 26.10.17.
//
#define GLOBAL_H_FILE
#include "charm_globals.h"



int charm_package_id = -1;

const char *charm_bnd_types[] ={
        "BOUND_INLET",
        "BOUND_OUTLET",
        "BOUND_WALL_SLIP",
        "BOUND_WALL_NO_SLIP",
        NULL
};


double scalar_prod(double v1[CHARM_DIM], double v2[CHARM_DIM])
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}


double vect_length(double v[CHARM_DIM])
{
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

void vect_prod(double v1[CHARM_DIM], double v2[CHARM_DIM], double res[CHARM_DIM])
{
    res[0] =  v1[1]*v2[2]-v1[2]*v2[1];
    res[1] = -v1[0]*v2[2]+v1[2]*v2[0];
    res[2] =  v1[0]*v2[1]-v1[1]*v2[0];
}


double charm_face_get_area(charm_data_t *d, int8_t face)
{
    return d->par.g.area[face];
}

double charm_face_get_normal(charm_data_t *d, int8_t face, double* n)
{
    memcpy(n, d->par.g.n[face], CHARM_DIM*sizeof(double));
    return d->par.g.area[face];
}

void charm_quad_get_center(charm_data_t *d, double* c)
{
    memcpy(c, d->par.g.c, 3*sizeof(double));
}

void charm_face_get_center(charm_data_t *d, int8_t face, double* c)
{
    memcpy(c, d->par.g.fc[face], CHARM_DIM*sizeof(double));
}

double charm_quad_get_volume(charm_data_t *d)
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

double gR = 8.314472;
void charm_mat_eos(p4est_t * p4est, charm_prim_t * p, int flag)
{
//    charm_ctx_t *ctx = (charm_ctx_t *)p4est->user_pointer;
//    charm_reg_t *reg = charm_reg_find_by_id(ctx, p->reg_id);
//    charm_mat_t *mat = charm_mat_find_by_id(ctx, reg->mat_id);
//    double Cp = mat->cp;
//    double M  = mat->m;
//    double Cv = Cp-gR/M;
//    double gam = Cp/Cv;
//    p->gam = gam;
//    p->cp = Cp;
//    p->cv = Cv;
//    switch (flag)
//    {
//        case 0:		// p=p(r,e)
//            p->p = p->r*p->e*(gam-1);
//            p->cz = sqrt(gam*p->p/p->r);
//            break;
//
//        case 1:		// e=e(r,p)
//            p->e = p->p/(p->r*(gam-1));
//            p->t = p->e/Cv;
//            break;
//
//        case 2:		// r=r(T,p)
//            p->r = p->p*M/(p->t*gR);
//            p->cz = sqrt(gam*p->p/p->r);
//            break;
//
//        case 3:
//            p->r  = p->p*M/(p->t*gR);
//            p->cz = sqrt(gam*p->p/p->r);
//            p->e  = p->p/(p->r*(gam-1));
//            break;
//
//        case 4:
//            p->p  = p->r*p->e*(gam-1);
//            p->cz = sqrt(gam*p->p/p->r);
//            p->t  = p->e/Cv;
//            break;
//
//        default:
//            CHARM_ASSERT(flag < 3);
//    }
//
}


void charm_param_cons_to_prim(p4est_t * p4est, charm_prim_t * p, charm_cons_t * c)
{
    p->reg_id = c->reg_id;
    p->r      = c->ro;
    p->u      = c->ru/c->ro;
    p->v      = c->rv/c->ro;
    p->w      = c->rw/c->ro;
    p->e_tot  = c->re/c->ro;
    p->e      = p->e_tot-0.5*(p->u*p->u+p->v*p->v+p->w*p->w);

    p->components_count = c->components_count;
    for (size_t i = 0; i < c->components_count; i++) {
        p->c[i] = c->rc[i] / c->ro;
    }

    charm_mat_eos(p4est, p, 4);  // {p,cz, t}=EOS(r,e)

    if (p->p < 0.) {
        int iii=0;
    }
}


void charm_param_prim_to_cons(p4est_t * p4est, charm_cons_t * c, charm_prim_t * p) {
    c->reg_id = p->reg_id;
    c->ro = p->r;
    c->ru = p->r * p->u;
    c->rv = p->r * p->v;
    c->rw = p->r * p->w;
    c->re = p->r * (p->e + 0.5 * (p->u * p->u + p->v * p->v + p->w * p->w));
    c->components_count = p->components_count;
    for (size_t i = 0; i < c->components_count; i++) {
        c->rc[i] = p->r * p->c[i];
    }
}

void charm_prim_cpy(charm_prim_t * dest, charm_prim_t * src)
{
    dest->reg_id = src->reg_id;
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

}


double charm_matr3_det(double a[3][3])
{
    double det_a = 0.;

    det_a += a[0][0]*a[1][1]*a[2][2];
    det_a += a[0][2]*a[1][0]*a[2][1];
    det_a += a[2][0]*a[0][1]*a[1][2];
    det_a -= a[0][2]*a[1][1]*a[2][0];
    det_a -= a[0][0]*a[1][2]*a[2][1];
    det_a -= a[0][1]*a[1][0]*a[2][2];

    return det_a;
}

void charm_matr3_inv(double a[3][3], double a_inv[3][3])
{
    double a_[3][3];
    int i, j;
    double det_a = charm_matr3_det(a);

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

void charm_matr_inv(double a_src[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT], double am[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT])
{
    int	       *mask;
    double	    fmaxval;
    int		    maxind;
    int		    tmpi;
    double	    tmp;
    double	  **a;
    int         N = CHARM_BASE_FN_COUNT;

    mask = (int*)malloc(N*sizeof(int));//   new int[N];
    a    = (double**)malloc(N*sizeof(double*)); //new double*[N];
    for (int i = 0; i < N; i++)
    {
        a[i] = (double*)malloc(N*sizeof(double)); //new double[N];
        for (int j = 0; j < N; j++)
        {
            a[i][j] = a_src[i][j];
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
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
    for (int i = 0; i < N; i++)
    {
        mask[i] = i;
    }
    for (int i = 0; i < N; i++)
    {
        maxind = i;
        fmaxval = fabs(a[i][i]);
        for (int ni = i + 1; ni < N; ni++)
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
            for (int nj = 0; nj < N; nj++)
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
        double aii = a[i][i];
        for (int j = 0; j < N; j++)
        {
            a[i][j] = a[i][j] / aii;
            am[i][j] = am[i][j] / aii;
        }
        for (int ni = 0; ni < N; ni++)
        {
            if (ni != i)
            {
                double fconst = a[ni][i];
                for (int nj = 0; nj < N; nj++)
                {
                    a[ni][nj] = a[ni][nj] - fconst *  a[i][nj];
                    am[ni][nj] = am[ni][nj] - fconst * am[i][nj];
                }
            }
        }
    }
    /**/
    for (int i = 0; i < N; i++)
    {
        if (mask[i] != i)
        {
            for (int j = 0; j < N; j++)
            {
                tmp				= a[i][j];
                a[i][j]			= a[mask[i]][j];
                a[mask[i]][j]	= tmp;
            }
        }
    }
    /**/
    for (int i = 0; i < N; i++)
    {
        free(a[i]);
    }
    free(a);
    free(mask);

}


void charm_matr_vect_mult(double a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT], double b[CHARM_BASE_FN_COUNT], double res[CHARM_BASE_FN_COUNT])
{
    int i, j;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        res[i] = 0.;
        for (j = 0; j < CHARM_BASE_FN_COUNT; j++) {
            res[i] += a[i][j]*b[j];
        }
    }
}


void charm_matr_add(double a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT], double b[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT])
{
    int i, j;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        for (j = 0; j < CHARM_BASE_FN_COUNT; j++) {
            a[i][j] += b[i][j];
        }
    }

}


void charm_matr_zero(double a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT])
{
    int i, j;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        for (j = 0; j < CHARM_BASE_FN_COUNT; j++) {
            a[i][j] = 0.;
        }
    }

}


void charm_vect_add(double a[CHARM_BASE_FN_COUNT], double b[CHARM_BASE_FN_COUNT])
{
    int i, j;
    for (i = 0; i < CHARM_BASE_FN_COUNT; i++) {
        a[i] += b[i];
    }

}


void charm_vect_zero(double a[CHARM_BASE_FN_COUNT])
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


void charm_abort(int err_code)
{
    int mpiret;

    sc_finalize ();
    MPI_Finalize ();
    exit(1);
}

charm_ctx_t* charm_get_ctx(p4est_t* p4est)
{
    return (charm_ctx_t*)(p4est->user_pointer);
}