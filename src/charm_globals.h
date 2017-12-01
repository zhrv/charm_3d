//
// Created by zhrv on 19.10.17.
//

#ifndef CHAMR_3D_CHARM_GLOBALS_H
#define CHAMR_3D_CHARM_GLOBALS_H

#include <p4est_to_p8est.h>

#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>

#define CHARM_DIM P4EST_DIM

#define CHARM_DEBUG

//#define SECOND_ORDER
#define CHARM_STRING "charm"

#define CHARM_RIM_NEWTON_STEPS 5000
#define CHARM_RIM_EPS 1.e-5

#define CHARM_EPS 1.e-11

#define ROOT_LEN   0.04


#ifdef CHARM_DEBUG

#define DBG_CH(R) {printf("Rank: %d. File: %s. Line: %d\n", (R), __FILE__, __LINE__);fflush(stdout);}

#else

#define DBG_CH(R) ((void)0)

#endif

#define CHARM_GET_H(LEVEL) ROOT_LEN*((double) P4EST_QUADRANT_LEN ((LEVEL)) / (double) P4EST_ROOT_LEN)

#define FLD_COUNT 5

#define GAM         1.4

#define _MAX_(X,Y) ((X)>(Y) ? (X) : (Y))
#define _MIN_(X,Y) ((X)<(Y) ? (X) : (Y))

#define CHARM_FACE_TYPE_INNER 0
#define CHARM_BND_MAX 128



typedef struct charm_param
{
    struct
    {
        double          ro;             /**< the state variable */
        double          ru;             /**< the state variable */
        double          rv;             /**< the state variable */
        double          rw;             /**< the state variable */
        double          re;             /**< the state variable */
    } c;

    struct
    {
        double          r;             /**< density */
        double          u;             /**< velosity */
        double          v;             /**< velosity */
        double          w;             /**< velosity */
        double          e;             /**< energy */
        double          e_tot;         /**< total energy */
        double          p;             /**< pressure */
        double          t;             /**< temperature */
        double          cz;            /**< sound velosity */
        double          gam;
        double          cp;
        double          cv;
    } p;

    struct
    {
        double          r[CHARM_DIM];             /**< density */
        double          u[CHARM_DIM];             /**< velosity */
        double          v[CHARM_DIM];             /**< velosity */
        double          w[CHARM_DIM];             /**< velosity */
        double          p[CHARM_DIM];             /**< pressure */
    } grad;

    struct geom
    {
        double          n[P4EST_FACES][CHARM_DIM];
        double          area[P4EST_FACES];
        double          volume;
        double          c[CHARM_DIM];
        double          fc[P4EST_FACES][CHARM_DIM];
    } g;
} charm_param_t;


typedef struct charm_data
{
    charm_param_t       par;
    double              drodt;          /**< the time derivative */
    double              drudt;          /**< the time derivative */
    double              drvdt;          /**< the time derivative */
    double              drwdt;          /**< the time derivative */
    double              dredt;          /**< the time derivative */

    int                 ref_flag;
} charm_data_t;



typedef struct charm_mat
{
    char name[64];
    int id;
    double m;
    double cp;
    double ml;
    double lambda;
    double k;
} charm_mat_t;


typedef struct charm_reg
{
    char name[64];
    int id;
    int cell_type;
    charm_mat_t *mat;
    double v[CHARM_DIM];
    double t;
    double p;
} charm_reg_t;


typedef void (*charm_bnd_cond_fn_t)(charm_param_t *par_in, charm_param_t *par_out, int8_t face, double* param);



#ifndef GLOBALS_H_FILE
extern const char *charm_bnd_types[];
#endif

typedef enum {
    BOUND_INLET,
    BOUND_OUTLET,
    BOUND_WALL_SLIP,
    BOUND_WALL_NO_SLIP
} bnd_types_t;

typedef struct charm_bnd
{
    char name[64];
    bnd_types_t type;
    int face_type;
    double *params;
    charm_bnd_cond_fn_t bnd_fn;

} charm_bnd_t;

typedef enum {
    CHARM_MESH_UNKNOWN,
    CHARM_MESH_GMSH_MSH,
    CHARM_MESH_GMSH_INP,
    CHARM_MESH_GMSH_UNV,
    CHARM_MESH_SALOME_UNV,
    CHARM_MESH_TETGEN
} charm_mesh_type_t;

typedef struct charm_mesh_info
{
    charm_mesh_type_t   type;
    char                filename[128];
} charm_mesh_info_t;

typedef struct charm_ctx
{
    double              max_err;            /**< maximum allowed global interpolation error */
    int                 refine_period;      /**< the number of time steps between mesh refinement */
    int                 repartition_period; /**< the number of time steps between repartitioning */
    int                 write_period;       /**< the number of time steps between writing vtk files */
    int                 min_level;          /**< the minimal level */
    int                 max_level;          /**< the allowed level */
    double              CFL;                /**< the CFL */
    double              dt;
    double              time;               /**< the max time */

    sc_array_t         *bnd;
    sc_array_t         *mat;
    sc_array_t         *reg;

    charm_mesh_info_t  *msh;
} charm_ctx_t;

typedef struct charm_tree_attr
{
    charm_bnd_t        *bnd[P4EST_FACES];
    charm_reg_t        *reg;
} charm_tree_attr_t;


double scalar_prod(double v1[CHARM_DIM], double v2[CHARM_DIM]);
double vect_length(double v[CHARM_DIM]);
void vect_prod(double v1[CHARM_DIM], double v2[CHARM_DIM], double res[CHARM_DIM]);


double charm_face_get_area(charm_data_t *d, int8_t face);
double charm_face_get_normal(charm_data_t *d, int8_t face, double* n);
void charm_quad_get_center(charm_data_t *d, double* c);
void charm_face_get_center(charm_data_t *d, int8_t face, double* c);
double charm_quad_get_volume(charm_data_t *d);


charm_mat_t * charm_mat_find_by_id(charm_ctx_t *ctx, int id);
charm_bnd_t * charm_bnd_find_by_face_type(charm_ctx_t *ctx, int type);

charm_mesh_type_t charm_mesh_get_type_by_str(char*);

charm_tree_attr_t * charm_get_tree_attr(p4est_t * p4est, p4est_topidx_t which_tree);

void charm_mat_eos(charm_mat_t * mat, charm_param_t * p, int variant);
void charm_param_cons_to_prim(charm_mat_t * mat, charm_param_t * p);
void charm_param_prim_to_cons(charm_mat_t * mat, charm_param_t * p);

void charm_prim_cpy(charm_param_t * dest, charm_param_t * src);


void dbg_print_param(charm_param_t *);

#endif //CHAMR_3D_CHARM_GLOBALS_H
