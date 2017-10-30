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

//#define CHARM_DEBUG

//#define SECOND_ORDER

//#define FLUX_RIM
#define FLUX_LF

#define RIM_EPS 1.e-5

#define ROOT_LEN   0.04


#ifdef CHARM_DEBUG

#define DBG_CH() printf("Line: %d\n", __LINE__)

#else

#define DBG_CH() ((void)0)

#endif

#define CHARM_GET_H(LEVEL) ROOT_LEN*((double) P4EST_QUADRANT_LEN ((LEVEL)) / (double) P4EST_ROOT_LEN)

#define FLD_COUNT 5

#define GAM         1.4

#define _MAX_(X,Y) ((X)>(Y) ? (X) : (Y))
#define _MIN_(X,Y) ((X)<(Y) ? (X) : (Y))

#define CHARM_FACE_TYPE_INNER 0
#define CHARM_BND_MAX 128

typedef void (*charm_bnd_cond_fn_t)(double ro, double ru, double rv, double rw, double re,
                                    double* ro_, double* ru_, double* rv_, double* rw_, double* re_,
                                    double* n, double* param);



typedef struct param
{
    struct cons
    {
        double          ro;             /**< the state variable */
        double          ru;             /**< the state variable */
        double          rv;             /**< the state variable */
        double          rw;             /**< the state variable */
        double          re;             /**< the state variable */
    } c;

    struct prim
    {
        double          r;             /**< density */
        double          u;             /**< velosity */
        double          v;             /**< velosity */
        double          w;             /**< velosity */
        double          e;             /**< energy */
        double          e_tot;         /**< total energy */
        double          p;             /**< pressure */
        double          t;             /**< temperature */
    } p;

    struct geom
    {
        double          n[P4EST_FACES][P4EST_DIM];
        double          area[P4EST_FACES];
        double          volume;
        double          c[P4EST_DIM];
        double          fc[P4EST_FACES][P4EST_DIM];
    } g;
} param_t;


typedef struct charm_data
{
    param_t             par;
    double              drodt;          /**< the time derivative */
    double              drudt;          /**< the time derivative */
    double              drvdt;          /**< the time derivative */
    double              drwdt;          /**< the time derivative */
    double              dredt;          /**< the time derivative */

//    double              dro[P4EST_DIM];
//    double              dru[P4EST_DIM];
//    double              drv[P4EST_DIM];
//    double              drw[P4EST_DIM];
//    double              dre[P4EST_DIM];

    int                 ref_flag;
} charm_data_t;


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
    int type;
    int id;
    double *params;
    charm_bnd_cond_fn_t bnd_fn;

} charm_bnd_t;

#define FACE_TYPE_COUNT 128

typedef struct charm_ctx
{
    double              max_err;            /**< maximum allowed global interpolation error */
    int                 refine_period;      /**< the number of time steps between mesh refinement */
    int                 repartition_period; /**< the number of time steps between repartitioning */
    int                 write_period;       /**< the number of time steps between writing vtk files */
    int                 allowed_level;      /**< the allowed level */
    int                 min_level;          /**< the minimal level */
    double              CFL;                /**< the CFL */
    double              dt;
    double              time;               /**< the max time */

    charm_bnd_t        *bnd[FACE_TYPE_COUNT];
} charm_ctx_t;

typedef struct charm_tree_attr
{
    int                 bnd_type[P4EST_FACES];
    int                 region;
} charm_tree_attr_t;


double scalar_prod(double v1[3], double v2[3]);
double vect_length(double v[3]);
void vect_prod(double v1[3], double v2[3], double res[3]);


double charm_face_get_normal(p4est_quadrant_t* q, int8_t face, double* n);
void charm_quad_get_center(p4est_quadrant_t* q, double* c);
void charm_face_get_center(p4est_quadrant_t* q, int8_t face, double* c);
double charm_quad_get_volume(p4est_quadrant_t* q);

#endif //CHAMR_3D_CHARM_GLOBALS_H
