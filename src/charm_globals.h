//
// Created by zhrv on 19.10.17.
//

#ifndef CHAMR_3D_CHARM_GLOBALS_H
#define CHAMR_3D_CHARM_GLOBALS_H

//#define CHARM_DEBUG

#define SECOND_ORDER

#define FLUX_RIM
//#define FLUX_LF

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

typedef struct param
{
    struct prim
    {
        double          ro;             /**< the state variable */
        double          ru;             /**< the state variable */
        double          rv;             /**< the state variable */
        double          rw;             /**< the state variable */
        double          re;             /**< the state variable */
    } p;

    struct cons
    {
        double          r;             /**< density */
        double          u;             /**< velosity */
        double          v;             /**< velosity */
        double          w;             /**< velosity */
        double          e;             /**< energy */
        double          e_tot;         /**< total energy */
        double          p;             /**< pressure */
        double          t;             /**< temperature */
    } c;
} param_t;


typedef struct charm_data
{
    param_t             par;
    double              drodt;          /**< the time derivative */
    double              drudt;          /**< the time derivative */
    double              drvdt;          /**< the time derivative */
    double              drwdt;          /**< the time derivative */
    double              dredt;          /**< the time derivative */

    double              dro[P4EST_DIM];
    double              dru[P4EST_DIM];
    double              drv[P4EST_DIM];
    double              drw[P4EST_DIM];
    double              dre[P4EST_DIM];

    int                 ref_flag;
} charm_data_t;

typedef struct charm_ctx
{
    double              max_err;            /**< maximum allowed global interpolation error */
    int                 refine_period;      /**< the number of time steps between mesh refinement */
    int                 repartition_period; /**< the number of time steps between repartitioning */
    int                 write_period;       /**< the number of time steps between writing vtk files */
    int                 allowed_level;      /**< the allowed level */
    int                 min_level;          /**< the minimal level */
    double              v_ref;              /**< the reference velosity */
    double              CFL;                /**< the CFL */
    double              dt;
    double              time;               /**< the max time */
} charm_ctx_t;

typedef struct charm_tree_attr
{
    int                 bnd_type[P4EST_FACES];
    int                 region;
} charm_tree_attr_t;




#endif //CHAMR_3D_CHARM_GLOBALS_H
