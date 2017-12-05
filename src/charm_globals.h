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

#ifdef CHARM_DEBUG
#define CHARM_LOG_LEVEL SC_LP_ESSENTIAL
#else
#define CHARM_LOG_LEVEL SC_LP_ESSENTIAL
#endif

/* log helper macros */
#define CHARM_GLOBAL_LOG(p,s)                           \
  SC_GEN_LOG (charm_package_id, SC_LC_GLOBAL, (p), (s))
#define CHARM_LOG(p,s)                                  \
  SC_GEN_LOG (charm_package_id, SC_LC_NORMAL, (p), (s))
void                CHARM_GLOBAL_LOGF (int priority, const char *fmt, ...)
__attribute__ ((format (printf, 2, 3)));
void                CHARM_LOGF (int priority, const char *fmt, ...)
__attribute__ ((format (printf, 2, 3)));
#ifndef __cplusplus
#define CHARM_GLOBAL_LOGF(p,f,...)                                      \
  SC_GEN_LOGF (charm_package_id, SC_LC_GLOBAL, (p), (f), __VA_ARGS__)
#define CHARM_LOGF(p,f,...)                                             \
  SC_GEN_LOGF (charm_package_id, SC_LC_NORMAL, (p), (f), __VA_ARGS__)
#endif

/* convenience global log macros will only print if identifier <= 0 */
#define CHARM_GLOBAL_TRACE(s) CHARM_GLOBAL_LOG (SC_LP_TRACE, (s))
#define CHARM_GLOBAL_LDEBUG(s) CHARM_GLOBAL_LOG (SC_LP_DEBUG, (s))
#define CHARM_GLOBAL_VERBOSE(s) CHARM_GLOBAL_LOG (SC_LP_VERBOSE, (s))
#define CHARM_GLOBAL_INFO(s) CHARM_GLOBAL_LOG (SC_LP_INFO, (s))
#define CHARM_GLOBAL_STATISTICS(s) CHARM_GLOBAL_LOG (SC_LP_STATISTICS, (s))
#define CHARM_GLOBAL_PRODUCTION(s) CHARM_GLOBAL_LOG (SC_LP_PRODUCTION, (s))
#define CHARM_GLOBAL_ESSENTIAL(s) CHARM_GLOBAL_LOG (SC_LP_ESSENTIAL, (s))
#define CHARM_GLOBAL_LERROR(s) CHARM_GLOBAL_LOG (SC_LP_ERROR, (s))
void                CHARM_GLOBAL_TRACEF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_LDEBUGF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_VERBOSEF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_INFOF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_STATISTICSF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_PRODUCTIONF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_ESSENTIALF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_GLOBAL_LERRORF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
#ifndef __cplusplus
#define CHARM_GLOBAL_TRACEF(f,...)                      \
  CHARM_GLOBAL_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define CHARM_GLOBAL_LDEBUGF(f,...)                     \
  CHARM_GLOBAL_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define CHARM_GLOBAL_VERBOSEF(f,...)                    \
  CHARM_GLOBAL_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define CHARM_GLOBAL_INFOF(f,...)                       \
  CHARM_GLOBAL_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define CHARM_GLOBAL_STATISTICSF(f,...)                         \
  CHARM_GLOBAL_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define CHARM_GLOBAL_PRODUCTIONF(f,...)                         \
  CHARM_GLOBAL_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define CHARM_GLOBAL_ESSENTIALF(f,...)                          \
  CHARM_GLOBAL_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define CHARM_GLOBAL_LERRORF(f,...)                     \
  CHARM_GLOBAL_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
#endif
#define CHARM_GLOBAL_NOTICE     CHARM_GLOBAL_STATISTICS
#define CHARM_GLOBAL_NOTICEF    CHARM_GLOBAL_STATISTICSF

/* convenience log macros that are active on every processor */
#define CHARM_TRACE(s) CHARM_LOG (SC_LP_TRACE, (s))
#define CHARM_LDEBUG(s) CHARM_LOG (SC_LP_DEBUG, (s))
#define CHARM_VERBOSE(s) CHARM_LOG (SC_LP_VERBOSE, (s))
#define CHARM_INFO(s) CHARM_LOG (SC_LP_INFO, (s))
#define CHARM_STATISTICS(s) CHARM_LOG (SC_LP_STATISTICS, (s))
#define CHARM_PRODUCTION(s) CHARM_LOG (SC_LP_PRODUCTION, (s))
#define CHARM_ESSENTIAL(s) CHARM_LOG (SC_LP_ESSENTIAL, (s))
#define CHARM_LERROR(s) CHARM_LOG (SC_LP_ERROR, (s))
void                CHARM_TRACEF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_LDEBUGF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_VERBOSEF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_INFOF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_STATISTICSF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_PRODUCTIONF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_ESSENTIALF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
void                CHARM_LERRORF (const char *fmt, ...)
__attribute__ ((format (printf, 1, 2)));
#ifndef __cplusplus
#define CHARM_TRACEF(f,...)                     \
  CHARM_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define CHARM_LDEBUGF(f,...)                    \
  CHARM_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define CHARM_VERBOSEF(f,...)                   \
  CHARM_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define CHARM_INFOF(f,...)                      \
  CHARM_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define CHARM_STATISTICSF(f,...)                        \
  CHARM_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define CHARM_PRODUCTIONF(f,...)                        \
  CHARM_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define CHARM_ESSENTIALF(f,...)                         \
  CHARM_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define CHARM_LERRORF(f,...)                    \
  CHARM_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
#endif
#define CHARM_NOTICE            CHARM_STATISTICS
#define CHARM_NOTICEF           CHARM_STATISTICSF



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


typedef void (*charm_bnd_cond_fn_t)(charm_param_t *par_in, charm_param_t *par_out, int8_t face, double* param, double* n);



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
    int                 log_period;         /**< the number of time steps between writing log */
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


extern int charm_package_id;

#endif //CHAMR_3D_CHARM_GLOBALS_H
