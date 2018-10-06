//
// Created by zhrv on 19.10.17.
//

#ifndef CHAMR_3D_CHARM_GLOBALS_H
#define CHAMR_3D_CHARM_GLOBALS_H

#include <p4est_to_p8est.h>

#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>

#define CHARM_DIM           P4EST_DIM
#define CHARM_ALLOC         P4EST_ALLOC
#define CHARM_FREE          P4EST_FREE
#define CHARM_HALF          P4EST_HALF
#define CHARM_CHILDREN      P4EST_CHILDREN
#define CHARM_FACES         P4EST_FACES
#define CHARM_CONNECT_FULL  P4EST_CONNECT_FULL
#define CHARM_CONNECT_FACE  P4EST_CONNECT_FACE
#define CHARM_REALLOC       P4EST_REALLOC
#define CHARM_QUADRANT_LEN  P4EST_QUADRANT_LEN



#ifdef CHARM_DEBUG

#define CHARM_LOG_LEVEL SC_LP_ESSENTIAL
#define DBG_CH(R) {printf("Rank: %d. File: %s. Line: %d\n", (R), __FILE__, __LINE__);fflush(stdout);}
#define CHARM_ASSERT P4EST_ASSERT

#else

#define CHARM_LOG_LEVEL SC_LP_ESSENTIAL
#define DBG_CH(R) ((void)0)
#define CHARM_ASSERT(R) ((void)0)

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

#define CHARM_STRING "charm_dg"

#define CHARM_RIM_NEWTON_STEPS 5000
#define CHARM_RIM_EPS 1.e-5

#define CHARM_EPS 1.e-12

#define FLD_COUNT 5

#define GAM         1.4

#define _MAX_(X,Y) ((X)>(Y) ? (X) : (Y))
#define _MIN_(X,Y) ((X)<(Y) ? (X) : (Y))
#define _SQR_(X) ((X)*(X))
#define _MAG_(X,Y,Z) (_SQR_(X)+_SQR_(Y)+_SQR_(Z))
#define _NORM_(X) ( (fabs(X) <= CHARM_EPS) ? 0. : (X) )

#define CHARM_FACE_TYPE_INNER 0
#define CHARM_BND_MAX 128

#define CHARM_BASE_FN_COUNT 4
#define CHARM_FACE_GP_COUNT 6
#define CHARM_QUAD_GP_COUNT 8

#define CHARM_MAX_COMPONETS_COUNT 128

#define CHARM_ARR_SET_ZERO(A) {int i; for (i = 0; i < CHARM_BASE_FN_COUNT; i++) A[i] = 0.; }

typedef struct charm_prim
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
    int             mat_id;
    double          c[CHARM_MAX_COMPONETS_COUNT]; // concentrations
} charm_prim_t;


typedef struct charm_cons
{
//    double          ro;
    double          ru;
    double          rv;
    double          rw;
    double          re;
    double          rc[CHARM_MAX_COMPONETS_COUNT];
    int             mat_id;
} charm_cons_t;


typedef enum {
    COMP_CONST,
    COMP_POLYNOM
} charm_comp_cp_type_t;

typedef void (*charm_eos_fn_t) (p4est_t * p4est, charm_prim_t * p, int flag);

typedef struct charm_comp
{
    char    name[64];
    int     id;
    double  m;
    double  ml;
    double  lambda;
    double  k;
    charm_comp_cp_type_t cp_type;
    sc_array_t *cp;
} charm_comp_t;

typedef struct charm_mat
{
    int                 id;
    char                name[64];
    charm_eos_fn_t      eos_fn;
//    sc_array_t         *comp_idx; /**< component's indexes in array ctx->comp; type: size_t   */
} charm_mat_t;

typedef struct charm_reg
{
    char            name[64];
    int             id;
    int             mat_id;
    double          v[CHARM_DIM];
    double          t;
    double          p;
    double          c[CHARM_MAX_COMPONETS_COUNT];
} charm_reg_t;


typedef struct charm_param
{
    struct
    {
//        double          ro[CHARM_BASE_FN_COUNT];             /**< the state variable */
        double          ru[CHARM_BASE_FN_COUNT];             /**< the state variable */
        double          rv[CHARM_BASE_FN_COUNT];             /**< the state variable */
        double          rw[CHARM_BASE_FN_COUNT];             /**< the state variable */
        double          re[CHARM_BASE_FN_COUNT];             /**< the state variable */
        double          rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];             /**< the state variable */
    } c;

    struct geom
    {
        double          n[CHARM_FACES][CHARM_DIM];
        double          face_gp[CHARM_FACES][CHARM_FACE_GP_COUNT][CHARM_DIM];
        double          face_gw[CHARM_FACES][CHARM_FACE_GP_COUNT];
        double          face_gj[CHARM_FACES][CHARM_FACE_GP_COUNT];
        double          quad_gp[CHARM_QUAD_GP_COUNT][CHARM_DIM];
        double          quad_gw[CHARM_QUAD_GP_COUNT];
        double          quad_gj[CHARM_QUAD_GP_COUNT];
        double          area[CHARM_FACES];
        double          volume;
        double          c[CHARM_DIM];
        double          fc[CHARM_FACES][CHARM_DIM];
        double          dh[CHARM_DIM];
        double          a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT];
        double          a_inv[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT];
    } g;

    int mat_id;

    struct lim
    {
        int         count;
        double      ro[CHARM_FACES+1];
        double      ru[CHARM_FACES+1];
        double      rv[CHARM_FACES+1];
        double      rw[CHARM_FACES+1];
        double      re[CHARM_FACES+1];
    } l;

    struct amr {
        double      grad_u[CHARM_DIM];
    } a;
} charm_param_t;


typedef struct charm_data
{
    charm_param_t       par;
//    double              int_ro[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    double              int_ru[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    double              int_rv[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    double              int_rw[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    double              int_re[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    double              int_rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];              /**< the time derivative */

    int                 ref_flag;
} charm_data_t;



typedef void (*charm_bnd_cond_fn_t)(charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, double* param, double* n);

typedef void (*charm_flux_fn_t)(charm_prim_t prim[2], double* qr, double* qu, double* qv, double* qw, double* qe, double n[3]);

#ifndef GLOBALS_H_FILE
extern const char *charm_bnd_types[];
#endif

typedef enum {
    BOUND_INLET,
    BOUND_OUTLET,
    BOUND_WALL_SLIP,
    BOUND_WALL_NO_SLIP,
    BOUND_UNKNOWN
} charm_bnd_types_t;

typedef struct charm_bnd
{
    char name[64];
    charm_bnd_types_t type;
    double *params;
    charm_bnd_cond_fn_t bnd_fn;

} charm_bnd_t;

typedef enum {
    CHARM_MESH_UNKNOWN,
    CHARM_MESH_GMSH_MSH,
    CHARM_MESH_GMSH_INP,
    CHARM_MESH_GMSH_UNV,
    CHARM_MESH_OPENFOAM,
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
    sc_array_t         *mat;  /**< materials */
    sc_array_t         *reg;  /**< regions */
    sc_array_t         *comp; /**< components */

    charm_mesh_info_t  *msh;
    charm_flux_fn_t     flux_fn;
} charm_ctx_t;

typedef struct charm_tree_attr
{
    charm_bnd_t        *bnd[CHARM_FACES];
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

charm_comp_t *  charm_get_comp(p4est_t * p4est, int i);
size_t          charm_get_comp_count(p4est_t* p4est);
charm_comp_t *  charm_comp_find_by_id(charm_ctx_t *ctx, int id);
int             charm_comp_index_find_by_id(charm_ctx_t *ctx, int id, size_t *index);
charm_mat_t  *  charm_mat_find_by_id(charm_ctx_t *ctx, int id);
int             charm_mat_index_find_by_id(charm_ctx_t *ctx, int id, size_t *index);
charm_reg_t  *  charm_reg_find_by_id(charm_ctx_t *ctx, int id);
charm_bnd_t  *  charm_bnd_find_by_face_type(charm_ctx_t *ctx, int type);

charm_mesh_type_t charm_mesh_get_type_by_str(char*);

charm_tree_attr_t * charm_get_tree_attr(p4est_t * p4est, p4est_topidx_t which_tree);

charm_data_t * charm_get_quad_data(p4est_quadrant_t *q);


void charm_param_cons_to_prim(p4est_t * p4est, charm_prim_t * p, charm_cons_t * c);
void charm_param_prim_to_cons(p4est_t * p4est, charm_cons_t * c, charm_prim_t * p);

void charm_prim_cpy(charm_prim_t * dest, charm_prim_t * src);

double charm_matr3_det(double a[3][3]);
void   charm_matr3_inv(double a[3][3], double a_inv[3][3]);
void   charm_matr_inv(double a_src[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT], double am[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT]);
void   charm_matr_vect_mult(double a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT], double b[CHARM_BASE_FN_COUNT], double res[CHARM_BASE_FN_COUNT]);
void   charm_matr_add(double a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT], double b[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT]);
void   charm_vect_add(double a[CHARM_BASE_FN_COUNT], double b[CHARM_BASE_FN_COUNT]);
void   charm_matr_zero(double a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT]);
void   charm_vect_zero(double a[CHARM_BASE_FN_COUNT]);

charm_ctx_t* charm_get_ctx(p4est_t* p4est);
void charm_abort(int err_code);

void dbg_print_param(charm_param_t *);


extern int charm_package_id;

#endif //CHAMR_3D_CHARM_GLOBALS_H
