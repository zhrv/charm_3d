//
// Created by zhrv on 19.10.17.
//

#ifndef CHAMR_3D_CHARM_GLOBALS_H
#define CHAMR_3D_CHARM_GLOBALS_H

#include <p4est_to_p8est.h>

#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>

#include "charm_def.h"

typedef double              charm_real_t;
typedef int                 charm_int_t;
typedef unsigned int        charm_uint_t;

typedef enum {
    COMP_CP_CONST,
    COMP_CP_POLYNOM
} charm_comp_cp_type_t;

typedef enum {
    COMP_KP_CONST,
    COMP_KP_SATHERLAND,
    COMP_KP_NV
} charm_comp_kp_type_t;

typedef enum {
    COMP_ML_CONST,
    COMP_ML_SATHERLAND,
    COMP_ML_NV
} charm_comp_ml_type_t;


typedef struct charm_comp
{
    char    name[64];
    int     id;
    charm_real_t  m;
    charm_real_t  ml0; //<! динамическая вязкость вещества при температуре T0
    charm_real_t  kp0; //<! теплопроводность вещества при температуре T0
    charm_real_t  t0;  //<! константа для формулы Сазерленда
    charm_real_t  ts;  //<! константа для формулы Сазерленда
    charm_real_t  sig; //<! параметры Леннарда-Джонса
    charm_real_t  ek;  //<! параметры Леннарда-Джонса
    charm_real_t  h0;  //<! энтальпия образования вещества
    charm_comp_cp_type_t cp_type;
    charm_comp_kp_type_t kp_type;
    charm_comp_ml_type_t ml_type;
    sc_array_t *cp;
} charm_comp_t;




typedef struct charm_prim
{
    charm_real_t          r;             /**< density        */
    charm_real_t          u;             /**< velosity       */
    charm_real_t          v;             /**< velosity       */
    charm_real_t          w;             /**< velosity       */
    charm_real_t          e;             /**< energy         */
    charm_real_t          e_tot;         /**< total energy   */
    charm_real_t          p;             /**< pressure       */
    charm_real_t          t;             /**< temperature    */
    charm_real_t          cz;            /**< sound velosity */
    charm_real_t          gam;
    charm_real_t          cp;
    charm_real_t          cv;
    charm_real_t          m;
    int                   mat_id;
    charm_real_t          c[CHARM_MAX_COMPONETS_COUNT]; // concentrations
} charm_prim_t;


typedef struct charm_cons
{
    charm_real_t          ru;
    charm_real_t          rv;
    charm_real_t          rw;
    charm_real_t          re;
    charm_real_t          rc[CHARM_MAX_COMPONETS_COUNT];
    int             mat_id;
} charm_cons_t;

typedef struct charm_tensor
{
    charm_real_t xx;
    charm_real_t yy;
    charm_real_t zz;
    charm_real_t xy;
    charm_real_t xz;
    charm_real_t yz;
} charm_tensor_t;

typedef struct charm_tensor_c
{
    charm_real_t xx[CHARM_BASE_FN_COUNT];
    charm_real_t yy[CHARM_BASE_FN_COUNT];
    charm_real_t zz[CHARM_BASE_FN_COUNT];
    charm_real_t xy[CHARM_BASE_FN_COUNT];
    charm_real_t xz[CHARM_BASE_FN_COUNT];
    charm_real_t yz[CHARM_BASE_FN_COUNT];
} charm_tensor_c_t;

typedef struct charm_vec {
    charm_real_t x[CHARM_BASE_FN_COUNT];
    charm_real_t y[CHARM_BASE_FN_COUNT];
    charm_real_t z[CHARM_BASE_FN_COUNT];
} charm_vec_t;

typedef struct charm_vec_c {
    charm_real_t x[CHARM_BASE_FN_COUNT];
    charm_real_t y[CHARM_BASE_FN_COUNT];
    charm_real_t z[CHARM_BASE_FN_COUNT];
} charm_vec_c_t;


typedef void (*charm_eos_fn_t) (p4est_t * p4est, charm_prim_t * p, int flag);

typedef struct charm_mat
{
    int                 id;
    char                name[64];
    charm_eos_fn_t      eos_fn;
} charm_mat_t;

typedef struct charm_reg
{
    char            name[64];
    int             id;
    int             mat_id;
    charm_real_t          v[CHARM_DIM];
    charm_real_t          t;
    charm_real_t          p;
    charm_real_t          c[CHARM_MAX_COMPONETS_COUNT];
    charm_real_t          grav[CHARM_DIM];
} charm_reg_t;

typedef struct charm_param
{
    struct
    {
        charm_real_t          ru[CHARM_BASE_FN_COUNT];             /**< the state variable */
        charm_real_t          rv[CHARM_BASE_FN_COUNT];             /**< the state variable */
        charm_real_t          rw[CHARM_BASE_FN_COUNT];             /**< the state variable */
        charm_real_t          re[CHARM_BASE_FN_COUNT];             /**< the state variable */
        charm_real_t          rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];             /**< the state variable */
    } c;

    struct
    {
        charm_real_t          ru[CHARM_BASE_FN_COUNT];             /**< the state variable */
        charm_real_t          rv[CHARM_BASE_FN_COUNT];             /**< the state variable */
        charm_real_t          rw[CHARM_BASE_FN_COUNT];             /**< the state variable */
        charm_real_t          re[CHARM_BASE_FN_COUNT];             /**< the state variable */
        charm_real_t          rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];             /**< the state variable */
    } c_old;

    union {
        struct {

        } euler;
        struct {
            charm_tensor_c_t tau;
            charm_vec_c_t q;
            charm_real_t d[CHARM_MAX_COMPONETS_COUNT];
            charm_real_t chem_rhs;
            struct {
                charm_real_t k;
                charm_real_t w;
                charm_real_t mu_t;

            } turb;
        } ns;
    } model;

    struct geom
    {
        charm_real_t          n[CHARM_FACES][CHARM_DIM];
        charm_real_t          face_gp[CHARM_FACES][CHARM_FACE_GP_COUNT][CHARM_DIM];
        charm_real_t          face_gw[CHARM_FACES][CHARM_FACE_GP_COUNT];
        charm_real_t          face_gj[CHARM_FACES][CHARM_FACE_GP_COUNT];
        charm_real_t          quad_gp[CHARM_QUAD_GP_COUNT][CHARM_DIM];
        charm_real_t          quad_gw[CHARM_QUAD_GP_COUNT];
        charm_real_t          quad_gj[CHARM_QUAD_GP_COUNT];
        charm_real_t          area[CHARM_FACES];
        charm_real_t          volume;
        charm_real_t          c[CHARM_DIM];
        charm_real_t          fc[CHARM_FACES][CHARM_DIM];
        charm_real_t          dh[CHARM_DIM];
        charm_real_t          a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT];
        charm_real_t          a_inv[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT];
        charm_real_t          y;
    } g;

    int         mat_id;
    charm_real_t      grav[CHARM_DIM];


    struct lim
    {
        int         count;
        charm_real_t      ru[CHARM_FACES+1];
        charm_real_t      rv[CHARM_FACES+1];
        charm_real_t      rw[CHARM_FACES+1];
        charm_real_t      re[CHARM_FACES+1];
        charm_real_t      rc[CHARM_MAX_COMPONETS_COUNT][CHARM_FACES+1];             /**< the state variable */
    } l;

    struct amr {
        charm_real_t      grad_u[CHARM_DIM];
    } a;

} charm_param_t;


typedef struct charm_data
{
    charm_param_t       par;
    charm_real_t              int_ru[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    charm_real_t              int_rv[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    charm_real_t              int_rw[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    charm_real_t              int_re[CHARM_BASE_FN_COUNT];          /**< the time derivative */
    charm_real_t              int_rc[CHARM_MAX_COMPONETS_COUNT][CHARM_BASE_FN_COUNT];              /**< the time derivative */

    charm_real_t              int_q_x[CHARM_BASE_FN_COUNT];
    charm_real_t              int_q_y[CHARM_BASE_FN_COUNT];
    charm_real_t              int_q_z[CHARM_BASE_FN_COUNT];

    charm_real_t              int_tau_xx[CHARM_BASE_FN_COUNT];
    charm_real_t              int_tau_yy[CHARM_BASE_FN_COUNT];
    charm_real_t              int_tau_zz[CHARM_BASE_FN_COUNT];
    charm_real_t              int_tau_xy[CHARM_BASE_FN_COUNT];
    charm_real_t              int_tau_xz[CHARM_BASE_FN_COUNT];
    charm_real_t              int_tau_yz[CHARM_BASE_FN_COUNT];
    int                 ref_flag;
} charm_data_t;

typedef void    (*charm_limiter_fn_t)               (p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data);
typedef void    (*charm_bnd_cond_fn_t)              (charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, charm_real_t* param, charm_real_t* n);
typedef void    (*charm_flux_fn_t)                  (p4est_t *p4est, charm_prim_t prim[2], charm_real_t* qu, charm_real_t* qv, charm_real_t* qw, charm_real_t* qe, charm_real_t qc[], charm_real_t n[3]);
typedef void    (*charm_timestep_single_fn_t)       (p4est_t * p4est, charm_real_t *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data);
typedef charm_real_t  (*charm_get_timestep_fn_t)    (p4est_t * p4est);
typedef void    (*charm_turb_model_fn_t)            (p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);

#ifndef GLOBALS_H_FILE
extern const char *charm_bnd_types[];
extern const char *charm_turb_models[];
#endif

typedef enum {
    TURB_MODEL_SST,
    TURB_MODEL_UNKNOWN
} charm_turb_models_t;

typedef enum {
    BOUND_INLET,
    BOUND_OUTLET,
    BOUND_WALL_SLIP,
    BOUND_WALL_NO_SLIP,
    BOUND_MASS_FLOW,
    BOUND_SYMMETRY,
    BOUND_UNKNOWN
} charm_bnd_types_t;

typedef struct charm_bnd
{
    char name[64];
    charm_bnd_types_t type;
    charm_real_t *params;
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

typedef struct charm_reaction
{
    int     left_comps[3];       /**< номера компонент реакции справа (int) */
    int     right_comps[3];      /**< номера компонент реакции слева  (int)*/
    charm_real_t  a;                   /**< предэкспоненциальный множитель */
    charm_real_t  e;                   /**< энергия активации */
    charm_real_t  n;                   /**< степень температуры */
} charm_reaction_t;

typedef struct charm_ctx
{
    charm_real_t              max_err;            /**< maximum allowed global interpolation error */
    int                 refine_period;      /**< the number of time steps between mesh refinement */
    int                 repartition_period; /**< the number of time steps between repartitioning */
    int                 write_period;       /**< the number of time steps between writing vtk files */
    int                 log_period;         /**< the number of time steps between writing log */
    int                 min_level;          /**< the minimal level */
    int                 max_level;          /**< the allowed level */
    charm_real_t              CFL;                /**< the CFL */
    charm_real_t              dt;
    charm_real_t              t;                  /**< the current time */
    charm_real_t              time;               /**< the max time */
    int                 timestep;

    union {
        struct {

        } euler;
        
        struct {
            int                         use_visc;
            int                         use_diff;
            charm_real_t                t_ref;
            struct {
                charm_turb_model_fn_t       model_fn;
                charm_turb_models_t         model_type;
                union {
                    struct {
                        charm_real_t a1;
                        charm_real_t sigma_k1;
                        charm_real_t sigma_k2;
                        charm_real_t sigma_w1;
                        charm_real_t sigma_w2;
                        charm_real_t beta_star;
                        charm_real_t beta_1;
                        charm_real_t beta_2;
                        charm_real_t kappa;
                    } sst;

                    struct {
                        double a1;
                    } sa;
                } param;
            } turb;

        } ns;
    } model;
//    charm_real_t              visc_m;
//    charm_real_t              visc_l;

    sc_array_t         *bnd;
    sc_array_t         *mat;       /**< materials */
    sc_array_t         *reg;       /**< regions */
    sc_array_t         *comp;      /**< components */
    sc_array_t         *reactions; /**< reactions */

    charm_mesh_info_t          *msh;
    charm_timestep_single_fn_t  timestep_single_fn;
    charm_get_timestep_fn_t     get_dt_fn;
    charm_flux_fn_t             flux_fn;
    charm_limiter_fn_t          lim_fn;
    
} charm_ctx_t;

typedef struct charm_tree_attr
{
    charm_bnd_t        *bnd[CHARM_FACES];
    charm_reg_t        *reg;
} charm_tree_attr_t;

#ifdef __cplusplus
    extern "C" {
#endif

charm_real_t scalar_prod(charm_real_t v1[CHARM_DIM], charm_real_t v2[CHARM_DIM]);

charm_real_t vect_length(charm_real_t v[CHARM_DIM]);

void vect_prod(charm_real_t v1[CHARM_DIM], charm_real_t v2[CHARM_DIM], charm_real_t res[CHARM_DIM]);
void vect_sub(charm_real_t v1[CHARM_DIM], charm_real_t v2[CHARM_DIM], charm_real_t res[CHARM_DIM]);
charm_real_t vect_dist(charm_real_t v1[CHARM_DIM], charm_real_t v2[CHARM_DIM]);


charm_real_t charm_face_get_area(charm_data_t *d, int8_t face);

charm_real_t charm_face_get_normal(charm_data_t *d, int8_t face, charm_real_t *n);

void charm_quad_get_center(charm_data_t *d, charm_real_t *c);

void charm_face_get_center(charm_data_t *d, int8_t face, charm_real_t *c);

charm_real_t charm_quad_get_volume(charm_data_t *d);

charm_comp_t *charm_get_comp(p4est_t *p4est, int i);

size_t charm_get_comp_count(p4est_t *p4est);

charm_reaction_t *charm_get_reaction(p4est_t *p4est, int i);

size_t charm_get_reactions_count(p4est_t *p4est);

charm_comp_t *charm_comp_find_by_id(charm_ctx_t *ctx, int id);

int charm_comp_index_find_by_id(charm_ctx_t *ctx, int id, size_t *index);

charm_mat_t *charm_mat_find_by_id(charm_ctx_t *ctx, int id);

int charm_mat_index_find_by_id(charm_ctx_t *ctx, int id, size_t *index);

charm_reg_t *charm_reg_find_by_id(charm_ctx_t *ctx, int id);

charm_bnd_t *charm_bnd_find_by_face_type(charm_ctx_t *ctx, int type);

charm_mesh_type_t charm_mesh_get_type_by_str(char *);

charm_tree_attr_t *charm_get_tree_attr(p4est_t *p4est, p4est_topidx_t which_tree);

charm_data_t *charm_get_quad_data(p4est_quadrant_t *q);


void charm_param_cons_to_prim(p4est_t *p4est, charm_prim_t *p, charm_cons_t *c);

void charm_param_prim_to_cons(p4est_t *p4est, charm_cons_t *c, charm_prim_t *p);

void charm_prim_cpy(charm_prim_t *dest, charm_prim_t *src);

charm_real_t charm_matr3_det(charm_real_t a[3][3]);

void charm_matr3_inv(charm_real_t a[3][3], charm_real_t a_inv[3][3]);

void charm_matr_inv(charm_real_t a_src[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT],
                    charm_real_t am[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT]);

void charm_matr_vect_mult(charm_real_t a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT], charm_real_t b[CHARM_BASE_FN_COUNT],
                          charm_real_t res[CHARM_BASE_FN_COUNT]);

void
charm_matr_add(charm_real_t a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT], charm_real_t b[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT]);

void charm_vect_add(charm_real_t a[CHARM_BASE_FN_COUNT], charm_real_t b[CHARM_BASE_FN_COUNT]);

void charm_matr_zero(charm_real_t a[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT]);

void charm_vect_zero(charm_real_t a[CHARM_BASE_FN_COUNT]);

charm_ctx_t *charm_get_ctx(p4est_t *p4est);

void charm_abort(p4est_t *p4est, int err_code);

void dbg_print_param(charm_param_t *);


extern int charm_package_id;

void charm_timesteps(p4est_t *p4est);

void charm_init_initial_condition(p4est_t *p4est, p4est_topidx_t which_tree,
                                  p4est_quadrant_t *q);

void charm_init_context(charm_ctx_t *ctx);


void charm_write_solution(p4est_t *p4est);

void charm_log_statistics(p4est_t *p4est, int timestep, charm_real_t time, charm_real_t dt, charm_real_t calc_time);


void charm_quad_get_vertices(p4est_t *p4est, p4est_quadrant_t *q, p4est_topidx_t treeid, charm_real_t v[8][CHARM_DIM]);

void charm_geom_quad_calc(p4est_t *p4est, p4est_quadrant_t *q, p4est_topidx_t treeid);


p4est_connectivity_t *charm_conn_create(charm_ctx_t *ctx);

charm_real_t charm_get_heat_k(p4est_t *p4est, charm_real_t *x, charm_data_t *data);

charm_real_t charm_get_visc_lambda(p4est_t *p4est, charm_data_t *data);

charm_real_t charm_get_visc_mu(p4est_t *p4est, charm_real_t *x, charm_data_t *data);

void charm_tensor_zero(charm_tensor_t *t);

void charm_tensor_add(charm_tensor_t *dest, charm_tensor_t *src);

void charm_tensor_sum(charm_tensor_t *t1, charm_tensor_t *t2, charm_tensor_t *result);

void charm_tensor_mul_scalar(charm_tensor_t *dest, charm_real_t x);

p4est_t *charm_get_p4est();

void charm_set_p4est(p4est_t *);


charm_real_t charm_comp_calc_cp(charm_comp_t *comp, charm_real_t t);

charm_real_t charm_comp_calc_cp_dt(charm_comp_t *comp, charm_real_t t);

charm_real_t charm_comp_calc_ml(charm_comp_t *comp, charm_real_t t);

charm_real_t charm_comp_calc_kp(charm_comp_t *comp, charm_real_t t);

charm_real_t charm_comp_calc_enthalpy(charm_comp_t *comp, charm_real_t t);


#ifdef __cplusplus
}
#endif


#endif //CHAMR_3D_CHARM_GLOBALS_H

