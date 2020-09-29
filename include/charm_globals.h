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
typedef double              charm_point_t[CHARM_DIM];
typedef double              charm_vector_t[CHARM_DIM];

typedef double              charm_vect_t[CHARM_BASE_FN_COUNT];
typedef double              charm_matr_t[CHARM_BASE_FN_COUNT][CHARM_BASE_FN_COUNT];

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
    charm_vect_t xx;
    charm_vect_t yy;
    charm_vect_t zz;
    charm_vect_t xy;
    charm_vect_t xz;
    charm_vect_t yz;
} charm_tensor_c_t;

typedef struct charm_vec {
    charm_real_t x;
    charm_real_t y;
    charm_real_t z;
} charm_vec_t;

typedef struct charm_vec_c {
    charm_vect_t x;
    charm_vect_t y;
    charm_vect_t z;
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
    charm_vector_t        v;
    charm_real_t          t;
    charm_real_t          p;
    charm_real_t          c[CHARM_MAX_COMPONETS_COUNT];
    charm_vector_t        grav;
} charm_reg_t;

typedef struct charm_param
{
    struct
    {
        charm_vect_t          ru;
        charm_vect_t          rv;
        charm_vect_t          rw;
        charm_vect_t          re;
        charm_vect_t          rc[CHARM_MAX_COMPONETS_COUNT];
    } c;

    struct
    {
        charm_vect_t          ru;
        charm_vect_t          rv;
        charm_vect_t          rw;
        charm_vect_t          re;
        charm_vect_t          rc[CHARM_MAX_COMPONETS_COUNT];
    } c_old;

    union {
        struct {
            charm_real_t v[3];
        } adv;
        struct {

        } euler;
        struct {
            charm_tensor_c_t tau;
            charm_vec_c_t q;
            charm_real_t d[CHARM_MAX_COMPONETS_COUNT];
            charm_real_t chem_rhs;
            struct {
                charm_real_t mu_t;
                union {
                    struct {
                        charm_real_t    nu_; // nu with tilde
                        charm_vector_t  grad_nu_;
                        charm_real_t    nu; // kinematic viscosity
                        charm_real_t    int_nu_;
                        charm_vector_t  grad_u[CHARM_DIM];
                    } sa;
                    struct {
                        charm_real_t k;
                        charm_real_t w;
                    } sst;
                } model;
            } turb;
        } ns;
    } model;

    struct geom
    {
        charm_vector_t        n[CHARM_FACES];
        charm_point_t         face_gp[CHARM_FACES][CHARM_FACE_GP_COUNT];
        charm_real_t          face_gw[CHARM_FACES][CHARM_FACE_GP_COUNT];
        charm_real_t          face_gj[CHARM_FACES][CHARM_FACE_GP_COUNT];
        charm_point_t         quad_gp[CHARM_QUAD_GP_COUNT];
        charm_real_t          quad_gw[CHARM_QUAD_GP_COUNT];
        charm_real_t          quad_gj[CHARM_QUAD_GP_COUNT];
        charm_real_t          area[CHARM_FACES];
        charm_real_t          volume;
        charm_point_t         c;
        charm_point_t         fc[CHARM_FACES];
        charm_vector_t        dh;
        charm_matr_t          a;
        charm_matr_t          a_inv;
        charm_real_t          y;
    } g;

    int         mat_id;
    charm_vector_t      grav;


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
        charm_vector_t      grad_u;
    } a;

} charm_param_t;


typedef struct charm_data
{
    charm_param_t       par;

    charm_vect_t        int_ru;
    charm_vect_t        int_rv;
    charm_vect_t        int_rw;
    charm_vect_t        int_re;
    charm_vect_t        int_rc[CHARM_MAX_COMPONETS_COUNT];

    charm_vect_t        int_q_x;
    charm_vect_t        int_q_y;
    charm_vect_t        int_q_z;

    charm_vect_t        int_tau_xx;
    charm_vect_t        int_tau_yy;
    charm_vect_t        int_tau_zz;
    charm_vect_t        int_tau_xy;
    charm_vect_t        int_tau_xz;
    charm_vect_t        int_tau_yz;
    int                 ref_flag;
} charm_data_t;

typedef void            (*charm_limiter_fn_t)               (p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data);
typedef void            (*charm_bnd_cond_fn_t)              (charm_prim_t *par_in, charm_prim_t *par_out, int8_t face, charm_real_t* param, charm_vector_t n);
typedef void            (*charm_flux_fn_t)                  (p4est_t *p4est, charm_prim_t prim[2], charm_real_t* qu, charm_real_t* qv, charm_real_t* qw, charm_real_t* qe, charm_real_t qc[], charm_vector_t n);
typedef void            (*charm_timestep_single_fn_t)       (p4est_t * p4est, charm_real_t *dt, p4est_ghost_t ** _ghost, charm_data_t ** _ghost_data);
typedef charm_real_t    (*charm_get_timestep_fn_t)          (p4est_t * p4est);
typedef void            (*charm_turb_model_fn_t)            (p4est_t * p4est, p4est_ghost_t * ghost, charm_data_t * ghost_data);
typedef void            (*charm_amr_init_fn_t)              (p4est_t *p4est);
typedef void            (*charm_amr_fn_t)                   (p4est_t *p4est, p4est_ghost_t *ghost, charm_data_t *ghost_data);

#ifndef GLOBALS_H_FILE
extern const char *charm_bnd_types[];
extern const char *charm_turb_models[];
#endif

typedef enum {
    TURB_MODEL_SA,
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
    BOUND_FREE_STREAM,
    BOUND_PRESSURE,
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
                        charm_real_t sigma;
                        charm_real_t kappa;
                        charm_real_t cb1;
                        charm_real_t cb2;
                        charm_real_t cw1;
                        charm_real_t cw2;
                        charm_real_t cw3;
                        charm_real_t cv1;
                        charm_real_t ct1;
                        charm_real_t ct2;
                        charm_real_t ct3;
                        charm_real_t ct4;
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
    charm_amr_init_fn_t         amr_init_fn;
    charm_amr_fn_t              amr_fn;
} charm_ctx_t;

typedef struct charm_tree_attr
{
    charm_bnd_t        *bnd[CHARM_FACES];
    charm_reg_t        *reg;
} charm_tree_attr_t;

#ifdef __cplusplus
    extern "C" {
#endif

charm_real_t scalar_prod(charm_vector_t v1, charm_vector_t v2);

charm_real_t vector_length(charm_vector_t v);

void vector_prod(charm_vector_t v1, charm_vector_t v2, charm_vector_t res);
void vector_sub(charm_vector_t v1, charm_vector_t v2, charm_vector_t res);
charm_real_t vector_dist(charm_vector_t v1, charm_vector_t v2);


charm_real_t charm_face_get_area(charm_data_t *d, int8_t face);

charm_real_t charm_face_get_normal(charm_data_t *d, int8_t face, charm_vector_t n);

void charm_quad_get_center(charm_data_t *d, charm_point_t c);

void charm_face_get_center(charm_data_t *d, int8_t face, charm_point_t c);

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

charm_real_t charm_prim_vel_mag(charm_prim_t * prim);

charm_real_t charm_matr3_det(charm_real_t a[3][3]);

void charm_matr3_inv(charm_real_t a[3][3], charm_real_t a_inv[3][3]);

void charm_matr_inv(charm_matr_t a_src, charm_matr_t am);

void charm_matr_vect_mult(charm_matr_t a, charm_vect_t b, charm_vect_t res);

void
charm_matr_add(charm_matr_t a, charm_matr_t b);

void charm_vect_add(charm_vect_t a, charm_vect_t b);

void charm_matr_zero(charm_matr_t a);

void charm_vect_zero(charm_vect_t a);

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


void charm_quad_get_vertices(p4est_t *p4est, p4est_quadrant_t *q, p4est_topidx_t treeid, charm_point_t v[8]);

void charm_geom_quad_calc(p4est_t *p4est, p4est_quadrant_t *q, p4est_topidx_t treeid);


p4est_connectivity_t *charm_conn_create(charm_ctx_t *ctx);

charm_real_t charm_get_heat_k(p4est_t *p4est, charm_real_t *x, charm_data_t *data);

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

