//
// Created by zhrv on 19.10.18.
//

#ifndef CHARM_DG_CHARM_LINSOLVER_H
#define CHARM_DG_CHARM_LINSOLVER_H

#include <petscksp.h>
#include "charm_globals.h"

typedef struct charm_lsolver {
    Vec x;          /* approx solution          */
    Vec b;          /* RHS                      */
    Mat A;          /* linear system matrix     */
    KSP ksp;        /* linear solver context    */
    int n_loc;      /* blocks count             */
    int n_glob;     /* blocks count             */
    int block_n;    /* block size               */
} charm_lsolver_t;

charm_lsolver_t * charm_lsolver_init(int n_glob, int n_loc, int block_n);
void charm_lsolver_ins_block(charm_lsolver_t *solver, int i, int j, double** matr);
void charm_lsolver_add_block(charm_lsolver_t *solver, int i, int j, double** matr);

void charm_lsolver_ins_rhs(charm_lsolver_t *solver, int i, double* vec);
void charm_lsolver_assembly(charm_lsolver_t *solver);
void charm_lsolver_solve(charm_lsolver_t *solver, double eps);

void charm_lsolver_get_sol(charm_lsolver_t *solver, int i, double* vec);


int ___petsc_test(int *argc, char*** argv);
void _petsc_test(p4est_t *p4est);

#endif //CHARM_DG_CHARM_LINSOLVER_H
