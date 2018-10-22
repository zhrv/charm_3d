//
// Created by zhrv on 19.10.18.
//

#include "charm_lsolver.h"
#include "charm_globals.h"

static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
-random_exact_sol : use a random exact solution vector\n\
-view_exact_sol : write exact solution vector to stdout\n\
-m <mesh_x> : number of mesh points in x-direction\n\
-n <mesh_n> : number of mesh points in y-direction\n\n";



charm_lsolver_t * charm_lsolver_init(int n_glob, int n_loc, int block_n)
{
    PetscErrorCode ierr;

    charm_lsolver_t * solver = CHARM_ALLOC(charm_lsolver_t, 1);
    solver->n_loc   = n_loc;
    solver->n_glob  = n_glob;
    solver->block_n = block_n;

    ierr = PetscInitialize(NULL, NULL, (char*)0, help);                                     CHKERRCONTINUE(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &(solver->A));                                       CHKERRCONTINUE(ierr);
    ierr = MatSetSizes(solver->A, PETSC_DECIDE, PETSC_DECIDE, n_glob*block_n, n_glob*block_n);        CHKERRCONTINUE(ierr);
    ierr = MatSetFromOptions           (solver->A);                                         CHKERRCONTINUE(ierr);
    ierr = MatMPIAIJSetPreallocation   (solver->A, 5, NULL, 5, NULL);                       CHKERRCONTINUE(ierr);
    ierr = MatSeqAIJSetPreallocation   (solver->A, 5, NULL);                                CHKERRCONTINUE(ierr);
    ierr = MatSeqSBAIJSetPreallocation (solver->A, 1, 5, NULL);                             CHKERRCONTINUE(ierr);
    ierr = MatMPISBAIJSetPreallocation (solver->A, 1, 5, NULL, 5, NULL);                    CHKERRCONTINUE(ierr);
    ierr = MatMPISELLSetPreallocation  (solver->A, 5, NULL, 5, NULL);                       CHKERRCONTINUE(ierr);
    ierr = MatSeqSELLSetPreallocation  (solver->A, 5, NULL);                                CHKERRCONTINUE(ierr);
    //MatSetOption(solver->A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    ierr = VecCreate(PETSC_COMM_WORLD, &(solver->b));                                       CHKERRCONTINUE(ierr);
    ierr = VecSetSizes(solver->b, n_loc*block_n, PETSC_DETERMINE);                              CHKERRCONTINUE(ierr);
    ierr = VecSetFromOptions(solver->b);                                                    CHKERRCONTINUE(ierr);
    ierr = VecDuplicate(solver->b, &(solver->x));                                           CHKERRCONTINUE(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD,&(solver->ksp));                                      CHKERRCONTINUE(ierr);

    return solver;
}


void charm_lsolver_add_block(charm_lsolver_t *solver, int i, int j, double** matr)
{
    PetscErrorCode ierr;
    const int nb = solver->block_n;
    int ib, jb, i_glob, j_glob;
    for (ib = 0, i_glob = i*nb; ib < nb; ib++, i_glob++) {
        for (jb = 0, j_glob = j*nb; jb < nb; jb++, j_glob++) {
            ierr = MatSetValues(solver->A, 1, &i_glob, 1, &j_glob, &(matr[ib][jb]),ADD_VALUES); CHKERRCONTINUE(ierr);
        }
    }
}


void charm_lsolver_ins_block(charm_lsolver_t *solver, int i, int j, double** matr)
{
    PetscErrorCode ierr;
    const int nb = solver->block_n;
    int ib, jb, i_glob, j_glob;
    double val;
    for (ib = 0, i_glob = i*nb; ib < nb; ib++, i_glob++) {
        for (jb = 0, j_glob = j*nb; jb < nb; jb++, j_glob++) {
            val = matr[ib][jb];
            ierr = MatSetValues(solver->A, 1, &i_glob, 1, &j_glob, &val,INSERT_VALUES); CHKERRCONTINUE(ierr);
        }
    }
}


void charm_lsolver_ins_rhs(charm_lsolver_t *solver, int i, double* vec)
{
    PetscErrorCode ierr;
    const int nb = solver->block_n;
    int ib, i_glob;
    for (ib = 0, i_glob = i*nb; ib < nb; ib++, i_glob++) {
        ierr = VecSetValues(solver->b, 1, &i_glob, &(vec[ib]), INSERT_VALUES); CHKERRCONTINUE(ierr);
    }
}


void charm_lsolver_add_rhs(charm_lsolver_t *solver, int i, double* vec)
{
    PetscErrorCode ierr;
    const int nb = solver->block_n;
    int ib, i_glob;
    for (ib = 0, i_glob = i*nb; ib < nb; ib++, i_glob++) {
        ierr = VecSetValues(solver->b, 1, &i_glob, &(vec[ib]), INSERT_VALUES); CHKERRCONTINUE(ierr);
    }
}


void charm_lsolver_solve(charm_lsolver_t *solver, double eps)
{
    PetscErrorCode ierr;
    ierr = MatAssemblyBegin(solver->A,MAT_FINAL_ASSEMBLY);                              CHKERRCONTINUE(ierr);
    ierr = MatAssemblyEnd(solver->A,MAT_FINAL_ASSEMBLY);                                CHKERRCONTINUE(ierr);
    ierr = VecAssemblyBegin(solver->b);                              CHKERRCONTINUE(ierr);
    ierr = VecAssemblyEnd(solver->b);                                CHKERRCONTINUE(ierr);
    ierr = KSPSetOperators(solver->ksp, solver->A, solver->A);                          CHKERRCONTINUE(ierr);
    ierr = KSPSetTolerances(solver->ksp, eps, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);    CHKERRCONTINUE(ierr);
    ierr = KSPSolve(solver->ksp, solver->b, solver->x);                                 CHKERRCONTINUE(ierr);
}


void charm_lsolver_get_sol(charm_lsolver_t *solver, int i, double* vec)
{
    PetscErrorCode ierr;
    const int nb = solver->block_n;
    int ib, i_glob;
    for (ib = 0, i_glob = i*nb; ib < nb; ib++, i_glob++) {
        ierr = VecGetValues(solver->x, 1, &i_glob, &(vec[ib])); CHKERRCONTINUE(ierr);
    }
}


void charm_lsolver_destroy(charm_lsolver_t *solver)
{
    PetscErrorCode ierr;
    ierr = KSPDestroy(&(solver->ksp));  CHKERRCONTINUE(ierr);
    ierr = VecDestroy(&(solver->x));    CHKERRCONTINUE(ierr);
    ierr = VecDestroy(&(solver->b));    CHKERRCONTINUE(ierr);
    ierr = MatDestroy(&(solver->A));    CHKERRCONTINUE(ierr);
    CHARM_FREE(solver);
}


void _petsc_test(p4est_t *p4est)
{
    int rank = p4est->mpirank;
    int size = p4est->mpisize;
    double **matr1, **matr2;
    int i, j, p;
//    double vec[12] = {15,27,21,36,60,44,60,96,68,57,90,63};
    double vec[3] = {4,6,4};

    PetscInt start, end;

    matr1 = CHARM_ALLOC(double*, 3);
    matr2 = CHARM_ALLOC(double*, 3);
    for (i = 0; i < 3; i++) {
        matr1[i] = CHARM_ALLOC(double, 3);
        matr2[i] = CHARM_ALLOC(double, 3);
        for (j = 0; j < 3; j++) {
            matr1[i][j] = (abs(i-j) < 2) ? 1. : 0.;
            matr2[i][j] = (abs(i-j) < 2) ? 2. : 0.;
        }
    }

    charm_lsolver_t *s = charm_lsolver_init(size, 1, 3);
    charm_lsolver_ins_block(s, rank, rank, matr2);
    charm_lsolver_ins_rhs(s, rank, &(vec[rank*3]));

    if (size > 1) {
        switch (rank) {
            case 0:
                charm_lsolver_ins_block(s, rank, rank + 1, matr1);
                break;
            case 1:
            case 2:
                charm_lsolver_ins_block(s, rank, rank - 1, matr1);
                charm_lsolver_ins_block(s, rank, rank + 1, matr1);
                break;
            case 3:
                charm_lsolver_ins_block(s, rank, rank - 1, matr1);
                break;
        }
    }


    charm_lsolver_solve(s, 1.e-5);



    int r = 0;
    for (p = 0; p < size; p++) {
        if (p == rank) {
            VecGetOwnershipRange(s->b, &start, &end);
            printf("*****  %d [%d, %d]  *****\n", rank, start, end);

                memset(vec, 0, sizeof(double)*3);
                charm_lsolver_get_sol(s, rank, (double*)vec);
                for (j = 0; j < 3; j++) {
                    printf("%d:  %f\n", r++, vec[j]);
                }

        }
        else {
            r += 1*3;
        }
        sc_MPI_Barrier(MPI_COMM_WORLD);
    }
//    VecView(s->b,PETSC_VIEWER_STDOUT_WORLD);

    charm_lsolver_destroy(s);
}

int ___petsc_test(int *argc, char*** argv)
{
    Vec x,b,u; /* approx solution, RHS, exact solution */
    Mat A; /* linear system matrix */
    KSP ksp; /* linear solver context */
    PetscRandom rctx; /* random number generator context */
    PetscReal norm; /* norm of solution error */
    PetscInt i,j,Ii,J,Istart,Iend,m = 8,n = 7,its;
    PetscErrorCode ierr;
    PetscBool flg = PETSC_FALSE;
    PetscScalar v;
#if defined(PETSC_USE_LOG)
    PetscLogStage stage;
#endif
    ierr = PetscInitialize(argc,argv,(char*)0,help);if (ierr) return ierr;
    ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Compute the matrix and right-hand-side vector that define
the linear system, Ax = b.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*
Create parallel matrix, specifying only its global dimensions.
When using MatCreate(), the matrix format can be specified at
runtime. Also, the parallel partitioning of the matrix is
determined by PETSc at runtime.
Performance tuning note: For problems of substantial size,preallocation of matrix memory is crucial for attaining good
performance. See the matrix chapter of the users manual for details.
*/
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,5,NULL);CHKERRQ(ierr);
    ierr = MatSeqSBAIJSetPreallocation(A,1,5,NULL);CHKERRQ(ierr);
    ierr = MatMPISBAIJSetPreallocation(A,1,5,NULL,5,NULL);CHKERRQ(ierr);
    ierr = MatMPISELLSetPreallocation(A,5,NULL,5,NULL);CHKERRQ(ierr);
    ierr = MatSeqSELLSetPreallocation(A,5,NULL);CHKERRQ(ierr);
/*
Currently, all PETSc parallel matrix formats are partitioned by
contiguous chunks of rows across the processors. Determine which
rows of the matrix are locally owned.
*/
    ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
/*
Set matrix elements for the 2-D, five-point stencil in parallel.
- Each processor needs to insert only elements that it owns
locally (but any non-local elements will be sent to the
appropriate processor during matrix assembly).
- Always specify global rows and columns of matrix entries.
Note: this uses the less common natural ordering that orders first
all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
instead of J = I +- m as you might expect. The more standard ordering
would first do all variables for y = h, then y = 2h etc.
*/
    ierr = PetscLogStageRegister("Assembly", &stage);CHKERRQ(ierr);
    ierr = PetscLogStagePush(stage);CHKERRQ(ierr);
    for (Ii=Istart; Ii<Iend; Ii++) {
        v = -1.0; i = Ii/n; j = Ii - i*n;
        if (i>0) {J = Ii - n; ierr =
                                      MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (i<m-1) {J = Ii + n; ierr =
                                        MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (j>0) {J = Ii - 1; ierr =
                                      MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        if (j<n-1) {J = Ii + 1; ierr =
                                        MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
        v = 4.0; ierr = MatSetValues(A,1,&Ii,1,&Ii,&v,ADD_VALUES);CHKERRQ(ierr);
    }
/*Assemble matrix, using the 2-step process:
MatAssemblyBegin(), MatAssemblyEnd()
Computations can be done while messages are in transition
by placing code between these two statements.
*/
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);
/* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
    ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
/*
Create parallel vectors.
- We form 1 vector from scratch and then duplicate as needed.
- When using VecCreate(), VecSetSizes and VecSetFromOptions()
in this example, we specify only the
vector’s global dimension; the parallel partitioning is determined
at runtime.
- When solving a linear system, the vectors and matrices MUST
be partitioned accordingly. PETSc automatically generates
appropriately partitioned matrices and vectors when MatCreate()
and VecCreate() are used with the same communicator.
- The user can alternatively specify the local vector and matrix
dimensions when more sophisticated partitioning is needed
(replacing the PETSC DECIDE argument in the VecSetSizes() statement
below).
*/
    ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
    ierr = VecSetSizes(u,PETSC_DECIDE,m*n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(u);CHKERRQ(ierr);
    ierr = VecDuplicate(u,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
/*
Set exact solution; then compute right-hand-side vector.
By default we use an exact solution of a vector with all
elements of 1.0; Alternatively, using the runtime option
-random_sol forms a solution vector with random components.
*/
    ierr =
            PetscOptionsGetBool(NULL,NULL,"-random_exact_sol",&flg,NULL);CHKERRQ(ierr);
    if (flg) {
        ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
        ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
        ierr = VecSetRandom(u,rctx);CHKERRQ(ierr);
        ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);
    } else {
        ierr = VecSet(u,1.0);CHKERRQ(ierr);}
    ierr = MatMult(A,u,b);CHKERRQ(ierr);
/*
View the exact solution vector if desired
*/
    flg = PETSC_FALSE;
    ierr =
            PetscOptionsGetBool(NULL,NULL,"-view_exact_sol",&flg,NULL);CHKERRQ(ierr);
    if (flg) {ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Create the linear solver and set various options
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*
Create linear solver context
*/
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
/*
Set operators. Here the matrix that defines the linear system
also serves as the preconditioning matrix.
*/
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
/*
Set linear solver defaults for this problem (optional).
- By extracting the KSP and PC contexts from the KSP context,
we can then directly call any KSP and PC routines to set
various options.
- The following two statements are optional; all of these
parameters could alternatively be specified at runtime via
KSPSetFromOptions(). All of these defaults can be
overridden at runtime, as indicated below.
*/
    ierr = KSPSetTolerances(ksp,1.e-2/((m+1)*(n+1)),1.e-50,PETSC_DEFAULT,
            PETSC_DEFAULT);CHKERRQ(ierr);
/*
Set runtime options, e.g.,
-ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
These options will override those specified above as long as
KSPSetFromOptions() is called _after_ any other customization
routines.
*/
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
/* - - - - - - - - - - - - - - - - - - - - - - - - - -Solve the linear system
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Check solution and clean up
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*
Check the error
*/
    ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
    ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
/*
Print convergence information. PetscPrintf() produces a single
print statement from all processes that share a communicator.
An alternative is PetscFPrintf(), which prints to a file.
*/
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations %D\n",(double)norm,its);CHKERRQ(ierr);
/*
Free work space. All PETSc objects should be destroyed when they
are no longer needed.
*/
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&u);CHKERRQ(ierr); ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
/*
Always call PetscFinalize() before exiting a program. This routine
- finalizes the PETSc libraries as well as MPI
- provides summary and diagnostic information if certain runtime
options are chosen (e.g., -log_view).
*/
    ierr = PetscFinalize();
    return ierr;
}



