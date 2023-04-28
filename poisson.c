#include "poisson.h"
// #include <mpi.h>

/*Called by poisson_solver at each time step*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*    -Impose zero mass flow here by changing value of U_star*/
/*    -Fill vector rhs*/
void computeRHS(double *rhs, PetscInt rowStart, PetscInt rowEnd, int m, int n, double inv_delta_tx, matrix *U, matrix *V)
{

    //YOU MUST IMPOSE A ZERO-MASS FLOW HERE ...
    for (int i = 0; i < n; i++){
        for (int j = 0; j< m; j++ ){
            // inv_delta_tx=dx/dt= (dx^2)/(dx*dt)
            rhs[i*m+j] = inv_delta_tx*(U->a[i+1][j+1]-U->a[i][j+1] + V->a[i+1][j+1]-V->a[i+1][j]); 
            //notes pour trop tards: size of U = (m, n+1) and V = (m+1, n)

            /*WRITE HERE (nabla dot u_star)/dt at each mesh point r*/
            /*Do not forget that the solution for the Poisson equation is defined within a constant.
            One point from Phi must then be set to an abritrary constant.*/
        }
    }

}

/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*    - Change the call to computeRHS as you have to modify its prototype too*/
/*    - Copy solution of the equation into your vector PHI*/
void poisson_solver(Poisson_data *data, double inv_delta_tx,int m, int n, matrix *U, matrix *V, matrix *phi)
{

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    int its;
    PetscInt rowStart, rowEnd;
    PetscScalar *rhs, *sol;

    KSP sles = data->sles;
    Vec b = data->b;
    Vec x = data->x;

    /* Fill the right-hand-side vector : b */
    VecGetOwnershipRange(b, &rowStart, &rowEnd);
    VecGetArray(b, &rhs);
    computeRHS(rhs, rowStart, rowEnd, m, n, inv_delta_tx, U, V); /*MODIFY THE PROTOTYPE HERE*/
    VecRestoreArray(b, &rhs);


    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

    VecGetArray(x, &sol);


    for(int r = rowStart; r<rowEnd; r++){
        /*YOUR VECTOR PHI[...]*/ // = sol[r];
        phi->data[r] = sol[r];
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
/*In its current state, it inserts unity on the main diagonal.*/
/*More than probably, you should need to add arguments to the prototype ... .*/
/*Modification to do in this function : */
/*   -Insert the correct factor in matrix A*/
void computeLaplacianMatrix(Mat A, int rowStart, int rowEnd, int m, int n)
{
    int r;
    // if bord 
    //coin inférieur gauche 
    MatSetValue(A, 0, 0 , -1.0, INSERT_VALUES);
    MatSetValue(A, 0, 1 , 0.5, INSERT_VALUES);
    MatSetValue(A, 0, m , 0.5, INSERT_VALUES);
    //coin superieur gauche 
    MatSetValue(A, m-1, m-1, -1.0, INSERT_VALUES);
    MatSetValue(A, m-1, m-2 , 0.5, INSERT_VALUES);
    MatSetValue(A, m-1, 2*m-1 , 0.5, INSERT_VALUES);
    //coin inférieur droit 
    MatSetValue(A, m*(n-1), m*(n-1), -1.0, INSERT_VALUES);
    MatSetValue(A, m*(n-1), m*(n-1)+1 , 0.5, INSERT_VALUES);
    MatSetValue(A, m*(n-1), m*(n-2) , 0.5, INSERT_VALUES);
    //coin supérieur droit 
    MatSetValue(A, m*n-1, m*n-1, -1.0, INSERT_VALUES);
    MatSetValue(A, m*n-1, m*n-2 , 0.5, INSERT_VALUES);
    MatSetValue(A, m*n-1, m*(n-1)-1 , 0.5, INSERT_VALUES);
    //bord gauche 
    for(r = 1; r < m-1; r++){
        MatSetValue(A, r, r , -1.0, INSERT_VALUES);
        MatSetValue(A, r, r-1 , 1.0/3.0, INSERT_VALUES);
        MatSetValue(A, r, r+1 , 1.0/3.0, INSERT_VALUES);
        MatSetValue(A, r, r+m , 1.0/3.0, INSERT_VALUES);
    }
    //bord droit 
    for(r = m*(n-1)+1; r < n*m-1; r++){
        MatSetValue(A, r, r , -1.0, INSERT_VALUES);
        MatSetValue(A, r, r-1 , 1.0/3.0, INSERT_VALUES);
        MatSetValue(A, r, r+1 , 1.0/3.0, INSERT_VALUES);
        MatSetValue(A, r, r-m , 1.0/3.0, INSERT_VALUES);
    }
    //bord inférieur 
    for(r = m; r < m*(n-1); r+=m){
        MatSetValue(A, r, r , -1.0, INSERT_VALUES);
        MatSetValue(A, r, r+1 , 1.0/3.0, INSERT_VALUES);
        MatSetValue(A, r, r-m , 1.0/3.0, INSERT_VALUES);
        MatSetValue(A, r, r+m , 1.0/3.0, INSERT_VALUES);
    }
    //bord supérieur
    for(r = 2*m-1; r < m*n-1; r+=m){
        MatSetValue(A, r, r , -1.0, INSERT_VALUES);
        MatSetValue(A, r, r-1 , 1.0/3.0, INSERT_VALUES);
        MatSetValue(A, r, r-n , 1.0/3.0, INSERT_VALUES);
        MatSetValue(A, r, r+n , 1.0/3.0, INSERT_VALUES);
    }
    for (int i = 1; i < n-1; i++){
        for (int j = 1; j< m-1; j++ ){
        r = i*m+j;
        MatSetValue(A, r, r , -1.0, INSERT_VALUES);
        MatSetValue(A, r, r-1 , 0.25, INSERT_VALUES);
        MatSetValue(A, r, r+1 , 0.25, INSERT_VALUES);
        MatSetValue(A, r, r+m , 0.25, INSERT_VALUES);
        MatSetValue(A, r, r-m , 0.25, INSERT_VALUES);
        /*USING MATSETVALUE FUNCTION, INSERT THE GOOD FACTOR AT THE GOOD PLACE*/
        /*Be careful; the solution from the system solved is defined within a constant.
        One point from Phi must then be set to an abritrary constant.*/
        }
    }
}

/*To call during the initialization of your solver, before the begin of the time loop*/
/*Maybe you should need to add an argument to specify the number of unknows*/
/*Modification to do in this function :*/
/*   -Specify the number of unknows*/
/*   -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver(Poisson_data* data, int m, int n)
{
    PetscInt rowStart = 0; /*rowStart = 0*/
    PetscInt rowEnd = m*n; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;
    int nphi = m*n; /*WRITE HERE THE NUMBER OF UNKNOWS*/

    /* Create the right-hand-side vector : b */
    VecCreate(PETSC_COMM_WORLD, &(data->b));
    VecSetSizes(data->b, PETSC_DECIDE, nphi);
    VecSetType(data->b,VECSTANDARD);

    /* Create the solution vector : x */
    VecCreate(PETSC_COMM_WORLD, &(data->x));
    VecSetSizes(data->x, PETSC_DECIDE, nphi);
    VecSetType(data->x,VECSTANDARD);

    /* Create and assemble the Laplacian matrix : A  */
    MatCreate(PETSC_COMM_WORLD, &(data->A));
    MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi , nphi);
    MatSetType(data->A, MATAIJ);
    MatSeqAIJSetPreallocation(data->A,5, NULL); // /*SET HERE THE NUMBER OF NON-ZERO DIAGONALS*/
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    computeLaplacianMatrix(data->A, rowStart, rowEnd, m, n);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    KSPSetType(data->sles,KSPGMRES); //KSPGMRES seems the best, it will not be used if PC LU.
    PC prec;
    KSPGetPC(data->sles, &prec);
    PCSetType(prec,PCLU);
    //KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles,1,4);
    KSPGMRESSetPreAllocateVectors(data->sles);

    PetscPrintf(PETSC_COMM_WORLD, "Assembly of Mattrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
/*Modification to do : nothing */
void free_poisson_solver(Poisson_data* data){
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
