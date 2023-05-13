// #include "poisson.h"
// #include <mpi.h>

// /*Called by poisson_solver at each time step*/
// /*More than probably, you should need to add arguments to the prototype ... */
// /*Modification to do :*/
// /*    -Impose zero mass flow here by changing value of U_star*/
// /*    -Fill vector rhs*/
// void computeRHS(double *rhs, PetscInt rowStart, PetscInt rowEnd, int m, int n, double dt, matrix *Uestim, matrix *Vestim)
// {
//     for (int i = 0; i < n; i++){
//         for (int j = 0; j< m; j++ ){
//             // inv_delta_tx=dx/dt= (dx^2)/(dx*dt)
//             rhs[i*m+j] = (Uestim->a[i+1][j+1]-Uestim->a[i][j+1] + Vestim->a[i+1][j+1]-Vestim->a[i+1][j])* 1.0/dt; 

//             /*WRITE HERE (nabla dot u_star)/dt at each mesh point r*/
//             /*Do not forget that the solution for the Poisson equation is defined within a constant.
//             One point from Phi must then be set to an  .*/
//         }
//     }
//     // for (int i = 0; i <(n-1)*m+(m-1);i++){
//     //     printf("%.6e ", rhs[i]);
//     //     printf("\n");
//     // }
    
//     rhs[m*n/2 + m/2] = 1e-10;
//     //rhs[0] = 100.0;
// }
// void computeRHS2(double *rhs, PetscInt rowStart, PetscInt rowEnd, int m, int n, double dt, matrix *Uestim, matrix *Vestim)
// {
//     int i,j;
//     for (j = 0; j< m; j++ ){
//         for (i = 0; i < n; i++){
//         // inv_delta_tx=dx/dt= (dx^2)/(dx*dt)
//         rhs[i*m+j] = Uestim->a[i][j]; 

//         /*WRITE HERE (nabla dot u_star)/dt at each mesh point r*/
//         /*Do not forget that the solution for the Poisson equation is defined within a constant.
//         One point from Phi must then be set to an  .*/
//         }
//     }
//     // rhs[m*n/2 + m/2] = 0.0;//69.0;
//     rhs[0] = 3 ;   //13.0/4.0;
// }

// void poisson_solver(Poisson_data *data, double deltax, double dt, int m, int n, matrix *Uestim, matrix *Vestim, matrix *phi)
// {

//     /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
//     int its;
//     PetscInt rowStart, rowEnd;
//     PetscScalar *rhs, *sol;

//     KSP sles = data->sles;
//     Vec b = data->b;
//     Vec x = data->x;

//     /* Fill the right-hand-side vector : b */
//     VecGetOwnershipRange(b, &rowStart, &rowEnd);
//     VecGetArray(b, &rhs);
//     computeRHS(rhs, rowStart, rowEnd, m, n, dt, Uestim, Vestim); /*MODIFY THE PROTOTYPE HERE*/
//     VecRestoreArray(b, &rhs);


//     /*Solve the linear system of equations */
//     KSPSolve(sles, b, x);
//     KSPGetIterationNumber(sles, &its);
//     PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

//     VecGetArray(x, &sol);


//     for(int r = rowStart; r<rowEnd; r++){
//         /*YOUR VECTOR PHI[...]*/ // = sol[r];
//         // int i = (int) r/m;
//         // int j = r%m;
//         // phi->a[i][j] = sol[r];
//         phi->data[r] = sol[r];
//     }

//     VecRestoreArray(x, &sol);
//     free(sol);
//     free(rhs);
// }

// void poisson_solver2(Poisson_data *data, double deltax, double dt, int m, int n, matrix *Uestim, matrix *Vestim, matrix *phi)
// {

//     /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
//     int its;
//     PetscInt rowStart, rowEnd;
//     PetscScalar *rhs, *sol;

//     KSP sles = data->sles;
//     Vec b = data->b;
//     Vec x = data->x;

//     /* Fill the right-hand-side vector : b */
//     VecGetOwnershipRange(b, &rowStart, &rowEnd);
//     VecGetArray(b, &rhs);
//     // computeRHS(rhs, rowStart, rowEnd, m, n, dt, Uestim, Vestim); /*MODIFY THE PROTOTYPE HERE*/
//     computeRHS2(rhs, rowStart, rowEnd, m, n, dt, Uestim, Vestim); /*MODIFY THE PROTOTYPE HERE*/
//     VecRestoreArray(b, &rhs);



//     /*Solve the linear system of equations */
//     KSPSolve(sles, b, x);
//     KSPGetIterationNumber(sles, &its);
//     PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

//     VecGetArray(x, &sol);


//     for(int r = rowStart; r<rowEnd; r++){
//         /*YOUR VECTOR PHI[...]*/ // = sol[r];
//         phi->data[r] = sol[r];
//     }


//     VecRestoreArray(x, &sol);

    
//     free(sol);
//     free(rhs);
// }

// /*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
// /*In its current state, it inserts unity on the main diagonal.*/
// /*More than probably, you should need to add arguments to the prototype ... .*/
// /*Modification to do in this function : */
// /*   -Insert the correct factor in matrix A*/
// void computeLaplacianMatrix(Mat A, int rowStart, int rowEnd, int m, int n, double deltax)
// {   
//     int r;
//     //coin inférieur gauche 
//     MatSetValue(A, 0, 0 , -2.0/deltax, INSERT_VALUES);
//     MatSetValue(A, 0, 1 , 1.0/deltax, INSERT_VALUES);
//     MatSetValue(A, 0, m , 1.0/deltax, INSERT_VALUES);
//     //coin superieur gauche 
//     MatSetValue(A, m-1, m-1, -2.0/deltax, INSERT_VALUES);
//     MatSetValue(A, m-1, m-2 , 1.0/deltax, INSERT_VALUES);
//     MatSetValue(A, m-1, 2*m-1 , 1.0/deltax, INSERT_VALUES);
//     //coin inférieur droit 
//     MatSetValue(A, m*(n-1), m*(n-1), -2.0/deltax, INSERT_VALUES);
//     MatSetValue(A, m*(n-1), m*(n-1)+1 , 1.0/deltax, INSERT_VALUES);
//     MatSetValue(A, m*(n-1), m*(n-2) , 1.0/deltax, INSERT_VALUES);
//     //coin supérieur droit 
//     MatSetValue(A, m*n-1, m*n-1, -2.0/deltax, INSERT_VALUES);
//     MatSetValue(A, m*n-1, m*n-2 , 1.0/deltax, INSERT_VALUES);
//     MatSetValue(A, m*n-1, m*(n-1)-1 , 1.0/deltax, INSERT_VALUES);
//     //bord gauche 
//     for(r = 1; r < m-1; r++){
//         MatSetValue(A, r, r , -3.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r-1 , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r+1 , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r+m , 1.0/deltax, INSERT_VALUES);
//     }
//     //bord droit 
//     for(r = m*(n-1)+1; r < n*m-1; r++){
//         MatSetValue(A, r, r , -3.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r-1 , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r+1 , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r-m , 1.0/deltax, INSERT_VALUES);
//     }
//     //bord inférieur 
//     for(r = m; r < m*(n-1); r+=m){
//         MatSetValue(A, r, r , -3.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r+1 , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r-m , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r+m , 1.0/deltax, INSERT_VALUES);
//     }
//     //bord supérieur
//     for(r = 2*m-1; r < m*n-1; r+=m){
//         MatSetValue(A, r, r , -3.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r-1 , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r-m , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r+m , 1.0/deltax, INSERT_VALUES);
//     }
//     for (int i = 1; i < n-1; i++){
//         for (int j = 1; j< m-1; j++ ){
//         r = i*m+j;
//         MatSetValue(A, r, r , -4.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r-1 , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r+1 , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r+m , 1.0/deltax, INSERT_VALUES);
//         MatSetValue(A, r, r-m , 1.0/deltax, INSERT_VALUES);
//         /*USING MATSETVALUE FUNCTION, INSERT THE GOOD FACTOR AT THE GOOD PLACE*/
//         /*Be careful; the solution from the system solved is defined within a constant.
//         One point from Phi must then be set to an arbitrary constant.*/
//         }
//     }
//     // MatSetValue(A, m+1, m+1 , 1.0, INSERT_VALUES);
//     // MatSetValue(A, m+1, 1 , 0.0, INSERT_VALUES);
//     // MatSetValue(A, m+1, 2*m+1 , 0.0, INSERT_VALUES);
//     // MatSetValue(A, m+1, m+2 , 0.0, INSERT_VALUES);
//     // MatSetValue(A, m+1, m , 0.0, INSERT_VALUES);

//     // MatSetValue(A, m/2 , m/2, 1.0, INSERT_VALUES);
//     // MatSetValue(A, m/2 , m/2 -1, 0.0, INSERT_VALUES);
//     // MatSetValue(A, m/2 , m/2 +m , 0.0, INSERT_VALUES);
//     // MatSetValue(A, m/2 , m/2 +1 , 0.0, INSERT_VALUES); // ca marche ici 

//     // MatSetValue(A, 0 , 0, 1.0, INSERT_VALUES);
//     // MatSetValue(A, 0 , 0 + 1, 0.0, INSERT_VALUES);
//     // MatSetValue(A, 0 , 0 + m , 0.0, INSERT_VALUES); 

//     MatSetValue(A, m*n/2 + m/2 , m*n/2 + m/2, 1.0, INSERT_VALUES);
//     MatSetValue(A, m*n/2 + m/2 , m*n/2 + m/2 -1, 0.0, INSERT_VALUES);
//     MatSetValue(A, m*n/2 + m/2 , m*n/2 + m/2 +m , 0.0, INSERT_VALUES);
//     MatSetValue(A, m*n/2 + m/2 , m*n/2 + m/2 +1 , 0.0, INSERT_VALUES);
//     MatSetValue(A, m*n/2 + m/2 , m*n/2 + m/2 -m , 0.0, INSERT_VALUES);
// }

// /*To call during the initialization of your solver, before the begin of the time loop*/
// /*Maybe you should need to add an argument to specify the number of unknows*/
// /*Modification to do in this function :*/
// /*   -Specify the number of unknows*/
// /*   -Specify the number of non-zero diagonals in the sparse matrix*/
// PetscErrorCode initialize_poisson_solver(Poisson_data* data, int m, int n, double deltax)
// {
//     // printf("we are initializing poisson \n");
//     PetscInt rowStart; /*rowStart = 0*/
//     PetscInt rowEnd; /*rowEnd = the number of unknows*/
//     PetscErrorCode ierr;
//     int nphi = m*n; /*WRITE HERE THE NUMBER OF UNKNOWS*/

//     /* Create the right-hand-side vector : b */
//     VecCreate(PETSC_COMM_WORLD, &(data->b));
//     VecSetSizes(data->b, PETSC_DECIDE, nphi);
//     VecSetType(data->b,VECSTANDARD);

//     /* Create the solution vector : x */
//     VecCreate(PETSC_COMM_WORLD, &(data->x));
//     VecSetSizes(data->x, PETSC_DECIDE, nphi);
//     VecSetType(data->x,VECSTANDARD);

//     /* Create and assemble the Laplacian matrix : A  */
//     MatCreate(PETSC_COMM_WORLD, &(data->A));
//     MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi , nphi);
//     MatSetType(data->A, MATAIJ);
//     MatSeqAIJSetPreallocation(data->A,5, NULL); // /*SET HERE THE NUMBER OF NON-ZERO DIAGONALS*/
//     MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

//     computeLaplacianMatrix(data->A, rowStart, rowEnd, m, n, deltax);
//     ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
//     CHKERRQ(ierr);
//     ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
//     CHKERRQ(ierr);

//     /* Create the Krylov context */
//     KSPCreate(PETSC_COMM_WORLD, &(data->sles));
//     KSPSetOperators(data->sles, data->A, data->A);
//     KSPSetType(data->sles,KSPGMRES); //KSPGMRES seems the best, it will not be used if PC LU.
//     PC prec;
//     KSPGetPC(data->sles, &prec);
//     PCSetType(prec,PCLU);
//     KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
//     KSPSetTolerances(data->sles, 1.e-12, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT);
//     KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
//     KSPSetUseFischerGuess(data->sles,1,4);
//     KSPGMRESSetPreAllocateVectors(data->sles);
    
    
//     return ierr;
// }

// // void print_petsc_matrix(Mat A){
// //     MatView(A, PETSCVIEWERASCII, PETSC_VIEWER_ASCII_INFO);
// //     //PETSC_VIEWER_ASCII_INFO 
// // }

// /*To call after the simulation to free the vectors needed for Poisson equation*/
// /*Modification to do : nothing */
// void free_poisson_solver(Poisson_data* data){
//     MatDestroy(&(data->A));
//     VecDestroy(&(data->b));
//     VecDestroy(&(data->x));
//     KSPDestroy(&(data->sles));
// }


// void test_poisson(Poisson_data *data, double deltax, double dt, int m, int n, matrix *U, matrix *V, matrix *phi) {
//     int i, j;
//     poisson_solver(data, deltax, dt, m, n, U, V, phi);
//     FILE *ptr = fopen("./results/test_poisson.txt", "w");
//     fprintf(ptr, "%d\n", n);
//     for (j = 0; j < m; j++) {
//         for (i = 0; i < n; i++) {
//             fprintf(ptr, "%.4le ", phi->data[i*m+j]);
//         }
//         fprintf(ptr, "\n");
//     } 
//     fclose(ptr);
// }


#include <mpi.h>
#include "poisson.h"

int PHI_IND(int i, int j, int Ny){
    return (Ny-1) * i + j;
}

/*Called by poisson_solver at each time step*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*  OKKKK  -Impose zero mass flow here by changing value of U_star*/
/*  OKKKK  -Fill vector rhs*/
void computeRHS_otis(double *rhs, PetscInt rowStart, PetscInt rowEnd, int nx, int ny, double **u_star, double **v_star, double dt)
{
    int i, j;
    // for(i = 1; i < nx - 2; i++){
    //     v_star[i][0] = 0.0;
    //     v_star[i][ny - 1] = 0.0;
    // }
    
    // for(i = 1; i < ny - 2; i++){
    //     u_star[0][i] = 0.0;
    //     u_star[nx - 1][i] = 0.0;
    // }

    int r = rowStart;
    for(r=rowStart; r<rowEnd ; r++){
        i = (int) r/(ny - 1);
        j = r%(ny - 1);
		rhs[r] =  1.0/dt * (u_star[i+1][j+1] - u_star[i][j+1] + v_star[i+1][j+1] - v_star[i+1][j]); /*WRITE HERE (nabla dot u_star)/dt at each mesh point r*/
        /*Do not forget that the solution for the Poisson equation is defined within a constant.
        One point from Phi must then be set to an abritrary constant.*/
    }
    rhs[(ny - 1)/2] = 0.0;
    // rhs[0] = 0.0;
    
}

/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
/*More than probably, you should need to add arguments to the prototype ... */
/*Modification to do :*/
/*   OKKKKK - Change the call to computeRHS as you have to modify its prototype too*/
/*   OKKKKK - Copy solution of the equation into your vector PHI*/
void poisson_solver_otis(Poisson_data *data, int nx, int ny, double **u_star, double **v_star, double dt, double *phi)
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
    computeRHS_otis(rhs, rowStart, rowEnd, nx, ny, u_star, v_star, dt); /*MODIFY THE PROTOTYPE HERE*/
    VecRestoreArray(b, &rhs);

    // VecView(b, PETSC_VIEWER_STDOUT_WORLD);


    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

    VecGetArray(x, &sol);

    int r;
    for(r=rowStart; r<rowEnd; r++){
        phi[r] = sol[r];
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
/*In its current state, it inserts unity on the main diagonal.*/
/*More than probably, you should need to add arguments to the prototype ... .*/
/*Modification to do in this function : */
/*  OKKKK -Insert the correct factor in matrix A*/
/*

Compute the discretized laplacian of phi
we need to set a value of phi, ex, phi[0] = 0

*/
void computeLaplacianMatrix_otis(Mat A, int rowStart, int rowEnd, int Nx, int Ny, double h)
{
    int i,j;
    int r;
    for(j = 1; j < Ny - 2; j++){
        r = rowStart + PHI_IND(0, j, Ny);
        MatSetValue(A, r, r , -3.0/h, INSERT_VALUES);
        MatSetValue(A, r, r+1 , 1.0/h, INSERT_VALUES);
        MatSetValue(A, r, r-1 , 1.0/h, INSERT_VALUES);
        MatSetValue(A, r, r+(Ny-1) , 1.0/h, INSERT_VALUES);
        r = rowStart + PHI_IND(Nx-2, j, Ny);
        MatSetValue(A, r, r , -3.0/h, INSERT_VALUES);
        MatSetValue(A, r, r+1 , 1.0/h, INSERT_VALUES);
        MatSetValue(A, r, r-1 , 1.0/h, INSERT_VALUES);
        MatSetValue(A, r, r-(Ny-1), 1.0/h, INSERT_VALUES);
    }
    for(i = 1; i < Nx - 2; i++){
        r = rowStart + PHI_IND(i, 0, Ny);
        MatSetValue(A, r, r , -3.0/h, INSERT_VALUES);
        MatSetValue(A, r, r+1 , 1.0/h, INSERT_VALUES);
        MatSetValue(A, r, r+(Ny-1) , 1.0/h, INSERT_VALUES);
        MatSetValue(A, r, r-(Ny-1) , 1.0/h, INSERT_VALUES);
        r = rowStart + PHI_IND(i, Ny - 2, Ny);
        MatSetValue(A, r, r , -3.0/h, INSERT_VALUES);
        MatSetValue(A, r, r-1 , 1.0/h, INSERT_VALUES);
        MatSetValue(A, r, r+(Ny-1) , 1.0/h, INSERT_VALUES);
        MatSetValue(A, r, r-(Ny-1) , 1.0/h, INSERT_VALUES);
        for(j = 1; j < Ny - 2; j++){
            r = rowStart + PHI_IND(i, j, Ny);
            MatSetValue(A, r, r , -4.0/h, INSERT_VALUES);
            MatSetValue(A, r, r+1 , 1.0/h, INSERT_VALUES);
            MatSetValue(A, r, r-1 , 1.0/h, INSERT_VALUES);
            MatSetValue(A, r, r+(Ny-1) , 1.0/h, INSERT_VALUES);
            MatSetValue(A, r, r-(Ny-1) , 1.0/h, INSERT_VALUES);
        }
    }

    r = rowStart + PHI_IND(0, 0, Ny);
    MatSetValue(A, r, r , -2.0/h, INSERT_VALUES);
    MatSetValue(A, r, r+1 , 1.0/h, INSERT_VALUES);
    MatSetValue(A, r, r+(Ny-1) , 1.0/h, INSERT_VALUES);

    r = rowStart + PHI_IND(0, Ny - 2, Ny);
    MatSetValue(A, r, r , -2.0/h, INSERT_VALUES);
    MatSetValue(A, r, r-1 , 1.0/h, INSERT_VALUES);
    MatSetValue(A, r, r+(Ny-1) , 1.0/h, INSERT_VALUES);

    r = rowStart + PHI_IND(Nx - 2, 0, Ny);
    MatSetValue(A, r, r , -2.0/h, INSERT_VALUES);
    MatSetValue(A, r, r+1 , 1.0/h, INSERT_VALUES);
    MatSetValue(A, r, r-(Ny-1) , 1.0/h, INSERT_VALUES);

    r = rowStart + PHI_IND(Nx - 2, Ny-2, Ny);
    MatSetValue(A, r, r , -2.0/h, INSERT_VALUES);
    MatSetValue(A, r, r-1 , 1.0/h, INSERT_VALUES);
    MatSetValue(A, r, r-(Ny-1) , 1.0/h, INSERT_VALUES);

    r = (Ny-1)/2;
    MatSetValue(A, r, r, 1.0, INSERT_VALUES);
    MatSetValue(A, r, r+1, 0.0, INSERT_VALUES);
    MatSetValue(A, r, r-1, 0.0, INSERT_VALUES);
    MatSetValue(A, r, r+(Ny-1), 0.0, INSERT_VALUES);
}

/*To call during the initialization of your solver, before the begin of the time loop*/
/*Maybe you should need to add an argument to specify the number of unknows*/
/*Modification to do in this function :*/
/*  OKKKKK -Specify the number of unknows*/
/*  OKKKKK -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver_otis(Poisson_data* data, int Nx, int Ny, double h)
{
    PetscInt rowStart; /*rowStart = 0*/
    PetscInt rowEnd; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

	int nphi = (Nx - 1) * (Ny - 1); /*WRITE HERE THE NUMBER OF UNKNOWS*/

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

    computeLaplacianMatrix_otis(data->A, rowStart, rowEnd, Nx, Ny, h);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    // MatView(data->A, PETSC_VIEWER_STDOUT_WORLD);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    KSPSetType(data->sles,KSPGMRES); //KSPGMRES seems the best, it will not be used if PC LU.
    PC prec;
    KSPGetPC(data->sles, &prec);
    PCSetType(prec,PCLU);
    KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles,1,4);
    KSPGMRESSetPreAllocateVectors(data->sles);

    // PetscPrintf(PETSC_COMM_WORLD, "Assembly of Mattrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
/*Modification to do : nothing */
void free_poisson_solver_otis(Poisson_data* data){
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
