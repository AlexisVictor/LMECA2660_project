#ifndef _POISSON_H_
#define _POISSON_H_

/*To include in the file in which you will call initialize_poisson_solver and poisson_solver*/

#include <petsc.h>
#include <petscsys.h>

//Structure storing petsc vectors

typedef struct {
	Vec b;
	Vec x;
	Mat A;
	KSP sles;
} Poisson_data;

typedef struct {
	int m, n; // dimensions de la matrice
	double *data; // tableau 1D de taille m*n contenant les entrees de la matrice
	double **a; // tableau 1D de m pointeurs vers chaque ligne, pour pouvoir appeler a[i][j]
} matrix;

matrix* callocate_matrix(int m, int n);

PetscErrorCode initialize_poisson_solver(Poisson_data* data, int m, int n);
void poisson_solver(Poisson_data *data, double inv_delta_tx,int m, int n, matrix *U, matrix *V, matrix *phi);
void free_poisson_solver(Poisson_data* data);

#endif

