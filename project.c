#include "poisson.h"

matrix* allocate_matrix(int m, int n) {
	matrix *mat = (matrix*) malloc(sizeof(matrix));
	mat->m = m, mat->n = n;
	mat->data = (double*)malloc(m*n*sizeof(double));
    mat->LU = 0;
	if(mat->data == NULL) return NULL;
	mat->a = (double**)malloc(m*sizeof(double*));
	if (mat->a == NULL) return NULL;
	for (int i = 0; i < m; i++)
		mat->a[i] = mat->data+i*n;
    mat->p = (int *) malloc(sizeof(int)*n);
    for (int i = 0; i<n; i++) {
        mat->p[i] = i;
    }
	return mat;
}

void free_matrix(matrix *mat) {
	if(mat == NULL) return;
	free(mat->a);
	free(mat->data); 
    free(mat->p);
    free(mat);   
}


int write_matrix_to_file(matrix *mat, char const *fileName, int nb_equation){
    FILE *file =fopen(fileName, "w");
    if (file == NULL)return -1;
    for (int i = 0; i <nb_equation;i++){
        for (int j=0; j < nb_equation;j++){
            fprintf(file, "%.6e ", mat->a[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    return 0;
}

int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);
    int m = 16; 
    int n = 16; 
    /*WRITE YOUR PROJECT ...*/
    matrix V;
    matrix U; 
    matrix V;
    matrix P; 
    matrix phi; 

    PetscFinalize();

}
