#include "poisson.h"

matrix* callocate_matrix(int m, int n) {
    // size_t nitems, size_t size
	matrix *mat = (matrix*) malloc(sizeof(matrix));
	mat->m = m, mat->n = n;
	mat->data = (double*) calloc(m*n, sizeof(double));
	if(mat->data == NULL) return NULL;
	mat->a = (double**)malloc(m*sizeof(double*));
	if (mat->a == NULL) return NULL;
	for (int i = 0; i < m; i++)
		mat->a[i] = mat->data+i*n;
	return mat;
}

void free_matrix(matrix *mat) {
	if(mat == NULL) return;
	free(mat->a);
	free(mat->data); 
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
    double delta_x ; //to specify 
    double delta_y ; //to specify 
    
    /*WRITE YOUR PROJECT ...*/
    matrix *V = callocate_matrix(m+1, n+1);
    matrix *U = callocate_matrix(m+1, n+1); 
    matrix *P = callocate_matrix(m, n); 
    matrix *phi = callocate_matrix(m, n); 
    matrix *H_y = callocate_matrix(m+1, n+1);
    matrix *H_x = callocate_matrix(m+1, n+1);

    // Boucle temporel 

    //calcul of the convective term away from the side of the domain 
    for (int i = 1; i < n; i++){
        for (int j = 1; j< m; j++ ){
            H_x->a[i][j] = 1.0/(4.0*delta_x) * (U->a[i+1][j]*U->a[i+1][j] - U->a[i-1][j]*U->a[i-1][j]) 
                            + 1/(4*delta_y) * ((V->a[i][j+1] + V->a[i-1][j+1]) * (U->a[i][j+1] - U->a[i][j])
                            + (V->a[i-1][j] + V->a[i-1][j-1]) * (U->a[i][j] - U->a[i][j-1]));
            H_y->a[i][j] = 1.0/(4.0*delta_x) * ((U->a[i+1][j] + U->a[i+1][j-1]) * (V->a[i+1][j] - V->a[i][j]) 
                            + (U->a[i][j] + U->[i][j-1]) * (V->a[i][j] - V->a[i-1][j]))
                            + 1/(4*delta_y) * (V->a[i][j+1]*V->a[i][j+1] - V->a[i][j-1]*V->a[i][j-1]);
        }
    }

    PetscFinalize();
    free_matrix(V);
    free_matrix(U);
    free_matrix(P);
    free_matrix(phi);

}
