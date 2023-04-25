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
int convective_velocity(matrix *V, matrix *U, matrix *H_x, matrix *H_y, double delta_x, double delta_y, int m, int n){
    int i,j;
    //calcul of the convective term away from the side of the domain 
    for (j = 0; j<m+1; j++){
        V->a[0][j] = -1/5 *(V->a[3][j] - 5*V->a[2][i] + 15*V->a[1][j]);
        V->a[n+1][j] = -1/5 *(V->a[n+1-3][j] - 5*V->a[n+1-2][i] + 15*V->a[n+1-1][j]);
    }
    for (i = 0; i<n+1; i++){
        U->a[i][0] = -1/5 *(U->a[3][i] - 5*U->a[2][i] + 15*U->a[1][i]);
        U->a[i][m+1] = U->[i][m];//-1/5 *(U->a[m+1-3][i] - 5*U->a[m+1-2][i] + 15*U->a[m+1-1][i])
    }
    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){
            H_x->a[i][j] = 1.0/(4.0*delta_x) * (U->a[i+1][j]*U->a[i+1][j] - U->a[i-1][j]*U->a[i-1][j]) 
                            + 1/(4*delta_y) * ((V->a[i][j+1] + V->a[i-1][j+1]) * (U->a[i][j+1] - U->a[i][j])
                            + (V->a[i-1][j] + V->a[i-1][j-1]) * (U->a[i][j] - U->a[i][j-1]));
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            H_y->a[i][j] = 1.0/(4.0*delta_x) * ((U->a[i+1][j] + U->a[i+1][j-1]) * (V->a[i+1][j] - V->a[i][j]) 
                            + (U->a[i][j] + U->[i][j-1]) * (V->a[i][j] - V->a[i-1][j]))
                            + 1/(4*delta_y) * (V->a[i][j+1]*V->a[i][j+1] - V->a[i][j-1]*V->a[i][j-1]);
        }
    }
    return 0;

}
int pressure_gradure(matrix *grad_Px, matrix *grad_Py, matrix *P, double delta_x){
    for (int i=0; i<P->n-1; i++){
        for (int j=0;j<P->m; j++ ){
            grad_Px->a[i][j] = 1/delta_x (P->a[i+1][j] - P->a[i][j]); 
        }
    }
    // for (int j=0;j<P->m; j++ ){
    //         grad_Px->a[n][j] = 1/delta_x ( - P->a[n][j]) 
    //     }
    for (int i=0; i<P->n; i++){
        for (int j=0;j<P->m-1; j++ ){
            grad_Px->a[i][j] = 1/delta_x (P->a[i][j+1] - P->a[i][j]); 
        }
        // grad_Px->a[i][m] = 1/delta_x (pressbords - P->a[i][m])
    }
    return 0; 
}

int laplacian_velocity(matrix *LapU, matrix *LapV, matrix *U, double delta_x, double delta_y, int m, int n){
    for (int i=1; i<n; i++){
        for (int j=1;j<m; j++ ){
            lapV->a[i][j] = 1/(delta_x*delta_x)*((V->a[i+1][j] - 2*V->a[i][j] + V->a[i-1][j])) 
                          + 1/(delta_y*delta_y)*(V->a[i][j+1] - 2*V->a[i][j] + V->a[i][j-1]); //A verifier on a compute ca trop vite

            LapU->a[i][j] = 1/(delta_x*delta_x)*((U->a[i+1][j] - 2*U->a[i][j] + U->a[i-1][j])) 
                          + 1/(delta_y*delta_y)*(U->a[i][j+1] - 2*U->a[i][j] + U->a[i][j-1]); //A verifier on a compute ca trop vite 
        }
    }
    return 0;
}

int laplacian_temperature(matrix *LapT, matrix *T, double delta_x, double delta_y, int m, int n){
    for (int i=1; i<n; i++){
        for (int j=1;j<m; j++){
            LapT->[i][j] = 1/(delta_x*delta_x)*((T->a[i+1][j] - 2*T->a[i][j] + T->a[i-1][j])) 
                      + 1/(delta_y*delta_y)*(T->a[i][j+1] - 2*T->a[i][j] + T->a[i][j-1]);
        }
    }
    return 0;
}

int convective_temperature(matrix Ht, matrix *T, double delta_x, double delta_y, int m, int n){
    
}
int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);
    int m = 16; 
    int n = 16; 
    double delta_x ; //to specify 
    double delta_y ; //to specify 
    
    /*WRITE YOUR PROJECT ...*/
    matrix *V = callocate_matrix(m+1, n+3);
    matrix *U = callocate_matrix(m+3, n+1); 
    matrix *P = callocate_matrix(m, n); 
    matrix *phi = callocate_matrix(m, n); 
    matrix *H_y = callocate_matrix(m+1, n+3);
    matrix *H_x = callocate_matrix(m+3, n+1);
    matrix *H_y_nplus1 = callocate_matrix(m+1, n+3);
    matrix *H_x_nplus1 = callocate_matrix(m+3, n+1);

    // Boucle temporel


    PetscFinalize();
    free_matrix(V);
    free_matrix(U);
    free_matrix(P);
    free_matrix(phi);

}
