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
int convective_velocity(matrix *V, matrix *U, matrix *H_x, matrix *H_y, double deltax, double deltay, int m, int n){
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
            H_x->a[i][j] = 1.0/(4.0*deltax) * (U->a[i+1][j]*U->a[i+1][j] - U->a[i-1][j]*U->a[i-1][j]) 
                            + 1/(4*deltay) * ((V->a[i][j+1] + V->a[i-1][j+1]) * (U->a[i][j+1] - U->a[i][j])
                            + (V->a[i-1][j] + V->a[i-1][j-1]) * (U->a[i][j] - U->a[i][j-1]));
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            H_y->a[i][j] = 1.0/(4.0*deltax) * ((U->a[i+1][j] + U->a[i+1][j-1]) * (V->a[i+1][j] - V->a[i][j]) 
                            + (U->a[i][j] + U->[i][j-1]) * (V->a[i][j] - V->a[i-1][j]))
                            + 1/(4*deltay) * (V->a[i][j+1]*V->a[i][j+1] - V->a[i][j-1]*V->a[i][j-1]);
        }
    }
    return 0;

}
int pressure_gradure(matrix *grad_Px, matrix *grad_Py, matrix *P, double deltax){
    for (int i=0; i<P->n-1; i++){
        for (int j=0;j<P->m; j++ ){
            grad_Px->a[i][j] = 1/deltax (P->a[i+1][j] - P->a[i][j]); 
        }
    }
    // for (int j=0;j<P->m; j++ ){
    //         grad_Px->a[n][j] = 1/deltax ( - P->a[n][j]) 
    //     }
    for (int i=0; i<P->n; i++){
        for (int j=0;j<P->m-1; j++ ){
            grad_Py->a[i][j] = 1/deltay (P->a[i][j+1] - P->a[i][j]); 
        }
        // grad_Px->a[i][m] = 1/deltax (pressbords - P->a[i][m])
    }
    return 0; 
}

int laplacian_velocity(matrix *LapU, matrix *LapV, matrix *U, matrix *V, double deltax, double deltay, int m, int n){
    for (int i=1; i<n; i++){
        for (int j=1;j<m; j++ ){
            lapV->a[i][j] = 1/(deltax*deltax)*((V->a[i+1][j] - 2*V->a[i][j] + V->a[i-1][j])) 
                          + 1/(deltay*deltay)*(V->a[i][j+1] - 2*V->a[i][j] + V->a[i][j-1]); //A verifier on a compute ca trop vite

            LapU->a[i][j] = 1/(deltax*deltax)*((U->a[i+1][j] - 2*U->a[i][j] + U->a[i-1][j])) 
                          + 1/(deltay*deltay)*(U->a[i][j+1] - 2*U->a[i][j] + U->a[i][j-1]); //A verifier on a compute ca trop vite 
        }
    }
    return 0;
}

int laplacian_temperature(matrix *LapT, matrix *T, double deltax, double deltay, int m, int n){
    for (int i=1; i<n; i++){
        for (int j=1;j<m; j++){
            LapT->[i][j] = 1/(deltax*deltax)*((T->a[i+1][j] - 2*T->a[i][j] + T->a[i-1][j])) 
                      + 1/(deltay*deltay)*(T->a[i][j+1] - 2*T->a[i][j] + T->a[i][j-1]);
        }
    }
    return 0;
}

int convective_temperature(matrix *Ht, matrix *T, double deltax, double deltay, int m, int n){
    for (int i=1; i<n+1; i++){
        for (int j=1;j<m+1; j++){
            Ht->a[i][j] = 0.5/deltax*(U->a[i][j]*(T->a[i][j]-T->a[i-1][j])) + (U->a[i+1][j]*(T->a[i+1][j]-T->a[i][j])) 
                        + 0.5/deltay*(V->a[i][j]*(T->a[i][j]-T->a[i][j-1])) + (V->a[i][j+1]*(T->a[i][j+1]-T->a[i][j]));
           
        }
        Ht->a[i][m+1] = -1/5*(T->a[i][m-2] - 5)
    }
    return 0;
}

int BorderTemp(matrix *T, double deltax, double deltay, int m, int n){

    /*Plus sur de l'ordre des index i et j à checker absolument!!!*/
    double l0 = 1e-3 ;
    for (int j=1;j<m; j++){
        //Left Surface : 
        T->a[0][j] = T->a[1][j];
        //Right Surface : 
        T->a[n+1][j] = T->a[n][j];
        }    
    double Tgamma;
    for (int i=1;i<n; i++){
        //Free Surface : 
        Tgamma = l0/(2*deltay+l0) * T->a[i][m];
        T->a[i][m+1] = -1/5* (T->a[i][m-2]-5*T->a[i][m-1]+15*T->a[i][m]-16*Tgamma);

        //Bottom Surface : 
        Tgamma = 2*deltay + T->a[i][1];
        T->a[i][0] = -1/5* (T->a[i][3]-5*T->a[i][2]+15*T->a[i][1]-16*Tgamma);
    }
    return 0;

}

int evalEstimVelocity(matrix *U, matrix *V, matrix *Hy, matrix *Hx, matrix *LapU, matrix *LapV,
                      matrix *Hyold, matrix *Hxold, matrix *grad_Px, matrix *grad_Py, matrix *T, matrix *Uestim, 
                      matrix *Vestim, double dt, double Gr){
    double temp; 

    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){
            temp = -1/2*(3*Hx->a[i][j] - Hxold->a[i][j]) - grad_Px->a[i][j] + Gr**(-1/2)*LapU->a[i][j];
            Uestim->a[i][j] = temp*dt + U->a[i][j];
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            temp = -1/2*(3*Hy->a[i][j] - Hyold->a[i][j]) - grad_Py->a[i][j] + Gr**(-1/2)*LapV->a[i][j] - T->a[i][j];
            Vestim->a[i][j] = temp*dt + V->a[i][j];
        }
    }
    return 0;
}   

int evalVelocity(matrix *phi, double dt, double deltax, double deltay, matrix *U, matrix *V, matrix *Uestim, matrix *Vestim){

    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){
            temp = (phi->a[i][j]-phi->a[i][j-1])/deltax;
            U->a[i][j] = temp*dt + Uestim->a[i][j];
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            temp = (phi->a[i][j]-phi->a[i-1][j])/deltay;
            V->a[i][j] = temp*dt + Vestim->a[i][j];
        }
    }
    return 0;
}

int evalTemperature(matrix *T, matrix *Ht, matrix *Htold, matrix *laplacian_temperature, double Pr, double Grashof){   
    double temp;
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m+1; j++){
            temp = -1/2*(3*Ht->a[i][j] - Htold->a[i][j]) + 1/(Gr**(1/2)*Pr)*LapT->a->[i][j] ;
            T->a[i][j] += temp*dt ;
        }
    }
    return 0;
}

/*
int initialPressure(){

    ICI FAIRE INTERVENIR LA DEPENDANCE EN Y DE LA PRESSION CINEMATIQUE adimensionnelle:

    DISCRETISER POUR AVOIR UNE VALEUR EN CHAQUE POINT : y = 1/m car on est en adimensionnel  

    ----> En fait je crois qu'on s'en balec pcq ce qui nous intéresse c'est tjs le gradient de pression et il est nul en t=0
    return 0;
}





*/

int main(int argc, char *argv[]){

    PetscInitialize(&argc, &argv, 0, 0);
    int n = 16; 
    int m = (int) 1.5*n; 
    double deltax ; //to specify 
    double deltay ; //to specify 
    double SimTime; //to specify
    double t = 0;
    double dt; //to specify
    
    /*WRITE YOUR PROJECT ...*/
    matrix *V = callocate_matrix(m+1, n+2);
    matrix *U = callocate_matrix(m+2, n+1); 
    matrix *Vestim = callocate_matrix(m+1, n+2);
    matrix *Uestim= = callocate_matrix(m+2, n+1);
    matrix *P = callocate_matrix(m, n); 
    matrix *T = callocate_matrix(m+2, n+2); 
    matrix *phi = callocate_matrix(m, n); 
    matrix *Hy = callocate_matrix(m+1, n+2);
    matrix *Hx = callocate_matrix(m+2, n+1);
    matrix *Ht = callocate_matrix(m+2, n+2);
    matrix *Hyold = callocate_matrix(m+1, n+2);
    matrix *Hxold = callocate_matrix(m+2, n+1);
    matrix *Htold = callocate_matrix(m+2, n+2);
    matrix *grad_Px = callocate_matrix(m, n);
    matrix *grad_Py = callocate_matrix(m, n);
    matrix *LapU = callocate_matrix(m+2, n+1);
    matrix *LapV = callocate_matrix(m+1, n+2);
    matrix *LapT = callocate_matrix(m+2, n+2);
    
    double Pr = 4;
    double Gr = 10e10

    // Boucle temporel
    int iter = 0;
    while (t<SimTime){

        //Actualisation
        memcpy( Hxold->a, Hx->a,  sizeof(double)*(m+2)*(n+1) ); 
        memcpy( Hyold->a, Hy->a,  sizeof(double)*(m+1)*(n+2) ); 
        memcpy( Htold->a, Ht->a,  sizeof(double)*(m+2)*(n+2) ); 
    
        convective_velocity(V, U, Hx, Hy, deltax, deltay, m, n);
        convective_temperature(Ht, T, deltax, deltay, m, n);
        laplacian_velocity(LapU, LapV, U, V, deltax, deltay, m, n);
        laplacian_temperature(LapT, T, deltax, deltay, m, n);
        pressure_gradure(grad_Px, grad_Py, P, deltax);

        evalEstimVelocity(U, V, Hy, Hx, LapU, LapV,
                        Hyold, Hxold, grad_Px, grad_Py, T, Uestim, 
                        Vestim, dt, Gr);
        /* Poisson */
        evalVelocity(phi,  dt,  deltax,  deltay, U, V, Uestim, Vestim);
        for (i = 0; i < n; i++){
            for (j = 0; j< m; j++){
                P->a[i][j] += phi->a[i][j];
            }
        }
        evalTemperature(T, Ht, Htold, laplacian_temperature, Pr, Grashof);
        BorderTemp(T, deltax,  deltay, m, n, qwallk, l0);

        iter++;
        t += dt;
    }


    PetscFinalize();
    free_matrix(V);
    free_matrix(U); 
    free_matrix(Vestim);
    free_matrix(Uestim);
    free_matrix(P); 
    free_matrix(T); 
    free_matrix(phi); 
    free_matrix(Hy);
    free_matrix(Hx);
    free_matrix(Ht);
    free_matrix(Hyold);
    free_matrix(Hxold);
    free_matrix(Htold);
    free_matrix(grad_Px);
    free_matrix(grad_Py);
    free_matrix(LapU);
    free_matrix(LapV);
    free_matrix(LapT);

}