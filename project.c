#include "poisson.h"
#include "math.h"
#include "stdio.h"


void copy_matrix(matrix *dest, matrix *src) {
    if (dest->m != src->m || dest->n != src->n) {
        fprintf(stderr, "Error: matrices have different dimensions\n");
        exit(1);
    }
    for (int i = 0; i < dest->n; i++) {
        for (int j = 0; j < dest->m; j++) {
            dest->a[i][j] = src->a[i][j];
        }
    }
}


matrix* callocate_matrix(int m, int n) {
    // size_t nitems, size_t size
	matrix *mat = (matrix*) malloc(sizeof(matrix));
	mat->m = m, mat->n = n;
	mat->data = (double*) calloc(m*n, sizeof(double));
	if(mat->data == NULL) return NULL;
	mat->a = (double**)malloc(n*sizeof(double*));
    // mat->a = (double**)malloc(m*sizeof(double*));
	if (mat->a == NULL) return NULL;
	// for (int i = 0; i < m; i++)
	// 	mat->a[i] = mat->data+i*n;
    for (int i = 0; i < n; i++)
		mat->a[i] = mat->data+i*m;
	return mat;
}

void MatrixWriteFile(const char *baseResultName, int iter, matrix *mat, double t)
{
    int i,j;
    const char *basename = "results/%s-%d.txt";
    char filename[100];
    sprintf(filename,basename,baseResultName,iter);
    FILE* file = fopen(filename,"w");
    if (file == NULL) {
        printf("Error : cannot create result file : did you create the output directory ? \n");
        exit(0); }
    
    for (i = 0; i < mat->n; ++i) {
    	for (j = 0; j < mat->m; ++j) {
            fprintf(file,"%d;%d;%le;\n",i,j,mat->a[i][j]);
            }
        }

    fclose(file);
}

void free_matrix(matrix *mat) {
	if(mat == NULL) return;
	free(mat->a);
	free(mat->data); 
    free(mat);   
}


int print_matrix(matrix *mat){
    for (int i = 0; i <mat->n;i++){
        for (int j=0; j < mat->m;j++){
            printf("%.6e ", mat->a[i][j]);
        }
        printf("\n");
    }
    return 0;
}
int convective_velocity(matrix *V, matrix *U, matrix *H_x, matrix *H_y, double deltax, double deltay, int m, int n){
    int i,j;
    //calcul of the convective term away from the side of the domain 
    for (j = 0; j<m+1; j++){
        V->a[0][j] = -1.0/5 *(V->a[3][j] - 5*V->a[2][j] + 15*V->a[1][j]);
        V->a[n+1][j] = -1.0/5 *(V->a[n+1-3][j] - 5*V->a[n+1-2][j] + 15*V->a[n+1-1][j]);
    }
    for (i = 0; i<n+1; i++){
        U->a[i][0] = -1.0/5 *(U->a[i][3] - 5*U->a[i][3] + 15*U->a[i][1]);
        U->a[i][m+1] = U->a[i][m];//-1/5 *(U->a[m+1-3][i] - 5*U->a[m+1-2][i] + 15*U->a[m+1-1][i])
    }
    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){ // decalle les indices des V de 1 vers la droite 
            // H_x->a[i][j] = 1.0/(4.0*deltax) * (U->a[i+1][j]*U->a[i+1][j] - U->a[i-1][j]*U->a[i-1][j]) 
            //                 + 1.0/(4*deltay) * ((V->a[i][j+1] + V->a[i-1][j+1]) * (U->a[i][j+1] - U->a[i][j])
            //                                   + (V->a[i-1][j] + V->a[i-1][j-1]) * (U->a[i][j] - U->a[i][j-1]));
            H_x->a[i][j] = 1.0/(4.0*deltax) * (U->a[i+1][j]*U->a[i+1][j] - U->a[i-1][j]*U->a[i-1][j])  //ligne pas encore verif 
                            + 1.0/(4*deltay) * ((V->a[i][j] + V->a[i+1][j]) * (U->a[i][j+1] - U->a[i][j])
                                              + (V->a[i][j-1] + V->a[i+1][j-1]) * (U->a[i][j] - U->a[i][j-1]));
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){ // decalle les indice des U de 1 vers le le haute
            H_y->a[i][j] = 1.0/(4.0*deltax) * ((U->a[i][j] + U->a[i][j-1]) * (V->a[i+1][j] - V->a[i][j]) 
                                             + (U->a[i-1][j] + U->a[i-1][j-1]) * (V->a[i][j] - V->a[i-1][j]))
                            + 1.0/(4*deltay) * (V->a[i][j+1]*V->a[i][j+1] - V->a[i][j-1]*V->a[i][j-1]);
        }
    }
    return 0;
}

int convective_velocity_div(matrix *V, matrix *U, matrix *H_x, matrix *H_y, double deltax, double deltay, int m, int n){
    int i,j;
    //calcul of the convective term away from the side of the domain 
    for (j = 0; j<m+1; j++){
        V->a[0][j] = -1.0/5 *(V->a[3][j] - 5*V->a[2][j] + 15*V->a[1][j]);
        V->a[n+1][j] = -1.0/5 *(V->a[n+1-3][j] - 5*V->a[n+1-2][j] + 15*V->a[n+1-1][j]);
    }
    for (i = 0; i<n+1; i++){
        U->a[i][0] = -1.0/5 *(U->a[i][3] - 5*U->a[i][3] + 15*U->a[i][1]);
        U->a[i][m+1] = U->a[i][m];//-1/5 *(U->a[m+1-3][i] - 5*U->a[m+1-2][i] + 15*U->a[m+1-1][i])
    }
    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){                                 
            // H_x->a[i][j] = 1.0/(4.0*deltax) * (U->a[i+1][j]*U->a[i+1][j] - 2*U->a[i+1][j]*U->a[i][j] - U->a[i-1][j]*U->a[i-1][j] - 2*U->a[i][j]*U->a[i-1][j])  
            //             + 1.0/(4*deltay) * (U->a[i][j+1]*V->a[i+1][j]) + U->a[i][j+1]*V->a[i][j-1] - U->a[i][j-1]*V->a[i+1][j-1] - U->a[i][j-1]*V->a[i][j-1];
            H_x->a[i][j] =  1.0/(4.0*deltax) * ((U->a[i+1][j] + U->a[i][j]) * (U->a[i+1][j] + U->a[i][j]) - (U->a[i][j] + U->a[i-1][j]) * (U->a[i][j] + U->a[i-1][j]))
                          + 1.0/(4.0*deltay) * ((U->a[i][j+1] + U->a[i][j]) * (V->a[i+1][j] + V->a[i-1][j]) - (U->a[i][j] + U->a[i][j-1]) * (V->a[i+1][j-1] + V->a[i][j-1]));

        }
    }
    for (i = 1; i < n+1; i++){
         for (j = 1; j< m; j++){ // decalle les indice des U de 1 vers le le haute
            H_y->a[i][j] = 1.0/(4.0*deltax) * ((U->a[i][j+1] + U->a[i][j]) * (V->a[i+1][j] + V->a[i][j]) - (U->a[i-1][j+1] + U->a[i-1][j]) * (V->a[i][j] + V->a[i-1][j]))
                         + 1.0/(4.0*deltay) * ((V->a[i][j+1] + V->a[i][j]) * (V->a[i][j+1] + V->a[i][j]) - (V->a[i][j] + V->a[i][j-1]) * (V->a[i][j] + V->a[i][j-1]));
    //         H_y->a[i][j] = 1.0/(4.0*deltax) * ((U->a[i][j] + U->a[i][j-1]) * (V->a[i+1][j] - V->a[i][j]) 
    //                                          + (U->a[i-1][j] + U->a[i-1][j-1]) * (V->a[i][j] - V->a[i-1][j]))
    //                         + 1.0/(4*deltay) * (V->a[i][j+1]*V->a[i][j+1] - V->a[i][j-1]*V->a[i][j-1]);
         }
    }
    return 0;

}


int pressure_gradure(matrix *grad_Px, matrix *grad_Py, matrix *P, double deltax, double deltay, int m, int n){
    for (int i=1; i<n; i++){
        for (int j=1;j<P->m+1; j++ ){
            grad_Px->a[i][j] = 1.0/deltax*(P->a[i][j-1] - P->a[i-1][j-1]); 
        }
    }
    // for (int j=0;j<P->m; j++ ){
    //         grad_Px->a[n][j] = 1/deltax ( - P->a[n][j]) 
    //     }
    for (int i=1; i<P->n+1; i++){
        for (int j=1;j<P->m; j++ ){
            grad_Py->a[i][j] = 1.0/deltay*(P->a[i-1][j] - P->a[i-1][j-1]); 
        }
        // grad_Px->a[i][m] = 1/deltax (pressbords - P->a[i][m])
    }
    return 0; 
} // checked 

int laplacian_velocity(matrix *LapU, matrix *LapV, matrix *U, matrix *V, double deltax, double deltay, int m, int n){
    for (int i=1; i<n; i++){
        for (int j=1;j<m+1; j++){
            LapU->a[i][j] = 1.0/(deltax*deltax)*((U->a[i+1][j] - 2*U->a[i][j] + U->a[i-1][j])) 
                          + 1.0/(deltay*deltay)*(U->a[i][j+1] - 2*U->a[i][j] + U->a[i][j-1]); //A verifier on a compute ca trop vite 
        }
    }
    for (int i=1; i<n+1; i++){
        for (int j=1;j<m; j++ ){
            LapV->a[i][j] = 1.0/(deltax*deltax)*((V->a[i+1][j] - 2*V->a[i][j] + V->a[i-1][j])) 
                          + 1.0/(deltay*deltay)*(V->a[i][j+1] - 2*V->a[i][j] + V->a[i][j-1]); //A verifier on a compute ca trop vite
            } 
        }
    return 0;
} // vite fait checked

int laplacian_temperature(matrix *LapT, matrix *T, double deltax, double deltay, int m, int n){
    for (int i=1; i<n+1; i++){
        for (int j=1;j<m+1; j++){
            LapT->a[i][j] = 1/(deltax*deltax)*((T->a[i+1][j] - 2*T->a[i][j] + T->a[i-1][j])) 
                      + 1/(deltay*deltay)*(T->a[i][j+1] - 2*T->a[i][j] + T->a[i][j-1]);
        }
    }
    return 0;
} // checked vite fait

int convective_temperature(matrix *Ht, matrix *T, matrix *U, matrix *V, double deltax, double deltay, int m, int n){
    for (int i=1; i<n+1; i++){
        for (int j=1;j<m+1; j++){
            Ht->a[i][j] = (1.0/2.0*deltax) * (U->a[i-1][j] * (T->a[i][j]-T->a[i-1][j])) + (U->a[i][j]*(T->a[i+1][j]-T->a[i][j])) 
                        + (1.0/2.0*deltay) * (V->a[i][j-1] *(T->a[i][j]-T->a[i][j-1])) + (V->a[i][j]*(T->a[i][j+1]-T->a[i][j]));
           
        }
    }
    return 0;
} // Checked, on des kog 



int evalEstimVelocity(matrix *U, matrix *V, matrix *Hy, matrix *Hx, matrix *LapU, matrix *LapV,
                      matrix *Hyold, matrix *Hxold, matrix *grad_Px, matrix *grad_Py, matrix *T, matrix *Uestim, 
                      matrix *Vestim, int m, int n, double dt, double Gr){
    double temp; 
    int i,j;
    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){
            temp = -1.0/2.0*(3.0*Hx->a[i][j] - Hxold->a[i][j]) - grad_Px->a[i][j] + pow(Gr,-1.0/2.0)*LapU->a[i][j];
            Uestim->a[i][j] = temp*dt + U->a[i][j];
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            temp = -1.0/2.0*(3.0*Hy->a[i][j] - Hyold->a[i][j]) - grad_Py->a[i][j] + pow(Gr,-1.0/2.0)*LapV->a[i][j] - (T->a[i][j]+T->a[i][j+1])/2.0;
            Vestim->a[i][j] = temp*dt + V->a[i][j];
        }
    }
    return 0;
}   //probleme au niveau de T (faire la moyenne pour avoir T au point souhaité)

int evalVelocity(matrix *phi, double dt, double deltax, double deltay, matrix *U, matrix *V, int m, int n, matrix *Uestim, matrix *Vestim){
    int i,j;
    double temp;
    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){
            temp = (phi->a[i][j-1]-phi->a[i-1][j-1])/deltax;
            U->a[i][j] = temp*dt + Uestim->a[i][j];
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            temp = (phi->a[i-1][j]-phi->a[i-1][j-1])/deltay;
            V->a[i][j] = temp*dt + Vestim->a[i][j];
        }
    }
    return 0;
}  //checked 

int evalTemperature(matrix *T, matrix *Ht, matrix *Htold, matrix *LapT, double deltax, double deltay, int m, int n, double Pr, double dt, double Gr){   
    double temp;
    int i,j;
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m+1; j++){
            temp = -1.0/2*(3*Ht->a[i][j] - Htold->a[i][j]) + 1.0/(pow(Gr,1.0/2.0)*Pr)*LapT->a[i][j] ;
            T->a[i][j] += temp*dt ;
        }
    }
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
        T->a[i][m+1] = -1.0/5 * (T->a[i][m-2]-5*T->a[i][m-1]+15*T->a[i][m]-16*Tgamma);

        //Bottom Surface : 
        Tgamma = 2*deltay + T->a[i][1];
        T->a[i][0] = -1.0/5 * (T->a[i][3]-5*T->a[i][2]+15*T->a[i][1]-16*Tgamma);
    }
    return 0;
} //potentiel erreur au niveau des conditions limites mais a priori ok 

/*
int initialPressure(){

    ICI FAIRE INTERVENIR LA DEPENDANCE EN Y DE LA PRESSION CINEMATIQUE adimensionnelle:

    DISCRETISER POUR AVOIR UNE VALEUR EN CHAQUE POINT : y = 1/m car on est en adimensionnel  

    ----> En fait je crois qu'on s'en balec pcq ce qui nous intéresse c'est tjs le gradient de pression et il est nul en t=0
    return 0;
}





*/

int main(int argc, char *argv[]){

    double L = 1.0;
    double H = 1.5*L;

    int n = 50 ; 
    int m = (int) (1.5*n);
   
    double deltax = L/n; //to specify 
    double deltay = H/m; //to specify 

    // printf("deltax = %.6e \t; deltay = %.6e \n", deltax, deltay);
    // printf("n = %d \t; m = %d \n", n, m);

    double SimTime = 2e-1; //to specify
    double t = 0.0;
    double dt= 2e-3; //to specify
    double Pr = 4.0;
    double Gr = 10e10;
    int iter = 0;
        
    /*WRITE YOUR PROJECT ...*/
    matrix *U = callocate_matrix(m+2, n+1); 
    matrix *V = callocate_matrix(m+1, n+2);
    matrix *LapU = callocate_matrix(m+2, n+1);
    matrix *LapV = callocate_matrix(m+1, n+2);
    matrix *Vestim = callocate_matrix(m+1, n+2);
    matrix *Uestim = callocate_matrix(m+2, n+1);
    matrix *grad_Px = callocate_matrix(m+2, n+1);
    matrix *grad_Py = callocate_matrix(m+1, n+2);
    matrix *Hx = callocate_matrix(m+2, n+1); 
    matrix *Hxold = callocate_matrix(m+2, n+1);
    matrix *Hyold = callocate_matrix(m+1, n+2);
    matrix *Hy = callocate_matrix(m+1, n+2);
    matrix *P = callocate_matrix(m, n); 
    matrix *T = callocate_matrix(m+2, n+2); 
    matrix *LapT = callocate_matrix(m+2, n+2);
    matrix *Ht = callocate_matrix(m+2, n+2);
    matrix *phi = callocate_matrix(m, n); 
    matrix *Htold = callocate_matrix(m+2, n+2);
    Poisson_data *data = (Poisson_data *)malloc(sizeof(Poisson_data));
    double inv_delta_tx = deltax/dt;

    MatrixWriteFile("U", iter, U, t);
    MatrixWriteFile("V", iter, V, t);
    MatrixWriteFile("P", iter, P, t);
    MatrixWriteFile("phi", iter, phi, t);
    MatrixWriteFile("T", iter, T, t);

    PetscInitialize(&argc, &argv, 0, 0);    
    while (t<SimTime){

        copy_matrix(Hxold, Hx);
        copy_matrix(Hyold, Hy);
        copy_matrix(Htold, Ht);
    
        // convective_velocity_div(V, U, Hx, Hy, deltax, deltay, m, n);
        convective_velocity(V, U, Hx, Hy, deltax, deltay, m, n);
        convective_temperature(Ht, T, U, V, deltax, deltay, m, n);
        laplacian_velocity(LapU, LapV, U, V, deltax, deltay, m, n); 
        laplacian_temperature(LapT, T, deltax, deltay, m, n); 
        pressure_gradure(grad_Px, grad_Py, P, deltax, deltay, m, n);
        evalEstimVelocity(U, V, Hy, Hx, LapU, LapV, Hyold, Hxold, grad_Px, grad_Py,
                          T, Uestim, Vestim, m, n, dt, Gr); 
        
        MatrixWriteFile("Uestim", iter, Uestim, t);
        
        /* Poisson */
        initialize_poisson_solver(data, m, n);
        poisson_solver( data,  inv_delta_tx, m,  n, Uestim, Vestim, phi);
        free_poisson_solver( data); 

        evalVelocity(phi,  dt,  deltax,  deltay, U, V, m, n, Uestim, Vestim);
        int i,j;
        for (i = 0; i < n; i++){
            for (j = 0; j< m; j++){
                P->a[i][j] += phi->a[i][j];
            }
        }
        evalTemperature(T, Ht, Htold, LapT, deltax, deltay, m, n, Pr, dt, Gr);

        iter++;
        t += dt;
        MatrixWriteFile("U", iter, U, t);
        MatrixWriteFile("V", iter, V, t);
        MatrixWriteFile("P", iter, P, t);
        MatrixWriteFile("T", iter, T, t);


        

        //print_matrix(U);
        //printf("Zizi Hihihi \n");
        
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