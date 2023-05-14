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
} //checked


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

void initDiagnost(const char *baseResultName)
{
    const char *basename = "results/%s.txt";
    char filename[100];
    sprintf(filename,basename,baseResultName);
}

void Diagnostwritefile(const char *filename, int iter, double value)
{
    FILE* file = fopen(filename,"w");
    if (file == NULL) {
        printf("Error : cannot create result file : did you create the output directory ? \n");
        exit(0); }
    fprintf(file,"%f;%d \n",value,iter);
}

void CloseFile(const char* filename)
{
    FILE* file = fopen(filename);
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
    // for (j = 1; j<m+1; j++){
    //     V->a[0][j] = -1.0/5.0 *(V->a[3][j] - 5*V->a[2][j] + 15*V->a[1][j]);
    //     V->a[n+1][j] = -1.0/5.0 *(V->a[n+1-3][j] - 5*V->a[n+1-2][j] + 15*V->a[n+1-1][j]);
    // }
    // for (i = 1; i<n+1; i++){
    //     U->a[i][0] = -1.0/5.0 *(U->a[i][3] - 5*U->a[i][2] + 15*U->a[i][1]);
    //     U->a[i][m+1] = -1.0/5.0 *(U->a[i][m+1-3] - 5*U->a[i][m+1-2] + 15*U->a[i][m+1-1]);
    // }
    //calcul of the convective term away from the side of the domain 
    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){ // decalle les indices des V de 1 vers la droite 
            // H_x->a[i][j] = 1.0/(4.0*deltax) * (U->a[i+1][j]*U->a[i+1][j] - U->a[i-1][j]*U->a[i-1][j]) 
            //                 + 1.0/(4*deltay) * ((V->a[i][j+1] + V->a[i-1][j+1]) * (U->a[i][j+1] - U->a[i][j])
            //                                   + (V->a[i-1][j] + V->a[i-1][j-1]) * (U->a[i][j] - U->a[i][j-1]));
            H_x->a[i][j] = 1.0/(4.0*deltax) * (U->a[i+1][j]*U->a[i+1][j] - U->a[i-1][j]*U->a[i-1][j])  //ligne pas encore verif 
                            + 1.0/(4.0*deltay) * ((V->a[i][j] + V->a[i+1][j]) * (U->a[i][j+1] - U->a[i][j])
                                              + (V->a[i][j-1] + V->a[i+1][j-1]) * (U->a[i][j] - U->a[i][j-1]));
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){ // decalle les indice des U de 1 vers le le haute
            H_y->a[i][j] = 1.0/(4.0*deltay) * ((U->a[i-1][j+1] + U->a[i-1][j]) * (V->a[i][j] - V->a[i-1][j]) 
                                             + (U->a[i][j+1] + U->a[i][j]) * (V->a[i+1][j] - V->a[i][j]))
                            + 1.0/(4.0*deltax) * (V->a[i][j+1]*V->a[i][j+1] - V->a[i][j-1]*V->a[i][j-1]);
        }
    }
    return 0;
}

int convective_velocity_div(matrix *V, matrix *U, matrix *H_x, matrix *H_y, double deltax, double deltay, int m, int n){
    int i,j;
    //calcul of the convective term away from the side of the domain 
    // for (j = 1; j<m+1; j++){
    //     V->a[0][j] = -1.0/5.0 *(V->a[3][j] - 5*V->a[2][j] + 15*V->a[1][j]);
    //     V->a[n+1][j] = -1.0/5.0 *(V->a[n+1-3][j] - 5*V->a[n+1-2][j] + 15*V->a[n+1-1][j]);
    // }
    // for (i = 1; i<n+1; i++){
    //     U->a[i][0] = -1.0/5.0 *(U->a[i][3] - 5*U->a[i][2] + 15*U->a[i][1]);
    //     U->a[i][m+1] = -1.0/5.0 *(U->a[i][m+1-3] - 5*U->a[i][m+1-2] + 15*U->a[i][m+1-1]);
    // }
    for (i = 1; i < n; i++){
        for (j = 1; j< m+1; j++ ){                                 
            // H_x->a[i][j] = 1.0/(4.0*deltax) * (U->a[i+1][j]*U->a[i+1][j] - 2*U->a[i+1][j]*U->a[i][j] - U->a[i-1][j]*U->a[i-1][j] - 2*U->a[i][j]*U->a[i-1][j])  
            //             + 1.0/(4*deltay) * (U->a[i][j+1]*V->a[i+1][j]) + U->a[i][j+1]*V->a[i][j-1] - U->a[i][j-1]*V->a[i+1][j-1] - U->a[i][j-1]*V->a[i][j-1];
            H_x->a[i][j] =  1.0/(4.0*deltax) * ((U->a[i+1][j] + U->a[i][j]) * (U->a[i+1][j] + U->a[i][j]) - (U->a[i][j] + U->a[i-1][j]) * (U->a[i][j] + U->a[i-1][j]))
                          + 1.0/(4.0*deltay) * ((U->a[i][j+1] + U->a[i][j]) * (V->a[i+1][j] + V->a[i][j]) - (U->a[i][j] + U->a[i][j-1]) * (V->a[i+1][j-1] + V->a[i][j-1]));

        }
    }
    for (i = 1; i < n+1; i++){
         for (j = 1; j< m; j++){ // decalle les indice des U de 1 vers le le haute
            H_y->a[i][j] = 1.0/(4.0*deltax) * (((U->a[i][j+1] + U->a[i][j]) * (V->a[i+1][j] + V->a[i][j])) - ((U->a[i-1][j+1] + U->a[i-1][j])*(V->a[i][j] + V->a[i-1][j])))
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
        for (int j=1;j<m+1; j++ ){
            grad_Px->a[i][j] = 1.0/deltax*(P->a[i][j-1] - P->a[i-1][j-1]); 
        }
    }
    // for (int j=0;j<P->m; j++ ){
    //         grad_Px->a[n][j] = 1/deltax ( - P->a[n][j]) 
    //     }
    for (int i=1; i<P->n+1; i++){
        for (int j=1;j<m; j++ ){
            grad_Py->a[i][j] = 1.0/deltay*(P->a[i-1][j] - P->a[i-1][j-1]); 
        }
        // grad_Px->a[i][m] = 1/deltax (pressbords - P->a[i][m])
    }
    return 0; 
} // checked 

int laplacian_velocity(matrix *LapU, matrix *LapV, matrix *U, matrix *V, double deltax, double deltay, int m, int n){
    for (int i=1; i<n; i++){
        for (int j=1;j<m+1; j++){
            LapU->a[i][j] = 1.0/(deltax*deltax)*((U->a[i+1][j] - 2.0*U->a[i][j] + U->a[i-1][j])) 
                          + 1.0/(deltay*deltay)*(U->a[i][j+1] - 2.0*U->a[i][j] + U->a[i][j-1]); //A verifier on a compute ca trop vite 
        }
    }
    for (int i=1; i<n+1; i++){
        for (int j=1;j<m; j++ ){
            LapV->a[i][j] = 1.0/(deltax*deltax)*((V->a[i+1][j] - 2.0*V->a[i][j] + V->a[i-1][j])) 
                          + 1.0/(deltay*deltay)*(V->a[i][j+1] - 2.0*V->a[i][j] + V->a[i][j-1]); //A verifier on a compute ca trop vite
            } 
        }
    return 0;
} // vite fait checked

int laplacian_temperature(matrix *LapT, matrix *T, double deltax, double deltay, int m, int n){
    for (int i=1; i<n+1; i++){
        for (int j=1;j<m+1; j++){
            LapT->a[i][j] = 1.0/(deltax*deltax)*((T->a[i+1][j] - 2.0*T->a[i][j] + T->a[i-1][j])) 
                      + 1.0/(deltay*deltay)*(T->a[i][j+1] - 2.0*T->a[i][j] + T->a[i][j-1]);
        }
    }
    return 0;
} // checked vite fait

int convective_temperature(matrix *Ht, matrix *T, matrix *U, matrix *V, double deltax, double deltay, int m, int n){
    for (int i=1; i<n+1; i++){
        for (int j=1;j<m+1; j++){
            Ht->a[i][j] = (1.0/(2*deltax)) * ((U->a[i-1][j] * (T->a[i][j]-T->a[i-1][j])) + (U->a[i][j]*(T->a[i+1][j]-T->a[i][j]))) 
                        + (1.0/(2*deltay)) * ((V->a[i][j-1] *(T->a[i][j]-T->a[i][j-1])) + (V->a[i][j]*(T->a[i][j+1]-T->a[i][j])));
           
        }
    }
    return 0;
} // Checked, on des kog 



// int evalEstimVelocity(matrix *U, matrix *V, matrix *Hy, matrix *Hx, matrix *LapU, matrix *LapV,
//                       matrix *Hyold, matrix *Hxold, matrix *grad_Px, matrix *grad_Py, matrix *T, matrix *Uestim, 
//                       matrix *Vestim, int m, int n, double dt, double Gr){
//     double temp; 
//     int i,j;
    
//     for (i = 1; i < n; i++){
//         for (j = 1; j< m+1; j++ ){
//             temp = -1.0/2.0*(3.0*Hx->a[i][j] - Hxold->a[i][j]) - grad_Px->a[i][j] + pow(Gr,-1.0/2.0)*LapU->a[i][j];
//             Uestim->a[i][j] = temp*dt + U->a[i][j];
//         }
//     }
//     for (i = 1; i < n+1; i++){
//         for (j = 1; j< m; j++){
//             temp = -1.0/2.0*(3.0*Hy->a[i][j] - Hyold->a[i][j]) - grad_Py->a[i][j] + pow(Gr,-1.0/2.0)*LapV->a[i][j] + (T->a[i][j]+T->a[i][j+1])/2.0;
//             Vestim->a[i][j] = temp*dt + V->a[i][j];
//             // Vestim->a[i][j] = V->a[i][j] + dt*(T->a[i][j]+T->a[i][j+1])/2.0;
//         }
//     }
    
//     return 0;
// }   // c'est ici qu'on suspecte la merde d'etre 

int evalEstimVelocity(matrix *U, matrix *V, matrix *Hy, matrix *Hx, matrix *LapU, matrix *LapV,
                      matrix *Hyold, matrix *Hxold, matrix *grad_Px, matrix *grad_Py, matrix *T, matrix *Uestim, 
                      matrix *Vestim, int m, int n, double dt, double Gr){
    double temp; 
    int i,j;

    for (j = 1; j< m+1; j++ ){
        for (i = 1; i < n; i++){
            temp = -0.5*(3.0*Hx->a[i][j] - Hxold->a[i][j]) - grad_Px->a[i][j] + pow(Gr,-0.5)*LapU->a[i][j];
            Uestim->a[i][j] = temp*dt + U->a[i][j];
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            temp = -0.5*(3.0*Hy->a[i][j] - Hyold->a[i][j]) - grad_Py->a[i][j] + pow(Gr,-0.5)*LapV->a[i][j] + (T->a[i][j]+T->a[i][j+1])/2.0;
            Vestim->a[i][j] = temp*dt + V->a[i][j];
        }
        Uestim->a[i][0] = -0.2 *(Uestim->a[i][3] - 5.0*Uestim->a[i][2] + 15.0*Uestim->a[i][1]);
        Uestim->a[i][m+1] = Uestim->a[i][m];  //-0.2 *(Uestim->a[i][m+1-3] - 5.0*Uestim->a[i][m+1-2] + 15.0*Uestim->a[i][m+1-1]);
    }

    for (j = 1; j< m+1; j++ ){
        Vestim->a[0][j] = -0.2 *(Vestim->a[3][j] - 5.0*Vestim->a[2][j] + 15.0*Vestim->a[1][j]);
        Vestim->a[n+1][j] = -0.2 *(Vestim->a[n+1-3][j] - 5.0*Vestim->a[n+1-2][j] + 15.0*Vestim->a[n+1-1][j]);
    }

    return 0;
}   //probleme au niveau de T (faire la moyenne pour avoir T au point souhaité)

int evalVelocity(matrix *phi, double dt, double deltax, double deltay, matrix *U, matrix *V, int m, int n, matrix *Uestim, matrix *Vestim){
    int i,j;
    double temp;
    for (i = 1; i < n; i++){
        for (j = 1; j < m+1; j++ ){
            temp = -(phi->a[i][j-1]-phi->a[i-1][j-1])/deltax;
            U->a[i][j] = temp*dt + Uestim->a[i][j];
        }
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            temp = -(phi->a[i-1][j]-phi->a[i-1][j-1])/deltay;
            V->a[i][j] = temp*dt + Vestim->a[i][j];
        }
    }
    // for (j = 1; j<m+1; j++){
    //      V->a[0][j] =  -1.0/5.0 *(V->a[3][j] - 5.0*V->a[2][j] + 15.0*V->a[1][j]);
    //      V->a[n+1][j] = -1.0/5.0 *(V->a[n+1-3][j] - 5.0*V->a[n+1-2][j] + 15.0*V->a[n+1-1][j]);
    // }
    // for (i = 1; i<n+1; i++){
    //      U->a[i][0] = -1.0/5.0 *(U->a[i][3] - 5.0*U->a[i][2] + 15.0*U->a[i][1]);
    //      U->a[i][m+1] =  -1.0/5.0 *(U->a[i][m+1-3] - 5*U->a[i][m+1-2] + 15*U->a[i][m+1-1]);
    //  }
    return 0;
}  //checked 

int evalTemperature(matrix *T, matrix *Ht, matrix *Htold, matrix *LapT, double deltax, double deltay, int m, int n, double Pr, double dt, double Gr){   
    double temp;
    int i,j;
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m+1; j++){
            temp = -1.0/2.0*(3*Ht->a[i][j] - Htold->a[i][j]) + 1.0/(pow(Gr,0.5)*Pr)*LapT->a[i][j] ;
            T->a[i][j] += temp*dt ;
        }
    }
    /*Plus sur de l'ordre des index i et j à checker absolument!!!*/
    double l0 = 1e-3*m*deltay ;
    for (int j=1;j<m+1; j++){
        //Left Surface : 
        T->a[0][j] = T->a[1][j];
        //Right Surface : 
        T->a[n+1][j] = T->a[n][j];
        }    
    double Tgamma;
    for (int i=1;i<n+1; i++){
        //Free Surface : 
        //T->a[i][m+1] = ((2*l0-deltay)/(2*l0+deltay))*T->a[i][m];
        Tgamma = ((4.0*l0-deltay)/(4.0*l0+deltay))*T->a[i][m];
        T->a[i][m+1] = -1.0/5.0 * (T->a[i][m-2]-5.0*T->a[i][m-1]+15.0*T->a[i][m]-16.0*Tgamma);

        //Bottom Surface : 
        T->a[i][0] = T->a[i][1]+deltay;
        // Tgamma = deltay/2 + T->a[i][1];
        // T->a[i][0] = -1.0/5 * (T->a[i][3]-5*T->a[i][2]+15*T->a[i][1]-16*Tgamma);
    }
    return 0;
} //potentiel erreur au niveau des conditions limites mais a priori ok 

int evalVelocity_norm( matrix* U, matrix* V, matrix* NormVelocity, int m, int n){
    double v,u;
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            v = (V->a[i+1][j] + V->a[i+1][j+1])/2.0;
            u = (U->a[i][j+1] + U->a[i+1][j+1])/2.0;
            NormVelocity->a[i][j] = sqrt(v*v + u*u);
        }
    }
    return 0;
}

int eval_vorticity(matrix* U, matrix* V, matrix* Vorticity, double deltax, double deltay, int m, int n){
    for (int i=0; i<n+1; i++){
        for(int j=0; j<m+1; j++){
            Vorticity->a[i][j] = (V->a[i+1][j] - V->a[i][j])/deltax + (U->a[i][j+1] - U->a[i][j])/deltay;
        }
    }   
    return 0;
}


int evalEstimVelocityMixer(matrix *U, matrix *V, matrix *Hy, matrix *Hx, matrix *LapU, matrix *LapV,
                      matrix *Hyold, matrix *Hxold, matrix *grad_Px, matrix *grad_Py, matrix *T, matrix *Uestim, 
                      matrix *Vestim, int m, int n, double dt, double Gr, double dtau, double L, double deltax, double deltay, int iter){
    double temp; 
    int i,j;
    double xu,yu,xv,yv;
    double xg = L /2.0;
    double yg = xg;
    double D = 3.0/5.0*L;
    double a = D/2 ;
    double halfd = D/10 ; 
    double omega = 0.1;
    double theta;
    double chiU, chiV;
    double norm;
    double Us;
    double Vs;

    for (j = 1; j< m+1; j++ ){
        for (i = 1; i < n; i++){

            xu = i*deltax;
            yu = j*deltay - deltay/2;
            theta = atan2(yu-yg, xu-xg);
            norm = sqrt((xu-xg)*(xu-xg)+(yu-yg)*(yu-yg));
            chiU = (norm<a*cos(3*(theta+iter*omega*dt))) || (norm<halfd) ? 1.0 : 0.0 ;
            Us = omega*(yu-yg);

            temp = -1.0/2.0*(3.0*Hx->a[i][j] - Hxold->a[i][j]) - grad_Px->a[i][j] + pow(Gr,-1.0/2.0)*LapU->a[i][j] + chiU*Us/dtau;
            Uestim->a[i][j] = (temp + U->a[i][j]/dt)/(1.0/dt + chiU/dtau);
        }
        Vestim->a[0][j] = -1.0/5.0 *(Vestim->a[3][j] - 5.0*Vestim->a[2][j] + 15.0*Vestim->a[1][j]);
        Vestim->a[n+1][j] = -1.0/5.0 *(Vestim->a[n+1-3][j] - 5.0*Vestim->a[n+1-2][j] + 15.0*Vestim->a[n+1-1][j]);
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            xv = i*deltax - deltax/2;
            yv = j*deltay;
            theta = atan2(yv-yg, xv-xg);
            norm = sqrt((xv-xg)*(xv-xg)+(yv-yg)*(yv-yg));
            chiV = (norm<a*cos(3*(theta+iter*omega*dt))) || (norm<halfd) ? 1.0 : 0.0 ;
            Vs = omega*(xv-xg) ;

            temp = -1.0/2.0*(3.0*Hy->a[i][j] - Hyold->a[i][j]) - grad_Py->a[i][j] + pow(Gr,-1.0/2.0)*LapV->a[i][j] + (T->a[i][j]+T->a[i][j+1])/2.0 + chiV*Vs/dtau ;
            Vestim->a[i][j] = (temp + V->a[i][j]/dt)/(1.0/dt + chiV/dtau);
        }
        Uestim->a[i][0] = -1.0/5 *(Uestim->a[i][3] - 5.0*Uestim->a[i][2] + 15.0*Uestim->a[i][1]);
        Uestim->a[i][m+1] =  Uestim->a[i][m]; //-1.0/5.0 *(U->a[i][m+1-3] - 5*U->a[i][m+1-2] + 15*U->a[i][m+1-1]-16*Uestim->a[i][m]);
    }

    return 0;
}

int evalTemperatureMixer(matrix *T, matrix *Ht, matrix *Htold, matrix *LapT, double deltax, double deltay, int m, int n, double Pr, double dt, double Gr, double dtau, double L, int iter){   

    double temp; 
    int i,j;
    double xT,yT;
    double xg = L /2.0;
    double yg = xg;
    double D = 3.0/5.0*L;
    double a = D/2 ;
    double halfd = D/10 ;
    double omega = 0.1; 
    double theta;
    double chiT;
    double norm;
    double Ts = 0;
    double nelem = 0;

    for (i = 1; i < n+1; i++){
        for (j = 1; j< n+1; j++){
            xT = i*deltax - deltax/2;
            yT = j*deltay - deltay/2;
            norm = sqrt((xT-xg)*(xT-xg)+(yT-yg)*(yT-yg));
            theta = atan2(xT-xg, xT-xg);
            nelem += norm<(D/2) ? 1.0 : 0.0 ;
            Ts += norm<(D/2) ? T->a[i][j] : 0.0 ;
        }
    }
    Ts /= nelem; 
    Diagnostwritefile("avTmixer", iter, Ts );
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m+1; j++){
            xT = i*deltax - deltax/2;
            yT = j*deltay - deltay/2;
            norm = sqrt((xT-xg)*(xT-xg)+(yT-yg)*(yT-yg));
            theta = atan2(yT-yg, xT-xg);
            chiT = (norm<a*cos(3*(theta+iter*omega*dt))) || (norm<halfd) ? 1.0 : 0.0;
            temp = -1.0/2.0*(3.0*Ht->a[i][j] - Htold->a[i][j]) +  pow(Gr,-1.0/2.0)*1.0/Pr*LapT->a[i][j] + chiT*Ts/dtau;
            T->a[i][j] = (temp + T->a[i][j]/dt)/(1.0/dt + chiT/dtau);
        }
    } 

    /*Plus sur de l'ordre des index i et j à checker absolument!!!*/
    double l0 = 1e-3*m*deltay ;
    for (int j=1;j<m+1; j++){
        //Left Surface : 
        T->a[0][j] = T->a[1][j]; //-1.0/5 * (T->a[3][j]-5*T->a[2][j]-T->a[1][j]); 
        //Right Surface : 
        T->a[n+1][j] = T->a[n][j]; //-1.0/5 * (T->a[n-2][j]-5*T->a[n-1][j]-T->a[n][j]);  
    }    
    double Tgamma;
    for (int i=1;i<n+1; i++){
        //Free Surface : 
        //T->a[i][m+1] = ((2*l0-deltay)/(2*l0+deltay))*T->a[i][m];
        Tgamma = ((4*l0-deltay)/(4*l0+deltay))*T->a[i][m];
        T->a[i][m+1] = -1.0/5 * (T->a[i][m-2]-5*T->a[i][m-1]+15*T->a[i][m]-16*Tgamma);

        //Bottom Surface : 
        T->a[i][0] = T->a[i][1]+deltay;
        // Tgamma = deltay/2 + T->a[i][1];
        // T->a[i][0] = -1.0/5 * (T->a[i][3]-5*T->a[i][2]+15*T->a[i][1]-16*Tgamma);
    }
    return 0;
} //potentiel erreur au niveau des conditions limites mais a priori ok 

int Diagnostic(matrix *T, int m, int n, double L, double H, int iter){
    //double Tcyl ; Calculé pour l'hélice
    double Tav = 0;
    double Trms = 0;
    double Omegaf = L*H;
    double qeav = 0;

    /*Fluid spatially-averaged temperature*/

    for (int i = 1; i < n+1; i++){
        qeav = (T->a[i][m+1]-T->a[i][m])/L;
        for (int j = 1; j< m+1; j++){
            Tav += T->a[i][j]/Omegaf;
            Trms += (T->a[i][j]-Tav)*(T->a[i][j]-Tav);
        }
    }Trms = sqrt(Trms/Omegaf);

    Diagnostwritefile("avFlux", iter, qeav);
    Diagnostwritefile("avTemp", iter, Tav);
    Diagnostwritefile("RmsTemp", iter, Trms);

    return 0;
}



int main(int argc, char *argv[]){
    int Sim=0;


    if(Sim==0){
        double H = 1.0;//1.5*L;
        double L = H/1.5;

        int m = 75 ; 
        int n = (int) (2.0/3*m);
    
        double deltax = L/n; //to specify 
        double deltay = H/m; //to specify 

        // sprintf("deltax = %.6e \t; deltay = %.6e \n", deltax, deltay);
        // printf("n = %d \t; m = %d \n", n, m);

        double SimTime = 10; //to specify
        double t = 0.0;
        // double dt = 1.0/300.0;
        double dt= 1.0/30;  //(60.0*1.5);//1.0/75.0; //to specify
        double dtau = 1e-3*dt;
        double Pr = 4.0;
        double Gr = 1e10;
        int iter = 0;
            
        int i,j;

        /*WRITE YOUR PROJECT ...*/
        matrix *U = callocate_matrix(m+2, n+1); 
        matrix *V = callocate_matrix(m+1, n+2);
        matrix *Vorticity = callocate_matrix(m+1, n+1);
        matrix *NormVelocity = callocate_matrix(m,n);
        matrix *LapU = callocate_matrix(m+2, n+1);
        matrix *LapV = callocate_matrix(m+1, n+2);
        matrix *Vestim = callocate_matrix(m+1, n+2);
        matrix *Uestim = callocate_matrix(m+2, n+1);
        matrix *grad_Px = callocate_matrix(m+2, n+1);
        matrix *grad_Py = callocate_matrix(m+1, n+2);
        matrix *Hx = callocate_matrix(m+2, n+1); 
        matrix *Hxold = callocate_matrix(m+2, n+1);
        matrix *Hy = callocate_matrix(m+1, n+2);
        matrix *Hyold = callocate_matrix(m+1, n+2);
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
        MatrixWriteFile("Vestim", iter, Vestim, t);
        MatrixWriteFile("P", iter, P, t);
        MatrixWriteFile("phi", iter, phi, t);
        MatrixWriteFile("Velocity", iter, NormVelocity, t);
        MatrixWriteFile("Vorticity", iter, NormVelocity, t);
        double l0 = 1e-3*H ;
        for (int i=1;i<n+1; i++){
        //Free Surface : 
        T->a[i][m+1] = ((2*l0-deltay)/(2*l0+deltay))*T->a[i][m];
        //double Tgamma = ((4*l0-deltay)/(4*l0+deltay))*T->a[i][m];
        //T->a[i][m+1] = -1.0/5.0 * (T->a[i][m-2]-5*T->a[i][m-1]+15*T->a[i][m]-16*Tgamma);

        //Bottom Surface : 
        T->a[i][0] = T->a[i][1]+deltay;
        // Tgamma = deltay/2 + T->a[i][1];
        // T->a[i][0] = -1.0/5 * (T->a[i][3]-5*T->a[i][2]+15*T->a[i][1]-16*Tgamma);
        }
        MatrixWriteFile("T", iter, T, t);

        initDiagnost("avFlux");
        initDiagnost("avTemp");
        initDiagnost("avTmixer");
        initDiagnost("RmsTemp");


        PetscInitialize(&argc, &argv, 0, 0);    
        initialize_poisson_solver(data, m, n);

        while (t<SimTime){

            copy_matrix(Hxold, Hx);
            copy_matrix(Hyold, Hy);
            copy_matrix(Htold, Ht);
        
            //convective_velocity_div(V, U, Hx, Hy, deltax, deltay, m, n);
            convective_velocity(V, U, Hx, Hy, deltax, deltay, m, n);
            convective_temperature(Ht, T, U, V, deltax, deltay, m, n);
            laplacian_velocity(LapU, LapV, U, V, deltax, deltay, m, n); 
            laplacian_temperature(LapT, T, deltax, deltay, m, n); 
            pressure_gradure(grad_Px, grad_Py, P, deltax, deltay, m, n);

            evalEstimVelocityMixer(U, V, Hy, Hx, LapU, LapV,
                        Hyold, Hxold, grad_Px, grad_Py, T, Uestim, Vestim, m, n, dt, Gr, dtau, L, deltax, deltay, iter);

            // evalEstimVelocity(U, V, Hy, Hx, LapU, LapV, Hyold, Hxold, grad_Px, grad_Py,
            //                 T, Uestim, Vestim, m, n, dt, Gr); 
            
            
            /* Poisson */
            
            poisson_solver(data, inv_delta_tx, m,  n, Uestim, Vestim, phi);
            
            
            // MatrixWriteFile("phi", iter, phi, t);
            evalVelocity(phi,  dt,  deltax,  deltay, U, V, m, n, Uestim, Vestim);

            for (i = 0; i < n; i++){
                for (j = 0; j< m; j++){
                    P->a[i][j] += phi->a[i][j];
                }
            }

            //evalTemperatureMixer(T, Ht, Htold, LapT, deltax, deltay, m, n, Pr, dt, Gr, dtau, L, iter);
            evalTemperature(T, Ht, Htold, LapT, deltax, deltay, m, n, Pr, dt, Gr);
            evalVelocity_norm(U,V, NormVelocity, m, n);
            eval_vorticity(U, V, Vorticity, deltax, deltay, m, n);

            iter++;
            t += dt;

            int niter = 30;
            if (iter%niter==0){
                MatrixWriteFile("U", iter/niter, U, t);
                MatrixWriteFile("P", iter/niter, P, t);
                MatrixWriteFile("T", iter/niter, T, t);
                MatrixWriteFile("V", iter/niter, V, t);
                MatrixWriteFile("Vestim", iter/niter, Vestim, t);
                MatrixWriteFile("Velocity", iter/niter, NormVelocity, t);
                MatrixWriteFile("Gradpy", iter/niter, grad_Py, t);
                MatrixWriteFile("LapV", iter/niter, LapV, t);
                MatrixWriteFile("Hy", iter/niter, Hy, t);
                MatrixWriteFile("phi", iter/niter, phi, t);
                MatrixWriteFile("Vorticity", iter/niter, phi, t);
                Diagnostic(T, m, n, L, H, iter/niter);
            }
            // MatrixWriteFile("U", iter, V, t);
            // MatrixWriteFile("P", iter, P, t);
            // MatrixWriteFile("T", iter, T, t);
            // MatrixWriteFile("V", iter, V, t);
            // MatrixWriteFile("Vestim", iter, Vestim, t);
            // MatrixWriteFile("Velocity", iter, NormVelocity, t);
            // MatrixWriteFile("Gradpy", iter, grad_Py, t);
            // MatrixWriteFile("LapV", iter, LapV, t);
            // MatrixWriteFile("Hy", iter, Hy, t);
        }

        CloseFile("avFlux");
        CloseFile("avTemp");
        CloseFile("avTmixer");
        CloseFile("RmsTemp");

        free_poisson_solver( data); 
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
    
    else{
        // int n = 3;
        // int m = 3;
        // matrix *phi = callocate_matrix(m, n); 
        // matrix *U = callocate_matrix(m+2, n+1);
        // matrix *Uold = callocate_matrix(m+2, n+1);
        // matrix *V = callocate_matrix(m+1, n+2); 
        // // matrix *Hx = callocate_matrix(m+2, n+1); 
        // // matrix *Hxold = callocate_matrix(m+2, n+1);
        // // matrix *Hyold = callocate_matrix(m+1, n+2);
        // matrix *Hy = callocate_matrix(m+1, n+2);
        // //matrix *gradV = callocate_matrix(m, n); 

        // matrix *phi_div = callocate_matrix(m, n); 
        // matrix *U_div = callocate_matrix(m+2, n+1);
        // matrix *V_div = callocate_matrix(m+1, n+2); 
        // matrix *Hx_div = callocate_matrix(m+2, n+1); 
        // // matrix *Hxold = callocate_matrix(m+2, n+1);
        // // matrix *Hyold = callocate_matrix(m+1, n+2);
        // matrix *Hy_div = callocate_matrix(m+1, n+2);
        // int i,j;
        // for (i = 1; i < n; i++){
        //     for (j = 1; j< m+1; j++){
        //         U->a[i][j] = i*(m+1)+j+1;
        //         U_div->a[i][j] = i*(m+1)+j+1;
        //     }
        // }

        // // for (i = 1; i < m+1; i++){
        // //     for (j = 1; j< n; j++){
        // //         V->a[i][j] = i-j;
        // //         V_div->a[i][j] = i-j;
        // //     }
        // // }
        // printf("this is matrix U \n");
        // print_matrix(U);
        // copy_matrix(Uold,U);
        // U->a[1][1]= 8.5;
        // printf("this is matrix U \n");
        // print_matrix(U);
        // printf("this is matrix Uold \n");
        // print_matrix(Uold);

        // // convective_velocity(V, U, Hx, Hy, deltax, deltay, m, n);
        // // convective_velocity_div(V_div, U_div, Hx_div, Hy_div, deltax, deltay, m, n);
        // // printf("This is div convective velocity Hy\n");
        // // print_matrix(Hy_div);
        // // printf("This is convective velocity Hy\n");
        // // print_matrix(Hy);
        // // printf("This is div convective velocity Hx\n");
        // // print_matrix(Hx_div);
        // // printf("This is convective velocity Hx\n");
        // // print_matrix(Hx);
        // // printf("This is U_div\n");
        // // print_matrix(U_div);
        // // printf("This is U\n");
        // // print_matrix(U);
        // // printf("This is V_div\n");
        // // print_matrix(V_div);
        // // printf("This is V\n");
        // // print_matrix(V);
        // // MatrixWriteFile("U", 0, U, 0);
        // // MatrixWriteFile("V", 0, V, 0);
        // // PetscInitialize(&argc, &argv, 0, 0);
        // // Poisson_data *data = (Poisson_data *)malloc(sizeof(Poisson_data));
        // // initialize_poisson_solver(data, m, n);
        // // test_poisson(data, dxdt, m, n, U, V, phi);
        // // free_poisson_solver( data);
        // free_matrix(phi); 
        // free_matrix(U); 
        // free_matrix(V); 
    }
}
