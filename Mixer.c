
int evalEstimVelocityMixer(matrix *U, matrix *V, matrix *Hy, matrix *Hx, matrix *LapU, matrix *LapV,
                      matrix *Hyold, matrix *Hxold, matrix *grad_Px, matrix *grad_Py, matrix *T, matrix *Uestim, 
                      matrix *Vestim, int m, int n, double dt, double Gr, double dtau, double L){
    double temp; 
    int i,j;
    double xu,yu,xv,yv;
    double xg = L /2.0;
    double yg = xg;
    double D = 3.0/5.0*L;
    double d = D/5 ; 
    double theta;
    double ChiU, ChiV;

    for (j = 1; j< m+1; j++ ){
        for (i = 1; i < n; i++){

            xu = i*dx;
            yu = j*dy - dy/2;
            theta = atan2(yu-yg, xu-xg);
            norm = sqrt((xu-xg)*(xu-xg)+(yu-yg)*(yu-yg));
            ChiU = (norm<a*cos(3*(theta+omega*iter*dt))) || (norm<d) ? 1.0 : 0.0 ;
            Us->a[i][j] = ChiU ? omega*norm*sin(theta) : 0.0;

            temp = -1.0/2.0*(3.0*Hx->a[i][j] - Hxold->a[i][j]) - grad_Px->a[i][j] + pow(Gr,-1.0/2.0)*LapU->a[i][j] + chiU*Us/dtau;
            Uestim->a[i][j] = (temp + U->a[i][j]/dt)/(1.0/dt + ChiU/dtau);
        }
        Vestim->a[0][j] = -1.0/5.0 *(Vestim->a[3][j] - 5.0*Vestim->a[2][j] + 15.0*Vestim->a[1][j]);
        Vestim->a[n+1][j] = -1.0/5.0 *(Vestim->a[n+1-3][j] - 5.0*Vestim->a[n+1-2][j] + 15.0*Vestim->a[n+1-1][j]);
    }
    for (i = 1; i < n+1; i++){
        for (j = 1; j< m; j++){
            xv = i*dx - dx/2;
            yv = j*dy;
            theta = atan2(yv-yg, xv-xg);
            norm = sqrt((xv-xg)*(xv-xg)+(yv-yg)*(yv-yg));
            chiV = (norm<a*cos(3*(theta+omega*iter/dt))) || (norm<d) ? 1.0 : 0.0 ;
            Vs = (norm<d) ? omega*norm*cos(theta) : 0.0;

            temp = -1.0/2.0*(3.0*Hy->a[i][j] - Hyold->a[i][j]) - grad_Py->a[i][j] + pow(Gr,-1.0/2.0)*LapV->a[i][j] + (T->a[i][j]+T->a[i][j+1])/2.0 + chiV*Vs/dtau ;
            Vestim->a[i][j] = (temp + V->a[i][j]/dt)/(1.0/dt + ChiU/dtau);
        }
        Uestim->a[i][0] = -1.0/5 *(Uestim->a[i][3] - 5.0*Uestim->a[i][2] + 15.0*Uestim->a[i][1]);
        Uestim->a[i][m+1] = Uestim->a[i][m]; //-1.0/5.0 *(U->a[i][m+1-3] - 5*U->a[i][m+1-2] + 15*U->a[i][m+1-1]);
    }
}



int evalTemperatureMixer(matrix *T, matrix *Ht, matrix *Htold, matrix *LapT, double deltax, double deltay, int m, int n, double Pr, double dt, double Gr, double dtau, double L){   

    double temp; 
    int i,j;
    double xu,yu,xv,yv;
    double xg = L /2.0;
    double yg = xg;
    double D = 3.0/5.0*L;
    double d = D/5 ; 
    double theta;
    double ChiU, ChiV;
    double Ts = 0
    double nelem = 0;

    for (i = 1; i < n+1; i++){
        for (j = 1; j< n+1; j++){
            nelem += (norm<a*cos(3*(theta+omega*iter/dt))) ? 1.0 : 0.0 ;
            Ts += (norm<a*cos(3*(theta+omega*iter/dt))) ? T[i][j] : 0.0 ;
        }
    }
    Ts /= nelem; 

    for (i = 1; i < n+1; i++){
        for (j = 1; j< m+1; j++){
            xT = i*dx - dx/2;
            yT = j*dy - dy/2;
            norm = sqrt((xT-xg)*(xT-xg)+(yT-yg)*(yT-yg));
            chiT = (norm<a*cos(3*(theta+omega*iter/dt))) ? 1.0 : 0.0;
            temp = -1.0/2.0*(3.0*Ht->a[i][j] - Htold->a[i][j]) +  pow(Gr,-1.0/2.0)*1.0/Pr*LapT->a[i][j] + chiT*Ts/dtau;
            Testim->a[i][j] = (temp + V->a[i][j]/dt)/(1.0/dt + ChiT/dtau);
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

int Diagnostic(matrix *T, int m, int n, double L, double H){
    //double Tcyl ; Calculé pour l'hélice
    double Tav = 0;
    double Trms = 0;
    double Omegaf = L*H;
    double qeav = 0;

    /*Fluid spatially-averaged temperature*/

    for (i = 1; i < n+1; i++){
        qeav = (T->a[i][m+1]-T->a[i][m])/L;
        for (j = 1; j< m+1; j++){
            Tav += T->a[i][j]/Omegaf;
            Trms += (T->a[i][j]-Tav)*(T->a[i][j]-Tav);
        }
    }Trms = sqrt(Trms/Omegaf);
}