void computeLaplacianMatrix(Mat A, int rowStart, int rowEnd, int Nx, int Ny, double h)
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

    r = rowStart + (Ny-1)/2;
    MatSetValue(A, r, r, 1.0, INSERT_VALUES);
    MatSetValue(A, r, r-1, 0.0, INSERT_VALUES);
    MatSetValue(A, r, r+1, 0.0, INSERT_VALUES);
    MatSetValue(A, r, r+(Ny-1), 0.0, INSERT_VALUES);
}