import numpy as np

A = np.zeros((9,9))
b = np.zeros(9)
U = np.zeros((4,5))
V = np.zeros((5,4))
m = 3
n = 3


for i in range(m+2):
    for j in range(m+1):
        U[j][i] = i*(m+1)+j+1
        V[i][j] = 20-(j*(n+2)+i)

#print(U, "\n", V)


for i in range(m):
     for j in range(n):
         b[i*m+j] = (U[i+1][j+1]-U[i][j+1]+V[i+1][j+1]-V[i+1][j])

print(b)


b[m+1] = 0

A[0][0] = -2
A[0][1] = 1
A[0][m] = 1

A[m-1][m-1] = -2
A[m-1][m-2] = 1
A[m-1][2*m-1] = 1

A[m*(n-1)][m*(n-1)] = -2
A[m*(n-1)][m*(n-1)+1] = 1
A[m*(n-1)][m*(n-2)] = 1

A[m*n-1][m*n-1] = -2
A[m*n-1][m*n-2] = 1
A[m*n-1][m*(n-1)-1] = 1

for r in range(1,m-1):
    A[r][r] = -3
    A[r][r-1] = 1
    A[r][r+1] = 1
    A[r][r+m] = 1

for r in range(m*(n-1)+1,m*n-1):
    A[r][r] = -3
    A[r][r-1] = 1
    A[r][r+1] = 1
    A[r][r-m] = 1

for r in range(m,m*(n-1),m):
    A[r][r] = -3
    A[r][r+1] = 1
    A[r][r-m] = 1
    A[r][r+m] = 1

for r in range(2*m-1,m*n-1,m):
    A[r][r] = -3
    A[r][r-1] = 1
    A[r][r-m] = 1
    A[r][r+m] = 1

for i in range(1,n-1):
    for j in range(1,m-1):
        r = i*m+j
        A[r][r] = -4
        A[r][r-1] = 1
        A[r][r+1] = 1
        A[r][r+m] = 1
        A[r][r-m] = 1
        
A[m+1][m+1] = 1.0
A[m+1][1] = 0.0
A[m+1][2*m+1] = 0.0   
A[m+1][m+2] = 0.0    
A[m+1][m] = 0.0


# print(A)

# print(np.linalg.det(A))

phi = np.linalg.solve(A,b)

print(phi)