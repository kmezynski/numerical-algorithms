%% Task 1: Solve system of linear equations
A = [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
b = [0; 0; 0; 5];
C = [A, b];

X = gaussNoPivot(C)
x = backSubstitution(X)

%% Task 2: Solve system of linear equations using the Gaussian elimination
A = [1 1 1; 1 1 2; 1 2 2];
b = [1; 2; 1];
C = [A, b];

% Gaussian elimination without pivoting
R0 = gaussPartialPivot(C);
x0 = backSubstitution(R0);

% Gaussian elimination with parital pivoting
R1 = gaussPartialPivot(C);
x1 = backSubstitution(R1);

% Gaussian elimination with complete pivoting
[R2, Q] = gaussCompletePivot(C);
y = backSubstitution(R2);
x2 = Q * y;

%% Task 3: Solve system of linear equations
A = [0.0001 1; 1 1];
b = [1; 2];
C = [A, b];



%% Task 4: Solve system of linear equations
A = [0.835 0667; 0.333 0.266];
b = [0.168; 0.067];
C = [A, b];

%% Task 5: Find inverse matrix of A to solve the system AX = I_3
A = [2 1 2; 1 2 3; 4 1 2];
inv_A = inv(A)
I = A * inv_A

%% Task 6: Solve system of linear equations using LU factorization, compute det(A)
A = [1 2 3 4; -1 1 2 1; 0 2 1 3; 0 0 1 1];
b = [1; 1; 1; 1];

[L, U] = LUFactor(A)
y = forwardSubstitution([L, b]);
x = backSubstitution([U, y])

det_A = det(L) * det(U)
det(A)

%% Task 7: Apply LU factorization to the matrix A, compute det(A)
A = hilb(5);

[L, U] = LUFactor(A)
det_A = det(L) * det(U)
det(A)

%% Task 8: 

