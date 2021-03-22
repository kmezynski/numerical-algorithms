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
A = [0.835 0.667; 0.333 0.266];
b = [0.168; 0.067];
C = [A, b];

R1 = gaussPartialPivot(C)
x1 = backSubstitution(R1)

b = [0.168; 0.066];
C = [A, b];

R1 = gaussPartialPivot(C)
x1 = backSubstitution(R1)

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

%% Task 8: IF matrix is symetric and positive-defiinite apply Cholesky factorization
%A = [1 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
A = [1,2;2,1]; %not positive-definite matrix

if(issymmetric(A))
    if(isPosDef(A))
        [L] = Cholesky_factor(A);
    else
        disp('matrix is not positive-definite')
    end
else
    disp('matrix is not symetric')
end

%% Task 9: A = pascal(100) if(posDef) -> Cholesky factorization. Interpret obtained factor
%Works for pascal up to 19, 
A=pascal(100);

    if(isPosDef(A))
        [L] = Cholesky_factor(A);
    else
        disp('matrix is not positive-definite')
    end

%% Task 13: Create random matrix with uniform distribution, Compare direct methods that solve the system. Different methods and different sizes of matrix
M = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]';
N = 50;
[m,~]=size(M);
MatrixSize =0;
Time = 0;
ResultsG = zeros(m,4);
ResultsGPP = zeros(m,4);
ResultsGCP = zeros(m,4);
ResultsGJ = zeros(m,4);
ResultsBIM = zeros(m,4);
ResultsLU = zeros(m,4);
ResultsCH = zeros(m,4);
rng('default');
for (i=1:m)

x = rand(M(i)*N,1);
C = rand(M(i),M(i));
A = kron(eye(N), C'*C);
b = A*x;
Xest = zeros(size(x));


disp('gauss');
%gaussNoPivot
tic;
Xe = gaussNoPivot([A,b]);
Xest = backSubstitution(Xe);
Time = toc;

M(i)
e = norm(Xest-x)/norm(x);
%e=(abs( Xest-x))./(abs(x));
%sum(e)/500

ResultsG(i,1) = M(i); % M value
ResultsG(i,2) = M(i)*N;  % Matrix size
ResultsG(i,3) = e;    % error
ResultsG(i,4) = Time; % Time



disp('gaussPP');
%gaussPartialPivot
tic;
Xe = gaussPartialPivot([A,b]);
Xest = backSubstitution(Xe);
Time = toc;

M(i)
e = norm(Xest-x)/norm(x);
ResultsGPP(i,1) = M(i); % M value
ResultsGPP(i,2) = M(i)*N;  % Matrix size
ResultsGPP(i,3) = e;    % error
ResultsGPP(i,4) = Time; % Time

disp('gaussCP');
%GaussCompletePivot
tic;
[R2, Q] = gaussCompletePivot([A,b]);
y = backSubstitution(R2);
Xest = Q * y;
Time = toc;

M(i)
e = norm(Xest-x)/norm(x);
ResultsGCP(i,1) = M(i); % M value
ResultsGCP(i,2) = M(i)*N;  % Matrix size
ResultsGCP(i,3) = e;    % error
ResultsGCP(i,4) = Time; % Time

disp('gaussJordan');
%Gauss Jordan
tic;

RR=gaussJordan([A,b]);
Xest = RR(:,501);
Time = toc;

M(i)
e = norm(Xest-x)/norm(x);
ResultsGJ(i,1) = M(i); % M value
ResultsGJ(i,2) = M(i)*N;  % Matrix size
ResultsGJ(i,3) = e;    % error
ResultsGJ(i,4) = Time; % Time

disp('BuiltIn');
%BuiltIn
tic;
Xest=A\b;
Time = toc;

M(i)
e = norm(Xest-x)/norm(x);
ResultsBIM(i,1) = M(i); % M value
ResultsBIM(i,2) = M(i)*N;  % Matrix size
ResultsBIM(i,3) = e;    % error
ResultsBIM(i,4) = Time; % Time



%LU
disp('LU');
tic;

[L, U] = LUFactor([A,b]);
y = forwardSubstitution([L, b]);
Xest = backSubstitution([U, y]);

Time = toc;

M(i)
e = norm(Xest-x)/norm(x);
ResultsLU(i,1) = M(i); % M value
ResultsLU(i,2) = M(i)*N;  % Matrix size
ResultsLU(i,3) = e;    % error
ResultsLU(i,4) = Time; % Time



%Cholesky
disp('Cholesky');
tic;

[L] = Cholesky_factor(A);
y = forwardSubstitution([L, b]);
Xest = backSubstitution([L', y]);

Time = toc;

M(i)
e = norm(Xest-x)/norm(x);
ResultsCH(i,1) = M(i); % M value
ResultsCH(i,2) = M(i)*N;  % Matrix size
ResultsCH(i,3) = e;    % error
ResultsCH(i,4) = Time; % Time

end