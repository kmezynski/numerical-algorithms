%% Task 1: Solve system of linear equations
A = [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
b = [0; 0; 0; 5];

X = gaussNoPivot([A, b])
x = backSubstitution(X)

%% Task 2: Solve system of linear equations using the Gaussian elimination
A = [1 1 1; 1 1 2; 1 2 2];
b = [1; 2; 1];

% Gaussian elimination without pivoting
R0 = gaussPartialPivot([A, b]);
x0 = backSubstitution(R0);

% Gaussian elimination with parital pivoting
R1 = gaussPartialPivot([A, b]);
x1 = backSubstitution(R1);

% Gaussian elimination with complete pivoting
[R2, Q] = gaussCompletePivot([A, b]);
y = backSubstitution(R2);
x2 = Q * y;

%% Task 3: Solve system of linear equations
A = [0.0001 1; 1 1];
b = [1; 2];

% Gaussian elimination without pivoting
R = gaussPartialPivot([A, b]);
x = backSubstitution(R)

% Gaussian elimination with complete pivoting
[R, Q] = gaussCompletePivot([A, b]);
y = backSubstitution(R);
x = Q * y

k = norm(A) * norm(inv(A))

%% Task 4: Solve system of linear equations
A = [0.835 0.667; 0.333 0.266];
b = [0.168; 0.067];

R = gaussPartialPivot([A, b]);
x = backSubstitution(R)

b = [0.168; 0.066];

R = gaussPartialPivot([A, b]);
x = backSubstitution(R)

k = norm(A) * norm(inv(A))

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
A = [1 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
%A = [1,2;2,1]; %not positive-definite matrix

if(issymmetric(A))
    if(isPosDef(A))
        [L] = Cholesky_factor(A);
    else
        disp('matrix is not positive-definite')
    end
else
    disp('matrix is not symetric')
end
inverse = (inv(L)')*(inv(L))
%% Task 9: A = pascal(100) if(posDef) -> Cholesky factorization. Interpret obtained factor
%Works for pascal up to 19, 
A=pascal(100);

    %if(isPosDef(A))
        [L] = Cholesky_factor(A);
    %else
    %    disp('matrix is not positive-definite')
    %end

%% Task 10: Transform matrix to RREF, determine rank(A)
A = [1 2 2 3 1;
     2 4 4 6 2;
     3 6 6 9 6;
     1 2 4 5 3];

X = RREF(A)
rank(A)

% x1, x3, x5 are basic variables (rank = 3)
% x2, x4 are free variable

%% Task 11: Compute the LU factorization, determine set of basic variables
    %find homogenous solution Ax=0
A = [1 3 3 2; 2 6 9 5; -1 -3 3 0];
b = [0; 0; 0];

[L, U] = lu(A)
A \ b

rref(A)
rank(A)

% x1, x3 are basic variables (rank = 2)
% x2, x4 are free variable

%% Task 13: Create random matrix with uniform distribution, Compare direct methods that solve the system. Different methods and different sizes of matrix
%Avalible methods: Gauss, GaussPartialPivot, GaussCompletePivot,
%GaussJordan, BuiltIn, LU, Cholesky
% If method is present in vactor Methods, it will be used
Methods = ["Gauss"; "GaussPartialPivot"; "GaussCompletePivot"; "GaussJordan"; "BuiltIn"; "LU"; "Cholesky"]; 
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
    
    if(ismember("Gauss",Methods))
        disp('gauss'); %gaussNoPivot
        disp(M(i));
        tic;
        Xe = gaussNoPivot([A,b]);
        Xest = backSubstitution(Xe);
        Time = toc;
        
        e = norm(Xest-x)/norm(x);
        
        ResultsG(i,1) = M(i); % M value
        ResultsG(i,2) = M(i)*N;  % Matrix size
        ResultsG(i,3) = e;    % error
        ResultsG(i,4) = Time; % Time
    end
    
    if(ismember("GaussPartialPivot",Methods))
        disp('gaussPP'); %gaussPartialPivot
        disp(M(i));
        tic;
        Xe = gaussPartialPivot([A,b]);
        Xest = backSubstitution(Xe);
        Time = toc;
        
        e = norm(Xest-x)/norm(x);
        ResultsGPP(i,1) = M(i); % M value
        ResultsGPP(i,2) = M(i)*N;  % Matrix size
        ResultsGPP(i,3) = e;    % error
        ResultsGPP(i,4) = Time; % Time
    end
    
    if(ismember("GaussCompletePivot",Methods))
        disp('gaussCP'); %GaussCompletePivot
        disp(M(i));
        tic;
        [R2, Q] = gaussCompletePivot([A,b]);
        y = backSubstitution(R2);
        Xest = Q * y;
        Time = toc;
        
        e = norm(Xest-x)/norm(x);
        ResultsGCP(i,1) = M(i); % M value
        ResultsGCP(i,2) = M(i)*N;  % Matrix size
        ResultsGCP(i,3) = e;    % error
        ResultsGCP(i,4) = Time; % Time
    end
    
    if(ismember("GaussJordan",Methods))
        disp('gaussJordan');%Gauss Jordan
        disp(M(i));
        tic;
        RR=gaussJordan([A,b]);
        Xest = RR(:,501);
        Time = toc;
        
        e = norm(Xest-x)/norm(x);
        ResultsGJ(i,1) = M(i); % M value
        ResultsGJ(i,2) = M(i)*N;  % Matrix size
        ResultsGJ(i,3) = e;    % error
        ResultsGJ(i,4) = Time; % Time
    end
    
    if(ismember("BuiltIn",Methods))
        disp('BuiltIn');%BuiltIn
        disp(M(i));
        
        tic;
        Xest=A\b;
        Time = toc;
        
        e = norm(Xest-x)/norm(x);
        ResultsBIM(i,1) = M(i); % M value
        ResultsBIM(i,2) = M(i)*N;  % Matrix size
        ResultsBIM(i,3) = e;    % error
        ResultsBIM(i,4) = Time; % Time
    end
    
    
    if(ismember("LU",Methods))
        disp('LU'); %LU
        disp(M(i));
        tic;
        
        [L, U] = LUFactor([A,b]);
        y = forwardSubstitution([L, b]);
        Xest = backSubstitution([U, y]);
        
        Time = toc;
        
        
        e = norm(Xest-x)/norm(x);
        ResultsLU(i,1) = M(i); % M value
        ResultsLU(i,2) = M(i)*N;  % Matrix size
        ResultsLU(i,3) = e;    % error
        ResultsLU(i,4) = Time; % Time
    end
    
    
    if(ismember("Cholesky",Methods))
        disp('Cholesky'); %Cholesky
        disp(M(i));
        tic;
        
        [L] = Cholesky_factor(A);
        y = forwardSubstitution([L, b]);
        Xest = backSubstitution([L', y]);
        
        Time = toc;
        
        e = norm(Xest-x)/norm(x);
        ResultsCH(i,1) = M(i); % M value
        ResultsCH(i,2) = M(i)*N;  % Matrix size
        ResultsCH(i,3) = e;    % error
        ResultsCH(i,4) = Time; % Time
    end
end

%% Task 14: Solve Ax=b, A=hilb(10), x=normrnd(0,1,[10,1])
A=hilb(10);
x=normrnd(0,1,[10,1]);
b = A*x;
k=norm(A)*norm(inv(A));

X_estimated=A\b;
E_estimated = norm(X_estimated-x)/norm(x);

A(5,5) = A(5,5)+0.001; %Small change of matrix A

X_perturbed=A\b;
E_perturbed = norm(X_perturbed-x)/norm(x);

%% Task 15: Find circuit currents, derive equations using KVL & KCL
A = [0 0 0 0 10 -10; 0 0 10 -10 0 0; 0 0 10 0 10 0; ...
    1 0 -1 -1 0 0; 0 1 1 0 -1 0;  0 1 0 -1 0 1];
b = [10; -10; 20; 0; 0; 0];

A \ b
