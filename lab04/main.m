%% Task 1: Find a condition on the numbers a, b, c such that the following
%          system of equations is consistent.
A = [1 3 1; -1 -2 1; 3 7 -1];
% (a)
x = [0; 3; 1];
% (b)
%x = [0; 0; 1];
% (c)
%x = [-2; 1; 0];
b = A * x;

x_focuss = FOCUSS(A, b, 1)
focuss_err = norm(b - A*x_focuss) / norm(b);

%% Task 2: Perform the forward projection of the exact solution x onto
%          the range space spanned by the columns in the matrix A
%A = [1 2 2 3 1; 2 4 4 6 2; 3 6 6 9 6; 1 2 4 5 3];
% a21 change
A = [1 2 2 3 1; 0 4 4 6 2; 3 6 6 9 6; 1 2 4 5 3];
x = [1; 0; 1; 1; 0];
b = A * x;

% FOCUSS algorithm
x_focuss = FOCUSS(A, b, 1)
focuss_err = norm(b - A*x_focuss) / norm(b)
focuss_sol_err = norm(x_focuss - x) / norm(x_focuss)

% SVD algorithm
x_svd = SVD(A, b)
svd_err = norm(b - A*x_svd) / norm(b)
svd_sol_err = norm(x_svd - x) / norm(x_svd)

% QR algorithm
x_qr = QR(A, b)
qr_err = norm(b - A*x_qr) / norm(b)
qr_sol_err = norm(x_qr - x) / norm(x_qr)

%% Task 3: Generate 5 sparse instantaneous signals such that at each time
%          sample at most 3 signals are non-zeros. Then perform the forward
%          projection of each time sample onto the range space spanned by
%          the columns of the matrix X
X = [0 1.1006 0 0.8886 0 1.5877 0.2157 2.5855 0 0;
     1.1174 1.5442 2.3505 0 0 0 0 0 0 0.3035; 
     0 0.0859 0 0 1.4193 0.6966 0 0.1873 0.8404 0;
     0.0326 0 0.7481 0 0.2916 0.8351 0.1049 0 0 0.4900;
     0.5525 0 0 0.4882 0.1978 0 0.7223 0 0.1001 0.7394];
 
A = [1 2 2 3 1; 
     0 4 4 6 2; 
     3 6 6 9 6; 
     1 2 4 5 3];
B = A*X;

% FOCUSS algorithm
x_focuss = [];
for i = 1:10 
    x_focuss(:, i) = FOCUSS(A, B(:, i), 1);
end
focuss_err = norm(B - A*x_focuss) / norm(B)
focuss_sol_err = norm(x_focuss - X) / norm(x_focuss)

% M-FOCUSS algorithm
x_mfocuss = MFOCUSS(A, B, 2)
mfocuss_err = norm(B - A*x_mfocuss) / norm(B)
mfocuss_sol_err = norm(x_mfocuss - X) / norm(x_mfocuss)

% plot & compare results
[m, n] = size(X);
x = 1:n;
for i = 1:m
    subplot(m, 1, i)
    plot(x, X(i, :), 'ko', x, x_focuss(i, :), 'rx', x, x_mfocuss(i, :), 'b+')
    legend('Exact solution', 'FOCUSS solution', 'M-FOCUSS solution')
    title(i + ". row values")
end

%% Task 4: Solve the problem, compare two cases, when p = 0 and p = 1, with
%          respect to the residual error
A = [2 3 -1 10 21 44 -9 1 -1; 1 2 2 8 15 35 8 -3 1; 3 1 1 6 16 53 -7 2 2];
b = [118; 77; 129];

x_focuss = FOCUSS(A, b, 0)
focuss_err = norm(b - A*x_focuss) / norm(b)

x_focuss = FOCUSS(A, b, 1)
focuss_err = norm(b - A*x_focuss) / norm(b)

%% Task 5: Solve the problem
A = [1 1 1 1; 1 3 1 1; 1 -1 3 1; 1 1 1 3; 1 1 1 -1];
b = [2; 1; 6; 3; 1;];
B = [1 1 1 -1; 1 -1 1 1; 1 1 -1 1];
d = [1; 3; -1];

% built-in MATLAB function
x_lsqlin = lsqlin(A, b, B, d)
lsqlin_err = norm(b - A*x_lsqlin) / norm(b)

% LS algorithm
x_LS = equalityLS(A, b, B, d)
LS_err = norm(b - A*x_LS) / norm(b)

%% Task 6: Solve the NNLS problem, which methods gives the best NNLS solution?
A = [73 71 52; 87 74 46; 72 2 7; 80 89 71];
b = [49; 67; 68; 20];

% interior-point algorithm
x_nnls = NNLS(A, b)
residual_nnls = norm(b - A*x_nnls)

% built-in MATLAB function
x_lsqnonneg = lsqnonneg(A, b)
lsqnonneg_err = norm(b - A*x_lsqnonneg)

%% Task 7: For the matrix A and b find the NNLS solution to Ax = b and
%          compare it with the ordinary LS solution.
A = [-4 -2 -4 -2; 2 -2 2 1; -4 1 -4 -2];
b = [-12; 3; -9];
% with additive white Gaussian noise of SNR = 20dB
%b = awgn(b, 20);

% built-in MATLAB function
x_lsqnonneg = lsqnonneg(A, b)
lsqnonneg_err = norm(b - A*x_lsqnonneg ) / norm(b)

% interior-point algorithm
x_nnls = NNLS(A, b)
residual_nnls = norm(b - A*x_nnls)

% FOCUSS algorithm
x_focuss = FOCUSS(A, b, 1);
focuss_err = norm(b - A*x_focuss) / norm(b)

% M-FOCUSS algorithm
x_mfocuss = MFOCUSS(A, b, 2)
mfocuss_err = norm(b - A*x_mfocuss) / norm(b)

% SVD algorithm
x_svd = SVD(A, b)
svd_err = norm(b - A*x_svd) / norm(b)

% QR algorithm
x_qr = QR(A, b)
qr_err = norm(b - A*x_qr) / norm(b)

%% Task 8: Applying the Galerkin discretization to the Fredholm integral
%          equation of the first kind, we get the matrix operator
A = []; b = []; s = []; bs = [];

for snr = [0, 10, 20, 30]
    for n = [1e1 1e2 1e3]
        h = 1/n;
        for i = 1:n
            for j = 1:n
                A(i,i) = h^2*(h*(i^2 - i + 1/4) - (i - 2/3));
                A(i,j) = h^2*(j - 1/2)*(h*(i - 1/2) - 1);
                s(j) = sin(4*pi*h*j);
                if s(j) > 0
                    bs(i, j) = A(i, j)*s(j);
                else
                    bs(i, j) = 0;
                end
            end
            b(i) = sum(bs(i, :));
            b = b';
        end

        sprintf('\nSNR: %d, n: %d', snr, n);

        if snr > 0
            b = awgn(b, snr);
        end
        % FOCUSS algorithm
        x_focuss = FOCUSS(A, b, 1);
        focuss_err = norm(b - A*x_focuss) / norm(b)
        
        % SVD algorithm
        x_svd = SVD(A, b)
        svd_err = norm(b - A*x_svd) / norm(b)
        
        % QR algorithm
        x_qr = QR(A, b)
        qr_err = norm(b - A*x_qr) / norm(b)
        
        % interior-point algorithm
        x_nnls = NNLS(A, b)
        residual_nnls = norm(b - A*x_nnls)
    end
end
