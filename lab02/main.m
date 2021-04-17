%% Task 1: Compute eigenpairs for the matrices
A1 = [1 0 0; 2 1 0; 0 0 3];
A2 = [0 -2 1; 1 3 -1; 0 0 1];
A3 = [4 1 0; 1 4 1; 0 1 4];
A4 = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];

x4 = 1;
V = [(-41+9*sqrt(41))/82*x4 3*sqrt(41)/41*x4 (41+3*sqrt(41))/82*x4 x4];
V/norm(V)

[X, L] = eig(A1)
% matrix A4 is singular (det(A4) = 0)
%% Task 2: Compute the smallest and the largest eigenvalue
A = [4 2 0 0; 1 4 1 0; 0 1 4 1; 0 0 2 4];
[max_X, max_L] = scalPowAlg(A)
[min_X, min_L] = invPowAlg(A)

[X, L] = eig(A)
%% Task 3: Solve the differential equation
%% Task 4: Find A^100 by diagonalizing matrix A
A = [4 3; 1 2];
[X, L] = eig(A);

A100 = X*L^(100)*inv(X)
A100simple = A^(100)
%% Task 5: Prove that matrix A isn't diagonalizable
A = [2 1 1; 2 1 -2; -1 0 -2];

[X, L] = eig(A);
A_test = X' * L * X;
% A != A_test thus matrix A isn't diagonalizable
%% Task 6: Compute eigenpairs for matrix A
% matrix must be symmetric!
A = [1 2 3 4; 1 2 2 3; 0 2 3 2; 0 0 3 4];

[ref_X, ref_L] = eig(A);
% basic QR iteration algorithm 
[X1, L1] = basicIterQR(A);
% single-shift QR iteration algorithm 
[X2, L2] = shiftIterQR(A);
%% Task 7: Draw the Gershgorincs disk, determine eigenvalues for the matrices
A = [-2 -1 0; 2 0 0; 0 0 2];
%A = [5 1 1; 0 6 1; 0 0 -5];
%A = [5.2 0.6 2.2; 0.6 6.4 0.5; 2.2 0.5 4.7];

cols = sum(abs(A));
rows = sum(abs(A), 2);
eig_vals = eig(A);
t = 0:pi / 100:2*pi;

% setup plot
grid on
hold on
axis equal
for i = 1:size(A,1)
    d = A(i,i);
    cols(i) = cols(i) - abs(d);
    rows(i) = rows(i) - abs(d);
    plot(real(eig_vals(i)), imag(eig_vals(i)), 'X', 'DisplayName', [num2str(i), '.eigenvalue']);
    plot(d + rows(i)*cos(t), rows(i)*sin(t), 'LineWidth',1, 'DisplayName', ['Gershgorins disc for ', num2str(i) '. row']);
    plot(d + cols(i)*cos(t), cols(i)*sin(t), 'LineWidth',3, 'DisplayName', ['Gershgorins disc for ', num2str(i) '. col'])
end
legend
%% Task 9: Perform the Schur decomposition to the matrices
A = [1 i; -i 1];
%A = [1 0 1+i; 0 2 0; 1-i 0 0];

[X, L] = eig(A)
[Q, R] = qr(X)              % applying QR factorization
U = inv(Q)*A*Q              % computing U
[U_ref, T_ref] = schur(A)   % comparing the results
%% Task 10: Compute the SVD of the matrices
%A = [1 0; 0 1; 1 0];
%A = [-1 2 2];
A = [2 2 2 2; 17/10 1/10 -17/10 -1/10; 3/5 9/5 -3/5 -9/5];

[S, U, X] = SVDDecomp(A);
A_check = U*S*X';
%% Task 11: Compute the SVD of 640x400 image matrix for various approximations
in_img = imread('task11_input.jpg');
A = double(im2gray(in_img)) / 255;

[S, U, X] = SVDDecomp(A);
[m, n] = size(A);
S2 = zeros(m, n);
k = 1;
for n = [10 20 40 60] %number of singular values to approximate with
    for i = 1:n
        S2(i, i) = S(i, i);   %creating new sigma matrix
    end
    out_img = U*S2*X';         % approximate image with new sigma
    figure(k);
    imshow(out_img);
    title(['Output for ', num2str(n), ' singular values', ]);
    k = k + 1;
end