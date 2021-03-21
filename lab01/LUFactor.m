% LU factorization using Crout's algorithm
%
% Input:
%   A - mxm square matrix
% Output:
%   L - mxm lower triangular matrix
%   U - mxm upper triangular matrix

function [L, U] = LUFactor(A)
    [m, ~] = size(A);
    L = zeros(m);
    U = eye(m);
    
    % perform Crout's algorithm
    for k = 1:m       
        for j = k:m
            U(k, j) = A(k, j) - dot(L(k, 1:(k-1)), U(1:(k-1), j));
        end
        
        for i = k:m
            L(i, k) = (A(i, k) - dot(L(i, 1:(k-1)), U(1:(k-1), k))) / U(k, k);
        end
    end
end