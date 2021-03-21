% Cholesky factorization
%
% Input:
%   A - mxm square matrix
% Output:
%   L - mxm lower triangular matrix


function [L] = Cholesky_factor(A)
    [rows, ~] = size(A);
    L = zeros(rows);
    
    
    
    for k = 1:rows
        L(k,k) = sqrt(A(k,k) - sum (L(k, 1:(k-1)).^2));
        for j=(k+1):rows
            L(j,k) = (1/L(k,k)) * (A(j,k) - sum((L(j, 1:k-1).*L(k, 1:k-1))));
        end
    end
end