% Back substitution algorithm
%
% Input:
%   A - mxn augmented matrix
% Output:
%   x - mx1 column vector of solutions

function [x] = forwardSubstitution(A)
    [m, n] = size(A);
    x = zeros(m, 1);

    for k = 1:m
        coefs_sum = A(k, 1:k) * x(1:k);
        x(k) = (A(k, n) - coefs_sum) / A(k, k);
    end
end