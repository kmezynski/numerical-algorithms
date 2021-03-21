% Back substitution algorithm
%
% Input:
%   A - mxn augmented matrix
% Output:
%   x - mx1 column vector of solutions

function [x] = backSubstitution(A)
    [m, n] = size(A);
    x = zeros(m, 1);

    for k = m:-1:1
        coefs_sum = A(k, k:m) * x(k:m);
        x(k) = (A(k, n) - coefs_sum) / A(k, k);
    end
end