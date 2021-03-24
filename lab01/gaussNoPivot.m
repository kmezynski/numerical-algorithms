% Gaussian elimination without pivoting
%
% Input:
%   A - mxn augmented matrix
% Output:
%   R - mxn row-echelon form matrix

function [R] = gaussNoPivot(A)
    [m, ~] = size(A);
    
    for k = 1:(m-1)
        if A(k, k) ~= 0
            i = k+1:m;
            A(i, k) = A(i, k) / A(k, k);
            A(i, :) = A(i, :) - A(i, k) * A(k, :);
        end
    end
    
    R = A;
end