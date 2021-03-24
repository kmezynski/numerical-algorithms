% Gaussian elimination with partial pivoting
%
% Input:
%   A - mxn augmented matrix
% Output:
%   R - mxn row-echelon form matrix

function [R] = gaussPartialPivot(A)
    [m, ~] = size(A);
    
    % record perumations of rows
    p = 1:m;
    
    for k = 1:(m-1)
        % find max column value in submatrix of A
        [~, max_row] = max(abs(A(p(k:m), k)));
        
        % transform relative row index to global row index of A
        max_row = max_row + k-1;
        
        % permute rows
        p([k max_row]) = p([max_row k]);
        
        % perform gaussian elimination with the greatest pivot in column
        if A(k, k) ~= 0
            i = (k+1):m;
            A(p(i), k) = A(p(i), k) / A(p(k), k);
            A(p(i), :) = A(p(i), :) - A(p(i), k) * A(p(k), :);
        end
    end
    
    R = A(p, :);
end