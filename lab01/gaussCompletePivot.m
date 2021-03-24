% Gaussian elimination with complete pivoting
%
% Input:
%   A - mxn augmented matrix
% Output:
%   R - mxn row-echelon form matrix
%   Q - columns permutations matrix

function [R, Q] = gaussCompletePivot(A)
    [m, n] = size(A);
    
    % record perumations of rows (p) and cols (q)
    p = 1:m; q = 1:n;
    
    for k = 1:(m-1)
        % find coords of max value (max_row, max_col) in submatrix of A
        [max_row_vals, row_indices] = max(abs(A(p(k:m), q(k:m))));
        [~, max_col] = max(max_row_vals);
        max_row = row_indices(max_col);
        
        % transform submatrix relative coords to global coords of A
        max_row = max_row + k-1;
        max_col = max_col + k-1;
        
        % permute rows and columns
        p([k max_row]) = p([max_row k]);
        q([k max_col]) = q([max_col k]);
        
        % perform gaussian elimination with the greatest pivot
        if A(k, k) ~= 0
            i = (k+1):m;
            A(p(i), q(k)) = A(p(i), q(k)) / A(p(k), q(k));
            A(p(i), :) = A(p(i), :) - A(p(i), q(k)) * A(p(k), :);
        end
    end
    
    R = A(p, q);
    Q = eye(m); Q = Q(:, q(1:m));
end