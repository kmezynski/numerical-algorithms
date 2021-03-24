% Reduced-row echelon form matrix computation
%
% Input:
%   A  - mxn augmented matrix
% Output:
%   RR - mxn reduced row-echelon form matrix

function [RR] = RREF(A)

[m, n] = size(A);
j = 0;

for k = 1:m
    j = j + 1;
    if j > n
        break
    end
    
    % make left-most coefficient 1 (pivot)
    row = A(k, :);
    if row(k) == 0
        j = j+1;
    end
    row = row/row(j);
    A(k, :) = row;
    
    for i = 1 : m
        if i ~= k
            A(i, :) = A(i, :)-(A(i, j)) * row;
        end
    end
    
    for i = k + 1 : m
       A(i:end, k+1:end);       % submatrix of A (in which we are looking for non-zero pivots)
       A(i:end, k+1);           % left-most column
       if ~any(A(i:end, k+1))   % if the left-most column has only zeros check the next one
           k = k+1;
       end
       A(i:end, k+1:end);
       if A(i, k+1) == 0
           non_zero_row = find(A(i:end, k+1), 1);
           if isempty(non_zero_row)
               continue
           end
           A([i, i+non_zero_row-1], :) = deal(A([i+non_zero_row-1, i], :));
       end
    end
end
RR = A;

end