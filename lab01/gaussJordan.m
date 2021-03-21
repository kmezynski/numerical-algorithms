% Gauss Jordan elimination
%
% Input:
%   A  - mxn augmented matrix
% Output:
%   RR - mxn reduced row-echelon form matrix

function [RR] = gaussJordan(A)
    [m, ~] = size(A);
    
    for k = 1:m
        % make k'th row pivot equal 1
        A(k, :) = A(k, :) / A(k, k);
        current_row = A(k, :);
        
        for i = 1:m
            if i ~= k
                A(i,:) = A(i, :) - A(i, k) * current_row;
            end
        end
    end
    RR = A;
end