% Cholesky factorization
%
% Input:
%   A - matrix
% Output:
%   posDef - bool: 1 - matrix is positive-definite , 0 - matrix is not positive-definite


function [posDef] = isPosDef(A)
    [rows, ~] = size(A);
    L = zeros(rows);
    Determinants= zeros(1,rows);
    posDef = 0;
    
    
    for i = 1:rows
            Determinants(i) = det(A(1:i,1:i)); %i upper left determinants
        end

        if (all(Determinants>0)) %chceck if all determinants are positive
            posDef=1;
        end 
end