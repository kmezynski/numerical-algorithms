function [TSVD] = truncSingValDecomp(A, b)
    % esential threshold
    tau = 0.1;
    [V, L] = eig(A'*A);
    
    % evaluate eigenpairs of matrix A
    L = sort(diag(sqrt(L)), 'descend'); % order values for Sigma
    S = zeros(length(L));
    % create and fill Sigma
    for i = 1:length(L)
        S(i, i) = L(i);
    end
    
    % evalute Sigma matrix
    [Sm, Sn] = size(S);
    tmp = zeros(Sm, Sn);
    for i = 1:min(Sm, Sn)
        tmp(i, i) = 1./S(i, i);
    end
    
    % evalute V factor
    V = fliplr(V);
    % evalute U factor
    U = A * V * tmp';
    
    TSVD = zeros(length(diag(S)), 1);
    for i = 1:length(diag(S))
        if abs(S(i, i)) > tau
            TSVD = TSVD + (U(:,i)'*b*V(:,i))/S(i,i);
        end
    end
end
