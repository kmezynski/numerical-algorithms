function [x] = SVD(A, b)
    [m, n] = size(A);
    S = zeros(m, n);        % prepare Sigma for further diagonal infill
    [V, L] = eig(A'*A);     % obtain eigenpairs of A
    L = sort(diag(sqrt(L)), 'descend'); % order values of future S matrix
    
    % create Sigma
    for i = 1:min(m, n) 
        S(i, i) = L(i);
    end
    
    [Sr, Sc] = size(S);
    tmp = zeros(Sr, Sc);
    for i = 1:min(Sr, Sc)
        tmp(i, i) = 1./S(i, i);
    end
    V = fliplr(V);  % edit V factor to fit
    U = A*V*tmp';   % compute U factor 
    
    r = rank(A);
    x = zeros(n, 1);
    for i = 1:r
        x = x + (U(:, i)' * b * V(:, i) / S(i, i));
    end
end
