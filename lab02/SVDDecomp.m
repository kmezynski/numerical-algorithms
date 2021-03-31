function [S, U, X] = SVDDecomp(A)
    [m, n] = size(A);
    S = zeros(m, n);        % prepare Sigma for further diagonal infill
    [X, L] = eig(A'*A);     % obtain eigenpairs of A
    L = sort(diag(sqrt(L)), 'descend'); % order values of future S matrix
    
    % create Sigma
    for i = 1:min(m, n) 
        S(i,i) = L(i);
    end
    
    [Sr, Sc] = size(S);
    tmp = zeros(Sr, Sc);
    for i = 1:min(Sr, Sc)
        tmp(i,i) = 1./S(i,i);
    end
    X = fliplr(X);  % edit V factor to fit
    U = A*X*tmp';   % compute U factor 
end