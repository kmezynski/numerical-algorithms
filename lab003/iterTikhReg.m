function [x_LS, error] = iterTikhReg(A, b)
    format long;
    [~, n] = size(A);
    x_LS = zeros(n, 1);
    L = eye(n);
    alpha = 1;
    error = 1;
    
    % evalute an error
    while(error > 0.00000000001)
        x0 = x_LS;
        x_LS = x_LS + inv(A'*A + alpha^2*L'*L) * A' * (b - A*x_LS);
        error = norm(b - A*x0) - norm(b - A*x_LS);
    end
end
