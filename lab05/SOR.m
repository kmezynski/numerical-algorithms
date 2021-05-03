function [x, i, res_err,sol_err] = SOR(A, b, x_exact)
    [~, n] = size(A);
    x = zeros(n, 1);
    i = 0;
    w = 1.4;        % relaxation parameter
    th = 0.001;     % initial threshold
    
    % evaluate spectral radius
    L = tril(A, -1);
    D = diag(diag(A));
    U = A - L - D;
    S = L + 1/w*D;
    T = -(U + (w - 1)/w*D);
    G = inv(S) * T;
    sr = max(abs(eigs(G)));
    if sr < 1
        div=w./diag(A);
        while(true & i < 1e6)
            i = i+1;
            x_prev = x;
            %x = inv(D + w*L) * (w*b - (w*U+(w-1)*D)*x_prev);
            for j=1:n
                x(j) = (1-w)*x(j) + div(j)*(b(j) - sum(A(j, 1:j-1)*x(1:j-1)) - sum(A(j, j+1:n)*x(j+1:n)));
            end
            res_err(i) = norm(b - A*x) / norm(b);
            sol_err(i) = norm(x - x_exact) / norm(x);
            
            if (norm(x - x_prev) < th)
                break;
            end
        end
    else
        i = 1;
        res_err = inf;
        sol_err = inf;
        disp('Spectral radius >= 1 for SOR method');
    end
end
