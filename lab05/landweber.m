function [x, i, res_err, sol_err] = landweber(A, b, x_exact)
    [m, n] = size(A);
    x = zeros(n, 1);
    i = 0;
    max_alpha = 2 * abs(max(eig(A'*A)))^(-1);
    alpha = max_alpha / 2;  % relaxation parameter
    th = 0.0001;            % initial threshold
    
    % evaluate spectral radius
    S = 1/alpha * eye(m);
    T = S - A'*A;
    G = inv(S) * T;
    % G = eye(m) - alpha*A'*A;
    sr = max(abs(eig(G)));
    
    if sr < 1
        B = alpha*A';
        while(true & i<1e3)
            i = i + 1;
            x_prev = x;
            x = x_prev + B*(b - A*x_prev);
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
        disp('Spectral radius >= 1 for Landweber method');
    end
end
