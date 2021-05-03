function  [x, i, res_err, sol_err]= gaussSeidel(A, b, x_exact)
    [~, n] = size(A);
    x = zeros(n, 1);
    i = 0;
    th = 0.0001; % initial threshold
    
    % evaluate spectral radius
    S = tril(A);
    T = -triu(A, 1);
    G = inv(S)*T;
    sr = max(abs(eigs(G)));
    
    if sr < 1
        div = 1./diag(A);
        while(true & i<1e6)
            i = i+1;
            x_prev = x;
            for j = 1:n
                 x(j) = div(j)*(b(j) - sum(A(j, 1:j-1)*x(1:j-1)) - sum(A(j, j+1:n)*x_prev(j+1:n)));
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
        disp('Spectral radius >= 1 for Gauss-Seidel method');
    end
end
