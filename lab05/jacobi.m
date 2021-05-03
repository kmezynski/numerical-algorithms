function  [x, i, res_err, sol_err]= jacobi(A, b, x_exact)
    [~, n] = size(A);
    x = zeros(n, 1);
    i = 0;
    th = 0.0001; % initial threshold
    
    % evaluate spectral radius
    s = diag(A);
    S = diag(s);
    T = -(A - S);
    G = inv(S) * T;
    sr = max(abs(eigs(G)));
    
    if sr < 1
        while(true & i < 1e6)
            i = i+1;
            x_prev = x;
            x = (T*x + b)./s;
            res_err(i) = norm(b - A*x) / norm(b);
            sol_err(i) = norm(x - x_exact) / norm(x);
            
            if(norm(x - x_prev) < th)
                break;
            end
        end
    else 
        i = 1;
        res_err = inf;
        sol_err = inf;
        disp('Spectral radius >= 1 for Jacobi method');
    end
end

