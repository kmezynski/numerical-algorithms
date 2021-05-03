function [x, i, res_err, sol_err] = SD(A, b, x_exact)
    [~, n] = size(A);
    x = zeros(n, 1);
    i = 0;
    th = 0.0001; % initial threshold
    
    while(true & i < 1e3)
        i = i+1;
        x_prev = x;
        g = A*x_prev - b;
        u = (g'*g)/(g'*A*g);
        x = x_prev - u*g;
        res_err(i) = norm(b - A*x) / norm(b);
        sol_err(i) = norm(x - x_exact) / norm(x);
        
        if (norm(x - x_prev) < th)
            break;
        end
    end
end

