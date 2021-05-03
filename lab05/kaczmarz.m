function [x, i, res_err, sol_err] = kaczmarz(A, b, x_exact)
    [m, n] = size(A);
    x = zeros(n, 1);
    i = 0;
    alpha = 1;
    th = 0.0001; % initial threshold
    
    At = A';
    for j = 1:m
       B(j) = alpha / norm(A(j, :))^2;
    end
    while(true & i < 1e6)
        i = i+1;
        x_prev = x;
        for j = 1:m
            x = x+(b(j)-A(j, :)*x)*At(:, j)*B(j);
        end
        res_err(i) = norm(b - A*x) / norm(b);
        sol_err(i) = norm(x - x_exact) / norm(x);
        
        if (norm(x - x_prev) < th)
            break;
        end
    end
end

