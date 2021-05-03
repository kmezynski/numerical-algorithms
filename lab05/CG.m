function [x, i, res_err, sol_err] = CG(A, b, x_exact)
    [~, n] = size(A);
    x = zeros(n, 1);
    i = 0;
    th = 0.0001; % initial threshold
    r = b - A*x;
    e = ones(n, 1);
    z = 0;
    beta = 1;
    possible = false;
    
    if(issymmetric(A) && all(eig(A)) > 0)   
        possible = true;
    end
    
    if possible
        while(true & i<1e3)
            x_prev = x;
            i = i + 1;
            rr = r'*r;
            beta = (rr)/(e'*e+eps);
            z = beta*z + r;
            Az = A*z;
            alpha = (rr)/(z'*A*z+eps);
            x = x + alpha*z;
            e = r;
            r = r - alpha*A*z;
            res_err(i) = norm(b - A*x)/norm(b);
            sol_err(i) = norm(x - x_exact)/norm(x);
            
            if (norm(x - x_prev) < th)
                break;
            end
        end
    else
        i = 1;
        res_err = inf;
        sol_err = inf;
        disp('The input matrix is not symmetric or not positive-definite');
    end
end
