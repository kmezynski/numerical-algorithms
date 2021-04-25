function [x, y] = NNLS(A, b)
    [~, n] = size(A);
    x0 = ones(n, 1);
    
    AtA = A'*A;
    Atb = A'*b;
    e = ones(n, 1);
    
    y0 = x0;
    k = 0;
    tolerance1 = n*eps;
    tolerance2 = n*sqrt(eps);
    xk = x0;
    yk = y0;
    loop = 0;
    while loop < 1000
        mk =(xk'*yk) ./ (n^2);
        Xk = diag(xk);
        Yk = diag(yk);
        Xk_inv = diag(1./xk);
        Xk_sqrt_inv = diag(1./sqrt(xk));
        Yk_sqrt = diag(sqrt(yk));
        Yk_sqrt_inv = diag(1./sqrt(yk));    
        C = [A; Xk_sqrt_inv*Yk_sqrt];
        d = [b-A*Xk*e; Xk_sqrt_inv*Yk_sqrt_inv*mk*e];
        uk = C \ d;
        vk = -Yk*e + Xk_inv*mk*e - Xk_inv*Yk*uk; 
        T1 = min(-xk(uk<0) ./ uk(uk<0));
        T2 = min(-yk(vk<0) ./ vk(vk<0));
        theta = 0.99995 * min(T1, T2);
        if isempty(theta)
            theta = 0;
            break;
        end
        xk = xk + theta*uk;
        yk = yk + theta*vk;
        if(xk'*yk < tolerance1) & (norm(AtA*xk - Atb-yk) < tolerance2)
            break;
        end
        loop = loop+1;
    end
    x = xk;
    y = yk;
end

