function [x] = FOCUSS(A, b, p)
    [m, n] = size(A);
    l = 1e-8;
    x = rand(n, 1);
    
    running = true; 
    th = 0.001;     % threshold
    
    while running
        prev_x = x;
        
        W = diag(abs(x).^(1 - p/2));
        tmp_A = A * W;
        tmp_x = tmp_A' * inv(tmp_A*tmp_A' + l*eye(m)) * b;
        x = W * tmp_x;
        if (norm(x - prev_x) < th)
            running = false;
        end
    end
end

