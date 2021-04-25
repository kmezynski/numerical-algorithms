function [X] = M_FOCUSS(A, B, p)
    [m, n] = size(A);
    l = 1e-8;
    X = rand(n,size(B, 2));
    
    running = true;
    th = 0.001;     % threshold
    i = 0;          % control iterations num
    
    while running
        i = i+1;
        prev_X = X;
        
        w = sqrt(sum(X.^2, 2));
        W = diag(w.^(1 - p/2));
        A = A * W;
        X = W * A'*inv(A*A' + l*eye(m)) * B;
        
        if (norm(X - prev_X) < th || i == 100)
            running = false;
        end
    end
end

