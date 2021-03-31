function [min_X, min_L] = invPowAlg(A)
    [m, ~] = size(A);
    
    min_X = rand(m, 1); %initial eigenvector
    prev_min_X = min_X;
    err = 1;            %initial error value
    alpha = 0.1;
    
    while err > 0.000001
        % update min value
        v = (A-alpha * eye(m)) \ min_X;
        min_X = v / norm(v);
        
        % compute difference between actual and previous
        err = abs(norm(min_X) - norm(prev_min_X));  
        prev_min_X = min_X;
    end
    
    % compute the lowest eigenvalue
    min_L = (min_X'*A*min_X)/(min_X'*min_X);
end