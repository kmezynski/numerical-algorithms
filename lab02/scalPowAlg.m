function [max_X, max_L] = scalPowAlg(A)
    [m, ~] = size(A);

    max_X = rand(m, 1); %initial eigenvector
    prev_max_X = max_X;
    err = 1;            %initial error value
    
    while err > 0.000001
        % update max value
        max_X = (A*max_X) / norm(A*max_X);
        
        % compute difference between actual and previous
        err = abs(norm(max_X) - norm(prev_max_X)) / norm(prev_max_X);  
        prev_max_X = max_X;
    end
    
    % compute the largest eigenvalue
    max_L = (max_X'*A*max_X) / (max_X'*max_X);
end