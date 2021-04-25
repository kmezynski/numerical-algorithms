function x = equalityLS(A, b, B, d)
    [Q, R] = qr(B');
    p = rank(B);
    n = size(Q, 1);
    
    R = R(1:p, 1:p);
    Q1 = Q(:, 1:p);
    Q2 = Q(:, p+1:n);
    
    x1 = Q1 * inv(R') * d;
    y2 = pinv(A*Q2) * (b-A*x1);
    
    x = x1 + Q2*y2; % unique solution
end

