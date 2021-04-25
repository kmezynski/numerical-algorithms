function [x] = QR(A, b)    
    [Q,R] = qr(A');
    r = rank(A);
    
    z = R(1:r, 1:r)' \ b(1:r);
    x = Q(:, 1:r) * z;
end

