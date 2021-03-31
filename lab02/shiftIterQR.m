function [X, L] = shiftIterQR(A)
    A_init = A;
    [m, ~] = size(A_init); %as we are dealing with square matrices
    sig = A(m,m);   %initial shift value
    Qr = eye(m);    %matrix prepared for the vector of eigenvectors
    for i = 1:1000
        [Q, R] = qr(A - sig * eye(m));  %iterative QR factoization
        A = R * Q + sig * eye(m);       %computing A matrix for the next iteration
        sig = A(m,m);   %updating shift value
        Qr = Qr * Q;    %computing vector of eigenvectors
    end
    
    % clean matrix with eigenvalues
    L = zeros(m);
    for i = 1:m
        L(i,i) = A(i,i);
    end
    X = Qr;
end