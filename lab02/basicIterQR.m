function [X, L] = basicIterQR(A)
    A_init = A;
    [m, ~] = size(A_init); %as we are dealing with square matrices
    Qr = eye(m);    %matrix prepared for the vector of eigenvectors
    for i = 1:1000
        [Q, R] = qr(A); %iterative QR factoization
        A = R * Q;      %computing A matrix for the next iteration
        Qr = Qr * Q;    %computing vector of eigenvectors
    end
    
    L = zeros(m);
    for i = 1:m         %cleaning up he eigenvalues       
        L(i,i) = A(i,i);
    end
    X = Qr;
end