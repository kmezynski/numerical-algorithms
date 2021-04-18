function [x_LS]=TikhonovRegularization(A,b,x0,err_th,alpha)
[~,col] = size(A);
x_LS = zeros(col,1)
L = eye(col);
err = 1;

while(err>err_th)
    x_old = x_LS; %previous solution
    x_LS = x_LS + inv(A'*A + alpha^2*L'*L)*A'*(b-A*x_LS); %calculating solutin
    err = norm(b-A*x0) - norm(b-A*x_LS); %calculating error
end

%COS nie dziala :(