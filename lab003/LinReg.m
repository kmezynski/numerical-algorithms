function [alpha,beta]=LinReg(X,Y)
beta = (Y*X' - size(Y,2)*mean(Y)*mean(X)) / (X*X' -size(Y,2)*mean(X)^2);
alpha = mean(Y) - beta*mean(X);
end