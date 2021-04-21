function [a, b] = linReg(X, Y)
    b = (Y*X' - size(Y, 2) * mean(Y) * mean(X)) / (X*X' - size(Y, 2) * mean(X)^2);
    a = mean(Y) - b * mean(X);
end