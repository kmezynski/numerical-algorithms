function [LS] = LSFit(A, b)
    LS = inv(A'*A) * A' * b;
end