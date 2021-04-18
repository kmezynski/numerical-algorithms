function [x]=LS_fit(A,b)
x = inv(A'*A)*A'*b;
end