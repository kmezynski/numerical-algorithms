A=hilb(10);
x=normrnd(0,1,[10,1]);
b = A*x;
k=norm(A)*norm(inv(A));

X_estimated=A\b;
E_estimated = norm(X_estimated-x)/norm(x);

A(5,5) = A(5,5)+0.001; %Small change of matrix A

X_perturbed=A\b;
E_perturbed = norm(X_perturbed-x)/norm(x);