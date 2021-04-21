%% Task 1: Find the solution that best approximates
    %the system of inconsistent linear equations
% (a)
A = [3 -1; 1 2; 2 1];                   b = [4; 0; 1];
% (b)
% A = [3 1 1; 2 3 -1; 2 -1 1; 3 -3 3];    b = [6; 1; 0; 8];
% (c)
% A = [1 1 -1; 2 -1 6; -1 4 1; 3 2 -1];   b = [5; 1; 0; 6];

A \ b
inv(A'*A)*A'*b

%% Task 2:  Find the least squares approximating function
    %of the form a_0 +a_1*x^2+a_2*sin((pi*x)/2)
clear all;
close all;
%(a)
x = [0 1 1 -1];
y = [3 0 -1 2];
%(b)
% x = [-1 0 2 3];
% y = [0.5 1 5 9];

A = zeros(length(y), 3);
for i = 1:length(y)
    A(i,:) = [ 1 x(i)^2 sin(x(i)*pi/2) ];
end
b = y';

res_bi = A \ b
% classical LS fitting
res_LSfit = inv(A'*A)*A'*b

figure(1);
plot(x, y, 'ro');
hold on;

x_plot = min(x):0.01:max(x);
y_plot = res_LSfit(1) + res_LSfit(2)*x_plot.^2 + res_LSfit(3)*sin(x_plot*pi/2);
plot(x_plot, y_plot, 'k');

%% Task 3: Find the best fit to the data in the table
    %with an equation in the form: a_0 +a_1*x_1 +a_2*x_2+a_3*x_3
clear all;
close all;
x1 = [50 40 35 40 30];
x2 = [18 20 14 12 16];
x3 = [10 16 10 12 14];
y = [28 30 21 23 23];

A = zeros(length(y), 4);
for i = 1:length(y)
    A(i,:) = [ 1 x1(i) x2(i) x3(i) ];
end
b = y';

res_bi = A \ b
% classical LS fitting
res_LSfit = inv(A'*A)*A'*b

%% Task 4: Using least squares find the "best" straight-line fit 
    %and the error estimates for the slope and intercept of that line
clear all;
close all;
x = [1 2 3 4 5 6 7 8];
y = [1.5 2 2.8 4.1 4.9 6.3 5.0 11.5];

A = zeros(length(y), 2);
for i = 1:length(y)
    A(i,:) = [ 1 x(i) ];
end
b = y';

res_bi = A \ b
% classical LS fitting
res_LSfit = inv(A'*A)*A'*b 

figure(1);
plot(x, y, 'ro');
hold on;

x_plot = 0:0.01:max(x);
y_plot = res_LSfit(1) + res_LSfit(2)*x_plot;
error = sqrt(sum((y - mean(y)).^2) / (length(y)-2)) / sqrt(sum((x-mean(x)).^2))
plot(x_plot, y_plot, 'k');

%% Task 5: Predict how far down range the missile will land
clear all;
close all;
x = [0 250 500 750 1000];
y = [0 8 15 19 20];

A = zeros(length(y), 3);
for i = 1:length(y)
    A(i,:) = [ 1 x(i) x(i).^2 ];
end
b = y';

res_bi = A \ b
% classical LS fitting
res_LSfit = inv(A'*A)*A'*b 

figure(1);
plot(x, y, 'ro');
hold on;

x_plot = min(x):0.01:(max(x)+1050);
y_plot = res_LSfit(1) + res_LSfit(2)*x_plot + res_LSfit(3)*x_plot.^2;
plot(x_plot, y_plot, 'k');

%% Task 6: Using least squares techniques, fit the following data
    %with a line y=a_0 + a_1*x and then fit the data with a quadratic
    %y =a_0+a_1*x+a_2*x^2. Determine which cuvre best fits data by
    %computing l_2 norm of errors.
clear all;
close all;
x = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
y = [2 7 9 12 13 14 14 13 10 8 4];

A = zeros(length(y), 2);
A2 = zeros(length(y), 3);
for i = 1:length(y)
    A(i,:) = [ 1 x(i) ];
    A2(i,:) = [ 1 x(i) x(i).^2 ];
end
b = y';

res_bi = A \ b
res_bi2 = A2 \ b
% classical LS fitting
res_LSfit = inv(A'*A)*A'*b
res_LSfit2 = inv(A2'*A2)*A2'*b

figure(1);
plot(x, y, 'ro');
hold on;

x_plot = min(x):0.01:max(x);
y_plot = res_LSfit(1) + res_LSfit(2)*x_plot;
error = zeros(length(y), 1);

y_plot2 = res_LSfit2(1) + res_LSfit2(2)*x_plot + res_LSfit2(3)*x_plot.^2;
error2 = zeros(length(y), 1);
for i = 1:length(y)
    error(i) = y(i) - (res_LSfit(1) + res_LSfit(2)*x(i));
    error2(i) = y(i) - (res_LSfit2(1) + res_LSfit2(2)*x(i) + res_LSfit2(3)*x(i)^2);
end
error = norm(error)
error2 = norm(error2)
plot(x_plot, y_plot, 'k', x_plot, y_plot2, 'b');

%% Task 7: Fit the polynomial g(x)to the function f(x) to minimalize the error
    %the error in the range [0, pi]
clear all;
close all;
x = 0:0.01:pi;
n = 10;
f = pi^2 - x.^2;

for i = 1:n
    A = [ cos(0*x') cos((1:i).*x') ];
    c = inv(A'*A)*A'*f';
    g = A*c;

    error(i) = norm(f' - g);
    
    figure(i);
    plot(x, f, 'k');
    hold on;
    plot(x, g, 'r');
    title(['n=', num2str(i)]);
    legend('f(x)=\pi^2-x^2','g(x)=\Sigma c_i cos(ix)');
end
error

%% Task 8: Determine the coeficients alpha_i of the cubic polynomial p(x)
    % that is as close to f(x) as possible
 x = 0:0.01:1;
 f = sin(pi*x);
 
 A = [x' (x.^2)' (x.^3)'];
 a = inv(A'*A)*A'*f';
 p = A*a;

 % estimate error
 error = trapz(abs((f'-p).^2));

 figure(1);
 plot(x, f, 'k');
 hold on;
 plot(x, p, 'r');
 legend('f(x)=sin(\pi x)', 'p(x)=\Sigma a_i x^i');

%% Task 9: For the matrix:
    %A) Find the solution Ax=b that has minimum Euclidean norm
    %B) estimate rank(A),cond(A), A^+
    %C) compute orthagonal projectors onto each of the four fundamental
        %subspaces associated with A
    %D) find the point in N(A)^⊥ that is closest to b
clear all;
close all;
A = [-4 -2 -4 -2; 2 -2 2 1; -4 1 -4 -2];
b = [-12; 3; -9];
cond(A) 
null(A)

format short
[U, sigma, V] = svd(A)

% LS solution using SVD
res_LSSVD = singValDecomp(A, b)
% LS solution using TSVD
res_LSTSVD = truncSingValDecomp(A, b)

% Iterative Tikhonov Regularization
[x_LS, error] = iterTikhReg(A, b)

% four subspaces
r = rank(A);
[m, n] = size(A);
u1 = U(:, 1:r);     % columnspace
u2 = U(:, r+1:m);   % nullspace for c-space
v1 = V(:, 1:r);     % rowspace
v2 = V(:, r+1:n);   % nullspace for r-space

% orthagonal projectors
P_r = u1*u1'
P_rh = v1*v1'
P_nh = u2*u2'
P_n = v2*v2'

%% Task 10: For matrix A and exact solution x_*:
    %A) determine the data vector for the model Ax=b, estimate the minimal
        %norm LS solution 
    %B) estimate rank(A), cond(A), A^+
    %C) find N(A) and the missing components that went onto N(A)
    %D) estimate solution error and residual error
A = [ 1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
xeq = [1; 2; 3];

% LS fitting with tikhonov regularization
b = A*xeq
[m,n] = size(A);
format long
if(m>n && rank(A) == n)
    x = inv(A'*A)*A'*b;
else 
    alpha = 0.0000001;
    display("Thikonov")
    x = inv(A'*A + alpha*eye(n)'*eye(n))*A'*b;
end

rank(A)
cond(A)
inv(A'*A) * A'
A' * inv(A*A') 
norm(x-xeq)
norm(b-A*x)
norm(A)

% Iterative Tikhonov Regularization
[x_LS, error] = iterTikhReg(A, b)

% LS solution using SVD
res_LSSVD = singValDecomp(A, b)

format short

%% Task 11: Fit the polynomial(^6) to the generated data in LS sense
    %compare the estimated coefficients to the true ones
    %Perturb the system with Gaussian noise
    % Ilustrate the fitting error (Euclidean norm) vs standard deviation
    %Estimate the maximal value of σ for which the fitting error 
    %exceeds 10^-6 assuming double precision arithmetic
clear all;
close all;

format long
x = 1:14;
a = [40; 10; 5; 3; 2; 1; 1];
b = (40 + 10*x + 5*x.^2 + 3*x.^3 + 2*x.^4 + x.^5 + x.^6)';

A = zeros(length(x), 7);
for i = 1:length(x)
    A(i,:) = [ 1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5 x(i)^6 ];
end

res_LSfit = inv(A'*A)*A'*b

sig = zeros(1, 200);        %creating array of standard deviation
euc_norm = zeros(1, 200);   %creating array of standard euclidean norm
sig(1) = 0.00000006;
b2 = b + randn(length(b), 1)*sig(1);
euc_norm(1) = norm(b - b2);
for i=1:199
    sig(i+1) = sig(i) + 0.000000001;
    b2 = b + randn(length(b), 1)*sig(i+1);  %perturbing b with the noise
    euc_norm(i+1) = norm(b - b2);   %calculating the fitting error
end
figure(1);
plot(sig, euc_norm, 'k');        %plotting standard deviation vs the error
hold on;
yline(10^(-6), 'r');

%% Task 12: For n=10,10^2,10^3,...:
    %A) estimate minimal norm LS solution
    %B) apply the iterative refinement to refine an approximate solution
    %C) estimate rank(A), cond(A),A^+
    %D) estimate solution error and residual error
    %E) draw decay of singular values of A and ilustrate the Picard
    %conditions

%% Task 13: Perform forward projection: A*x^*=b and:
    %A) Solve the system of linear equations with regularization methods
        %estimate the minimal norm LS solution. Which method gives best
        %approximation
    %B) draw the errors
    %C) estimate rank(A), cond(A)
    %D) draw a decay of singular values of A and ilustrate the Picard
        %conditions
    %E) find N(A) and the missing components that went onto N(A)
