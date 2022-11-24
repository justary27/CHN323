% Name: Aryan Ranjan
% Enrollment No: 20112026
% Quiz 2

% Clearing memory and terminal
clear, clc

% Start measuring program runtime.
tic;

% Preparing the data
x = [5; randn(4,1); 0.5];
V = randn(1);
k = 1;
Q = 200;
tol = 1e-8;
f = fn(x,V,Q,k);
iter = 0;

% Newton's method
while sum(f<tol) ~= size(f,1)
    X = [x(2:5); V];
    j = jacobian(x,V,Q,k);
    X_new = X - (j\f);
    x = [x(1); X_new(1:4); x(6)];
    V = X_new(5);
    f = fn(x,V,Q,k);
    iter = iter + 1;
end

% Printing to 4th decimal
fprintf("The volume is %.4f\n",V);

% Printing the total no. of iterations.
disp("No. of iterations = " + iter);

% Stop measuring runtime and printing it.
timetaken = toc;
disp("The runtime is " + timetaken + " sec.");

% Function to calculate the eqn matrix.
function y = fn(x,V,Q,k)
    y = zeros(5,1);

    for i=2:6
        y(i-1) = (Q * x(i-1)) - (Q * x(i)) - (V* k *x(i)^2);
    end

end

% Function to calculate the jacobian of eqn.
function j = jacobian(x,V,Q,k)

    j = zeros(5,5);

    for i=1:5
        if i ~= 1
            j(i,i-1) = Q;
        end
        if i~=5
            j(i,i) = -(Q + (V * k * x(i)));
        end
        j(i,5) = -(k * x(i)^2);
    end

end