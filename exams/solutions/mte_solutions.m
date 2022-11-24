% Clearing terminal & memory, measuring exec time.
clc;
clear;
tic;

% Switching off fsolve text
options = optimset('Display','off');

% Q1) Determine the values of lambda ...
% Sol)
% Beta - 0.1 - 1.1
% Zeta - 3000 - 4000

beta = 0.1 : 0.1 : 1.1;
zeta = 3000: 100 : 4000;

% (i)
% Take beta = 0.5, zeta = 3500
b = 0.5;
z = 3500;
lo = 1;
lamda = fsolve(@(l) lamdaR(l, b, z), lo, options);

disp("Q1");
disp("The value of lambda for beta = " + b + ...
    " and zeta = " + z + " is " + lamda);
disp(" ");

% (ii) & (iii)
betaR = beta(1 : 5);
zetaR = zeta(1 : 5);

lval = zeros(length(zetaR), length(beta));
for i = 1 : length(zetaR)
    for j = 1 : length(beta)
        lval(i, j) = fsolve(@(l) lamdaR(l, beta(j), zetaR(i)), lo, options);
    end
end

fig1 = figure();
plot(beta, lval);

lval2 = zeros(length(betaR), length(zeta));
for i = 1 : length(betaR)
    for j = 1 : length(zeta)
        lval2(i, j) = fsolve(@(l) lamdaR(l, betaR(i), zeta(j)), lo, options);
    end
end

fig2 = figure();
plot(beta, lval2);

%  Q2) A mixture of ....
% Sol2) 
% (ii)
% First solve for F3, F4, F5, F6. We'll be making coefficient matrix
% A and value matrix B to calculate the Feed flow rate matrix F

disp("Q2");

A = [ 0.07, 0.18, 0.15, 0.24;
      0.04, 0.24, 0.10, 0.65;
      0.54, 0.42, 0.54, 0.10;
      0.35, 0.16, 0.21, 0.01
    ];
B = [10.5; 17.5; 28.0; 14.0];

F = A\B;

disp(F);

function lvalue = lamdaR(l, b, z)
    lvalue = 1/sqrt(l) + 4*log(b + 5/(z * sqrt(l))) - 2;
end