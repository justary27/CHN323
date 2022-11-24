% Indian Institute of Technology Roorkee
% Autumn Semester 2022-23
% CHN 323 Assignment 5

% Aryan Ranjan
% 20112026

% Clearing terminal & memory, measuring exec time.
clc;
clear;
tic;

% Q3) Minimize the cost of waste-water treatment ...

% Sol3) We have the below relations:
%   Wi = xiPi                 in mg/day
%   ci = (1-xi)Pi/Qij           in mg/L, i, j = (1,3) & (2,3)
%   c3 = [R1(1-x1)P1 + R2(1-x2)P2 + (1-x3)P3/(Q34)] in mg/L
%   c4 = {R3[R1(1-x1)P1 + R2(1-x2)P2 + (1-x3)P3] + (1-x4)P4}/Q45  in mg/L
%   di = ci * Qij

%   Now we have to minimize the total cost.

% Given Constant Values
P = [1e9, 2e9, 4e9, 2.5e9];
Q = [1e7, 5e7, 11e7, 25e7];
d = [2e-6, 2e-6, 4e-6, 4e-6];
R = [0.5, 0.35, 0.6, 0];

% Initial Guess for Fractions Removed at Each City.
X0 = [0, 0, 0, 0];

% Storing the Waste Fractions Removed & the Minimum Cost.
[X, minTotalCost] = fmincon(@(X)totalCost(X, P, d), X0, [], [], [], [], zeros(1, 4), [1, 1, 1, 1], @(X)streamConc(X, P, R, Q));

disp("So, the minimum cost is $" + minTotalCost + "." );
disp("It occurs at waste fractions removed of " ...
    + X(1) + ", " + X(2)+ ", " + X(3)+ ", " + X(4) ...
    + " respectively for each city 1 to 4.");

timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.")

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

% % % % % % % % % % Functions Used % % % % % % % % % % % % % % % 
function [conc, conceq] = streamConc(X, P, R, Q)
    conc = zeros(1, 4);
    minConC = 20;
    
%   MinConC is being subtracted to satisy the 
%   inequality conc(X)<0 for fmincon

    conc(1) = (((1 - X(1)) .* P(1)) / Q(1)) - minConC;
    conc(2) = (((1 - X(2)) .* P(2)) / Q(2)) - minConC;
    conc(3) = (((R(1) .* (conc(1) + minConC) .*  Q(1)) + (R(2) .* (conc(2) + minConC) .* Q(2)) + ((1-X(3)) .* P(3))) / Q(3)) - minConC;
    conc(4) = (((R(3) .* (conc(3) + minConC) .* Q(3))  + ((1 - X(4)) .* P(4))) / Q(4)) - minConC;
    
    conceq = [];
end

function cost = totalCost(X, P, d)
    costVec = zeros(1, 4);
    
    costVec(1) = X(1) .* P(1) .* d(1);
    costVec(2) = X(2) .* P(2) .* d(2);
    costVec(3) = X(3) .* P(3) .* d(3);
    costVec(4) = X(4) .* P(4) .* d(4);

    cost = mean(costVec) * 4;
end