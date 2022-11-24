% Post MTE summary.

% Table of contents:
%  - Parameter estimation.
%  - Optimization techniques.
%  - Interpolation methods.
%  - Solving ODEs for IVPs and BVPs.

%% Parameter Estimation %%

% Steps:
% a) Assume a model for the given measured quantites and result variable
%    as a mathematical relation.
% b) For the given data, taking an initial guess, for model params, 
%    calculate model predicted value of result.
% c) Calculate RMSE for the model predicted and actual result value.
% d) Minimize RMSE till a desired value using fminsearch.

% Example) Consider the data set shown ... (Check Parameter_estimation.pdf)
% Sol) 

% Clearing terminal & memory, measuring exec time.
clc;
clear;
tic;

% Assuming Nu = a.Re^b.Pr^c.Mu^d.

% Nu values from experiment.
Nu_exp = [277, 348, 421, 223, 77, 114.8, 95.9, 68.4, 49.1, 56, 39.9, 47, 94.2, 99.9, 83.1, 35.9];

% Measured variable values from experiment.
Re_exp = [49e3, 68.6e3, 84.8e3, 34.2e3, 22.9e3, 1321, 931, 518, 346, 122.9, 54, 84.6, 1249, 1021, 465, 54.8];
Pr_exp = [2.3, 2.28, 2.27, 2.32, 2.36, 246, 247, 251, 273, 1518, 1590, 1521, 107.4, 186, 414, 1302];
Mu_exp = [0.947, 0.954, 0.959, 0.943, 0.936, 0.952, 0.583, 0.579, 0.290, 0.294, 0.279, 0.267, 0.724, 0.612, 0.512, 0.273];

% Initial Guess for a, b, c & d.
Params_0 = [1, 1, 1, 1];

% Calculating param values for min RMSE and printing the result.
[param, RMS] = fminsearch(@(param)RMSE(Nu_exp, Nu(Re_exp, Pr_exp, Mu_exp, param), 16), Params_0);

disp("So the predicted relation for Nusselt number is:");
disp("Nu = " ...
    + param(1) + " Re^(" + param(2) ...
    + ") Pr^(" + param(3) + ") Mu^(" ...
    + param(4) + ")");

% Plotting the Graph of Nu_exp and Nu_pred.
Nu_pred = Nu(Re_exp, Pr_exp, Mu_exp, param);
plot(Nu_exp);
hold on
plot(Nu_pred, '-o');
hold on

% Adding Metadata.
title("Comparison of Predicted & Experimental Values of Nu");
ylabel("Re Value");
xlabel("Trial No.");
legend("Nu Predicted","Nu Measured");
hold off

% Calulating exec time.
timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

%% Optimization Techniques %%

% Steps:
% a) Write down the equations concerned with the problem.
% b) Find the variable that's to be optimised and related conditions on it.
% c) Use fmincon to find the minimum value by specifying the conditions.
%    If you need to find max value take negative of that variable in 
%    fmincon.

% Clearing terminal & memory, measuring exec time.
clc;
clear;
tic;

% Example 1) An object is ... (Check Optimization.pdf)
z0 = 100;
m = 80;
g = 9.81;
c = 15;
v0 = 55;

% Plotting the graph to find the point of max height.
tspan = 0 : 15;
z = height(z0, m, g, c, v0, tspan);

fig1 = figure();
plot(tspan, z);
yline(0);
xline(0);
grid();

% Adding Metadata.
title("Height(z) vs time");
ylabel("Height");
xlabel("time");

% Initial Guess for time
t0 = 7;

% Constraints

%  t > 0
A = -1;
B = 0;

% Using fmincon to find time at max height.
[t, h] = fmincon(@(t) -1 * height(z0, m, g, c, v0, t), t0, A, B);

disp("So max height = " + -1 * h + " which occurs at t = " + t);

% Calulating exec time.
timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

% Clearing memory, measuring exec time.
clear;
tic;

% Example 2) A total charge Q ... (Check Optimization.pdf)
% Sol)

% Constants
q = 2e-5;
Q = 2e-5;
a = 0.9;
e0 = 8.85e-12;

% Initial Guess for x
x0 = 2;

% Constraints

%  t > 0
A = -1;
B = 0;

% Using fmincon to find time at max force.
[x, F] = fmincon(@(x) -1 * eforce(q, Q, e0, a, x), x0, A, B);

disp("So max Force = " + -1 * F + " which occurs at x = " + x);

% Calulating exec time.
timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

% Clearing memory, measuring exec time.
clear;
tic;

% Example 3) A train (series) of four ... (Check Optimization.pdf)
% Sol) We have lower bound on all volumes to be greater than 0 
% & upper bound of 20. Also we need to minimize c4 for 
% max product in outlet.

% Constants
q = 71/3600;
Vt = 40;
c0 = 20;
k = 6.25e-3;
n = 2.5;

% Initial Guess for Volume
V0 = [5, 5, 5, 5];

% Constraints

% Vi > 0 for i = 1 : 4
lb = [0, 0, 0, 0];
ub = [20, 20, 20, 20];

% Sum of all volume = 20
Aeq = [1, 1, 1, 1];
Beq = 20;

% Using fmincon to find time at max product conc.
[V, c4] = fmincon(@(V)prodConc(q, V, c0, k, n), V0, [], [], Aeq, Beq, lb, ub);  

disp("The Volume of reactors are " + V(1) + ", " + V(2) + ...
    ", " + V(3) + ", " + V(4) + " (in m^3) respectively at min " + ...
    "outlet concentration of " + c4);

% Calulating exec time.
timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

%% Interpolation Methods %%

% Steps
% a) Given a data between input and output variable, calculate (delta)^n, 
%    where n is the order of preciseness.
% b) for y = f(xi + ah)             h is step size of data.
%    

%% Solving ODEs for IVPs and BVPs. %%

% Steps
% 
% For Parameter estimation with ODEs
% a) For the given data, taking an initial guess, for model params, 
%    calculate model predicted value of result.
% b) The above will be done by solving the diff eqns by using ode45.
% c) Calculate RMSE for the model predicted and actual result value.
% d) Minimize RMSE till a desired value using fminsearch.
% 
% For Euler's method
% a) Let dy/dt = f(t, y), then yi+1 - yi/ delta(t) = f(ti, yi)
% a) f(ti+1) = f(ti) + f(ti, yi) * delta(t)
% b) Perform a for loop for a large number of points for
%    accurate graph.

% Clearing terminal & memory, measuring exec time.
clc;
clear;
tic;

% Example 1) A first order reaction occurs ... (Check ODE_IVP.pdf)
% Measured variable values from experiment.
t = 0 : 50 : 500;
C_exp=[0.01,0.0084,0.0068,0.0054,0.0042,0.0034,0.0029,0.0027,0.0025,0.0024,0.0024];
T_exp=[300,303.30,306.2,308.62,310.47,311.75,312.55,313.03,313.31,313.48,313.58];

% Initial Guess for Alpha and Beta
param_0 = [1, 1];

% Calculating param values for min RMSE and printing the result.
[param, RMS] = fminsearch(@(params) odeRMS(C_exp, T_exp, 11, params, t), param_0);

disp("The value of parameters alpha and beta are " + ...
    param(1) + " and " + param(2) + " respectively.")

% Calulating exec time.
timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % %

% Clearing memory, measuring exec time.
clear;
tic;

% Example) Solve ... (Check ODE_IVP.pdf)
% Sol)
% f' = 4e^0.8t - 0.5y
y0 = 2;

samplespace = zeros(10, 1);
samplespace(1) = y0;
t = 0 : 0.1 : 5;

for i = 2: 51
    samplespace(i) = f(samplespace(i - 1), t(i), 1);
end

% Printing the result
disp("      y         t");
disp([samplespace, t']);

% Plotting y vs t.
plot(t, samplespace)
title("y vs t")
xlabel("t")
ylabel("y")

% Calulating exec time.
timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

%% % % % % % % % % % Functions Used % % % % % % % % % % % % % % %%

% Parameter Estimation

function y = Nu(Re, Pr, Mu, Par)
    y = Par(1) .* Re.^ Par(2) .* Pr .^ Par(3) .* Mu .^ Par(4);
end

function r = RMSE(Nu_exp, Nu_pred, n)
    err = Nu_exp - Nu_pred;
    r = sqrt(sum(err.^2)/n);
end

% Optimization Techniques %

% Example 1
function z = height(z0, m, g, c, v0, t)
    z = z0 + m/c * (v0 + m*g/c) * (1 - exp(-c*t/m)) - m*g*t/c;
end

% Example 2
function f = eforce(q, Q, e, a, x)
    f = 1/(4*pi*e) * q*Q* x/(x.^2 + a.^2).^3/2;
end

% Example 3
function c = concentration(q, V, k, n, cprev, ccurr)
    c = q*cprev - q*ccurr - V*k * ccurr .^ n;
end

function c4 = prodConc(q, Vg, c0, k, n)
    c = zeros(4, 1);
    cg = 1;
    c(1) = fsolve(@(C)concentration(q, Vg(1), k, n, c0, C), cg);
    c(2) = fsolve(@(C)concentration(q, Vg(2), k, n, c(1), C), cg);
    c(3) = fsolve(@(C)concentration(q, Vg(3), k, n, c(2), C), cg);
    c(4) = fsolve(@(C)concentration(q, Vg(4), k, n, c(3), C), cg);
    
    c4 = c(4);
end

% Solving ODEs for IVPs and BVPs % 

% Example 1
function dSdt = concTempdiff(S, t, param)
    dSdt = zeros(2,1);
    dSdt(1) = 5e-5 - S(1) * (5e-3 + exp(param(1) - 11324/S(2)));
    dSdt(2) = 1.74 - 0.0057 * S(2) + S(1) * exp(param(2) - 11324/S(2));
end

function r = odeRMS(Cexp, Texp, n, params, tspan)
    
    S0 = [Cexp(1), Texp(1)];
    [~, Spred]  = ode45(@(t, S)concTempdiff(S, t, params), tspan, S0);
    
    Cpred = Spred(:, 1)';
    Tpred = Spred(:, 2)';
    
    cerr = Cexp - Cpred;
    terr = Texp - Tpred;
    
   
    Crms = sqrt(sum(cerr.^2)/n);
    Trms = sqrt(sum(terr.^2/n));

    r = Crms + Trms;
end

% Example
function y = f(yi, ti, delta)
    y = (4 * exp(0.8 * ti) - yi) * delta + yi;
end