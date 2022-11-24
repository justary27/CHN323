% CHN 323 ETE

% Aryan Ranjan 
% 20112026

% Clearing memory & terminal, measuring exec time.
clc;
clear;
tic;

% Q1) In an experiment, the following ...
% Sol)
% 
x = 35 : 50;
y = [405, 285, 207, 154, 118, 93, 74, 60, 49, 41, 35, 30, 26, 22, 20, 17];

param0 = [1, 1, 1];

[soln, err] = fminsearch(@(params)RMSE(x, y, 16, params), param0);

disp("Q1)");
disp("The values of a, b and c are " + soln(1) + soln(2) + soln(3) + " respectvely.");
disp(" ");

% Q2) At certain specified conditions ...
% Sol)
% We need to solve the sysytem of equations obtained, we can use ode45

% Initial concentrations of A, B, C & D respectively.
S0 = [50; 5; 0; 0];

% Taking a time span of 10 hrs, note that as all reactions
% are first order, so it will never end.
tspan = linspace(0, 10, 100);

% Using ode45 to find time variations of concentration.
[t, S] = ode45(@(t, S)diffConc(S, t), tspan, S0);

disp("Q2) ")
disp("As time varies till 10 hrs, the concentrations of A, B, C and D"); 
disp("varies as shown in columns from left to right"); 
disp("respectively as below:");
disp(" ");
disp("        A        B        C        D");
disp(S);

% To find out optimal residence time we need time at max B
% conc.

% Conc. of B
Sb = S(:, 2);

% Max Conc. of B
Sbmax = max(Sb);

% Finding optimal time.
index = 1;
for i = 1 : length(Sb)
    if Sb(i) == Sbmax
        index = i;
        break;
    end
end

disp("The optimal residence time is " + t(i) + " hours.")
disp(" ");

% Q4) Heat Transfer through a variable ...
% Sol)
% We get the following system of eqns:

A = [-130, 56, 0, 0; 56, -98, 40, 0; 0, 40, -66, 24; 0, 1, 3, -4];
B = [-72; 0; 0; 0];

% Temperatures T2, T3, T4, T5
X = A\B;

% Temperature variation in fin.
T = [1; X];

disp("Q3")
disp("The temperature varies as below:")
disp(T);
disp(" ")

% Q4) Solve the above problem using ...
% Sol)

xmesh = linspace(0, 1, 10);
solinit = bvpinit(xmesh, [1; 1]);
sol = bvp4c(@diffTemp, @boundTemp, solinit);
Temperatures = sol.y(1,:)';
disp("Q4")
disp("Temperature varies with x as follows: ")
disp(Temperatures)

% Comparing Solution from Q3 and Q4
fig1 = figure();
plot(sol.x, Temperatures, '-o');
hold on
plot(0:0.25:1, T, '-x');

% Adding metadata
grid()
legend("bvp4c", "Point difference");
xlabel("x");
ylabel("Temperature");
title("Temperature vs x");

% Q5) The growth estimation of populations ...
pop = zeros(21, 1);
pop(1) = 2560;
for i = 2 : 21
    pop(i) = rk4(pop(i - 1), i, @f);
end

fig2 = figure();
plot(1950 : 5 : 2050, pop, '-o');

% Adding metadata
xlabel("Year");
ylabel("Population (in millions)");
title("Population vs Year");
grid();

data = pop(1: 11)';
actual_data = [2560, 2780, 3040, 3350, 3710, 4090, 4450, 4850, 5280, 5690, 6080];

rmse = sqrt(sum((actual_data - data).^2)/11);
meanabs = sum(abs(actual_data - data))/11;

disp("Q5")
disp("The RMSE and MAE are " + rmse + " and " + meanabs + " (in millions) respectively.");

% Printing out the exec time.
timeTaken = toc;
disp(" ");
disp("The exec time is " + timeTaken + " seconds.");

function rms = RMSE(x, y, n, params)
    ypred = exp(params(1) + params(2)./(x - params(3)));
    err = y - ypred;
    
    rms = sqrt(sum(err.^2)/n);
end

function dSdt = diffConc(S, ~)
    dSdt = zeros(4, 1);
    dSdt(1) = S(2) - 2 * S(1);
    dSdt(2) = 2 * S(1) - 1.8 * S(2);
    dSdt(3) = 0.2 * S(2);
    dSdt(4) = 0.6 * S(2);
end

function dTdx = diffTemp(x, y)
    dTdx = [y(2)
        1/(5-4*x) * (4 * y(2) + 2 * y(1))];
end

function b = boundTemp(ya, yb)
    b = [ya(1) - 1
        yb(2)];
end

function val = f(~, p)
    val = 0.0178 * p;
end
function y = rk4(yprev, xprev, f)
    h = 5;
    k1 = f(xprev, yprev);
    k2 = f(xprev + h/2, yprev + k1/2);
    k3 = f(xprev + h/2, yprev + k2/2);
    k4 = f(xprev + h, yprev + k3);
    y = yprev + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
end