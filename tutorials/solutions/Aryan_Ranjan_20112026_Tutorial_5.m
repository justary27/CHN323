% Indian Institute of Technology Roorkee
% Autumn Semester 2022-23
% CHN 323 Assignment 5

% Aryan Ranjan
% 20112026

% Clearing terminal & memory, measuring exec time.
clc;
clear;
tic;

% Q1) A first order reaction occurs in a ...

% Sol 1) 
% Time Range of Experiment.
tspan = 0:50:500;


% Experimental Result Values.
C_exp = [0.01, 0.0084, 0.0068, 0.0054, 0.0042, 0.0034, 0.0029, 0.0027, 0.0025, 0.0024, 0.0024];
T_exp = [300, 303.30, 306.2, 308.62, 310.47, 311.75, 312.55, 313.03, 313.31, 313.48, 313.58];


% Initial Guess for Parameters Alpha & Beta.
param0 = [10, 0];

% Minimizing the error using fminsearch & calculating params with ode45. 
[a, b] = fminsearch(@(param) RMSE(param, C_exp, T_exp, tspan), param0);
[~, y] = ode45(@(t,Y) eqn(t, Y, a), tspan, [0.0100, 300.00]);

% Plotting the Graphs.

% Plotting Concentration.
figure(1)
plot(tspan, C_exp , tspan, y(:,1)');
grid
xlabel("Time(s)")
ylabel("Concentration")
title("Comparison of Predicted & Experimental Values of Concentration", FontSize=9)
legend("C actual","C Estimated")

% Plotting Temperature
figure(2)
plot(tspan, T_exp, tspan, y(:,2)');
grid
xlabel("Time(s)")
ylabel("Temperature")
title("Comparison of Predicted & Experimental Values of Temperature", FontSize=9)
legend("T actual","T Estimated")

% Printing the result.
fprintf("Estimated Paramters are: \n Alpha = %f \n Beta = %f \n Error = %f \n",a(1),a(2),b)

timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

% % % % % % % % % % Functions Used % % % % % % % % % % % % % % %

function y_diff = eqn(~, Y, param)
    y1 = 0.00005 - Y(1) * (0.005 + exp(param(1) - 11324/Y(2)));
    y2 = 1.74 - 0.0057 * Y(2) + Y(1) * exp(param(2) - 11324/Y(2));
    y_diff = [y1;y2];
end

function rmse = RMSE(param, C_exp, T_exp, tspan)
    [~, y] = ode45(@(t,Y)eqn(t, Y, param), tspan, [0.0100, 300.00]);
    C_pred = y(:,1)';
    T_pred = (y(:,2)');

    % Calculating Error Values.
    abs1 = ((C_exp - C_pred));
    abs2 = ((T_exp - T_pred));

    % Calculating the RMSE errors in C & T.
    e1 = sqrt(mean(abs1.^2));
    e2 = sqrt(mean(abs2.^2));

    % Taking the total RMSE to be mean of both RMSEs.
    rmse = mean([e1,e2]);
end

