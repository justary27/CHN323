% Indian Institute of Technology Roorkee
% Autumn Semester 2022-23
% CHN 323 Assignment 5

% Aryan Ranjan
% 20112026

% Clearing terminal & memory, measuring exec time.
clc;
clear;
tic;

% Q2) The following model is proposed to ...

% Given Data
% Experiment Values
x = 0 : 1 : 10;

% Experimental Result Values.
Y = [12, 11.0937, 10.3516, 9.7441, 9.2466, 8.8394, 8.506, 8.233, 8.0095, 7.8265, 7.6767];

% Initial Guess For Parameters.
X0 = [1, 1, 1];

% Storing the Parameter Values & Final RMSE Value.
[param, RMS] = fminsearch(@(X)RMSE(Y, x, X), X0);
disp("So the predicted relation for y is:");
disp("y = " ...
    + param(1) + " e^(" + param(2) ...
    + "x) + " + param(3));

% Final Theoretically Predicted Nu Values
YFinal = YRes(x, param);

% Plotting the graphs
p = plot(YFinal);
p(1).Marker = "o";
hold on;
p = plot(Y, 'r--x');


title("Comparison of Predicted & Experimental Values of y");
ylabel("y value");
xlabel("Trial No.");
legend("y Predicted","y Measured");

timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

% % % % % % % % % % Functions Used % % % % % % % % % % % % % % % 

function y = YRes(x, X)
    y = X(1) .* exp(x .* X(2)) + X(3);
end

function rms = RMSE(y, x, X)
    rms = sqrt(mean((YRes(x, X) - y).^2));
end