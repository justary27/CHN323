% Indian Institute of Technology Roorkee
% Autumn Semester 2022-23
% CHN 323 Assignment 5

% Aryan Ranjan
% 20112026

% Clearing terminal & memory, measuring exec time.
clc;
clear;
tic;

% Q1) Consider the data set shown in table ...

% Given Data
% Experiment Values.
Re = [49000, 68600, 84800, 34200, 22900, 1321, 931, 518, 346, 122.9, 54, 84.6, 1249, 1021, 465, 54.8];
Pr = [2.3 , 2.28 , 2.27 , 2.32, 2.36 , 246 , 247 , 251 , 273 , 1518 , 1590 , 1521, 107.4 , 186 , 414 , 1302];
Mu = [0.947 , 0.954 , 0.959 , 0.943 , 0.936 , 0.592 , 0.583 , 0.579 , 0.29 , 0.294 , 0.279 , 0.267 , 0.724 , 0.612 , 0.512 , 0.273];

% Experimental Result Values.
Nu = [277 , 348 , 421 , 223, 177, 114.8 , 95.9 , 68.3, 49.1 , 56 , 39.9 , 47 , 94.2 , 99.9 , 83.1 , 35.9];

% Initial Guess For Parameters.
X0 = [1, 1, 1, 1];

% Storing the Parameter Values & Final RMSE Value.
[param, RMS] = fminsearch(@(X)RMSE(Nu, Re, Pr, Mu, X), X0);
disp("So the predicted relation for Nusselt number is:");
disp("Nu = " ...
    + param(1) + " Re^(" + param(2) ...
    + ") Pr^(" + param(3) + ") Mu^(" ...
    + param(4) + ")");

% Final Theoretically Predicted Nu Values
NuFinal = NuRes(Re, Pr, Mu, param);

% Plotting the graphs
p = plot(NuFinal);
p(1).Marker = "o";
hold on;
p = plot(Nu);


title("Comparison of Predicted & Experimental Values of Nu");
ylabel("Re Value");
xlabel("Trial No.");
legend("Nu Predicted","Nu Measured");

timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

% % % % % % % % % % Functions Used % % % % % % % % % % % % % % % 

function Nu = NuRes(Re, Pr, Mu, X)
    Nu = X(1) .* Re.^X(2) .* Pr.^X(3) .* Mu.^X(4);
end

function rms = RMSE(Nu, Re, Pr, Mu, X)
    Nu_est=NuRes(Re, Pr, Mu,  X);
    rms = sqrt(mean((Nu - Nu_est).^2));
end
