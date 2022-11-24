% Aryan Ranjan
% 20112026

% Surprise Quiz 1

% Secant Method

% Initial Guesses
x0 = 0;
x1 = 1;

% Tolerance value
TOLERANCE = 1e-4;

% Function f(x)
fun = @(x) 3*x + sin(x) - exp(x);

% Algorithm starts here

if abs(fun(x0)) < abs(fun(x1))
        temp = x0;
        x0 = x1;
        x1 = temp;
end

x2 = x1 - (fun(x1) * (x1 - x0)/(fun(x1) - fun(x0)));
while(fun(x2)) > TOLERANCE
    x2 = x1 - (fun(x1) * (x1 - x0)/(fun(x1) - fun(x0)));
    x0 = x1;
    x1 = x2;
end

% Algorithm ends here

disp("Root of f(x) is x = " + x2);

Tc = 405.5;
Pc = 111.3;
R = 0.08206;

a = (27/64) * Pc* (R*Tc/Pc).^2;
b = (R*Tc)/(8*Pc);

fun = @(V, P, T, a, b, R) (P + (a/V.^2)) * (V-b) - R * T;

R = 0.08206;
P = 56;
T = 450;

fun2 = @(V) fun(V, P, T, a, b, R);
v = fzero(fun2, 0.1)
fun = @(v, g, m, Cd, t) v - sqrt((g * m)/Cd) * tanh(sqrt((g * Cd)/m) * t);

g = 9.81;
v = 36;
t = 4;
Cd = 0.25;

criticalMass = fzero(@(m) fun(v, g, m, Cd, t), 1)
x0 = 0;
x1 = 1;
TOLERANCE = 1e-4;
fun = @(x) 3*x + sin(x) - exp(x);
counter = 0
% Bisection method/ Newton Leibiniz Method
while abs(x1-x0) > 2 * TOLERANCE
    x2 = (x1 + x0) / 2;
    if fun(x2) * fun(x0) < 0
        x1 = x2;
    else
        x0 = x2;
    end
    counter = counter+ 1;
end
x0
x1
counter
TOLERANCE = 1e-4;
fun = @(x) 3*x + sin(x) - exp(x);
gun = @(x) 3 + cos(x) - exp(x);
counter = 0;
x0 = 1;
x1 = x0 - fun(x0)/gun(x0);

if fun(x0) ~= 0 && gun(x0) ~= 0
    while abs(x1-x0) > TOLERANCE || abs(fun(x1)) > TOLERANCE
        x1 = x0;
        x0 = x0 - fun(x0)/gun(x0);
        counter = counter + 1;
    end
end

x1
counter