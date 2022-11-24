% Indian Institute of Technology Roorkee
% Autumn Semester 2022-23
% CHN 323 Assignment 6

% Project Members:
% a) Aryan Ranjan [20112026]
% b) Madhumita [20112061]
% c) Ishika [20112052]

% Clearing terminal & memory, measuring exec time, setting print options.
clc;
clear;
tic;
options = optimset('Display','off');

% Q) Let us consider a simple isothermal reaction A -> B ...

% Sol) 
% Before t = 0s, Caf = 5 and as system is at steady state,
% so dCa/dt = 0. Therefore, F(Caf - Ca) - k.Ca^2.V = 0

F = 9;
k = 2;
Caf = 5;
V = 1;

dCadt = @(Ca) F * (Caf - Ca) - k * Ca.^2 * V;

% Finding initial value of Ca before and at t = 0.
Cai = fsolve(dCadt, 2, options);

% Now value Caf changes to 6.
Caf = 6;
tspan = [0, 10];
[t, Ca] = ode45(@(t, Ca)F * (Caf - Ca) - k * Ca.^2 * V, tspan, Cai);

Ca_final_intial = Ca(size(Ca));
Ca_final = Ca_final_intial(1);

% a) Variation of concentration of A wrt time.
figure(1)
plot(t, Ca, '-o');
grid;
xlabel("Time(min)");
ylabel("Concentration of A (mol/m^3)");
title("Variation of concentration of A with time");

disp("a) The final steady state concentration of A " + ...
    "in the MFR is "+ Ca_final + " mol/m^3.");


disp(" ");


% b) Now considering a train of MFRs to reduce 
% final Ca conc. to 1 unit value. In this case,
% the input for successing MFRs will be the
% steady state value of Ca in the previous 
% reactor. Also dCa/dt = 0 for all reactors.
targetConc = 1;
reactorCount = NoOfReactors(Cai, targetConc, options);

disp("b) Number of reactors required (to approach a " + ...
    "steady state output of Ca = 1 (approx)) " + ...
    "is " + reactorCount + ".");


disp(" ");


% c) Plotting the concentration variation of A 
% in the nth reactor, given that n is within 
% realistic bounds. Here we take example of 
% n = 5.
NthReactorConc(5, reactorCount, options);


timeTaken = toc;
disp(" ");
disp("The execution time is: " + timeTaken + " seconds.");

% % % % % % % % % % End of Program % % % % % % % % % % % % % % % 

% % % % % % % % % % Functions Used % % % % % % % % % % % % % % %
function Nr = NoOfReactors(initalConc, targetConc, options)
    Nr = 0;
    F = 9;
    k = 2;
    V = 1;
    Caf = initalConc;
    
    while (initalConc > targetConc)
        dCadt = @(Ca) F * (Caf - Ca) - k * Ca.^2 * V;
        initalConc = fsolve(dCadt, 2, options);
        Caf = initalConc;
        Nr = Nr + 1;
    end

    % Adding 1 more reactor as we have started 
    % from first reactor's output as initial
    % concentration.
    Nr = Nr + 1;
end

function NthReactorConc(n, maxReactors, options)
    if n > maxReactors
        disp("c) Invalid input as n value is greater " + ...
            "that maximum possible reactors!")
    elseif n <=0
        disp("c) N can't be negative or 0!");
    else
        F = 9;
        k = 2;
        V = 1;
        Caf = 5;
        tspan = [0, 10];

        C_AI = zeros(n, 1);
        C_AF = zeros(n, 1);
        
        % Calculating the steady state feed concentration of the nth 
        % reactor, Caf at steady state, initially when Caf of 1st 
        % reactor is 5, so; dCa/dt = 0.
        for i = 1 : n
            dCadt = @(Ca) F * (Caf - Ca) - k * Ca.^2 * V;
            Caf = fsolve(dCadt, 2, options);
            
            C_AI(i) = Caf;
        end
        
        % Calculating the steady state feed concentration of the nth 
        % reactor, Caf at steady state, initially when Caf of 1st 
        % reactor is 6, so; dCa/dt = 0.
        Caf = 6;
        for i = 1 : n
            dCadt = @(Ca) F * (Caf - Ca) - k * Ca.^2 * V;
            Caf = fsolve(dCadt, 2, options);
            
            C_AF(i) = Caf;
        end

        % Given the initial feed concentration of the n reactors, 
        % solving for the time variation of Ca
        % and plotting it.
        [tn, Cn] = ode45(@(t, Ca) EqnGen(t, Ca, n), tspan, C_AI);
        
        figure(2);
        plot(tn, Cn(:, n), '-o');
        grid;
        xlabel("Time(min)");
        ylabel("Concentration of A (mol/m^3)");
        title("Variation of concentration of A with " + ...
            "time in reactor " + n + ".");

        disp("c) For initial feed concentration Caf = "+ C_AI(n) + ",")
        disp("the steady state concentration of A, i.e Ca " + ...
            "in reactor " + n + " is " + ...
            + C_AF(n) + ".");

    end
end

% Simultaneous Differential Equation Generator
% to feed into ode45.
function E = EqnGen(~, Ca, n)
    F = 9;
    k = 2;
    V = 1;
    E = zeros(n, 1);
    E(1) = F * (6 - Ca(1)) - k * Ca(1).^2 * V;

    if n >= 2
        for i = 2 : n
            E(i) = F * (Ca(i-1) - Ca(i)) - k * Ca(i).^2 * V;
        end
    end
    
end