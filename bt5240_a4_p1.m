% BT5240 - Assignment 4 - Problem 1
% Anirudh Rao (BE21B004)

clc, clearvars;

%% Problem 1 (a) - Differential equations for population

% This function defines the system of ODEs governing the population dynamics
% of species A, B, C and signalling molecule X
function derivatives = population_dynamics(t, V, params)

    A = V(1);
    B = V(2);
    C = V(3);
    X = V(4);

    % Extract parameters from input structure
    g_A = params.g_A;
    g_B = params.g_B;
    g_C = params.g_C;
    d_A = params.d_A;
    d_B = params.d_B;
    d_C = params.d_C;
    alpha = params.alpha;
    beta = params.beta;
    delta = params.delta;
    epsilon = params.epsilon;
    eta = params.eta;

    % Define the differential equations
    dAdt = (g_A * A) - (d_A * A) + (delta * A * B) - (epsilon * A * C);
    dBdt = (g_B * B) - (d_B * B) + (beta * eta * X * B);
    dCdt = (g_C * C) - (d_C * C) + (epsilon * A * C);
    dXdt = (alpha * A) - (beta * X * B);

    % Return derivatives as a column vector
    derivatives = [dAdt; dBdt; dCdt; dXdt]; 

end

%% Problem 1 (b) - Simulation of population dynamics for three cases

% Initial conditions: A, B, C, X
V0 = [5; 5; 5; 0];
tspan = [0 5];  % Simulate for 5 hours
options = odeset('NonNegative', 1:4);  % Ensure no negative values

% --- Case 1 ---
params1 = struct('g_A', 0.2, 'd_A', 0.1, ...
                 'g_B', 0.2, 'd_B', 0.1, ...
                 'g_C', 0.2, 'd_C', 0.1, ...
                 'alpha', 0.6, 'beta', 0.85, 'delta', 0.80, ...
                 'epsilon', 0.30, 'eta', 0.6);

[t, V] = ode23(@(t, V) population_dynamics(t, V, params1), tspan, V0, options);
figure;
plot(t, V(:,1), 'r', 'DisplayName', 'A'); hold on;
plot(t, V(:,2), 'g', 'DisplayName', 'B');
plot(t, V(:,3), 'b', 'DisplayName', 'C');
xlabel('Time (h)');
ylabel('Biomass (g)');
title('Case 1');
legend('Location', 'northwest');

% --- Case 2 ---
params2 = struct('g_A', 0.65, 'd_A', 0.02, ...
                 'g_B', 0.7, 'd_B', 0.02, ...
                 'g_C', 0.1, 'd_C', 0.01, ...
                 'alpha', 0.6, 'beta', 0.85, 'delta', 0.80, ...
                 'epsilon', 0.30, 'eta', 0.6);

[t, V] = ode23(@(t, V) population_dynamics(t, V, params2), tspan, V0, options);
figure;
plot(t, V(:,1), 'r', 'DisplayName', 'A'); hold on;
plot(t, V(:,2), 'g', 'DisplayName', 'B');
plot(t, V(:,3), 'b', 'DisplayName', 'C');
xlabel('Time (h)');
ylabel('Biomass (g)');
title('Case 2');
legend('Location', 'northwest');

% --- Case 3 ---
params3 = struct('g_A', 0.3, 'd_A', 0.05, ...
                 'g_B', 0.25, 'd_B', 0.01, ...
                 'g_C', 0, 'd_C', 5, ...
                 'alpha', 0.6, 'beta', 0.85, 'delta', 0.80, ...
                 'epsilon', 0.30, 'eta', 0.6);

[t, V] = ode23(@(t, V) population_dynamics(t, V, params3), tspan, V0, options);
figure;
plot(t, V(:,1), 'r', 'DisplayName', 'A'); hold on;
plot(t, V(:,2), 'g', 'DisplayName', 'B');
plot(t, V(:,3), 'b', 'DisplayName', 'C');
xlabel('Time (h)');
ylabel('Biomass (g)');
title('Case 3');
legend('Location', 'northwest');

%% Problem 1 (c) - Introduction of new species

% New system of ODEs with an additional species D that interacts with C
function new_derivatives = new_population_dynamics(t, V, params)

    A = V(1);
    B = V(2);
    C = V(3);
    X = V(4);
    D = V(5);

    % Extract parameters
    g_A = params.g_A;
    g_B = params.g_B;
    g_C = params.g_C;
    d_A = params.d_A;
    d_B = params.d_B;
    d_C = params.d_C;
    alpha = params.alpha;
    beta = params.beta;
    delta = params.delta;
    epsilon = params.epsilon;
    eta = params.eta;
    g_D = params.g_D;
    d_D = params.d_D;
    gamma = params.gamma;

    % Updated differential equations including D's effect on C
    dAdt = (g_A * A) - (d_A * A) + (delta * A * B) - (epsilon * A * C);
    dBdt = (g_B * B) - (d_B * B) + (beta * eta * X * B);
    dCdt = (g_C * C) - (d_C * C) + (epsilon * A * C) - (gamma * C * D);
    dXdt = (alpha * A) - (beta * X * B);
    dDdt = (g_D * D) - (d_D * D);

    new_derivatives = [dAdt; dBdt; dCdt; dXdt; dDdt]; 

end

% Initial conditions including new species D
V0 = [5; 5; 5; 0; 5];
tspan = [0 5];
options = odeset('NonNegative', 1:5);

% Parameters with interaction term gamma for D-C interaction
params = struct('g_A', 0.2, 'd_A', 0.1, ...
                'g_B', 0.2, 'd_B', 0.1, ...
                'g_C', 0.2, 'd_C', 0.1, ...
                'alpha', 0.6, 'beta', 0.85, 'delta', 0.80, ...
                'epsilon', 0.30, 'eta', 0.6, ...
                'g_D', 0.2, 'd_D', 0.1, 'gamma', 0.8);

[t, V] = ode23(@(t, V) new_population_dynamics(t, V, params), tspan, V0, options);
figure;
plot(t, V(:,1), 'r', 'DisplayName', 'A'); hold on;
plot(t, V(:,2), 'g', 'DisplayName', 'B');
plot(t, V(:,3), 'b', 'DisplayName', 'C');
plot(t, V(:,5), 'k', 'DisplayName', 'D');
xlabel('Time (h)');
ylabel('Biomass (g)');
title('Case 1 with species D');
legend('Location', 'northwest');
