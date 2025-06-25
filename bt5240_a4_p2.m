% BT5240 - Assignment 4 - Problem 2
% Anirudh Rao (BE21B004)

clc, clearvars;

%% Problem 2 (b) - Truth table

fprintf('<strong>Truth table</strong>\n\n');

% Generate all 16 binary combinations for 4 proteins
combinations = dec2bin(0:15) - '0';  % Each row: [P53, MYC, MDM2, RB]

% Preallocate output vectors for apoptosis and proliferation outcomes
A = zeros(16,1);
P = zeros(16,1);

% Evaluate each protein combination for apoptosis and proliferation
for i = 1:16
    P53 = combinations(i,1);
    MYC = combinations(i,2);
    MDM2 = combinations(i,3);
    RB  = combinations(i,4);

    A(i) = P53 && ~MDM2;    % Apoptosis if P53 is active and MDM2 is inactive
    P(i) = MYC && ~RB;      % Proliferation if MYC is active and RB is inactive
end

% Create and display the truth table
T = table((1:16)', combinations(:,1), combinations(:,2), combinations(:,3), combinations(:,4), ...
          A, P, ...
          'VariableNames', {'Index', 'P53', 'MYC', 'MDM2', 'RB', 'Apoptosis', 'Proliferation'});
disp(T);

%% Problem 2 (c) - Boolean model simulation

% Simulates Boolean network dynamics given an initial state
function simulateTumourBoolean(init, T, fixed)
    if nargin < 3
        fixed = struct();  % No fixed protein values by default
    end

    S = zeros(T, 4);        % Matrix to store states over time
    A = zeros(T, 1);        % Track apoptosis over time
    P = zeros(T, 1);        % Track proliferation over time
    past = strings(T, 1);   % Keep record of past states to detect cycles

    S(1, :) = init;               % Set initial state
    past(1) = mat2str(init);      % Store initial state as string
    stop = "";                    % Stopping reason

    % Evaluate initial state
    A(1) = init(1) && ~init(3);
    P(1) = init(2) && ~init(4);

    for t = 2:T
        prev = S(t-1, :);
        next = zeros(1, 4);

        % Apply Boolean update rules
        next(1) = prev(2) && ~prev(3);  % P53
        next(2) = ~prev(4);             % MYC
        next(3) = S(1, 3);              % MDM2 remains constant
        next(4) = ~prev(2);             % RB

        % Override values if drugs are fixing certain proteins
        if isfield(fixed, 'P53'),  next(1) = fixed.P53;  end
        if isfield(fixed, 'MYC'),  next(2) = fixed.MYC;  end
        if isfield(fixed, 'MDM2'), next(3) = fixed.MDM2; end
        if isfield(fixed, 'RB'),   next(4) = fixed.RB;   end

        S(t, :) = next;

        % Evaluate current state
        A(t) = next(1) && ~next(3);
        P(t) = next(2) && ~next(4);
        str = mat2str(next);

        % Check for stopping conditions
        if A(t) && ~P(t)
            stop = "Apoptosis"; break;
        elseif isequal(next, prev)
            stop = "Steady state"; break;
        elseif any(past(1:t-1) == str)
            stop = "Oscillation"; break;
        end

        past(t) = str;
    end

    % Trim output to actual number of iterations
    S = S(1:t, :);
    A = A(1:t);
    P = P(1:t);

    % Display simulation results
    T = table((1:t)', S(:,1), S(:,2), S(:,3), S(:,4), A, P, ...
              'VariableNames', {'Time', 'P53', 'MYC', 'MDM2', 'RB', 'Apoptosis', 'Proliferation'});
    disp(T)
    disp("Stopped because: " + stop)
end

% Run simulations for different initial states
T = 20;
initial1 = [0 0 0 0];  % No proteins active initially
initial2 = [1 1 1 1];  % All proteins active initially
initial3 = [0 0 0 1];  % Only RB is active

fprintf("\n<strong>Initial Condition 1</strong>\n\n")
simulateTumourBoolean(initial1, T);
fprintf("\n<strong>Initial Condition 2</strong>\n\n")
simulateTumourBoolean(initial2, T);
fprintf("\n<strong>Initial Condition 3</strong>\n\n")
simulateTumourBoolean(initial3, T);

%% Problem 2 (d) - Boolean model simulation with drugs

T = 20;
initial = [0 1 0 0];  % Starting state for drug simulations

% Apply each drug and simulate effects
fprintf("\n<strong>Drug A</strong> - inhibits MDM2\n\n")
simulateTumourBoolean(initial, T, struct('MDM2', 0));

fprintf("\n<strong>Drug B</strong> - inhibits MYC\n\n")
simulateTumourBoolean(initial, T, struct('MYC', 0));

fprintf("\n<strong>Drug C</strong> - activates RB\n\n")
simulateTumourBoolean(initial, T, struct('RB', 1));
