function IDD_score_list_and_param = AIF_signaling_RecordScore_RandomSampling_rescaled_stable(N,u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter random sampling for pure AIF model (manuscript ODE)
%%%
%%% States:
%%%   y  : regulated species
%%%   z1 : controller species 1
%%%   z2 : controller species 2
%%%
%%% Dynamics:
%%%   dy  = a1*u*(1-y)/(K1 + u + (1-y)) ...
%%%       - b*z1*y/(K2 + z1 + y) ...
%%%       - r1*P1*y/(Kp1 + P1 + y);
%%%
%%%   dz1 = mu          - eta*z1*z2 - r2*z1;
%%%   dz2 = theta*(1-y) - eta*z1*z2 - r3*z2;
%%%
%%% For each random parameter set, this script computes:
%%%   - adaptation_ratio, adaptation_flag
%%%       (2-step input: u = u1 for t < t_switch, u = u2 for t >= t_switch)
%%%   - IDD score, IDA score
%%%       (pulse input: u(t) = u for 0 <= t <= D, 0 otherwise, D = logspace(0,2,17))
%%%
%%% Output per row:
%%%   [adaptation_ratio, adaptation_flag, IDD, IDA, parameters...]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('default');

%% ------------------------------------------------------------------------
% Random sampling ONLY for required AIF parameters
%% ------------------------------------------------------------------------

a1    = 10 .^ (2*rand(N,1));      % 10^0 – 10^2
b     = 10 .^ (2*rand(N,1));      % 10^0 – 10^2

K1    = 10 .^ (3*rand(N,1));      % 10^-1 – 10^1
K2    = 10 .^ (3*rand(N,1));      % 10^-1 – 10^1
Kp1   = 10 .^ (3*rand(N,1));      % 10^-1 – 10^0

r1    = 10 .^ (2*rand(N,1));      % 10^0 – 10^2
r2    = 10 .^ (2*rand(N,1)-2);      % 10^0 – 10^2 (first-order decay of z1)
r3    = 10 .^ (2*rand(N,1)-2);      % 10^0 – 10^2 (first-order decay of z2)
P1    = 10 .^ (1*rand(N,1)-1);    % 10^-1 – 10^0

mu    = 10 .^ (2*rand(N,1)-1);    % 10^-1 – 10^1
theta = 10 .^ (2*rand(N,1)-1);    % 10^-1 – 10^1
eta   = 10 .^ (3*rand(N,1)+1);    % 10^-1 – 10^2

% Parameter layout (per row):
% [a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta]
param = [a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta];

%% ------------------------------------------------------------------------
% Evaluate (adaptation, IDD, IDA) for each parameter set
%% ------------------------------------------------------------------------
% Columns: [adaptation_ratio, adaptation_flag, IDD, IDA]
score_list = zeros(N,4);

for idx = 1:N
    score_list(idx,:) = main_all_scores(param(idx,:),u); %#ok<PFBNS>
end

% Final output:
% [adaptation_ratio, adaptation_flag, IDD, IDA, a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta]
IDD_score_list_and_param = [score_list, param];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main wrapper: compute adaptation_ratio, adaptation_flag, IDD, IDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scores = main_all_scores(param_row,u)

    % ---- 1) Adaptation (2-step input) ----
    [adaptation_ratio, adaptation_flag] = evaluate_adaptation(param_row,u);

    % ---- 2) IDD, IDA (pulse input) ----
    [IDD, IDA] = evaluate_IDD_IDA(param_row,u);

    scores = [adaptation_ratio, adaptation_flag, IDD, IDA];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptation evaluation: 2-step input u1 -> u2 at t_switch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [adaptation_ratio, adaptation_flag] = evaluate_adaptation(param_row,u)

    % Extract parameters
    a1    = param_row(1);
    b     = param_row(2);
    K1    = param_row(3);
    K2    = param_row(4);
    Kp1   = param_row(5);
    r1    = param_row(6);
    r2    = param_row(7);
    r3    = param_row(8);
    P1    = param_row(9);
    mu    = param_row(10);
    theta = param_row(11);
    eta   = param_row(12);

    % 2-step input parameters (aligned with AIF_simulate_adaptation)
    u1       = u*0.1; % 0.01;   % input amplitude for t < t_switch
    u2       = u;     % 0.10;   % input amplitude for t >= t_switch
    t_switch = 600;    % [h]

    % ODE parameter vector
    % [a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta, u1, u2, t_switch]
    param_AIF = [a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta, u1, u2, t_switch];

    % Simulation settings
    T_final = 1200;                      % [h]
    tRange  = linspace(0, T_final, 1200);
    Y0      = zeros(3,1);                % [y, z1, z2]
    opts    = odeset('AbsTol',1e-8, 'RelTol',1e-6, 'MaxStep',1);

    % Solve ODE
    [tSol, YSol] = ode23s(@(t,Y) model_AIF_twoStep_adapt(t, Y, param_AIF), ...
                          tRange, Y0, opts);

    y = YSol(:,1);

    % t_ref = value at t_switch
    idx_ref = find(tSol >= t_switch, 1, 'first');
    if isempty(idx_ref)
        adaptation_ratio = NaN;
        adaptation_flag  = 0;
        return;
    end
    y_ref = y(idx_ref);

    % y_peak: maximum after t_switch
    [y_peak, idx_peak_rel] = max(y(idx_ref:end));
    idx_peak = idx_ref + idx_peak_rel - 1;

    % y_inf: mean of last 10% of trajectory, after the peak
    n_points   = numel(y);
    tail_start = max(idx_peak + 1, round(0.9 * n_points));
    if tail_start > n_points
        adaptation_ratio = NaN;
        adaptation_flag  = 0;
        return;
    end
    y_inf = mean(y(tail_start:end));

    % Adaptation ratio: (y_peak - y_inf)/(y_peak - y_ref)
    denom = (y_peak - y_ref);
    if denom <= 0
        adaptation_ratio = NaN;
        adaptation_flag  = 0;
    else
        adaptation_ratio = (y_peak - y_inf) / denom;
        adaptation_flag  = (adaptation_ratio > 0.9) && (y_peak > y_inf);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDD/IDA evaluation: pulse input 0 -> u -> 0, D = logspace(0,2,17)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IDD, IDA] = evaluate_IDD_IDA(param_row,u)

    % Extract parameters
    a1    = param_row(1);
    b     = param_row(2);
    K1    = param_row(3);
    K2    = param_row(4);
    Kp1   = param_row(5);
    r1    = param_row(6);
    r2    = param_row(7);
    r3    = param_row(8);
    P1    = param_row(9);
    mu    = param_row(10);
    theta = param_row(11);
    eta   = param_row(12);

    % Pulse input amplitude
    %u = 0.01;

    % Pulse durations (temporal dose)
    D_list = logspace(0, 2, 17);

    % Simulation settings
    T_final = 1000;
    tRange  = linspace(0, T_final, 2000);
    Y0      = zeros(3,1);  % [y, z1, z2]
    
    AUC_list   = zeros(size(D_list));
    tau_list   = zeros(size(D_list));
    theta_list = zeros(size(D_list));
    Amp_list   = zeros(size(D_list)); %#ok<NASGU> % keep for possible debugging/analysis

    for j = 1:numel(D_list)
        D = D_list(j);

        % Parameter for pulse ODE:
        % [a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta, u, D]
        param_AIF = [a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta, u, D];
        opts    = odeset('AbsTol',1e-8, 'RelTol',1e-6, 'MaxStep',min(D_list(j),10));

        [tSol, YSol] = ode23s(@(t,Y) model_AIF_pulse_IDD(t, Y, param_AIF), ...
                              tRange, Y0, opts);

        y = YSol(:,1);

        I  = calculate_AUC(tSol, y);
        Tm = calculate_T(tSol, y);
        Qm = calculate_Q(tSol, y);

        AUC_list(j) = I;
        if I > 0
            tau_list(j)   = Tm / I;
            theta_list(j) = sqrt(Qm / I - (Tm / I)^2);
        else
            tau_list(j)   = NaN;
            theta_list(j) = NaN;
        end
        Amp_list(j) = max(y);
    end

    if ~all(isreal(theta_list))
        IDD = NaN;
        IDA = NaN;
        return;
    end

    IDD = evaluate_IDD(D_list, theta_list);
    IDA = evaluate_IDA(D_list, AUC_list);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE: 2-step input (for adaptation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dYdt = model_AIF_twoStep_adapt(t, Y, param)
    % param = [a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta, u1, u2, t_switch]

    a1       = param(1);
    b        = param(2);
    K1       = param(3);
    K2       = param(4);
    Kp1      = param(5);
    r1       = param(6);
    r2       = param(7);
    r3       = param(8);
    P1       = param(9);
    mu       = param(10);
    theta    = param(11);
    eta      = param(12);
    u1       = param(13);
    u2       = param(14);
    t_switch = param(15);

    if t < t_switch
        U = u1;
    else
        U = u2;
    end

    y  = Y(1);
    z1 = Y(2);
    z2 = Y(3);

    dy  = a1*U*(1-y) / (K1 + U + (1-y)) ...
        - b*z1*y    / (K2 + z1 + y) ...
        - r1*P1*y   / (Kp1 + P1 + y);

    dz1 = mu            - eta*z1*z2 - r2*z1;
    dz2 = theta*(1-y)   - eta*z1*z2 - r3*z2;

    dYdt = [dy; dz1; dz2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE: pulse input (for IDD/IDA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dYdt = model_AIF_pulse_IDD(t, Y, param)
    % param = [a1, b, K1, K2, Kp1, r1, r2, r3, P1, mu, theta, eta, u, D]

    a1    = param(1);
    b     = param(2);
    K1    = param(3);
    K2    = param(4);
    Kp1   = param(5);
    r1    = param(6);
    r2    = param(7);
    r3    = param(8);
    P1    = param(9);
    mu    = param(10);
    theta = param(11);
    eta   = param(12);
    u     = param(13);
    D     = param(14);

    if t <= D
        U = u;
    else
        U = 0;
    end

    y  = Y(1);
    z1 = Y(2);
    z2 = Y(3);

    dy  = a1*U*(1-y) / (K1 + U + (1-y)) ...
        - b*z1*y    / (K2 + z1 + y) ...
        - r1*P1*y   / (Kp1 + P1 + y);

    dz1 = mu            - eta*z1*z2 - r2*z1;
    dz2 = theta*(1-y)   - eta*z1*z2 - r3*z2;

    dYdt = [dy; dz1; dz2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-moment integrals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = calculate_AUC(t,y)
    I = 0;
    for i = 2:length(t)
        I = I + (t(i)-t(i-1)) * y(i);
    end
end

function Tm = calculate_T(t,y)
    Tm = 0;
    for i = 2:length(t)
        Tm = Tm + (t(i)-t(i-1)) * y(i) * (t(i)+t(i-1))/2;
    end
end

function Qm = calculate_Q(t,y)
    Qm = 0;
    for i = 2:length(t)
        ti_mid = (t(i)+t(i-1))/2;
        Qm = Qm + (t(i)-t(i-1)) * y(i) * (ti_mid^2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDD and IDA evaluation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IDD_score = evaluate_IDD(D_list,theta)
    idx_max = 1;
    idx_min = 1;
    max_min_list = [];
    count = 0;
    for i = 2:length(D_list)
        if theta(i-1) <= theta(i)
            if idx_max ~= i-1
                count = count + 1;
                max_min_list(count,:) = [theta(idx_max), theta(idx_min)];
                idx_max = i;
                idx_min = i;
            else
                idx_max = i;
                idx_min = i;
            end
        elseif theta(i-1) > theta(i)
            if i == length(D_list)
                count = count + 1;
                idx_min = i;
                max_min_list(count,:) = [theta(idx_max), theta(idx_min)];
            end
            idx_min = i;
        end
    end
    IDD = [];
    if count == 0
        IDD = 0;
    else
        for j = 1:count
            IDD(count) = max_min_list(j,1) / min(max_min_list(j:end,2)) - 1;
        end
    end
    IDD_score = max(IDD);
end

function IDA_score = evaluate_IDA(D_list,AUC)
    idx_max = 1;
    idx_min = 1;
    max_min_list = [];
    count = 0;
    for i = 2:length(D_list)
        if AUC(i-1) <= AUC(i)
            if idx_max ~= i-1
                count = count + 1;
                max_min_list(count,:) = [AUC(idx_max), AUC(idx_min)];
                idx_max = i;
                idx_min = i;
            else
                idx_max = i;
                idx_min = i;
            end
        elseif AUC(i-1) > AUC(i)
            if i == length(D_list)
                count = count + 1;
                idx_min = i;
                max_min_list(count,:) = [AUC(idx_max), AUC(idx_min)];
            end
            idx_min = i;
        end
    end
    IDA = [];
    if count == 0
        IDA = 0;
    else
        for j = 1:count
            IDA(count) = max_min_list(j,1) / min(max_min_list(j:end,2)) - 1;
        end
    end
    IDA_score = max(IDA);
end
