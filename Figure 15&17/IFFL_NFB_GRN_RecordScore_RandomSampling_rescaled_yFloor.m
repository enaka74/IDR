function IDD_score_list_and_param = IFFL_NFB_GRN_RecordScore_RandomSampling_rescaled_yFloor(N,u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3: m-n plane %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constants
%N=100; % number of trials

rng('default');

% Convert numbers from 0 to2^n-1 into binaries, which are used as powers for parameters
% Parameter regimes are adopted from a previous work by Qiao et al.(https://elifesciences.org/articles/76188/figures#content)

b = 10 .^ (2*rand(N,2)); 
r = 10 .^ (2*rand(N,2)-2); 
K = 10 .^ (4*rand(N,4)-2); 
h = randi(4,[N,4]);

param = [b,r,K,h];

IDD_score_list = nan(N,3);   % [IDD, IDA, failcode]

parfor idx=1:N
    [IDD, IDA, failcode] = main_IDD(param(idx,:),u);
    IDD_score_list(idx,:) = [IDD, IDA, failcode];
end

IDD_score_list_and_param = [IDD_score_list, b, r, K, h];

end

function [IDD, IDA, failcode] = main_IDD(param,u)

D_list = logspace(-1,1,17);

failcode = 0;
Y0 = [0;0];

AUC = nan(size(D_list));
theta = nan(size(D_list));

Tend = 40;
tolNegAbs = 1e-10;

% ===== ADDED: screening threshold (absolute) =====
y_floor = 1e-4;        % tentative
y_max_allD = 0;        % max over all D and time
% ================================================

for j = 1:numel(D_list)
    r = param(3:4); tau = sum(1 ./ r);
    Tfinal = Tend*tau;
    param_ = [param,u,D_list(j)*tau];
    tRange = linspace(0,Tfinal,2000);
    opts = odeset('RelTol', 1e-6,'AbsTol',1e-8*tau,'MaxStep',min(0.1*tau,5),'NonNegative',[1,2]);

    [tSol,YSol] = ode15s(@(t,Y) model_IFFL_NFL_GRN(t,Y,param_), tRange, Y0, opts);

    if any(~isfinite(YSol(:)))
        failcode = max(failcode, 1);
        theta(j) = NaN; AUC(j) = NaN;
        continue;
    end

    y = YSol(:,end);

    % ===== ADDED: update max y across D =====
    y_max_allD = max(y_max_allD, max(y));
    % =======================================

    % --- your integration (rect_int or trapz) ---
    I = calculate_AUC(tSol,y);
    Tm = calculate_T(tSol,y);
    Qm = calculate_Q(tSol,y);

    AUC(j) = I;

    if ~isfinite(I) || I <= 1e-12 || ~isfinite(Tm) || ~isfinite(Qm)
        failcode = max(failcode, 2);
        theta(j) = NaN;
    else
        var_t = Qm/I - (Tm/I)^2;
        if var_t < 0
            if var_t > -tolNegAbs
                var_t = 0;
            else
                failcode = max(failcode, 3);
                theta(j) = NaN;
                continue;
            end
        end
        theta(j) = sqrt(var_t);
    end
end

% ===== ADDED: screening based on y_max =====
% If system is essentially OFF for all D, treat as non-decoding:
if isfinite(y_max_allD) && (y_max_allD < y_floor)
    IDD = 0;
    IDA = 0;
    failcode = max(failcode, 5);
    return;
end
% ===========================================

if any(isnan(theta)) || ~all(isreal(theta))
    IDD = NaN; IDA = NaN;
    failcode = max(failcode, 4);
    return;
end

IDD = evaluate_IDD(D_list, theta);
IDA = evaluate_AUC(D_list, AUC);

if ~isfinite(IDD) || ~isfinite(IDA)
    failcode = max(failcode, 4);
end
end


function dYdt = model_IFFL_NFL_GRN(t,Y,param)
    b = param(1:2);
    r = param(3:4);
    K = param(5:8);
    h = param(9:12);
    u = param(13); D = param(14);
    
    u = IFFL_input(t,u,D);
    x = Y(1); y = Y(2);
    
    dx = b(1)*(u/K(1))^h(1)/(1+(u/K(1))^h(1)+(y/K(2))^h(2)) - x*r(1); 
    dy = b(2)*(u/K(3))^h(3)/(1+(u/K(3))^h(3)+(x/K(4))^h(4)) - y*r(2);
    dYdt = [dx;dy];
end

function U = IFFL_input(t,u,D)
    if t <= D
        U = u;
    elseif t > D 
        U = 0;
    end
end

function I = calculate_AUC(t,y)
    I = 0;
    for i=2:length(t)
        I = I + (t(i)-t(i-1))*(y(i));
    end
end

function T = calculate_T(t,y)
    T = 0;
    for i=2:length(t)
        T = T + (t(i)-t(i-1))*(y(i))*(t(i)+t(i-1))/2;
    end
end

function Q = calculate_Q(t,y)
    Q = 0;
    for i=2:length(t)
        Q = Q + (t(i)-t(i-1))*(y(i))*((t(i)+t(i-1))/2)^2;
    end
end

function IDD_score = evaluate_IDD(D_list,theta)
    idx_max = 1;
    idx_min = 1;
    max_min_list = [];
    count = 0;
    for i=2:length(D_list)
        if theta(i-1) <= theta(i)
            if idx_max ~= i-1
                count = count +1;
                %idx_list(count,:) = [idx_max,idx_min];
                max_min_list(count,:) = [theta(idx_max),theta(idx_min)];
                idx_max = i;
                idx_min = i;
            else
                idx_max = i;
                idx_min = i;
            end
        elseif theta(i-1) > theta(i)
            if i==length(D_list)
                count = count + 1;
                idx_min = i;
                max_min_list(count,:) = [theta(idx_max),theta(idx_min)];
            end
            idx_min = i;
        end
    end
    IDD = [];
    if count == 0
        IDD = 0;
    else
        for j=1:count
            IDD(j) = max_min_list(j,1)/min(max_min_list(j:end,2))-1;
        end
    end
    IDD_score = max(IDD);
end

function IDD_score = evaluate_AUC(D_list,AUC)
    idx_max = 1;
    idx_min = 1;
    max_min_list = [];
    count = 0;
    for i=2:length(D_list)
        if AUC(i-1) <= AUC(i)
            if idx_max ~= i-1
                count = count +1;
                %idx_list(count,:) = [idx_max,idx_min];
                max_min_list(count,:) = [AUC(idx_max),AUC(idx_min)];
                idx_max = i;
                idx_min = i;
            else
                idx_max = i;
                idx_min = i;
            end
        elseif AUC(i-1) > AUC(i)
            if i==length(D_list)
                count = count + 1;
                idx_min = i;
                max_min_list(count,:) = [AUC(idx_max),AUC(idx_min)];
            end
            idx_min = i;
        end
    end
    IDD = [];
    if count == 0
        IDD = 0;
    else
        for j=1:count
            IDD(j) = max_min_list(j,1)/min(max_min_list(j:end,2))-1;
        end
    end
    IDD_score = max(IDD);
end

function I = calculate_AUC_trapz(t,y)
I = trapz(t, y);
end

function Tm = calculate_T_trapz(t,y)
Tm = trapz(t, t.*y);
end

function Qm = calculate_Q_trapz(t,y)
Qm = trapz(t, (t.^2).*y);
end