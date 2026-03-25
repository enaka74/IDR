function IDD_score_list_and_param = IFFL_signaling_RecordScore_RandomSampling_rescaled_stable(N,u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3: m-n plane %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constants
%N=100; % number of trials

% Convert numbers from 0 to2^n-1 into binaries, which are used as powers for parameters
% Parameter regimes are adopted from a previous work by Qiao et al.(https://elifesciences.org/articles/76188/figures#content)

rng('default');

a = 10 .^ (3*rand(N,2)-1); % random number from 10^0 to 10^2
b = 10 .^ (3*rand(N,1)-1); % random number from 10^0 to 10^2
K = 10 .^ (4*rand(N,3)-2); % random number from 10^0 to 10^3
Kp = 10 .^ (4*rand(N,2)-2); % random number from 10^0 to 10^3
r = 10 .^ (3*rand(N,2)-1); % random number from 10^0 to 10^2
P = 10 .^ (1*rand(N,2)-1); % random number from 10^-1 to 10^0

param = [a,b,K,Kp,r,P];

%{
% ======== ADDED: effective saturation diagnostics (minimal change) ========
u0   = 0.01;   % must match main_IDD() input amplitude
x0   = 0.0;    y0 = 0.0;   % initial state proxy
xmid = 0.5;    ymid = 0.5; % mid state proxy (state-dependent saturation)

K1  = K(:,1);  K2  = K(:,2);  K3  = K(:,3);
Kp1 = Kp(:,1); Kp2 = Kp(:,2);
P1  = P(:,1);  P2  = P(:,2);

% Activation "saturation fraction" wrt u (and (1-x) term in denominator)
sat_u_x_init = u0 ./ (K1 + u0 + (1 - x0));
sat_u_y_init = u0 ./ (K2 + u0 + (1 - y0));
sat_u_x_mid  = u0 ./ (K1 + u0 + (1 - xmid));
sat_u_y_mid  = u0 ./ (K2 + u0 + (1 - ymid));

% Simple scale ratios (useful sanity checks)
ratio_u_K1 = u0 ./ K1;
ratio_u_K2 = u0 ./ K2;

% Deactivation saturation wrt substrate (x or y) at mid-state
sat_dx_mid = xmid ./ (Kp1 + P1 + xmid);
sat_dy_mid = ymid ./ (Kp2 + P2 + ymid);

% Inhibition term effective "gain" at mid-state: xy/(K3+x+y)
inh_eff_mid = (xmid * ymid) ./ (K3 + xmid + ymid);

% Optionally: inhibition saturation fraction (how close x+y is to dominating K3)
sat_inh_mid = (xmid + ymid) ./ (K3 + xmid + ymid);

% Pack diagnostics (columns)
diag_metrics = [ ...
    ratio_u_K1, ratio_u_K2, ...
    sat_u_x_init, sat_u_y_init, sat_u_x_mid, sat_u_y_mid, ...
    sat_dx_mid, sat_dy_mid, ...
    inh_eff_mid, sat_inh_mid ...
];
% =========================================================================
%}

parfor idx=1:N
    IDD_score_list(idx,:) = main_IDD(param(idx,:),u);
end

%IDD_score_list;
IDD_score_list_and_param = [IDD_score_list,a, b, K, Kp, r, P];

end

function IDD_score = main_IDD(param,u)
    % Constants
    %u = 1;
    % Variables
    D_list = logspace(0,2,17);
    
    for j=1:length(D_list)
        Y0 = zeros(1,2);
        %r = param(9:10); tau = sum(1 ./ r);
        T = 1000;
        param_ = [param,u,D_list(j)];
        tRange = linspace(0,T,2000);
        % opts = odeset('AbsTol',1e-6,'MaxStep',min(D_list(j),10));
        opts = odeset('RelTol', 1e-6,'AbsTol',1e-8,'MaxStep',min(D_list(j),10));
        
        [tSol,YSol] = ode23s(@(t,Y)model_IFFL(t,Y,param_),tRange,Y0,opts);

        I = calculate_AUC(tSol,YSol(:,end));
        T = calculate_T(tSol,YSol(:,end));
        Q = calculate_Q(tSol,YSol(:,end));
        AUC(j) = I;
        signaling_time(j) = T/I; % tau
        Amp(j) = max(YSol(:,end));

        if I <= 1e-12 || ~isfinite(I) || ~isfinite(T) || ~isfinite(Q)
            theta = NaN;
        else
            var_t = Q/I - (T/I)^2;
            if var_t < 0
                % 数値誤差の範囲だけを潰す（閾値は調整可）
                if var_t > -1e-10
                    var_t = 0;
                else
                    var_t = NaN; % “本当におかしい” ものは NaN のまま
                end
            end
            theta = sqrt(var_t);
        end
        signal_duration(j) = theta; % theta
        
    end    
    if ~(all(isreal(signal_duration)))
        result_IDD = NaN;
        result_AUC = NaN;
        
    else
        result_IDD = evaluate_IDD(D_list,signal_duration);
        result_AUC = evaluate_AUC(D_list,AUC);
        
    end
    
    IDD_score = [result_IDD,result_AUC];
end

function dYdt = model_IFFL(t,Y,param)
    a = param(1:2); b = param(3); 
    K = param(4:6); Kp = param(7:8);
    r = param(9:10); P = param(11:12);
    u = param(13); D = param(14);
    
    u = IFFL_input(t, u, D);
    x = Y(1); y = Y(2);
    
    dx = a(1)*u*(1-x) / (K(1)+u+(1-x)) - r(1)*P(1)*x / (Kp(1)+P(1)+x);
    dy = a(2)*u*(1-y) / (K(2)+u+(1-y)) - b*x*y / (K(3)+x+y) - r(2)*P(2)*y / (Kp(2)+P(2)+y);

    dYdt = [dx; dy];
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
            %IDD(count) = max_min_list(j,1)/min(max_min_list(j:end,2))-1;
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
            %IDD(count) = max_min_list(j,1)/min(max_min_list(j:end,2))-1;
            IDD(j) = max_min_list(j,1)/min(max_min_list(j:end,2))-1;
        end
    end
    IDD_score = max(IDD);
end
