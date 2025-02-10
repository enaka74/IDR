function IDD_score_list_and_param = NFL_signaling_RecordScore_RandomSampling_minscale(N)

% Model 1

%
% Constants

% Convert numbers from 0 to2^n-1 into binaries, which are used as powers for parameters
% Parameter regimes are adopted from a previous work by Qiao et al.(https://elifesciences.org/articles/76188/figures#content)

a = 10 .^ (2*rand(N,2)); % random number from 10^-1 to 10^1
b = 10 .^ (2*rand(N,1)); % random number from 10^-1 to 10^1
K = 10 .^ (3*rand(N,3)); % random number from 10^-3 to 10^2
Kpx = 10 .^ (3*rand(N,1)); % random number from 10^-3 to 10^2
Kpy = 10 .^ (3*rand(N,1)); % random number from 10^-3 to 10^2
r = 10 .^ (2*rand(N,2)); % random number from 10^-1 to 10^1
Px = 10 .^ (1*rand(N,1)-1); % random number from 10^-1 to 10^0
Py = 10 .^ (1*rand(N,1)-1); % random number from 10^-1 to 10^0

param = [a,b,K,Kpx,Kpy,r,Px,Py];

parfor idx=1:N
    IDD_score_list(idx,:) = main_IDD_1(param(idx,:));
end

%IDD_score_list;
IDD_score_list_and_param = [IDD_score_list,param];

%{
nexttile(2);
imagesc(m_list,n_list,result_peak);
xlabel('m') 
ylabel('n');
xticks([1 2 3 4 5 6])
yticks([1 2 3 4 5 6])
title('AUC peak')
set(gca,'YDir','normal','FontSize',20)
colorbar;
%}

end

function IDD_score = main_IDD_1(param)
    % Constants
    D=1; u = 0.01;

    % Variables
    %D_list = logspace(-1,1,17);
    D_list = logspace(0,2,17);
    %u_list = logspace(-3,1,21);
    %Dose_list = D_list;
    
    for j=1:length(D_list)
        %tic
        %param_ = [param,u,Dose_list(j)];
        T=1000;
        Y0 = zeros(1,2);
        param_ = [param,u,D_list(j)];
        tRange = linspace(0,T,2000);
        opts = odeset('AbsTol',1e-6,'MaxStep',min(D_list(j),10));
        [tSol,YSol] = ode23s(@(t,Y)model_NFL_signaling(t,Y,param_),tRange,Y0,opts);
        
        I = calculate_AUC(tSol,YSol(:,end));
        T = calculate_T(tSol,YSol(:,end));
        Q = calculate_Q(tSol,YSol(:,end));
        AUC(j) = I;
        signaling_time(j) = T/I; % tau
        signal_duration(j) = sqrt(Q/I-(T/I)^2); % theta
        Amp(j) = max(YSol(:,end));
        
    end    
    if ~all(isreal(signal_duration))
        result_IDD = 0; %NaN;
        result_AUC = 0; %NaN;
    else
        result_IDD = evaluate_IDD(D_list,signal_duration);
        result_AUC = evaluate_AUC(D_list,AUC);
    end
    
    IDD_score = [result_IDD,result_AUC];
end

function dYdt = model_NFL_signaling(t,Y,param)
    a = param(1:2)'; b = param(3); 
    K = param(4:6)'; Kpx = param(7); Kpy = param(8);
    r = param(9:10)'; Px = param(11); Py = param(12);
    u = param(13); D = param(14);
    
    u = IFFL_input(t, u, D);
    x = Y(1); y = Y(2);
    
    dx = a(1).* y.*(1-x)./(K(1)+y+(1-x)) - r(1) .* Px .* x ./ (Kpx+Px+x);
    dy = a(2) .* u .* (1-y) ./ (K(2)+u+(1-y)) - b .* x .* y ./ (K(3)+x+y) - r(2) .* Py .* y ./ (Kpy+Py+y);

    dYdt = [dx'; dy'];
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
            IDD(count) = max_min_list(j,1)/min(max_min_list(j:end,2))-1;
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
            IDD(count) = max_min_list(j,1)/min(max_min_list(j:end,2))-1;
        end
    end
    IDD_score = max(IDD);
end
