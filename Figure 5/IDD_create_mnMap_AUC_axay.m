function IDD_create_mnMap_AUC_axay(param_idx,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3: m-n plane %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constants

m_range = 1:1:6 ; 
n_range = 1:1:6 ;
[m_list,n_list] = meshgrid(m_range,n_range);
IDD_matrix = [];

parfor idx=1:numel(m_list)
    current_m=m_list(idx);
    current_n=n_list(idx);
    IDD_matrix(idx) = main_IDD(current_m,current_n,param)
end

IDD_matrix = reshape(IDD_matrix,[6,6]);

nexttile(param_idx);
imagesc(m_range,n_range,IDD_matrix);
xlabel('m') 
ylabel('n');
xticks([1 2 3 4 5 6])
yticks([1 2 3 4 5 6])
title('\alpha_x=' + string(param(1)) + ', \alpha_y=' + string(param(2)))
set(gca,'YDir','normal','FontSize',20)
colorbar;
clim([0,2.5]);

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

function result_IDD = main_IDD(m,n,param)
    % Constants
    ax = param(1); ay = param(2); r = 0.1; u=0.1; D=1;
    r_x = repmat(1,1,m-1);
    r_y = repmat(0.1,1,n-1);
    r_ = [1,r_x, 0.1, r_y, r];
    
    % Variables
    D_list = logspace(-1,1,21);
    %u_list = logspace(-3,0,16)*5;
    Dose_list = D_list;
    
    parfor j=1:length(Dose_list)
        param = [ax,ay,m,n,r_,u,D_list(j)];
        %param = [a,b,m,n,r_,0.5,D_list(j)];
        Y0 = zeros(1,m+n+1);
        %tRange = linspace(0,200,1000);
        %opts = odeset('AbsTol',1e-5,'MaxStep',0.01);
        tRange = linspace(0,200,201);
        opts = odeset('AbsTol',1e-4,'MaxStep',min(Dose_list(j),1));
        [tSol,YSol] = ode23s(@(t,Y)model_1_m_n(t,Y,param),tRange,Y0,opts);
        I = calculate_AUC(tSol,YSol(:,end));
        T = calculate_T(tSol,YSol(:,end));
        Q = calculate_Q(tSol,YSol(:,end));
        AUC(j) = I;
        signaling_time(j) = T/I; % tau
        signal_duration(j) = sqrt(Q/I-signaling_time(j)^2); % theta
    end    
    result_IDD = evaluate_IDD(Dose_list,AUC);
    %result_peak(k,j) = detect_peak(D_list,AUC);
end

function dYdt=NDmodel_NL(t,Y,param)
    x=Y(1); y=Y(2); z=Y(3);
    %a=100; b=100; r=param(1); D=param(2);
    a=100; b=100; r=param(1); D=param(2);
    u=IFFL_input(t,1,D);
    
    dxdt = u*(1-x)-x;
    dydt = (x*(1-y)-y)*r;
    dzdt = (a*x*(1-z)-b*y*z-z)*r;
    dYdt=[dxdt;dydt;dzdt];
end

function dYdt = model_1_m_n(t, Y, param)
    ax = param(1); ay = param(2); m = param(3); n = param(4); r = param(5:5+m+n)'; u = param(end-1); D = param(end);
    Kux = 1; Kuy = 1; Kx = 1; Ky = 1; Kpx = 1; Kpy = 1; Kpz = 1; Kxz = 1; Kyz = 1; P = 1;
    a=100; b=10;
    u = IFFL_input(t, u, D);
    x = Y(1:m); y = Y(m+1:m+n); z = Y(m+n+1);
    
    dx(1) = u .* (1-x(1)) ./ (Kux+u+(1-x(1))) - r(1) * P .* x(1) ./ (Kpx+P+x(1));
    dy(1) = u .* (1-y(1)) ./ (Kuy+u+(1-y(1))) - r(m+1) * P .* y(1) ./ (Kpx+P+y(1));

    if m > 1
        dx(2:m) = ax .* x(1:m-1) .* (1-x(2:m)) ./ (Kx+x(1:m-1)+(1-x(2:m))) - r(2:m) .* P .* x(2:m) ./ (Kpx+P+x(2:m));
    end

    if n > 1
        dy(2:n) = ay .* y(1:n-1) .* (1-y(2:n)) ./ (Ky+y(1:n-1)+(1-y(2:n))) - r(m+2:m+n) .* P .* y(2:n) ./ (Kpy+P+y(2:n));
    end

    dz = a .* x(m) .* (1-z) ./ (Kxz+x(m)+(1-z)) - b .* y(n) .* z ./ (Kyz+y(n)+z) - r(end) .* 0.2 .* z ./ (Kpz+0.2+z);

    dYdt = [dx'; dy'; dz];
end

function t_10 = calculate_t_10(t,y)
    y_max = max(y);
    idx = y>(y_max*0.10);
    t_high = t(idx);
    %t_10 = t_high(2)-t_high(1);
    t_10 = max(t_high)-min(t_high);
end

function Y = value_at_T(T,t,y)
    idx_l = T>t;
    idx_h = T<=t;
    Y_l = y(idx_l);
    Y_h = y(idx_h);
    Y = (Y_l(end)+Y_h(1))/2; 
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

function IDD = test_IDD(D_list,theta)
    IDD = theta(1) > theta(end);
    %IDD = 1;
    %for i=1:length(D_list)-1
    %    IDD = IDD * (theta(i+1)-theta(i)<0);
    %end
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

function peak = detect_peak(D_list,theta)
    [theta_max,ind] = max(theta);
    peak = 1;
    peak = peak * ((ind~=1)&(ind~=length(D_list)));
    for i=1:ind-1
        peak = peak*(theta(i+1)-theta(i)>0);
    end
    for i=ind:length(D_list)-1
        peak = peak*(theta(i+1)-theta(i)<0);
    end    
end