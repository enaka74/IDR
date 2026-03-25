%{
Define constant parameters
a=; b=; r=; kx= ; ky=; kz=; u=; T=;  
%}

var_num = 3;
color_num = 2;
COLOR(:,:,1) = [245/256 183/256 177/256; 203/256 67/256 53/256];
COLOR(:,:,2) = [174/256 214/256 241/256; 46/256 134/256 193/256];
Color_RGBs=[];
for i=1:color_num
    R_v = linspace(COLOR(1,1,i),COLOR(2,1,i),var_num)';
    B_v = linspace(COLOR(1,2,i),COLOR(2,2,i),var_num)';
    G_v = linspace(COLOR(1,3,i),COLOR(2,3,i),var_num)';
    Color_RGBs(:,:,i) = [R_v B_v G_v];
end

R = linspace(1,203/256,256); G = linspace(1,67/256,256); B = linspace(1,53/256,256);
map = [R' G' B'];

height = 430;
width = 1350;

f = figure;
f.Position= [0, 0, width, height];
t=tiledlayout(1,3);
t.TileSpacing = 'compact';

%%% parameters %%%
u=1; a_0=2; b_0=20; c_0=0.1; p_0=2;
step=0.05; t_end=20;
%%% Fig a,b %%%

p_list = [1,2,3];

for p=p_list

AUC = [];
signaling_time = []; % tau
signal_duration = []; % theta
Amp = [];      

nexttile();
D_list = linspace(1,20,20)*step;
u_list = linspace(0.5,10,20);
for i=1:length(D_list)
    D = D_list(i);
    for j=1:length(u_list)
        u = u_list(j);
        param = [a_0, b_0, c_0, p, u, D];
        Y0 = [0,0];
        %tRange = [0,28.8]; matching for dimensional model on 20230121
        tRange = linspace(0,t_end,8000);
        opts = odeset('AbsTol',1e-6,'MaxStep',0.01);
        [tSol,YSol] = ode23s(@(t,Y)simple_NF_nd(t,Y,param),tRange,Y0,opts);

        I = calculate_AUC(tSol,YSol(:,end));
        T = calculate_T(tSol,YSol(:,end));
        Q = calculate_Q(tSol,YSol(:,end));
        AUC(j,i) = I;
        signaling_time(j,i) = T/I; % tau
        signal_duration(j,i) = sqrt(Q/I-(T/I)^2); % theta
        Amp(j,i) = max(YSol(:,end));    
    end
    %title('r = '+string(r_list(i)));
    lgd = legend;
    lgd.NumColumns = 1;
    legend('FontSize',13)
    
end
imagesc(D_list,u_list,AUC);
%xlabel('Pulse duration') 
xlabel('D_{in}') 
%yticks([0.1 1 10]*a_0)
title('\theta')
set(gca,'YDir','normal','FontSize',20)
%set(gca,'YScale','log')
ylabel('u','FontSize',30);
ylim([u_list(1),u_list(end)]);
colormap(map)
colorbar;

end



function dYdt = simple_NF_nd(t,Y,param)
    a=param(1); b=param(2); c=param(3); p=param(4); u=param(5); D=param(6);
    u = IFFL_input(t,u,D);
    x=Y(1); y=Y(2);
    
    dxdt = y - x;
    dydt = a*u - b*(x^(p))*y - c*y;

    dYdt = [dxdt;dydt];
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