%{
1. Vary several parameters and calculate signal duration
    - For simplicity, use the same parameter value along each cascade
    - Parameters at the terminal node (alpha and gamma) are varied

2. Create following plots
    2.1 Plots for all the nodes
    2.2 Several plots for different input duration
    2.3 n (length of y cascade) vs signal duration
    2.4 m (length of x cascade) vs signal duration (with constant n)
    2.5 Signal duration on n-m plane
    2.6 Parameter value (alpha, gamma) vs AUC (color map)
    
%}

%%%

red = [120, 40, 31; 148, 49, 38; 203, 67, 53; 231, 76, 60; 236, 112, 99; 241, 148, 138]/256;
blue = [36, 113, 163; 41, 128, 185; 84, 153, 199; 127, 179, 213; 169, 204, 227; 212, 230, 241]/256;
green = [17, 122, 101; 19, 141, 117; 22, 160, 133; 69, 179, 157; 115, 198, 182; 162, 217, 206]/256;
yellow = [185, 119, 14; 214, 137, 16; 243, 156, 18; 245, 176, 65; 248, 196, 113; 250, 215, 160]/256;
gray = [97, 106, 107; 112, 123, 124; 127, 140, 141; 153, 163, 164; 178, 186, 187; 204, 209, 209]/256;


%%% Figure setting %%%

var_num = 20;
COLOR = [41/256 128/256 185/256; 241/256 196/256 15/256];
R_v = linspace(COLOR(1,1),COLOR(2,1),var_num)';
B_v = linspace(COLOR(1,2),COLOR(2,2),var_num)';
G_v = linspace(COLOR(1,3),COLOR(2,3),var_num)';
Color_RGBs = [R_v B_v G_v];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 1: n (length of y cascade) vs signal duration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
m = 6; n= 1; 

a = [2,2,2,2,2,2,   2,    100]*1; 
b = 10; 
K = [1,1,1,1,1,1,   1,    1,1];  % K values for y controls threshold of dose
Kp = 1;
Kpz = 1;
%r = [0.1,0.1,0.1,0.1,0.1,0.1,   0.1,0.1,0.1,0.1,0.1,0.1,    0.1];
r = [1,1,1,1,1,1,   0.1,  0.1];
P = 1	;	
Pz = 0.2 ;

u=0.1; D=1;

% Variables
%D_list = linspace(1,20,20)*0.15;
%u_list = linspace(2,40,20)*0.010;
D_list = logspace(-1,1,21);
%u_list = logspace(-4,1,26)*4;
Dose_list = D_list;

% Graph Setting 
COLOR = [41/256 128/256 185/256; 241/256 196/256 15/256];
R_v = linspace(COLOR(1,1),COLOR(2,1),m+n+1)';
B_v = linspace(COLOR(1,2),COLOR(2,2),m+n+1)';
G_v = linspace(COLOR(1,3),COLOR(2,3),m+n+1)';
Color_RGBs = [R_v B_v G_v];

f = figure(1);
f.Position = [0,0,800,600];
tiledlayout(3,2);

x_max = 0;
y_max = 0;
z_max = 0;

%idx_list = [1,17];
idx_list = [1,13];
for i=1:length(idx_list)
    %param = [m,n,a,b,K,Kp,Kpz,r,P,Pz,Dose_list(idx_list(i)),D];
    param = [m,n,a,b,K,Kp,Kpz,r,P,Pz,u,Dose_list(idx_list(i))];
    Y0 = zeros(1,m+n+1);
    %Y0 = repmat(0,1,m+n+1);
    tRange = [0,100];
    opts = odeset('AbsTol',1e-5,'MaxStep',0.1);
    [tSol,YSol] = ode23s(@(t,Y)model_1_m_n(t,Y,param),tRange,Y0,opts);
    x_max = max([x_max,max(YSol(:,1:m),[],"all")]);
    y_max = max([y_max,max(YSol(:,m+1:m+n),[],"all")]);
    z_max = max([z_max,max(YSol(:,m+n+1),[],"all")]);

    % plot x
    nexttile(i);
    hold on
    for j=1:m
        plot(tSol,YSol(:,j),'LineWidth',4,'LineStyle','-','Color',red(1,:)+(red(end,:)-red(1,:))*j/m) %v
    end
    hold off
    box on

    % plot y
    nexttile(2+i);
    hold on
    for j=1:n
        plot(tSol,YSol(:,m+j),'LineWidth',4,'LineStyle','-','Color',blue(1,:)+(blue(end,:)-blue(1,:))*j/n) %v
    end
    hold off
    box on

    % plot z
    nexttile(4+i);
    hold on
    plot(tSol,YSol(:,end),'LineWidth',4,'LineStyle','-','Color',gray(1,:)) %v
    hold off
    box on

end

length_list = [m,n,1];
label_list = ['x','y','z'];
lim_list = [x_max,y_max,z_max];
for j=1:3
    nexttile(2*(j-1)+1);
    legend(label_list(j)+string(linspace(1,length_list(j),length_list(j))));
    xlabel('Time') 
    ylabel('Concentration'); %ylim([0,1]);
    xlim(tRange)
    ylim([0,lim_list(j)])
    %ylim([0,inf])
    if j==1
        %title('D_{in} = '+string(Dose_list(1)*D))
    end
    set(gca,'FontSize',20);
    lgd = legend; lgd.FontSize = 15;
    legend

    nexttile(2*(j-1)+2);
    legend(label_list(j)+string(linspace(1,length_list(j),length_list(j))));
    xlabel('Time') 
    ylabel('Concentration'); %ylim([0,1]);
    xlim(tRange)
    ylim([0,lim_list(j)])
    %ylim([0,inf])
    if j==1
        %title('D_{in} = '+string(Dose_list(end)*D))
    end
    set(gca,'FontSize',20);
    lgd = legend; lgd.FontSize = 15;
    legend
end

set(gcf, 'Renderer', 'painters');
print('filename.eps', '-depsc', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 2: Several plots for different input duration %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants

% Variables

% Graph Setting 
var_num = length(Dose_list);
COLOR = [41/256 128/256 185/256; 241/256 196/256 15/256];
R_v = linspace(COLOR(1,1),COLOR(2,1),var_num)';
B_v = linspace(COLOR(1,2),COLOR(2,2),var_num)';
G_v = linspace(COLOR(1,3),COLOR(2,3),var_num)';
Color_RGBs = [R_v B_v G_v];

f = figure(2);
tiledlayout(1,4,"TileSpacing","loose");

nexttile(1);
hold on
parfor i=1:length(Dose_list)
    param = [m,n,a,b,K,Kp,Kpz,r,P,Pz,u,Dose_list(i)];
    %param = [m,n,a,b,K,Kp,Kpz,r,P,Pz,Dose_list(i),D];
    Y0 = zeros(1,m+n+1);
    %Y0 = repmat(0,1,m+n+1);
    tRange = linspace(0,200,201);
    opts = odeset('AbsTol',1e-6,'MaxStep',0.1);
    tic
    [tSol,YSol] = ode23s(@(t,Y)model_1_m_n(t,Y,param),tRange,Y0,opts);
    calc_time(i) = toc
    %{
    I_x = calculate_AUC(tSol,YSol(:,end));
    T_x = calculate_T(tSol,YSol(:,end));
    Q_x = calculate_Q(tSol,YSol(:,end));
    signaling_time_x(i) = T/I; % tau
    signal_duration_x(i) = sqrt(Q/I-signaling_time(i)^2); % theta
    
    I_y = calculate_AUC(tSol,YSol(:,end));
    T_y = calculate_T(tSol,YSol(:,end));
    Q_y = calculate_Q(tSol,YSol(:,end));
    signaling_time_y(i) = T/I; % tau
    signal_duration_y(i) = sqrt(Q/I-signaling_time(i)^2); % theta
    %}
    I = calculate_AUC(tSol,YSol(:,end));
    T = calculate_T(tSol,YSol(:,end));
    Q = calculate_Q(tSol,YSol(:,end));
    AUC(i) = I;
    signaling_time(i) = T/I; % tau
    signal_duration(i) = sqrt(Q/I-(T/I)^2); % theta
    t_10(i) = calculate_t_10(tSol,YSol(:,end));
    plot(tSol,YSol(:,end),'LineWidth',2,'LineStyle','-','Color',Color_RGBs(i,:)) %v
    %plot(tSol,YSol(:,end)/max(YSol(:,end)),'LineWidth',2,'LineStyle','-','Color',Color_RGBs(i,:)) %v
end    
hold off

result_IDD = evaluate_IDD(Dose_list,signal_duration)
result_peak = detect_peak(Dose_list,AUC)
legend('D='+string(flip(Dose_list)));
xlabel('Time') 
ylabel('Concentration'); %ylim([0,1]);
set(gca,'FontSize',20);
lgd = legend; lgd.NumColumns = 2; lgd.FontSize = 15;
legend

nexttile(2);
plot(Dose_list,signaling_time,'LineWidth',3,'LineStyle','-','Color','black')
%xlabel('Input duration')
xlabel('Dose')
ylabel('Signaling time')
title('Signaling time')
ylim([0,inf]);
set(gca,'FontSize',20);

nexttile(3);
hold on
%plot(Dose_list,signal_duration,'LineWidth',3,'LineStyle','-','Color','black')
%plot(Dose_list,t_10,'LineWidth',3,'LineStyle',':','Color','black')
semilogx(Dose_list,signal_duration,'LineWidth',3,'LineStyle','-','Color','black')
semilogx(Dose_list,t_10,'LineWidth',3,'LineStyle',':','Color','black')
%xlabel('Input duration')
xlabel('Dose')
ylabel('Output duration')
title('Output duration')
legend('\theta','t_{10}')
ylim([0,inf]);
set(gca,'FontSize',20);
grid on
hold off

nexttile(4);
plot(Dose_list,AUC,'LineWidth',3,'LineStyle','-','Color','black')
%semilogx(log10(Dose_list),AUC,'LineWidth',3,'LineStyle','-','Color','black')
%xlabel('Input duration')
xlabel('Dose')
ylabel('AUC')
title('AUC')
ylim([0,inf]);
set(gca,'FontSize',20);


function dYdt=NDmodel_L(t,Y,param)
    x=Y(1); y=Y(2); z=Y(3);
    %a=100; b=100; r=param(1); D=param(2);
    a=100; b=100; r=param(1); D=param(2);
    u=IFFL_input(t,1,D);

    dxdt = u-x;
    dydt = (x-y)*r;
    dzdt = (a*x-b*y-z)*r;
    dYdt=[dxdt;dydt;dzdt];
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

function dYdt = model_1_m_n(t,Y,param)
    m = param(1); n = param(2);
    
    a = param(3:3+m+n)'; b = param(4+m+n); 
    K = param(5+m+n:6+2*m+2*n)'; Kp = param(7+2*m+2*n); Kpz = param(8+2*m+2*n);
    r = param(9+2*m+2*n:9+3*m+3*n)'; P = param(10+3*m+3*n); Pz = param(11+3*m+3*n);
    u = param(12+3*m+3*n); D = param(13+3*m+3*n);
    
    u = IFFL_input(t, u, D);
    x = Y(1:m); y = Y(m+1:m+n); z = Y(m+n+1);
    
    dx(1) = a(1) * u .* (1-x(1)) ./ (K(1)+u+(1-x(1))) - r(1) * P .* x(1) ./ (Kp+P+x(1));
    dy(1) = a(m+1) * u .* (1-y(1)) ./ (K(m+1)+u+(1-y(1))) - r(m+1) * P .* y(1) ./ (Kp+P+y(1));

    if m > 1
        dx(2:m) = a(2:m) .* x(1:m-1) .* (1-x(2:m)) ./ (K(2:m)+x(1:m-1)+(1-x(2:m))) - r(2:m) .* P .* x(2:m) ./ (Kp+P+x(2:m));
    end

    if n > 1
        dy(2:n) = a(m+2:m+n) .* y(1:n-1) .* (1-y(2:n)) ./ (K(m+2:m+n)+y(1:n-1)+(1-y(2:n))) - r(m+2:m+n) .* P .* y(2:n) ./ (Kp+P+y(2:n));
    end
    
    dz = a(m+n+1) .* x(m) .* (1-z) ./ (K(m+n+1)+x(m)+(1-z)) - b .* y(n) .* z ./ (K(m+n+2)+y(n)+z) - r(end) .* Pz .* z ./ (Kpz+Pz+z);

    dYdt = [dx'; dy'; dz];
end

function dYdt=Dmodel(t,Y)
    x=Y(1); y=Y(2); z=Y(3);
    a=param(1); b=param(2); r=param(3); kx=param(4); ky=param(5);
    u=param(6); T=param(7);
    u = IFFL_input(t,u,T);
    
    dxdt=a*u-kx*x;
    dydt=b*x-ky*y;
    dzdt=r*x*(1-z)-r*y*z-kz*z;
    
    dYdt=[dxdt;dydt;dzdt];
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
            IDD(count) = max_min_list(j,1)/min(max_min_list(j:end,2))-1;
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