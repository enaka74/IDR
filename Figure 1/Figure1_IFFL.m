%{
Define constant parameters
a=; b=; r=; kx= ; ky=; kz=; u=; T=;  
%}

var_num = 5;
COLOR = [41/256 128/256 185/256; 241/256 196/256 15/256];
R_v = linspace(COLOR(1,1),COLOR(2,1),var_num)';
B_v = linspace(COLOR(1,2),COLOR(2,2),var_num)';
G_v = linspace(COLOR(1,3),COLOR(2,3),var_num)';
Color_RGBs = [R_v B_v G_v];

%b_list = linspace(1,10,10);
%D_list = linspace(1,20,var_num)*0.05;
%D_list = linspace(1,10,var_num)*0.1;
D_list = [0.2, 0.4, 0.6, 0.8, 1];
%D_list = linspace(1,2*var_num-1,var_num)*600*1e-4; %matchig for dimensional model on 20230121
%r_list = linspace(4,20,5);
%T_sense = 20;

figure()
hold on   
for i=1:length(D_list)
    D = D_list(i);
    u = 1;
    %param = [1.185091394	93.71243808	0.108495406	2 D];
    param = [2	20	0.1	2 u D];
    Y0 = [0,0];
    %tRange = [0,28.8]; matching for dimensional model on 20230121
    tRange = linspace(0,80,1000);
    opts = odeset('AbsTol',1e-5,'MaxStep',0.01);
    [tSol(:,i),YSol(:,:,i)] = ode23s(@(t,Y)simple_IFFL_nd(t,Y,param),tRange,Y0,opts);
    plot(tSol(:,i),YSol(:,2,i),'LineWidth',3,'LineStyle','-','Color',Color_RGBs(i,:),'DisplayName',['D_{in}='+string(D)]) %v
    xlabel('Time') 
    ylabel('Concentration');
    set(gca,'FontSize',20);
    box on
    input=[];

    I = calculate_AUC(tSol,YSol(:,end));
    T = calculate_T(tSol,YSol(:,end));
    Q = calculate_Q(tSol,YSol(:,end));
    AUC(i) = I;
    signaling_time(i) = T/I; % tau
    signal_duration(i) = sqrt(Q/I-(T/I)^2); % theta
    Amp(i) = max(YSol(:,end));

    lgd = legend;
    lgd.NumColumns = 1;
    lgd.FontSize = 10;
    legend
    
end

xlim([0,20])
hold off

fig = figure();
hold on
yyaxis left
plot(D_list,signal_duration,"Color",[0,0,0],"LineStyle","-","LineWidth",2)
xlabel('D_{in}')
ylabel('Output duration D_{out}')
%ylim([3,9])
ylim([3,10])

yyaxis right
plot(D_list,AUC,"Color",[0,0,0],"LineStyle","--","LineWidth",2)
xlabel('D_{in}')
ylabel('Output AUC')
ylim([0.5,3])
set(gca,'FontSize',20);
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
box on
     
legend('D_{out}','AUC')
hold off

time_point = linspace(0,10,21);
for i=[1,length(D_list)]
    figure()
    tiledlayout(4,5);
    %title("D = "num2str(D_list(i))
    for j=1:length(time_point)-1
        nexttile
        hold on    
        [x,y]=locate_point(tSol(:,i),YSol(:,:,i),time_point(j));
        x_dotzero = x_dot_zero(linspace(0,1,length(YSol(:,2,i))),param);
        y_dotzero = y_dot_zero(linspace(0,1,length(YSol(:,1,i))),param);
        plot(YSol(:,1,i),YSol(:,2,i),'LineWidth',3,'LineStyle','-','DisplayName',['D='+string(D_list(i))],'Color',[0.7,0.7,0.7]) %v
        plot(linspace(0,1,length(YSol(:,1,i))),y_dotzero,'LineWidth',3,'LineStyle','--','DisplayName','y\_dot=0');
        plot(x_dotzero,linspace(0,1,length(YSol(:,2,i))),'LineWidth',3,'LineStyle',':','DisplayName','x\_dot=0');
        scatter(x,y,50,[0,0,0],"filled");
        xlabel('x');
        ylabel('y');
        %xlim([0,max([YSol(:,1,i);YSol(:,2,i)])]);
        %ylim([0,max([YSol(:,1,i);YSol(:,2,i)])]);
        xlim([0,1])
        ylim([0,1]);
        set(gca,'FontSize',20);
        box on
        
        lgd = legend;
        lgd.NumColumns = 1;
        lgd.FontSize = 10;
        %legend
        
        hold off
    end
    
end

function dYdt = simple_IFFL_nd(t,Y,param)
    a=param(1); b=param(2); c=param(3); p=param(4); u=param(5); D=param(6);
    u = IFFL_input(t,u,D);
    x=Y(1); y=Y(2);
    
    dxdt = u - x;
    dydt = a*u - b*(x^(p))*y - c*y;

    dYdt = [dxdt;dydt];
end

function [Px,Py] = locate_point(t,v,t0)
    x = v(:,1); y = v(:,2);
    [~,t_idx] = min(abs(t-t0));
    Px = x(t_idx);
    Py = y(t_idx);
    %P = [x(t_idx),y(t_idx)];
end

function y_eq = y_dot_zero(x,param)
    a=param(1); b=param(2); c=param(3); p=param(4); u=param(5); D=param(6);
    y_eq = a*u./(c+b*x.^p);
end

function x_eq = x_dot_zero(y,param)
    a=param(1); b=param(2); c=param(3); p=param(4); u=param(5); D=param(6);
    x_eq = repmat(u,1,length(y));
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