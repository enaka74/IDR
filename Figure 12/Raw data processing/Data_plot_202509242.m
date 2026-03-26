function [t_data,ON_data_1_,ON_data_2_,ON_data_3_]  = Data_plot_202509242(num_total,Sample_number,t_ctr_min,F_ctr_min,t_ctr_baseline,F_ctr_baseline)

%Color_RGBs = ["#0072BD";"#D95319";"#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F"];    
var_num = num_total;
Blu = [26, 82, 118; 234, 242, 248]/256;
R_v = linspace(Blu(1,1),Blu(2,1),var_num)'; B_v = linspace(Blu(1,2),Blu(2,2),var_num)'; G_v = linspace(Blu(1,3),Blu(2,3),var_num)';
Color_RGBs_1 = [R_v B_v G_v]; % blue

Red = [148, 20, 20; 253, 200, 200]/256;
R_v = linspace(Red(1,1),Red(2,1),var_num)'; B_v = linspace(Red(1,2),Red(2,2),var_num)'; G_v = linspace(Red(1,3),Red(2,3),var_num)';
Color_RGBs_2 = [R_v B_v G_v]; % red

Yel = [230, 126, 34; 252, 243, 207]/256;
R_v = linspace(Yel(1,1),Yel(2,1),var_num)'; B_v = linspace(Yel(1,2),Yel(2,2),var_num)'; G_v = linspace(Yel(1,3),Yel(2,3),var_num)';
Color_RGBs_3 = [R_v B_v G_v]; % yellow

sz = 30;
lw = 1.5;

% Inputs from user
number_of_plots = 1;
number_of_channels = 3;
S = double.empty(number_of_plots,0);
date = "25-09-24_2";

group1 = "G8C1"; group2 = "G3R1"; group3 = "G1S1"; 

dV = [2.200, 2.4, 3.5;
      2.200, 2.4, 3.5;
      2.200, 2.4, 3.5];


V_final = 60;

time_pts = [1,10,25,30;
            1,10,25,30;
            1,10,25,30]; % each time point indicates the first time point after adding each reagent

% Visualize each plot

% Import fluorometer data (Reporter node)

S = readmatrix(date + ".txt");
t_data_1 = []; t_data_2 = []; t_data_3 = [];
F_data_1 = []; F_data_2 = []; F_data_3 = [];
for i=1:number_of_plots
    t_data_2(:,i) = S(3:end,18*Sample_number(i)-4)/60/1000/60; % time [min]
    t_data_1(:,i) = S(3:end,18*Sample_number(i)-2)/60/1000/60; % time [min]
    t_data_3(:,i) = S(3:end,18*Sample_number(i))/60/1000/60; % time [min]
    F_data_2(:,i) = S(3:end,18*Sample_number(i)-3); % Fluorescence signal
    F_data_1(:,i) = S(3:end,18*Sample_number(i)-1); % Fluorescence signal
    F_data_3(:,i) = S(3:end,18*Sample_number(i)+1); % Fluorescence signal
end 

time_pts = cat(2,time_pts, repmat(length(t_data_1(:,1)),size(time_pts,1),1));

%cat(1,event_time_pts,repmat(size(F_data,1),1,size(event_time_pts,2)));
%F_data_1_temp = F_data_1;

% Adjust volume change
for i=1:number_of_plots
    for j=1:number_of_channels
        for k=1:length(time_pts)-1
            eval(['F_data_' num2str(j) '(time_pts(Sample_number(i),k)+1:time_pts(Sample_number(i),k+1),i) = (V_final-sum(dV(Sample_number(i),k:end)))/V_final * F_data_' num2str(j) '(time_pts(Sample_number(i),k)+1:time_pts(Sample_number(i),k+1),i);'])
        end

        % Cut initial phase before adding dA
        eval(['F_data_' num2str(j) ' = F_data_' num2str(j) '(time_pts(Sample_number(i),2):time_pts(Sample_number(i),2)+240,i);'])
    
    end
end

%time_pts_new = time_pts(:,2:end) - (time_pts(:,2) - 1) + [0,0,1,0;0,0,1,0;0,0,1,0];
time_pts_new = time_pts(:,2:end) - (time_pts(:,2) - 1) + repmat([0,0,1,0],1,number_of_plots);
%ref_time_pts = [1,3,3];
ref_time_pts = [1,1;
                3,3;
                1,3];

%%% Convert fluorescence intensity to Fraction ON
for i=1:number_of_plots
    for j=1:number_of_channels
        eval(['F_baseline(:,j) = F_ctr_baseline(1:length(F_data_' num2str(j) '(:,i)),j);'])
        eval(['F_min(:,j) = F_ctr_min(1:length(F_data_' num2str(j) '(:,1)),j);'])
        eval(['ON_data_' num2str(j) '(:,i) = normalize_decaying_baseline(F_data_' num2str(j) '(:,i),F_baseline(:,j),F_min(:,j),time_pts_new(i,ref_time_pts(j,:)));'])
        %eval(['ON_data_' num2str(j) '(:,i) = (F_baseline(:,j) - F_data_' num2str(j) '(:,i))./(F_baseline(:,j) - F_min(:,j));'])
    end
end

%%% Normalize all the data from 0 min (T7 addition) to 240 min %%%

for i=1:number_of_plots
    t_data_temp(:,i) = t_data_3(time_pts_new(i,3):end,i)-t_data_3(time_pts_new(i,3),i);
    t_data(:,i) = t_data_temp(1:121,i);
    ON_data_1_(:,i) = ON_data_1(time_pts_new(i,3):time_pts_new(i,3)+120,i);
    ON_data_2_(:,i) = ON_data_2(time_pts_new(i,3):time_pts_new(i,3)+120,i);
    ON_data_3_(:,i) = ON_data_3(time_pts_new(i,3):time_pts_new(i,3)+120,i);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% End of main function %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y_normed = normalize_flat_baseline(y,y_min)
    y_max = y(1);
    y_normed = (y_max - y)/(y_max - y_min);
end

function y_normed = normalize_decaying_baseline(y,y_baseline,y_min,time_pt_ref)
    y_min = y(time_pt_ref(1))/y_baseline(time_pt_ref(1)) * y_min;
    y_max = y(time_pt_ref(2))/y_baseline(time_pt_ref(2)) * y_baseline;

    y_normed = (y_max - y)./(y_max - y_min);
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