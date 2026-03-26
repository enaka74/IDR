clear all
close all

Colors_red = ["#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026"];

Colors_green = ["#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45","#005a32"];

Colors_blue= ["#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594"];

Colors_RGB = ["#d73027","#fc8d59","#fee08b","#d9ef8b","#91cf60","#1a9850"];

Colors_div = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628"];

dates = [20250918,202509222,202509243; % 10min
         202509192,202509241,202510061;% 20min 202509251
         202509221,202509242,202509252; % 30min   % 202509221 will be replaced
         202509261,20250930,20251001; % 40min 
         202510062,20251007, 20260126; % 60min
         20260126, 20260126, 20260126]; % 240min 

durations = [10,20,30,40,60,240]; % input durations

Sample_numbers = [2,2,2;
                  2,2,2;
                  2,2,2;
                  2,2,2;
                  2,2,4;
                  1,2,3];

num_total=6; 

Sample_labels = ["D_{in}=10 min","D_{in}=20 min","D_{in}=30 min","D_{in}=40 min","D_{in}=60 min","D_{in}=240 min"];
Sample_labels_three = ["D_{in}=10 min","D_{in}=30 min","D_{in}=60 min"];

replicates = ["Replicate #1","Replicate #2","Replicate #3"];

% Define control sample

ctr_date = "25-11-13";
ctr_min_sample = 1;
ctr_baseline_sample = 2;
number_of_plots = 2;
number_of_channels = 3;
Sample_number = [1,2];

time_pts = [1,10,25,30;
            1,10,25,30;
            1,10,25,30];

dV = [1.626, 2.4, 0;
      0, 2.4, 0;
      0, 0, 0];

V_final = 60;

S = readmatrix(ctr_date + ".txt");
t_data_1 = []; t_data_2 = []; t_data_3 = [];
F_data_1 = []; F_data_2 = []; F_data_3 = [];
for i=1:number_of_plots
    t_data_2(:,i) = S(3:end,18*Sample_number(i)-4)/60/1000; % time [min]
    t_data_1(:,i) = S(3:end,18*Sample_number(i)-2)/60/1000; % time [min]
    t_data_3(:,i) = S(3:end,18*Sample_number(i))/60/1000; % time [min]
    F_data_2(:,i) = S(3:end,18*Sample_number(i)-3); % Fluorescence signal
    F_data_1(:,i) = S(3:end,18*Sample_number(i)-1); % Fluorescence signal
    F_data_3(:,i) = S(3:end,18*Sample_number(i)+1); % Fluorescence signal
end

time_pts = cat(2,time_pts, repmat(length(t_data_1(:,1)),number_of_channels,1));

for i=1:number_of_plots
    for j=1:3
        for k=1:length(time_pts)-1
            eval(['F_data_' num2str(j) '(time_pts(i,k):time_pts(i,k+1)-1,i) = (V_final-sum(dV(i,k:end)))/V_final*F_data_' num2str(j) '(time_pts(i,k):time_pts(i,k+1)-1,i);'])
        end
    end
end
for j=1:number_of_channels
    eval(['t_ctr_min(:,j) = t_data_' num2str(j) '(time_pts(ctr_min_sample,2):end,ctr_min_sample);'])
    eval(['F_ctr_min(:,j) = F_data_' num2str(j) '(time_pts(ctr_min_sample,2):end,ctr_min_sample);'])
    eval(['t_ctr_baseline(:,j) = t_data_' num2str(j) '(time_pts(ctr_baseline_sample,2):end,ctr_baseline_sample);'])
    eval(['F_ctr_baseline(:,j) = F_data_' num2str(j) '(time_pts(ctr_baseline_sample,2):end,ctr_baseline_sample);'])
    %
end

%


[m, n] = size(dates);
f=figure(1); tiledlayout(m,3); f.Position = [0,0,1200,800];

% === Added: explicit D_in values (minutes) to name CSV files ===
% This must align with the row order of `dates` (i = 1..m).
D_in_minutes = [10,20,30,40,60,240];

for i=1:m
    
    t_data = [];
    ON_data = [];
    
    figure(1);
    for j=1:n 
        d = dates(i,j);
        eval(['[t_data(:,j),ON_data_1(:,j),ON_data_2(:,j),ON_data_3(:,j)] = Data_plot_' num2str(d) '(n,Sample_numbers(i,j),t_ctr_min,F_ctr_min,t_ctr_baseline,F_ctr_baseline);'])
        AUC(i,j) = calculate_AUC(t_data(:,j),ON_data_3(:,j));
    end
    x = t_data(:,1);
    y_mean_1 = mean(ON_data_1, 2);       % mean across replicates for channel 1
    y_std_1  = std(ON_data_1, 0, 2);     % std across replicates for channel 1
    err_1 = [-y_std_1,y_std_1];
    
    y_mean_2 = mean(ON_data_2, 2);       % mean across replicates for channel 2
    y_std_2  = std(ON_data_2, 0, 2);     % std across replicates for channel 2
    err_2 = [-y_std_2,y_std_2];
    
    y_mean_3 = mean(ON_data_3, 2);       % mean across replicates for channel 3
    y_std_3  = std(ON_data_3, 0, 2);     % std across replicates for channel 3
    err_3 = [-y_std_3,y_std_3];
    
    % === Added: write per-D_in CSV (time in hours) ===
    % Writes: time_h, ON_mean_ch1, ON_mean_ch2, ON_mean_ch3
    % File name pattern: D_in_XXmin.csv (e.g., D_in_10min.csv)
    T_out = table(x, y_mean_1, y_mean_2, y_mean_3, ...
                  'VariableNames', {'time_h','ON_mean_ch1','ON_mean_ch2','ON_mean_ch3'});
    writetable(T_out, sprintf('Din_%dmin_means.csv', D_in_minutes(i)));

    % shadedErrorBar(x, y,y_std, {'-b','LineWidth',2}, 0.2);
    
    y = bsxfun(@plus,err_3,y_mean_3); % Add the sine wave to the random data and offset the data
    % Plot with shaded error bar
    
    nexttile(3*i-2)
    hold on
    for j=1:n
        plot(x, ON_data_1(:,j), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    end
    xlim([0,4])
    ylim([0,1])
    xlabel('Time [h]')
    ylabel('Active fraction')
    box on;
    set(gca,'FontSize',20);
    hold off
    
    nexttile(3*i-1)
    hold on;   
    for j=1:n
        plot(x, ON_data_2(:,j), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    end
    xlim([0,4])
    ylim([0,1])
    xlabel('Time [h]')
    ylabel('Active fraction')
    box on;
    set(gca,'FontSize',20);
    hold off;
    
    nexttile(3*i)
    hold on;   
    for j=1:n
        plot(x, ON_data_3(:,j), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    end
    xlim([0,4])
    ylim([0,1])
    xlabel('Time [h]')
    ylabel('Active fraction')
    box on;
    legend(replicates)
    set(gca,'FontSize',20);
    hold off;
    
    f=figure(7);
    f.Position = [0,0,700,300];
    hold on;   
    plot(x, y_mean_1, '-', 'LineWidth', 3, 'Color', Colors_blue(i));
    xlim([0,4])
    ylim([0,1])
    xlabel('Time [h]')
    ylabel('Active fraction')
    box on;
    set(gca,'FontSize',25);
    hold off;

    f=figure(8);
    f.Position = [0,0,700,300];
    hold on;   
    plot(x, y_mean_2, '-', 'LineWidth', 3, 'Color', Colors_green(i));
    xlim([0,4])
    ylim([0,1])
    xlabel('Time [h]')
    ylabel('Active fraction')
    box on;
    set(gca,'FontSize',25);
    hold off;

    f=figure(9);
    f.Position = [0,0,700,300];
    hold on;   
    plot(x, y_mean_3, '-', 'LineWidth', 3, 'Color', Colors_red(i));
    xlim([0,4])
    ylim([0,1])
    xlabel('Time [h]')
    ylabel('Active fraction')
    box on;
    set(gca,'FontSize',25);
    hold off;

end

figure(7)
legend(Sample_labels)
figure(8)
legend(Sample_labels)
figure(9)
legend(Sample_labels)

AUC_mean = mean(AUC,2);
AUC_std = std(AUC,0,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% End of main function %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = calculate_AUC(t,y)
    I = 0;
    for i=2:length(t)
        I = I + (t(i)-t(i-1))*(y(i));
    end
end

