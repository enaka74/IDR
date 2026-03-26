clear all
close all

Colors_red = ["#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026"];

Colors_green = ["#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45","#005a32"];

Colors_blue= ["#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594"];

dates = [20250305,20250306,20250307;  % 10 nM
         20250305,20250306,20250307;  % 25 nM
         20250305,20250306,20250307]; % 50 nM

input_concs = [10,25,50]; % input input_concs

concs = [10,25,50]; % input concs

Sample_numbers = [1,1,1;
                  2,2,2;
                  3,3,3];

num_total=6; 

Sample_labels = ["10 nM G8C1 / 2 nM G8C3","25 nM G8C1 / 5 nM G8C3","50 nM G8C1 / 10 nM G8C3"];


replicates = ["Replicate #1","Replicate #2","Replicate #3"];

% Define control sample

ctr_date = "24-07-11";
ctr_min_sample = 1;
ctr_baseline_sample = 4;
number_of_plots = 4;
number_of_channels = 3;
Sample_number = [1,2,3,4];

time_pts = [1,15,39,44;
            1,16,40,45;
            1,16,40,45;
            1,15,40,45];

dV = [2.762, 0, 0;
      1.908, 0, 0;
      1.908, 0, 0;
      0,     0, 0];

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

time_pts = cat(2,time_pts, repmat(length(t_data_1(:,1)),number_of_plots,1));

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

f=figure(1); tiledlayout(6,3); f.Position = [0,0,1000,1000];

[m, n] = size(dates);

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
    y_mean_1 = mean(ON_data_1, 2);       % mean for each column (1×3 row vector)
    y_std_1  = std(ON_data_1, 0, 2);     % standard deviation for each column (1×3 row vector)
    err_1 = [-y_std_1,y_std_1];
    
    y_mean_2 = mean(ON_data_2, 2);       % mean for each column (1×3 row vector)
    y_std_2  = std(ON_data_2, 0, 2);     % standard deviation for each column (1×3 row vector)
    err_2 = [-y_std_2,y_std_2];
    
    y_mean_3 = mean(ON_data_3, 2);       % mean for each column (1×3 row vector)
    y_std_3  = std(ON_data_3, 0, 2);     % standard deviation for each column (1×3 row vector)
    err_3 = [-y_std_3,y_std_3];
    
    % shadedErrorBar(x, y,y_std, {'-b','LineWidth',2}, 0.2);
    
    y = bsxfun(@plus,err_3,y_mean_3); % Add the sine wave to the random data and offset the data
    % Plot with shaded error bar
    
    nexttile(3*i-2)
    hold on
    plot(x, ON_data_1(:,1), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    plot(x, ON_data_1(:,2), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    plot(x, ON_data_1(:,3), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    xlim([0,4])
    %ylim([0,1])
    xlabel('Time [h]')
    ylabel('Active fraction')
    box on;
    set(gca,'FontSize',20);
    hold off
    
    nexttile(3*i-1)
    hold on;   
    plot(x, ON_data_2(:,1), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    plot(x, ON_data_2(:,2), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    plot(x, ON_data_2(:,3), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    xlim([0,4])
    ylim([0,1])
    xlabel('Time [h]')
    ylabel('Active fraction')
    box on;
    set(gca,'FontSize',20);
    hold off;
    
    nexttile(3*i)
    hold on;   
    plot(x, ON_data_3(:,1), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    plot(x, ON_data_3(:,2), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
    plot(x, ON_data_3(:,3), '-', 'LineWidth', 3)%, 'Color',[41, 128, 185]/256);
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