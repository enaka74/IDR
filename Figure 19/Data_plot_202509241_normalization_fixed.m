%% ============================================================
%  plot_202509242_3x3_correctedSpec.m
%
%  3x3 tiles: Rows=channels, Cols=
%   (1) Raw: sample (solid), control baseline (dashed), control min (dotted)
%   (2) Reference-point scaled: sample y(t) (solid), y_max(t)=alpha*baseline(t) (dashed), min(t) (dotted)
%   (3) Fraction ON after enzyme: ON(t) aligned at T7+1 (t0), 0..4h
%
%  NOTE: time_pts values/logic are kept identical to Data_plot_202509242.
%% ============================================================

clear all
close all

%% ===== Choose replicate (ONE only) =====
replicate_sample_number = 2;  % 1 or 2
rep_label = sprintf('Replicate #%d', replicate_sample_number);

Channels = ["G_{in}C_{out} (TYE563)", ...
            "G_{r}R_{out} (TEX615)", ...
            "G_{out}S (6-FAM)"];

root_dir = 'Input pathname to the directory you saved 25-11-13.txt and "25-09-24_1.txt';

%% ===== Control settings (same as your driver) =====
ctr_date = "25-11-13";
ctr_min_sample = 1;
ctr_baseline_sample = 2;
Sample_number_ctr = [1,2];

number_of_channels = 3;
V_final = 60;

time_pts_ctr = [1,10,25,30;
                1,10,25,30;
                1,10,25,30];

dV_ctr = [1.626, 2.4, 0;
          0,     2.4, 0;
          0,     0,   0];

%% ============================================================
%  (A) Load & preprocess control (produce baseline/min after dA-cut)
% ============================================================
Sctr = readmatrix(strcat(root_dir, ctr_date, '.txt'));

for i=1:2
    t_ctr_raw_2(:,i) = Sctr(3:end,18*Sample_number_ctr(i)-4)/60/1000; % time (as in your driver)
    t_ctr_raw_1(:,i) = Sctr(3:end,18*Sample_number_ctr(i)-2)/60/1000;
    t_ctr_raw_3(:,i) = Sctr(3:end,18*Sample_number_ctr(i))/60/1000;

    F_ctr_raw_2(:,i) = Sctr(3:end,18*Sample_number_ctr(i)-3);
    F_ctr_raw_1(:,i) = Sctr(3:end,18*Sample_number_ctr(i)-1);
    F_ctr_raw_3(:,i) = Sctr(3:end,18*Sample_number_ctr(i)+1);
end

time_pts_ctr = cat(2, time_pts_ctr, repmat(length(F_ctr_raw_1(:,1)), 3, 1));

% Volume correction (control)
for i=1:2
    for ch=1:3
        for k=1:(size(time_pts_ctr,2)-1)
            eval(['F_ctr_raw_' num2str(ch) '(time_pts_ctr(ch,k):time_pts_ctr(ch,k+1)-1,i) = ' ...
                  '(V_final-sum(dV_ctr(ch,k:end)))/V_final * F_ctr_raw_' num2str(ch) '(time_pts_ctr(ch,k):time_pts_ctr(ch,k+1)-1,i);'])
        end
    end
end

% ---- FULL control traces (from true t0) ----
F_ctr_baseline_full = [F_ctr_raw_1(:,ctr_baseline_sample), ...
                       F_ctr_raw_2(:,ctr_baseline_sample), ...
                       F_ctr_raw_3(:,ctr_baseline_sample)];

F_ctr_min_full      = [F_ctr_raw_1(:,ctr_min_sample), ...
                       F_ctr_raw_2(:,ctr_min_sample), ...
                       F_ctr_raw_3(:,ctr_min_sample)];

% ---- dA-cut control traces (keep for legacy / debugging) ----
F_ctr_baseline_dA = F_ctr_baseline_full(time_pts_ctr(ctr_baseline_sample,2):end, :);
F_ctr_min_dA      = F_ctr_min_full(     time_pts_ctr(ctr_min_sample,2):end, :);


%% ============================================================
%  (B) Extract sample raw (dA-cut), and compute y_max/ON with SAME time_pts logic
% ============================================================
[t_raw, y_raw, b_raw, m_raw, y_max_raw, y_min_raw, ON_T7, t_T7] = ...
    extract_202509242_correctedSpec(replicate_sample_number, ...
                                    F_ctr_baseline_full, F_ctr_min_full, root_dir);

%% ============================================================
%  (C) Plot 3x3
% ============================================================
f = figure(1); clf;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
f.Position = [0,0,1500,900];

col1 = "Raw";
col2 = "Ref-point scaled (y, y_{max}=\alpha b, m)";
col3 = "Fraction ON (T7+1 aligned)";

for ch = 1:3

    % ----- Column 1: raw -----
    nexttile(3*(ch-1)+1);
    hold on
    plot(t_raw, y_raw(:,ch), 'k-',  'LineWidth', 2.5);
    plot(t_raw, b_raw(:,ch), 'k--', 'LineWidth', 2.0);
    plot(t_raw, m_raw(:,ch), 'k:',  'LineWidth', 2.0);
    hold off
    box on; set(gca,'FontSize',15,'LineWidth',2);
    xlabel('Time [h]'); ylabel('Fluorescence (a.u.)');
    xlim([0,5]);
    title(sprintf('%s\n%s', Channels(ch), col1), 'Interpreter','tex');
    if ch==1
        legend({'sample','ctrl baseline','ctrl min'}, 'Location','best'); legend boxoff;
    end

    % ----- Column 2: ref-point scaled 3 trajectories -----
    nexttile(3*(ch-1)+2);
    hold on
    plot(t_raw, y_raw(:,ch),     'k-',  'LineWidth', 2.5);
    plot(t_raw, y_max_raw(:,ch), 'k--', 'LineWidth', 2.0);
    plot(t_raw, y_min_raw(:,ch),     'k:',  'LineWidth', 2.0);
    hold off
    box on; set(gca,'FontSize',15,'LineWidth',2);
    ax=gca; ax.XColor='black'; ax.YColor='black'; 
    xlabel('Time [h]'); ylabel('Fluorescence-scaled (a.u.)');
    xlim([0,5]);
    title(col2, 'Interpreter','none');

    % ----- Column 3: ON aligned at T7+1 -----
    nexttile(3*(ch-1)+3);
    plot(t_T7, ON_T7(:,ch), 'k-', 'LineWidth', 2.5);
    box on; set(gca,'FontSize',15,'LineWidth',2);
    xlabel('Time after T7 [h]'); ylabel('Active fraction');
    xlim([0,4]); ylim([0,1]);
    title(col3, 'Interpreter','none');

end

sgtitle(sprintf('202509242 | %s', rep_label), 'Interpreter','none');

%% ============================================================
%  Helper: extract sample raw + compute y_max and ON (time_pts logic unchanged)
%% ============================================================
function [t_raw, y_raw, b_raw, m_raw, y_max_raw, y_min_raw, ON_T7, t_T7] = ...
    extract_202509242_correctedSpec(Sample_number, F_ctr_baseline_full, F_ctr_min_full, root_dir)

    % --- constants identical to Data_plot_202509242 ---
    date = "25-09-24_2";
    number_of_plots = 1;
    number_of_channels = 3;

    dV = [2.200, 2.4, 3.5;
          2.200, 2.4, 3.5;
          2.200, 2.4, 3.5];

    V_final = 60;

    time_pts = [1,10,25,30;
                1,10,25,30;
                1,10,25,30];

    % --- read sample ---
    S = readmatrix(strcat(root_dir, date, '.txt'));

    t_data_3 = [];
    F_data_1 = []; F_data_2 = []; F_data_3 = [];

    for i=1:number_of_plots
        % keep exact divisors as original Data_plot_202509242 (hours)
        t_data_3(:,i) = S(3:end,18*Sample_number(i))/60/1000/60;

        F_data_2(:,i) = S(3:end,18*Sample_number(i)-3);
        F_data_1(:,i) = S(3:end,18*Sample_number(i)-1);
        F_data_3(:,i) = S(3:end,18*Sample_number(i)+1);
    end

    time_pts = cat(2, time_pts, repmat(length(F_data_1(:,1)), size(time_pts,1), 1));

    % --- volume correction (same indexing style as your original) ---
    for i=1:number_of_plots
        for ch=1:number_of_channels
            for k=1:length(time_pts)-1
                eval(['F_data_' num2str(ch) '(time_pts(Sample_number(i),k)+1:time_pts(Sample_number(i),k+1),i) = ' ...
                      '(V_final-sum(dV(Sample_number(i),k:end)))/V_final * F_data_' num2str(ch) '(time_pts(Sample_number(i),k)+1:time_pts(Sample_number(i),k+1),i);'])
            end
        end
    end

    % --- FULL raw traces (from true t0) ---
    y_raw = [F_data_1(:,1), F_data_2(:,1), F_data_3(:,1)];
    t_raw = t_data_3(:,1);
    t_raw = t_raw - t_raw(1); % true t0 alignment

    % control full traces (assumed same time grid)
    b_raw = F_ctr_baseline_full;
    m_raw = F_ctr_min_full;

    % --- define dA-cut region (241 points) for reference-point logic ---
    idx_dA0_raw = time_pts(Sample_number(1),2);
    idx_dA1_raw = idx_dA0_raw + 240;

    % --- time_pts_new and ref points (keep your original logic) ---
    number_of_plots = 1;
    time_pts_new = time_pts(:,2:end) - (time_pts(:,2) - 1) + repmat([0,0,1,0],1,number_of_plots);
    ref_time_pts = [1,1;
                    1,3;
                    1,3]; % in dA-index coordinates (starting at 1)

    % --- compute y_max_raw and ON_full over FULL length ---
    y_max_raw = nan(size(y_raw));
    ON_full   = nan(size(y_raw));

    for ch=1:3
        tp_ref_dA_min  = time_pts_new(1, ref_time_pts(ch,1));             % dA-index
        tp_ref_min = (idx_dA0_raw - 1) + tp_ref_dA_min;
        
        tp_ref_dA_max  = time_pts_new(1, ref_time_pts(ch,2));             % dA-index
        tp_ref_max = (idx_dA0_raw - 1) + tp_ref_dA_max;                 % convert to raw index

        alpha = y_raw(tp_ref_max, ch) / b_raw(tp_ref_max, ch);
        beta = y_raw(tp_ref_min, ch) / m_raw(tp_ref_min, ch);

        y_max_raw(:,ch) = alpha * b_raw(:,ch);
        y_min_raw(:,ch) = beta * m_raw(:,ch);
        ON_full(:,ch)   = (y_max_raw(:,ch) - y_raw(:,ch)) ./ (y_max_raw(:,ch) - y_min_raw(:,ch));
    end

    % --- align Col3 at T7+1 (dA-index start = time_pts_new(1,3)) ---
    idx_T7p1_dA  = time_pts_new(1,3);
    idx_T7p1_raw = (idx_dA0_raw - 1) + idx_T7p1_dA;

    ON_T7 = ON_full(idx_T7p1_raw:idx_T7p1_raw+120, :);

    t_T7 = t_raw(idx_T7p1_raw:idx_T7p1_raw+120) - t_raw(idx_T7p1_raw);
end

