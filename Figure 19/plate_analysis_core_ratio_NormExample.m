% =====================================================
% Common analysis pipeline for plate reader experiments
% Expects all experiment-specific variables to be set in
% the caller script before this file is run.
% =====================================================

% ---------------
% Plotting colors
% ---------------
Color_RGBs = [0 0.4470 0.7410;
              0.8500 0.3250 0.0980;
              0.9290 0.6940 0.1250;
              0.4940 0.1840 0.5560;
              0.4660 0.6740 0.1880;
              0.3010 0.7450 0.9330;
              0.6350 0.0780 0.1840];

Color_Yellow = [241, 196, 15;
                243, 156, 18; 
                230, 126, 34; 
                211, 84, 0]/256;

var_num = N;
COLOR = [41/256 128/256 185/256; 241/256 196/256 15/256];
R_v = linspace(COLOR(1,1),COLOR(2,1),var_num)';
B_v = linspace(COLOR(1,2),COLOR(2,2),var_num)';
G_v = linspace(COLOR(1,3),COLOR(2,3),var_num)';
Color_RGBs_grad = [R_v B_v G_v];

sz = 30;
lw = 1.5;
linewidth = [2,2,5];
axLineWidth = 2;  % unified axes/tick line width (frame & ticks)
% ==========================================
% Import fluorometer (sample) data
% ==========================================

% Input sanity checks
assert(numel(well_list) == N_s, ...
    'well_list must have N_s (%d) entries, but got %d.', N_s, numel(well_list));
assert(numel(row_list) == 2*N_n, ...
    'row_list must have 2*N_n (%d) entries, but got %d.', 2*N_n, numel(row_list));

% Convert well IDs to Excel column letters
col_list = well2col_converter(well_list);     % 1×N_s string

% Pre-allocate (channels × samples)
cell_list      = strings(N_n, N_s);
cell_list_time = strings(N_n, N_s);

for i = 1:N_n
    % Data ranges
    cell_list(i,:)      = append(col_list, row_list(2*(i-1)+1), ":", col_list, row_list(2*(i-1)+2));
    % Time ranges (column A)
    cell_list_time(i,:) = append("A",      row_list(2*(i-1)+1), ":", "A",      row_list(2*(i-1)+2));
end

% (Duplicated assignment kept for backward compatibility)
col_list = well2col_converter(well_list);
for i=1:N_n
    cell_list(i,:)      = append(col_list, row_list(2*(i-1)+1), ":", col_list, row_list(2*(i-1)+2));
    cell_list_time(i,:) = append("A",      row_list(2*(i-1)+1), ":", "A",      row_list(2*(i-1)+2));
end

% Load each range into a cell (allowing variable length)
T = cell(N_n, N_s);   % channel × sample
L = zeros(N_n, N_s);  % length per series

for i = 1:N_s
    for j = 1:N_n
        tmp = readmatrix(filename, "Sheet", sheetname, "range", cell_list(j,i));
        tmp = tmp(:).';          % row vector
        T{j,i} = tmp;
        L(j,i) = numel(tmp);
    end
end

% NaN-pad into matrices F_data_1..N_n, each N_s × Lmax
Lmax = max(L, [], 'all');
for j = 1:N_n
    eval(sprintf('F_data_%d = nan(N_s, Lmax);', j));
    for i = 1:N_s
        x = T{j,i};
        eval(sprintf('F_data_%d(i,1:numel(x)) = x;', j));
    end
end

% Time axes (minutes) for each channel
for j = 1:N_n
    eval(sprintf('t_data_%d = (0:Lmax-1) * interval;', j));  % minutes
end

% ===================================
% Import control (baseline/min) data
% ===================================

% Sanity checks for control ranges
assert(numel(ctr_well_list) == 2, 'ctr_well_list must have 2 entries.');
assert(numel(ctr_row_list)  == 2*N_n, ...
    'ctr_row_list must have 2*N_n (%d) entries, but got %d.', 2*N_n, numel(ctr_row_list));

ctr_col_list = well2col_converter(ctr_well_list);  % 1×2

% Pre-allocate (channels × 2 control wells)
ctr_cell_list      = strings(N_n, 2);
ctr_cell_list_time = strings(N_n, 2);

for i = 1:N_n
    ctr_cell_list(i,:)      = append(ctr_col_list, ctr_row_list(2*(i-1)+1), ":", ctr_col_list, ctr_row_list(2*(i-1)+2));
    ctr_cell_list_time(i,:) = append("A",          ctr_row_list(2*(i-1)+1), ":", "A",          ctr_row_list(2*(i-1)+2));
end

% (Duplicated assignment kept for backward compatibility)
ctr_col_list = well2col_converter(ctr_well_list);
for i=1:N_n
    ctr_cell_list(i,:)      = append(ctr_col_list, ctr_row_list(2*(i-1)+1), ":", ctr_col_list, ctr_row_list(2*(i-1)+2));
    ctr_cell_list_time(i,:) = append("A",          ctr_row_list(2*(i-1)+1), ":", "A",          ctr_row_list(2*(i-1)+2));
end

% Load controls into cells → NaN-pad into ctr_F_data_1..N_n (2 × Lmax_ctr)
CTR = cell(N_n, 2);
Lc  = zeros(N_n, 2);
for i = 1:2
    for j = 1:N_n
        tmp = readmatrix(ctr_filename, "Sheet", ctr_sheetname, "range", ctr_cell_list(j,i));
        tmp = tmp(:).';
        CTR{j,i} = tmp;
        Lc(j,i)  = numel(tmp);
    end
end

Lmax_ctr = max(Lc, [], 'all');

for j = 1:N_n
    eval(sprintf('ctr_F_data_%d = nan(2, Lmax_ctr);', j));
    for i = 1:2
        x = CTR{j,i};
        eval(sprintf('ctr_F_data_%d(i,1:numel(x)) = x;', j));
    end
end

% Control time axis (minutes)
ctr_t_data = (0:Lmax_ctr-1) * interval;

% ===================
% Figures 1–3: Raw + Smoothed (3×3 per channel)
% ===================
for j = 1:N_n
    figure(j); clf
    tiledlayout(3,3)
end

% Full-scale per channel (NaN-safe)
for j = 1:N_n
    eval(['F_max(j) = nanmax(F_data_' num2str(j) '(:));'])
    eval(['F_min(j) = nanmin(F_data_' num2str(j) '(:));'])
end

for i = 1:N_s
    for j = 1:N_n
        figure(j);
        nexttile(i);
        hold on

        eval(['t_vec_min = t_data_' num2str(j) ';'])   % minutes
        t_data_hr = t_vec_min / 60;                    % hours
        eval(['F_row = F_data_' num2str(j) '(i,:);'])

        % Smoothing only over valid (non-NaN) region
        valid = ~isnan(t_data_hr) & ~isnan(F_row);
        t_loc = t_data_hr(valid);
        y_loc = F_row(valid);

        % Segment indices (channel-dependent) in index space, then clipped
        % For ch1: [dA, dR] → [1,3]; for ch2/3: [dA, T7] → [1,2]
        tp_idx   = [1,3; 1,2; 1,2];  % ch1: dA→dR, ch2/3: dA→T7
        tps      = [time_pts(i, tp_idx(j,1)), time_pts(i, tp_idx(j,2))];
        tps_local = tps;
        tps_local(~isnan(tps_local)) = max(1, min(numel(y_loc), tps_local(~isnan(tps_local))));

        % Plot raw data
        scatter(t_loc, y_loc, sz, Color_RGBs(j,:), 'filled');

        % Piecewise smoothed data (stitched with smoothed joints)
        smoothed = separate_smooth_join(t_loc, y_loc, tps_local, method(j), window(j,:));
        plot(t_loc, smoothed, "LineWidth", lw, "Color", [0,0,0])

        % Vertical lines at reagent additions (absolute time)
        for k=1:size(time_pts,2)
            if ~isnan(time_pts(i,k))
                tline = interval/60 * ((2*time_pts(i,k)-1)/2);
                plot([tline tline], [0, F_max(j)], 'k--', 'LineWidth', 2)
            end
        end

        legend(Channels(j), 'Location', 'best'); legend boxoff
        box on
        xlabel('Time [hours]');
        ylabel('Count');
        xlim([0, inf]);
        ylim([0, F_max(j)]);
        title(Sample_label(i))
        set(gca,'FontSize',15,'LineWidth',axLineWidth);
        hold off
    end
end

% ======================================
% Prepare arrays for Fraction ON (121 pt, T7+1 aligned)
% ======================================
Fraction_ON_1 = nan(N_s, 121);
Fraction_ON_2 = nan(N_s, 121);
Fraction_ON_3 = nan(N_s, 121);

figure(4); clf
tiledlayout(3,3)

% ================================
% Compute Fraction ON (T7+1 aligned)
% ================================
for i=1:N_s
    nexttile(i);
    hold on

    for j=1:N_n
        % Full traces
        eval(['t_data_full_min = t_data_' num2str(j) ';']) % minutes
        eval(['F_data_full     = F_data_' num2str(j) '(i,:);'])
        t_data_full_hr = t_data_full_min / 60;

        % Control traces for min/baseline (smoothed)
        eval(['F_min_c = ctr_F_data_' num2str(j) '(ctr_min,:);'])
        F_min_smoothed = separate_smooth_join(ctr_t_data/60, F_min_c, ctr_time_pts(ctr_min,:), ...
                                              method(j), window(j,:));

        eval(['F_baseline_c = ctr_F_data_' num2str(j) '(ctr_baseline,:);'])
        % If a separate baseline-only window is used, replace window(j,:) below.
        F_baseline_smoothed = separate_smooth_join(ctr_t_data/60, F_baseline_c, ctr_time_pts(ctr_baseline,:), ...
                                                   method(j), window(j,:));

        % Sample piecewise smoothing (after removing NaNs)
        valid      = ~isnan(t_data_full_hr) & ~isnan(F_data_full);
        t_full_loc = t_data_full_hr(valid);
        y_full_loc = F_data_full(valid);

        % Segment indices for smoothing
        % For ch1: [dA, dR]; for ch2/3: [dA, T7]
        tp_idx = [1,3; 1,2; 1,2];
        time_pts_norm_local = [time_pts(i, tp_idx(j,1)), time_pts(i, tp_idx(j,2))];
        time_pts_norm_local(~isnan(time_pts_norm_local)) = ...
            max(1, min(numel(y_full_loc), time_pts_norm_local(~isnan(time_pts_norm_local))));
        y_smooth = separate_smooth_join(t_full_loc, y_full_loc, time_pts_norm_local, ...
                                        method(j), window(j,:));

        % Normalize (length alignment handled inside the function)
        frac_full = normalize_decaying_baseline(y_smooth, F_baseline_smoothed, F_min_smoothed, time_pts_ref(j));

        % ---- Extract 121 points with T7+1 as time zero ----
        % Take 121 points starting from the index just after the original T7 index
        t0_idx  = time_pts(i,2) + 1;          % skip the exact T7 sample
        t0_idx  = max(1, min(t0_idx, numel(frac_full)));
        idx_end = t0_idx + 120;

        if t0_idx > numel(frac_full)
            % T7+1 lies beyond the data length: full NaN window
            frac_win = nan(1,121);
        elseif idx_end > numel(frac_full)
            tail  = frac_full(t0_idx:end);
            padN  = idx_end - numel(frac_full);
            frac_win = [tail, nan(1, padN)];
        else
            frac_win = frac_full(t0_idx:idx_end);
        end

        % Save window
        switch j
            case 1, Fraction_ON_1(i,:) = frac_win;
            case 2, Fraction_ON_2(i,:) = frac_win;
            case 3, Fraction_ON_3(i,:) = frac_win;
        end

        % Plot (T7+1 as 0h)
        plot(t_T7_axis, frac_win, "LineWidth", linewidth(j), "Color", Color_RGBs(j,:));
    end

    % Vertical lines at reagent additions, relative to T7+1 (0..4h only)
    rel_hours = (time_pts(i,:) - time_pts(i,2) - 1) * interval / 60;
    for k=1:length(rel_hours)
        if ~isnan(rel_hours(k)) && rel_hours(k) >= 0 && rel_hours(k) <= 4
            plot([rel_hours(k) rel_hours(k)], [0,1], 'k--', 'LineWidth', 2)
        end
    end

    lgd = legend(Channels); set(lgd,'color','none'); legend boxoff
    title(Sample_label(i))
    box on; set(gca,'FontSize',15,'LineWidth',axLineWidth);
    xlim([0, inf]); ylim([0, 1]);
    hold off
end

% ======================
% Figure 5: condition × channel = 3×3 (Rep#1..3 overlaid)
% ======================
figure(5); clf
tiledlayout(3,3)  % rows = channels, columns = conditions

for j = 1:N_n    % row: channel 1..3
    for c = 1:3  % col: condition 1..3
        nexttile( (j-1)*3 + c ); hold on

        % Sample IDs belonging to this condition (clipped to 1..N)
        ids = [c, c+3, c+6];
        ids = ids(ids<=N & ids>=1);

        % Select Fraction_ON for channel j
        switch j
            case 1, FON = Fraction_ON_1;
            case 2, FON = Fraction_ON_2;
            case 3, FON = Fraction_ON_3;
        end

        for k = 1:numel(ids)
            id = ids(k);
            plot(t_T7_axis, FON(id,:), "LineWidth", linewidth(j), "Color", Color_RGBs_grad(id,:));
        end

        if ~isempty(ids)
            lgd = legend(Sample_label(ids)); set(lgd,'color','none'); legend boxoff
        end

        xlabel('Time [hours]'); ylabel('Fraction ON');
        ylim([0,1]); xlim([0,inf]); box on; set(gca,'FontSize',15,'LineWidth',axLineWidth);

        if j == 1, title(cond_labels(c)); end
        if c == 1, ylabel(sprintf('%s\nFraction ON', Channels(j))); end

        hold off
    end
end

% ======================
% Figure 6: mean ± std (shaded) for each (channel × condition)
% ======================
% Output arrays for later use
% Shape: (channel=1..3, condition=1..3, time=121)
FON_mean = nan(N_n, 3, 121);
FON_std  = nan(N_n, 3, 121);
FON_sem  = nan(N_n, 3, 121);
FON_nEff = zeros(N_n, 3, 121);   % effective N at each time point

% Replicate IDs per condition (1..9)
cond_ids = {
    [1, 4, 7];   % condition 1 (e.g., 5 nM)
    [2, 5, 8];   % condition 2 (e.g., 10 nM)
    [3, 6, 9]    % condition 3 (e.g., 25 nM)
};

figure(6); clf
tiledlayout(3,3)

for j = 1:N_n
    % Select Fraction_ON for channel j
    switch j
        case 1, FON = Fraction_ON_1;
        case 2, FON = Fraction_ON_2;
        case 3, FON = Fraction_ON_3;
    end

    for c = 1:3
        ids = cond_ids{c};
        ids = ids(ids>=1 & ids<=size(FON,1));

        % data: (rep ≤3) × 121
        data = FON(ids, :);

        % NaN-safe statistics at each time point
        mu   = nanmean(data, 1);
        sd   = nanstd(data, 0, 1);
        nEff = sum(~isnan(data), 1);
        sem  = sd ./ sqrt(max(nEff,1));

        % Store into 3D arrays
        FON_mean(j, c, :) = mu;
        FON_std (j, c, :) = sd;
        FON_sem (j, c, :) = sem;
        FON_nEff(j, c, :) = nEff;

        % Plot
        nexttile((j-1)*3 + c); hold on

        % Shaded region: mean ± std
        y_up = mu + sd;
        y_lo = mu - sd;

        X = [t_T7_axis, fliplr(t_T7_axis)];
        Y = [y_lo,       fliplr(y_up)      ];

        p = patch(X, Y, 1); %#ok<NASGU>
        set(p, 'FaceColor', Color_RGBs(j,:), ...
               'FaceAlpha', 0.20, ...
               'EdgeColor', 'none')

        % Mean line
        plot(t_T7_axis, mu, 'LineWidth', linewidth(j), 'Color', Color_RGBs(j,:));

        xlabel('Time [hours]');
        ylabel('Fraction ON');
        ylim([-0.2,1.2]); xlim([0,inf]); box on; set(gca,'FontSize',15,'LineWidth',axLineWidth);

        if j == 1, title(cond_labels(c)); end
        if c == 1, ylabel(sprintf('%s\nFraction ON', Channels(j))); end

        hold off
    end
end

% ======================================
% AUC (G1S1 / 6-FAM, T7+1 aligned 0–4h window)
% ======================================
% Assumptions:
%  - t_T7_axis: 0..4 [hours], 121 points
%  - Fraction_ON_3: (N_s x 121) for G1S1 (6-FAM), aligned to T7+1

AUC_6FAM = nan(N_s, 1);
for i = 1:N_s
    y = Fraction_ON_3(i, :);
    t = t_T7_axis;
    valid = isfinite(t) & isfinite(y);
    if nnz(valid) >= 2
        AUC_6FAM(i) = calculate_AUC(t(valid), y(valid));
    else
        AUC_6FAM(i) = NaN;
    end
end

% Replicate mean/std per condition ([1,4,7], [2,5,8], [3,6,9])
cond_ids = {
    [1, 4, 7];
    [2, 5, 8];
    [3, 6, 9]
};

AUC_6FAM_by_cond_mean = nan(1,3);
AUC_6FAM_by_cond_std  = nan(1,3);
AUC_6FAM_by_cond_sem  = nan(1,3);
AUC_6FAM_by_cond_n    = zeros(1,3);

for c = 1:3
    ids = cond_ids{c};
    ids = ids(ids>=1 & ids<=N_s);
    vals = AUC_6FAM(ids);
    AUC_6FAM_by_cond_mean(c) = nanmean(vals);
    AUC_6FAM_by_cond_std(c)  = nanstd(vals, 0);
    AUC_6FAM_by_cond_n(c)    = sum(isfinite(vals));
    AUC_6FAM_by_cond_sem(c)  = AUC_6FAM_by_cond_std(c) / max(AUC_6FAM_by_cond_n(c),1);
end

% Convenient struct for AUC summary
AUC_6FAM_stats.mean   = AUC_6FAM_by_cond_mean;
AUC_6FAM_stats.std    = AUC_6FAM_by_cond_std;
AUC_6FAM_stats.sem    = AUC_6FAM_by_cond_sem;
AUC_6FAM_stats.n      = AUC_6FAM_by_cond_n;
AUC_6FAM_stats.labels = {"[G3R1]=5 nM","[G3R1]=10 nM","[G3R1]=25 nM"};

% ============================================
% ============ Local helper functions =========
% ============================================

function col_list = well2col_converter(well_list)
    % Convert 96/384-well ID (e.g., 'A1') into Excel column letters with an offset.
    % 'alphabets' maps index to Excel letters. 'Z' is a sentinel so A=2, B=3, ...
    alphabets = ["Z" "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y"];
    col_list = {};
    for i=1:size(well_list,2)
        well = well_list{i};
        row  = find(alphabets==well(1)) - 1;  % row index (A->1, B->2, ...)
        col  = str2double(well(2:end));       % numeric column
        idx  = 24*(row-1) + col + 2;          % '+2' corresponds to an Excel offset
        digit = floor(log2(idx)/log2(26)) + 1;

        col_temp = "";
        if digit == 1
            col_temp = append(alphabets(idx+1), col_temp);
        else
            for j=1:digit
                if j==digit
                    col_temp = append(alphabets(floor((idx-1)/26^(j-1)) + 1), col_temp);
                else
                    col_temp = append(alphabets(floor(mod(idx,26^j)/26^(j-1)) + 1), col_temp);
                end
            end
        end
        col_list{i} = col_temp;
    end
    col_list = string(col_list);
end

function data_joined = separate_smooth_join(t,y,time_pts,method,window)
    % Piecewise smooth y(t) in segments defined by 'time_pts' (indices).
    % 'time_pts' may contain NaN to indicate "stop after this segment".
    % 'window' can be a 1xK array; if shorter than needed, it is padded.
    % t, y are row vectors of equal length; t is used as SamplePoints.
    %
    % Joints are returned as smoothed values (not raw values).

    if iscolumn(t), t = t.'; end
    if iscolumn(y), y = y.'; end

    assert(numel(t)==numel(y), 'separate_smooth_join: t and y must have equal length');

    % time_pts_ uses 0 as a sentinel so the first segment starts at index 1
    time_pts_ = [0, time_pts];
    data_joined = [];

    % Pad window length to the number of segments
    needed = numel(time_pts_);
    if numel(window) < needed
        window = [window, repmat(window(end), 1, needed - numel(window))];
    end

    for i=1:length(time_pts_)
        if i < length(time_pts_) && isnan(time_pts_(i+1))
            % Last segment: from current sentinel to the end
            seg_idx   = time_pts_(i)+1:numel(y);
            data_each = smoothdata(y(seg_idx), method, window(i), "SamplePoints", t(seg_idx));
            data_joined = [data_joined, data_each];
            break
        end
        if i == length(time_pts_)
            % Final segment
            seg_idx   = time_pts_(i)+1:numel(y);
            data_each = smoothdata(y(seg_idx), method, window(i), "SamplePoints", t(seg_idx));
            data_joined = [data_joined, data_each];
        else
            % Intermediate segment: [start+1 .. next_break]
            seg_idx   = time_pts_(i)+1:time_pts_(i+1);
            data_each = smoothdata(y(seg_idx), method, window(i), "SamplePoints", t(seg_idx));
            % No trimming here because the next segment starts at the next index
            data_joined = [data_joined, data_each];
        end
    end
end

function y_normed = normalize_decaying_baseline(y, y_baseline, y_min, time_pt_ref)
    % Length-safe normalization against a decaying baseline.
    % Scales y so that at 'time_pt_ref' it matches y_baseline(time_pt_ref),
    % then maps to Fraction ON via (y_max - y) / (y_max - y_min),
    % where y_max is a baseline-proportional target curve.

    if iscolumn(y),          y          = y.'; end
    if iscolumn(y_baseline), y_baseline = y_baseline.'; end
    if iscolumn(y_min),      y_min      = y_min.'; end

    % Align lengths to the shortest input
    L = min([numel(y), numel(y_baseline), numel(y_min)]);
    y          = y(1:L);
    y_baseline = y_baseline(1:L);
    y_min      = y_min(1:L);

    % Clamp reference index to [1, L]
    time_pt_ref = max(1, min(time_pt_ref, L));

    % Baseline-proportional target
    denom = y_baseline(time_pt_ref);
    if isnan(denom) || denom==0
        scale = 1;
    else
        scale = y(time_pt_ref) / denom;
    end
    y_max = scale * y_baseline;
    y_min = scale * y_min;
    
    % Fraction ON
    denom2 = (y_max - y_min);
    denom2(denom2==0) = NaN;
    y_normed = (y_max - y) ./ denom2;
end

function I = calculate_AUC(t,y)
    % Simple right Riemann sum AUC: sum over (t_i - t_{i-1}) * y_i
    I = 0;
    for i=2:length(t)
        I = I + (t(i)-t(i-1))*(y(i));
    end
end

% =========================================
% Normalization visualization (Figure 7)
%  - 3 channels × 3 panels (smoothed only)
%    Col 1: smoothed (pre-scaling) sample + controls
%    Col 2: smoothed (post-scaling) sample + controls (scaled at ref)
%    Col 3: 0–4 h window (T7+1 = 0), same as Fig6
% =========================================

% If axis line width is not defined in this script, set a default
if ~exist('axLineWidth','var')
    axLineWidth = 1.5;
end

% --------- CONFIG: which sample to visualize ---------
sample_idx_for_norm_viz = 1;

% --------- CONFIG: line colors and styles for this figure ---------
norm_color_raw_sample        = [0.0 0.0 0.0];
norm_color_raw_ctr_baseline  = [0.0 0.4470 0.7410];
norm_color_raw_ctr_min       = [0.8500 0.3250 0.0980];

norm_ls_raw_sample           = '-';
norm_ls_raw_ctr_baseline     = '--';
norm_ls_raw_ctr_min          = ':';

norm_lw_raw_sample           = 1.5;
norm_lw_raw_ctrl             = 1.2;

font_size                    = 16;

% For 0–4 h window (Col 3)
norm_color_frac_window = [0.1 0.6 0.1];
norm_ls_frac_window    = '-';
norm_lw_frac_window    = 2.0;

% --------- Segment index pattern per channel ---------
tp_idx_mat = [1,3; 1,2; 1,2];

% --------- Figure 7 layout ---------
figure(7); clf
tiledlayout(N_n, 3);   % rows = channels, cols = 3 panels (=> 9 tiles)

for j = 1:N_n
    % ==========================
    % Common data for this channel (smoothed)
    % ==========================

    % ---- Sample: smoothed (piecewise) ----
    eval(['t_sample_min  = t_data_' num2str(j) ';']);              % minutes
    eval(['y_sample_full_raw = F_data_' num2str(j) '(sample_idx_for_norm_viz,:);']);
    t_sample_hr_full = t_sample_min / 60;                           % hours

    valid_samp        = isfinite(t_sample_hr_full) & isfinite(y_sample_full_raw);
    t_sample_loc_hr   = t_sample_hr_full(valid_samp);
    y_sample_loc_raw  = y_sample_full_raw(valid_samp);

    tp_idx = tp_idx_mat(j,:);
    time_pts_norm_local = [time_pts(sample_idx_for_norm_viz, tp_idx(1)), ...
                           time_pts(sample_idx_for_norm_viz, tp_idx(2))];
    time_pts_norm_local(~isnan(time_pts_norm_local)) = ...
        max(1, min(numel(y_sample_loc_raw), time_pts_norm_local(~isnan(time_pts_norm_local))));

    y_smooth_full = separate_smooth_join(t_sample_loc_hr, y_sample_loc_raw, ...
                                         time_pts_norm_local, method(j), window(j,:));
    t_full_hr = t_sample_loc_hr;

    % ---- Controls: smoothed min & baseline ----
    eval(['y_ctr_min_raw      = ctr_F_data_' num2str(j) '(ctr_min,:);']);
    eval(['y_ctr_baseline_raw = ctr_F_data_' num2str(j) '(ctr_baseline,:);']);
    t_ctr_hr = ctr_t_data / 60;

    F_min_smoothed = separate_smooth_join(t_ctr_hr, y_ctr_min_raw, ...
                                          ctr_time_pts(ctr_min,:), method(j), window(j,:));

    win_base = window(j,:);
    if ~isempty(ctr_baseline_window)
        win_base = ctr_baseline_window(j,:);
    end
    F_baseline_smoothed = separate_smooth_join(t_ctr_hr, y_ctr_baseline_raw, ...
                                               ctr_time_pts(ctr_baseline,:), method(j), win_base);

    % ===== Common length for Col 1 & Col 2 =====
    L_norm = min([numel(y_smooth_full), numel(F_baseline_smoothed), numel(F_min_smoothed)]);
    y_samp_n   = y_smooth_full(1:L_norm);
    base_n     = F_baseline_smoothed(1:L_norm);
    min_n      = F_min_smoothed(1:L_norm);
    t_full_hr_norm = t_full_hr(1:L_norm);
    t_ctr_hr_norm  = t_ctr_hr(1:L_norm);

    % Reference index (same rule as main pipeline)
    ref_idx = max(1, min(time_pts_ref(j), L_norm));

    % Scaling so that sample(ref_idx) == baseline(ref_idx)
    denom = base_n(ref_idx);
    if isnan(denom) || denom == 0
        scale2 = 1;
    else
        scale2 = y_samp_n(ref_idx) / denom;
    end
    y_ctr_baseline_scaled = scale2 * base_n;
    y_ctr_min_scaled      = scale2 * min_n;

    % ========= NEW: draw dashed lines at t_dA and t_T7 for ALL channels =========
    % Use the same absolute-time convention as Figures 1–3:
    % tline = (idx - 0.5) * interval / 60
    idx_dA = time_pts(sample_idx_for_norm_viz, 1);
    idx_T7 = time_pts(sample_idx_for_norm_viz, 2);
    tline_dA = interval/60 * ((2*idx_dA - 1)/2);
    tline_T7 = interval/60 * ((2*idx_T7 - 1)/2);
    % ==========================================================================

    % ==========================
    % Col 1: smoothed pre-scaling
    % ==========================
    nexttile((j-1)*3 + 1); hold on

    plot(t_full_hr_norm, y_samp_n, ...
        'Color', norm_color_raw_sample, 'LineStyle', norm_ls_raw_sample, 'LineWidth', norm_lw_raw_sample);
    plot(t_ctr_hr_norm, base_n, ...
        'Color', norm_color_raw_ctr_baseline, 'LineStyle', norm_ls_raw_ctr_baseline, 'LineWidth', norm_lw_raw_ctrl);
    plot(t_ctr_hr_norm, min_n, ...
        'Color', norm_color_raw_ctr_min, 'LineStyle', norm_ls_raw_ctr_min, 'LineWidth', norm_lw_raw_ctrl);

    xline(tline_dA, 'k--', 'LineWidth', 1.0);
    xline(tline_T7, 'k--', 'LineWidth', 1.0);

    xlabel('Time [hours]'); xlim([0,5.5]);
    ylabel('Smoothed intensity (a.u.)');
    title(sprintf('%s – smoothed (pre-scaling)', Channels(j)));
    legend({'Sample (smoothed)','Control baseline (smoothed)','Control min (smoothed)'}, ...
           'Location','best','Box','off');
    box on; set(gca,'FontSize',font_size,'LineWidth',axLineWidth);
    hold off

    % ==========================
    % Col 2: smoothed post-scaling
    % ==========================
    nexttile((j-1)*3 + 2); hold on

    plot(t_full_hr_norm, y_samp_n, ...
        'Color', norm_color_raw_sample, 'LineStyle', norm_ls_raw_sample, 'LineWidth', norm_lw_raw_sample);
    plot(t_ctr_hr_norm, y_ctr_baseline_scaled, ...
        'Color', norm_color_raw_ctr_baseline, 'LineStyle', norm_ls_raw_ctr_baseline, 'LineWidth', norm_lw_raw_ctrl);
    plot(t_ctr_hr_norm, y_ctr_min_scaled, ...
        'Color', norm_color_raw_ctr_min, 'LineStyle', norm_ls_raw_ctr_min, 'LineWidth', norm_lw_raw_ctrl);

    % dashed at t_dA and t_T7 (requested)
    xline(tline_dA, 'k--', 'LineWidth', 1.0);
    xline(tline_T7, 'k--', 'LineWidth', 1.0);

    xlabel('Time [hours]'); xlim([0,5.5]);
    ylabel('Smoothed intensity (a.u.)');
    title(sprintf('%s – smoothed (scaled)', Channels(j)));
    legend({'Sample (smoothed)','Control baseline (scaled)','Control min (scaled)'}, ...
           'Location','best','Box','off');
    box on; set(gca,'FontSize',font_size,'LineWidth',axLineWidth);
    hold off

    % ==========================
    % Col 3: 0–4 h window (T7+1 as t = 0), same as Fig6
    % ==========================
    % Full-length Fraction ON (needed only to extract the window)
    y_max_n  = scale2 * base_n;
    y_min_n  = scale2 * min_n;
    frac_full = (y_max_n - y_samp_n) ./ (y_max_n - y_min_n);

    t0_idx  = time_pts(sample_idx_for_norm_viz, 2) + 1;   % T7+1 as t0
    t0_idx  = max(1, min(t0_idx, numel(frac_full)));
    idx_end = t0_idx + 120;

    if t0_idx > numel(frac_full)
        frac_win_T7 = nan(1,121);
    elseif idx_end > numel(frac_full)
        tail        = frac_full(t0_idx:end);
        padN        = idx_end - numel(frac_full);
        frac_win_T7 = [tail, nan(1, padN)];
    else
        frac_win_T7 = frac_full(t0_idx:idx_end);
    end

    nexttile((j-1)*3 + 3); hold on
    plot(t_T7_axis, frac_win_T7, ...
        'Color', norm_color_frac_window, 'LineStyle', norm_ls_frac_window, 'LineWidth', norm_lw_frac_window);

    % Requested dashed lines (relative to T7+1 = 0)
    rel_dA = (idx_dA - idx_T7 - 1) * interval / 60;
    rel_T7 = (idx_T7 - idx_T7 - 1) * interval / 60;  % = -interval/60
    xline(rel_dA, 'k--', 'LineWidth', 1.0);
    xline(rel_T7, 'k--', 'LineWidth', 1.0);

    xlabel('Time [hours]');
    ylabel('Fraction ON');
    xlim([0, 4]);                  % same convention as Fig6
    ylim([0, 1]);
    title(sprintf('%s – Fraction ON (after T7 RNAP added)', Channels(j)));
    box on; set(gca,'FontSize',font_size,'LineWidth',axLineWidth);
    hold off
end


