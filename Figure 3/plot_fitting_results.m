%% ============================================================
%  analyze_IFFL_AUC_from_results.m
%  - Independent script to compare Exp1/Exp2/Exp3 AUC
%    between experiment and simulation using IFFL_fit_results.
%  - Run AFTER run_iffl_fit_include_autoinhibition() so that
%    IFFL_fit_results exists in the base workspace.
%
%  Notes:
%  - This script is self-contained: it includes simulate_IFFL(),
%    General_genelet_model_segment(), initial_conditions(), etc.
%  - It uses IFFL_fit_results.final as the fitted parameter struct p.
%  - Exp2/Exp3 use p.tuned.exp2 / p.tuned.exp3 directly.
% ============================================================

% clear; close all; clc;

%% ======== 0. Load fitted results from base workspace ====================
try
    R = evalin('base','IFFL_fit_results');
catch
    error('IFFL_fit_results not found in base workspace. Run the fitting code first.');
end

if ~isfield(R,'final')
    error('IFFL_fit_results.final not found.');
end

p_fit = R.final;

% ======== (NEW) Exp0.5-only scaling factor from Step4 (fit alpha) ========
% This alpha must be applied ONLY to Fig1.5 Exp0.5 simulation (Tile 2).
% It must NOT modify p_fit (global fitted parameters) or affect other experiments.
alpha_exp05 = 1.0;

if isfield(R,'exp4') && isfield(R.exp4,'alpha_best') && isfinite(R.exp4.alpha_best)
    alpha_exp05 = R.exp4.alpha_best;
elseif isfield(R,'stages')
    % Fallback: try to find a stage that stored alpha_best
    for si = 1:numel(R.stages)
        if isfield(R.stages(si),'alpha_best') && isfinite(R.stages(si).alpha_best)
            alpha_exp05 = R.stages(si).alpha_best;
        end
    end
else
    warning('Exp0.5 alpha not found in IFFL_fit_results. Using alpha=1.0 (no scaling).');
end

fprintf('[Info] Exp0.5 alpha (Step4) = %.6g (will be used ONLY for Fig1.5 Exp0.5 sim)\n', alpha_exp05);


% Model size (should match fitting code)
n = 3;
m = 5;

%% ======== (NEW) Plot parameters (Fig1+Fig2) =============================
Plot = struct();

% ===== Figure 1 (TS) =====
Plot.TS_LW_data = 4.0;    % solid (data)
Plot.TS_LW_sim  = 1.5;    % dashed (sim)

% ===== Figure 1/2: frame/tick thickness =====
Plot.AxesLineWidth = 2.5; % outer box thickness (also affects tick thickness in MATLAB)
Plot.TickLineWidth = 1.5; % kept as a parameter (MATLAB ties tick thickness to Axes LineWidth)

% ===== Figure 2 (AUC tiles): marker size =====
Plot.AUC_MarkerSize_exp = 12;   % experimental data markers (filled)
Plot.AUC_MarkerSize_sim = 12;   % simulation markers

% ===== Figure 2 (AUC tiles): colors (nominal=black, variants=purple/pink) =====
Plot.ColorOrder_exp1  = [0 0 0; ...
                         hex2rgb_mat("#7B2CBF"); ... % purple
                         hex2rgb_mat("#FF4FA3")];    % bright pink
Plot.ColorOrder_exp23 = Plot.ColorOrder_exp1;

% ===== Existing-ish style params =====
Plot.ExpLineWidth   = 3;
Plot.SimLineWidth   = 3;
Plot.FontSize_Fig2  = 20;
Plot.AUC_ShadowAlpha = 0.18;  % <-- NEW: ±1σ shadow transparency

% ===== Figure 2 fixed axes =====
Plot.AUC_xlim   = [0 70];
Plot.AUC_ylim   = [-0.3 3.2];
Plot.AUC_yticks = [0 0.5 1 1.5 2 2.5 3];

%% ======== TS overlay plots (match old Fig7/8/9 styling) =================

dirpath = fileparts(mfilename('fullpath'));

% --- read mean time-series CSVs (same as fitting inputs) ---
filePattern  = 'Din_%dmin_means.csv';
timeCol      = 'time_h';
ch1Col       = 'ON_mean_ch1';
ch2Col       = 'ON_mean_ch2';
ch3Col       = 'ON_mean_ch3';
clipZero     = false;

T_list_min   = [];   % [] => auto-detect Din from filenames

dataTS = read_iffl_dir(dirpath, ...
    'Pattern',filePattern, ...
    'T_list_min',T_list_min, ...
    'TimeCol',timeCol, ...
    'Ch1Col',ch1Col, ...
    'Ch2Col',ch2Col, ...
    'Ch3Col',ch3Col, ...
    'ClipZero',clipZero);

% --- styling copied from your legacy code ---
Colors_red   = ["#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026"];
Colors_green = ["#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45","#005a32"];
Colors_blue  = ["#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594"];

% Auto-build from detected Din (minutes)
Din_detected = arrayfun(@(c) round(c.T_sec/60), dataTS.conds);
[Din_sorted, sortIdx] = sort(Din_detected);

% reorder conditions in ascending Din so colors/legend match nicely
dataTS.conds = dataTS.conds(sortIdx);
nCond = numel(dataTS.conds);

Sample_labels = strings(1,nCond);
for k = 1:nCond
    Sample_labels(k) = sprintf('D_{in}=%d min', round(dataTS.conds(k).T_sec/60));
end

% --- figure/tiled layout ---
figTS = figure('Name','Time-series (Data solid) vs (Sim dashed)');
figTS.Position = [0,0,400,600];
tlTS = tiledlayout(figTS, 3, 1, 'TileSpacing','compact');

% line widths (from Plot)
LW_data = Plot.TS_LW_data;
LW_sim  = Plot.TS_LW_sim;

ax1 = nexttile(tlTS); hold(ax1,'on'); box(ax1,'on');
ax2 = nexttile(tlTS); hold(ax2,'on'); box(ax2,'on');
ax3 = nexttile(tlTS); hold(ax3,'on'); box(ax3,'on');

for k = 1:nCond
    ck = dataTS.conds(k);

    % Pick colors by channel, index by Din order (k)
    c1 = Colors_blue(min(k, numel(Colors_blue)));
    c2 = Colors_green(min(k, numel(Colors_green)));
    c3 = Colors_red(min(k, numel(Colors_red)));

    % --- Data (solid) ---
    plot(ax1, ck.time_s/3600, ck.y_ch1, '-', 'LineWidth', LW_data, 'Color', c1);
    plot(ax2, ck.time_s/3600, ck.y_ch2, '-', 'LineWidth', LW_data, 'Color', c2);
    plot(ax3, ck.time_s/3600, ck.y_ch3, '-', 'LineWidth', LW_data, 'Color', c3);

    % --- Sim (dashed) ---
    p_ij = p_fit;
    p_ij.T = ck.T_sec;

    [tSim, YSim, ok] = simulate_IFFL(p_ij, n, m);
    if ok
        frac1 = YSim(:,1) ./ p_ij.G_tot(1);
        frac2 = YSim(:,3) ./ p_ij.G_tot(3);
        frac3 = YSim(:,5) ./ p_ij.G_tot(5);

        plot(ax1, tSim/3600, frac1, '--', 'LineWidth', LW_sim, 'Color', c1);
        plot(ax2, tSim/3600, frac2, '--', 'LineWidth', LW_sim, 'Color', c2);
        plot(ax3, tSim/3600, frac3, '--', 'LineWidth', LW_sim, 'Color', c3);
    end
end

% --- axes formatting (match your legacy Fig7/8/9) ---
for ax = [ax1 ax2 ax3]
    xlim(ax,[0,4]);
    ylim(ax,[0,1]);
    xlabel(ax,'Time [h]');
    ylabel(ax,'Active fraction');
    set(ax,'FontSize',25);

    grid(ax,'off');
    box(ax,'on');

    % frame/tick thickness
    ax.LineWidth = Plot.AxesLineWidth;
    ax.TickDir   = 'in';
end

title(ax1,'Ch1','Interpreter','none');
title(ax2,'Ch2','Interpreter','none');
title(ax3,'Ch3','Interpreter','none');

% legend(ax1, Sample_labels, 'Location','best');
% legend(ax2, Sample_labels, 'Location','best');
% legend(ax3, Sample_labels, 'Location','best');

title(tlTS, 'Time-series: data (solid) vs simulation (dashed)', 'Interpreter','none');


%% ======== 1. Embedded experimental AUC data =============================
AUC_Din_vals_min = [20 40 60];
AUC = build_embedded_AUC_data(AUC_Din_vals_min);

%% ======== 2. Analysis settings =========================================
Channel = 'ch3';  % 'ch1' | 'ch2' | 'ch3' or 1|2|3

% Time range for AUC integration [s]
Tfinal_s     = 4*3600;
AUC_tStart_s = 0;
AUC_tEnd_s   = [];     % [] means use Tfinal_s

[chIdx, denomIdx, chLabel] = map_channel(Channel);

%% ======== 3. Compute AUC tables for Exp1/2/3 ============================
DurMin = AUC.Din_vals(:).';

% --- Exp1: rows = G3R1 levels [5 10 25], cols = Din [20 40 60]
exp1_levels = AUC.exp1.G3R1_levels_nM(:).';
AUC_sim_exp1 = nan(numel(exp1_levels), numel(DurMin));

for i = 1:numel(exp1_levels)
    for j = 1:numel(DurMin)
        p_ij = override_total_G3R1_and_dB3(p_fit, exp1_levels(i));
        p_ij.T = DurMin(j)*60;

        [tSim,YSim,ok] = simulate_IFFL(p_ij, n, m);
        if ~ok, continue; end

        frac = YSim(:,chIdx) ./ p_ij.G_tot(denomIdx);

        [t_clip_h, f_clip] = clip_for_auc(tSim, frac, AUC_tStart_s, AUC_tEnd_s, Tfinal_s);
        if isempty(t_clip_h), continue; end

        AUC_sim_exp1(i,j) = calculate_AUC(t_clip_h, f_clip);
    end
end

% --- Exp2: rows = toehold [8 6 4], cols = Din [20 40 60]
exp2_toe = AUC.exp2.toehold_nt(:).';
AUC_sim_exp2 = nan(numel(exp2_toe), numel(DurMin));
for i = 1:numel(exp2_toe)
    idx = find(p_fit.tuned.exp2.toehold == exp2_toe(i), 1);
    if isempty(idx), error('Exp2: toehold %d not found in p_fit.tuned.exp2.toehold', exp2_toe(i)); end

    for j = 1:numel(DurMin)
        p_ij = p_fit;

        % Apply tuned component-3 rates directly (absolute values)
        p_ij.k_ar(3)  = p_fit.tuned.exp2.k_ar3(idx);
        p_ij.k_gar(3) = p_fit.tuned.exp2.k_gar3(idx);

        p_ij.T = DurMin(j)*60;

        [tSim,YSim,ok] = simulate_IFFL(p_ij, n, m);
        if ~ok, continue; end

        frac = YSim(:,chIdx) ./ p_ij.G_tot(denomIdx);

        [t_clip_h, f_clip] = clip_for_auc(tSim, frac, AUC_tStart_s, AUC_tEnd_s, Tfinal_s);
        if isempty(t_clip_h), continue; end

        AUC_sim_exp2(i,j) = calculate_AUC(t_clip_h, f_clip);
    end
end

% --- Exp3: rows = toehold [8 6 4], cols = Din [20 40 60]
exp3_toe = AUC.exp3.toehold_nt(:).';
AUC_sim_exp3 = nan(numel(exp3_toe), numel(DurMin));
for i = 1:numel(exp3_toe)
    idx = find(p_fit.tuned.exp3.toehold == exp3_toe(i), 1);
    if isempty(idx), error('Exp3: toehold %d not found in p_fit.tuned.exp3.toehold', exp3_toe(i)); end

    for j = 1:numel(DurMin)
        p_ij = p_fit;

        % Apply tuned component-3 rates directly (absolute values)
        p_ij.k_bc(3)  = p_fit.tuned.exp3.k_bc3(idx);
        p_ij.k_gbc(3) = p_fit.tuned.exp3.k_gbc3(idx);

        p_ij.T = DurMin(j)*60;

        [tSim,YSim,ok] = simulate_IFFL(p_ij, n, m);
        if ~ok, continue; end

        frac = YSim(:,chIdx) ./ p_ij.G_tot(denomIdx);

        [t_clip_h, f_clip] = clip_for_auc(tSim, frac, AUC_tStart_s, AUC_tEnd_s, Tfinal_s);
        if isempty(t_clip_h), continue; end

        AUC_sim_exp3(i,j) = calculate_AUC(t_clip_h, f_clip);
    end
end

%% ======== 4. Plot  ======================================================

% ============================================================
% Figure 1.5: (Tile 1) TS-derived AUC  + (Tile 2) Exp0.5 dose sweep
% ============================================================
fig15 = figure('Name','Fig1.5: TS-derived AUC + Exp0.5 dose sweep');
tl15  = tiledlayout(fig15, 1, 2, 'TileSpacing','compact');

% -------------------------
% (Tile 1) TS-based AUC  (moved from Fig2)
% -------------------------
ax = nexttile(tl15); hold(ax,'on'); box(ax,'on'); grid(ax,'off');

% --- compute TS-based AUC (DATA) from mean traces in dataTS (Ch3) ---
DinMin_TS = arrayfun(@(c) round(c.T_sec/60), dataTS.conds);
[DinMin_TS, ord] = sort(DinMin_TS);
conds_sorted = dataTS.conds(ord);

% --- (NEW) x-axis "break" mapping for Tile 1 only ---
x_raw = DinMin_TS;

x0    = 70;    % break starts (left visible end)
x1    = 220;   % break ends (right visible start; no data expected here)
gap   = 6;     % visual gap width (in x units)
scale = 0.15;  % compression factor for the right segment

x_plot = x_raw;
maskR  = (x_raw >= x1);
x_plot(maskR) = x0 + gap + (x_raw(maskR) - x1)*scale;

% optional: if any points fall inside (x0, x1), they will appear in the "gap";
% (you likely have none, so no special handling needed)

AUC_ts_data = nan(size(DinMin_TS));
AUC_ts_sim  = nan(size(DinMin_TS));

for k = 1:numel(conds_sorted)
    ck = conds_sorted(k);

    % DATA AUC (from mean trace)
    t_h_data = ck.time_s(:)/3600;
    y_data   = ck.y_ch3(:);
    AUC_ts_data(k) = calculate_AUC(t_h_data, y_data);

    % SIM AUC (Ch3 fraction ON)
    p_ij = p_fit;
    p_ij.T = ck.T_sec;
    [tSim,YSim,ok] = simulate_IFFL(p_ij, n, m);
    if ok
        frac3 = YSim(:,5) ./ p_ij.G_tot(5);
        AUC_ts_sim(k) = calculate_AUC(tSim(:)/3600, frac3(:));
    end
end

% --- colors: keep red ramp for this tile (nominally Ch3 summary) ---
C = make_hex_ramp(["#fed976","#b10026"], numel(DinMin_TS));  % light->dark

% Use provided Exp0 mean/std if sizes match
AUC_ts_data_mean = AUC_ts_data;
AUC_ts_data_std  = nan(size(AUC_ts_data));
if isfield(AUC,'exp0') && isfield(AUC.exp0,'mean_5nM') && isfield(AUC.exp0,'std_5nM') && ...
        numel(AUC.exp0.mean_5nM) == numel(DinMin_TS) && numel(AUC.exp0.std_5nM) == numel(DinMin_TS)
    AUC_ts_data_mean = AUC.exp0.mean_5nM(:).';
    AUC_ts_data_std  = AUC.exp0.std_5nM(:).';
end

% ±1σ shadow + mean line (black), then colored filled points
% --- (NEW) do NOT connect across 60->240 ---
idx_conn = (DinMin_TS <= 60);
idx_iso  = (DinMin_TS >= 240);

% Connected segment only (<=60): shadow + line
plot_auc_single_curve_with_std(ax, x_plot(idx_conn), AUC_ts_data_mean(idx_conn), AUC_ts_data_std(idx_conn), ...
    [0 0 0], Plot.ExpLineWidth, Plot.AUC_MarkerSize_exp, Plot.AUC_ShadowAlpha, ...
    'Data (5 nM nominal)');

% Isolated point (>=240): add a small vertical "std band" with finite x-width
if any(idx_iso)
    xi = x_plot(idx_iso);
    yi = AUC_ts_data_mean(idx_iso);
    si = AUC_ts_data_std(idx_iso);

    w  = 4;  % half-width of the band in x-units (tweak for aesthetics)
    xx = [xi-w, xi+w, xi+w, xi-w];
    yy = [yi+si, yi+si, yi-si, yi-si];

    fill(ax, xx, yy, [0 0 0], ...
        'FaceAlpha', Plot.AUC_ShadowAlpha, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    % then the point + errorbar (no connecting line)
    errorbar(ax, xi, yi, si, 'o', ...
        'LineWidth', Plot.ExpLineWidth, ...
        'MarkerSize', Plot.AUC_MarkerSize_exp, ...
        'Color', [0 0 0], ...
        'MarkerFaceColor', [0 0 0], ...
        'MarkerEdgeColor', [0 0 0], ...
        'HandleVisibility','off');
end


% Isolated point (>=240): show as errorbar only (no connecting line)
if any(idx_iso)
    errorbar(ax, x_plot(idx_iso), AUC_ts_data_mean(idx_iso), AUC_ts_data_std(idx_iso), ...
        'o', 'LineWidth', Plot.ExpLineWidth, 'MarkerSize', Plot.AUC_MarkerSize_exp, ...
        'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], ...
        'HandleVisibility','off');
end


for k = 1:numel(DinMin_TS)
    plot(ax, x_plot(k), AUC_ts_data_mean(k), 'o', ...
        'LineWidth', Plot.ExpLineWidth, ...
        'MarkerSize', Plot.AUC_MarkerSize_exp, ...
        'Color', C(k,:), ...
        'MarkerFaceColor', C(k,:), ...
        'MarkerEdgeColor', C(k,:), ...
        'HandleVisibility','off');
end

% SIM: unfilled + dashed connector (hidden from legend)
for k = 1:numel(DinMin_TS)
    plot(ax, x_plot(k), AUC_ts_sim(k), 'o', ...
        'LineWidth', Plot.SimLineWidth, ...
        'MarkerSize', Plot.AUC_MarkerSize_sim, ...
        'Color', C(k,:), ...
        'MarkerFaceColor','none', ...
        'MarkerEdgeColor', C(k,:));
end
% Sim connector: only connect <=60 (do not connect 60->240)
plot(ax, x_plot(idx_conn), AUC_ts_sim(idx_conn), '--', ...
    'LineWidth', Plot.SimLineWidth, ...
    'Color', [0 0 0], ...
    'HandleVisibility','off');

title(ax, 'TS-derived AUC (Ch3): mean trace vs sim', 'Interpreter','none');
xlabel(ax, '$D_{in}$ [min]', 'Interpreter','latex');
ylabel(ax, 'AUC (Ch3, fraction \cdot h)');
set(ax,'TickDir','in','FontSize',Plot.FontSize_Fig2);
% --- y-axis style: match Fig2, but x-axis is custom (broken) ---
ylim(ax, Plot.AUC_ylim);
yticks(ax, Plot.AUC_yticks);
box(ax,'on'); grid(ax,'off');
ax.LineWidth = Plot.AxesLineWidth;
ax.TickDir   = 'in';

% --- x-axis custom ticks/labels (show raw Din values) ---
x240_plot = x0 + gap + (240 - x1)*scale;
xticks(ax, [20 40 60 x240_plot]);
xticklabels(ax, {'20','40','60','240'});

xlim(ax, [0, x240_plot + 6]);  % margin on the right

% --- draw axis break marker ("//") near the bottom of axis ---
yB = Plot.AUC_ylim(1) + 0.06*range(Plot.AUC_ylim);  % a bit above bottom
xB = x0 + gap/2;
dx = 1.0;  dy = 0.10;  % marker size (tweak if desired)

plot(ax, [xB-dx xB-dx/3], [yB-dy yB+dy], 'k-', 'LineWidth', Plot.AxesLineWidth);
plot(ax, [xB+dx/3 xB+dx], [yB-dy yB+dy], 'k-', 'LineWidth', Plot.AxesLineWidth);


% -------------------------
% (Tile 2) Exp0.5: dose sweep (NEW)
% -------------------------
ax = nexttile(tl15); hold(ax,'on'); box(ax,'on'); grid(ax,'off');

% Experimental dose conditions (categorical x)
dose_labels = ["10 nM G8C1 / 2 nM G8C3", ...
               "25 nM G8C1 / 5 nM G8C3", ...
               "50 nM G8C1 / 10 nM G8C3"];
xDose = [1,2.5,5];

% DATA (mean ± std) — same helper, same style
plot_auc_single_curve_with_std(ax, xDose, AUC.exp0.mean_dose, AUC.exp0.std_dose, ...
    [0 0 0], Plot.ExpLineWidth, Plot.AUC_MarkerSize_exp, Plot.AUC_ShadowAlpha, ...
    'Data');

% SIM (dashed, unfilled markers) — choose Din for this dose-sweep simulation
% NOTE: Din for Exp0.5 was not specified in your message; here we assume Din = 20 min.
Din_for_doseSweep_min = 30; 
AUC_sim_dose = nan(1,3);

G2C1_list = [10 25 50];
G2C3_list = [ 2  5 10];

for k = 1:3
    p_ij = override_input_G2C1_G2C3(p_fit, G2C1_list(k), G2C3_list(k));
    p_ij.T = Din_for_doseSweep_min * 60;

    % Apply alpha ONLY for Exp0.5 simulation (local copy; no side effects)
    p_ij = apply_alpha_kp_scale_local(p_ij, alpha_exp05);

    [tSim,YSim,ok] = simulate_IFFL(p_ij, n, m);
    if ok
        frac = YSim(:,chIdx) ./ p_ij.G_tot(denomIdx);
        [t_clip_h, f_clip] = clip_for_auc(tSim, frac, AUC_tStart_s, AUC_tEnd_s, Tfinal_s);
        if ~isempty(t_clip_h)
            AUC_sim_dose(k) = calculate_AUC(t_clip_h, f_clip);
        end
    end
end


plot(ax, xDose, AUC_sim_dose, '--o', ...
    'LineWidth', Plot.SimLineWidth, ...
    'MarkerSize', Plot.AUC_MarkerSize_sim, ...
    'Color', [0 0 0], ...
    'MarkerFaceColor','none', ...
    'HandleVisibility','off');

% Axes formatting (y-axis style matches Fig2; x is categorical)
set(ax,'XLim',[0 6], 'XTick',xDose, 'XTickLabel',dose_labels);
ax.XTickLabelRotation = 20;

ylim(ax, Plot.AUC_ylim);
yticks(ax, Plot.AUC_yticks);
ax.LineWidth = Plot.AxesLineWidth;
ax.TickDir   = 'in';
set(ax,'FontSize',Plot.FontSize_Fig2);

title(ax, sprintf('Exp0.5: dose sweep (sim @ D_{in}=%d min)', Din_for_doseSweep_min), 'Interpreter','none');
xlabel(ax, 'Input node dose condition');
ylabel(ax, sprintf('AUC (%s, fraction \\cdot h)', chLabel));

legend(ax,'Location','best','Box','off');

title(tl15, 'Fig1.5: TS-derived AUC + Exp0.5 dose sweep', 'Interpreter','none');


% ============================================================
% Figure 2: Exp1/Exp2/Exp3 only (Tile 1 removed)
% ============================================================
fig = figure('Name','Figure 2: AUC vs Din (Exp1/Exp2/Exp3)');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing','compact');

% =========================
% (Tile 1) Exp1
% =========================
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'off');

labels_exp1 = ["5 nM (nominal)","10 nM","25 nM"];
plot_auc_experiment_matrix_with_std(ax, DurMin, AUC.exp1.Y, AUC.exp1.Ystd, Plot.ColorOrder_exp1, ...
    Plot.ExpLineWidth, Plot.AUC_MarkerSize_exp, Plot.AUC_ShadowAlpha, labels_exp1);
plot_auc_sim_matrix(ax, DurMin, AUC_sim_exp1, Plot.ColorOrder_exp1, ...
    Plot.SimLineWidth, Plot.AUC_MarkerSize_sim);

legend(ax,'Location','best','Box','off');

title(ax, 'Exp1: sweep G3R1 (5/10/25 nM)', 'Interpreter','none');
xlabel(ax, '$D_{in}$ [min]', 'Interpreter','latex');
ylabel(ax, sprintf('AUC (%s, fraction \\cdot h)', chLabel));
set(ax,'TickDir','in','FontSize',Plot.FontSize_Fig2);
apply_auc_axes_style(ax, Plot);

% =========================
% (Tile 2) Exp2
% =========================
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'off');

labels_exp23 = ["8 nt (nominal)","6 nt","4 nt"];
plot_auc_experiment_matrix_with_std(ax, DurMin, AUC.exp2.Y, AUC.exp2.Ystd, Plot.ColorOrder_exp23, ...
    Plot.ExpLineWidth, Plot.AUC_MarkerSize_exp, Plot.AUC_ShadowAlpha, labels_exp23);
plot_auc_sim_matrix(ax, DurMin, AUC_sim_exp2, Plot.ColorOrder_exp23, ...
    Plot.SimLineWidth, Plot.AUC_MarkerSize_sim);

legend(ax,'Location','best','Box','off');

title(ax, 'Exp2: toehold dA1 (8/6/4 nt)', 'Interpreter','none');
xlabel(ax, '$D_{in}$ [min]', 'Interpreter','latex');
ylabel(ax, sprintf('AUC (%s, fraction \\cdot h)', chLabel));
set(ax,'TickDir','in','FontSize',Plot.FontSize_Fig2);
apply_auc_axes_style(ax, Plot);

% =========================
% (Tile 3) Exp3
% =========================
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'off');

plot_auc_experiment_matrix_with_std(ax, DurMin, AUC.exp3.Y, AUC.exp3.Ystd, Plot.ColorOrder_exp23, ...
    Plot.ExpLineWidth, Plot.AUC_MarkerSize_exp, Plot.AUC_ShadowAlpha, labels_exp23);
plot_auc_sim_matrix(ax, DurMin, AUC_sim_exp3, Plot.ColorOrder_exp23, ...
    Plot.SimLineWidth, Plot.AUC_MarkerSize_sim);

legend(ax,'Location','best','Box','off');

title(ax, 'Exp3: toehold dB1 (8/6/4 nt)', 'Interpreter','none');
xlabel(ax, '$D_{in}$ [min]', 'Interpreter','latex');
ylabel(ax, sprintf('AUC (%s, fraction \\cdot h)', chLabel));
set(ax,'TickDir','in','FontSize',Plot.FontSize_Fig2);
apply_auc_axes_style(ax, Plot);

title(tl, 'Figure 2: AUC vs input duration (Exp1/Exp2/Exp3)', 'Interpreter','none');


%% ======== 5. Print quick tables ========================================
fprintf('\n===== AUC tables (SIM) =====\n');
fprintf('Exp1 rows=[G3R1 5/10/25], cols=[Din 20/40/60]\n'); disp(AUC_sim_exp1);
fprintf('Exp2 rows=[toe 8/6/4], cols=[Din 20/40/60]\n'); disp(AUC_sim_exp2);
fprintf('Exp3 rows=[toe 8/6/4], cols=[Din 20/40/60]\n'); disp(AUC_sim_exp3);

%% ==================== Local helper functions ============================

function apply_auc_axes_style(ax, Plot)
% Apply Figure-2 specific axis settings (no grid lines, fixed limits/ticks, and frame thickness)
    xlim(ax, Plot.AUC_xlim);
    ylim(ax, Plot.AUC_ylim);
    yticks(ax, Plot.AUC_yticks);

    box(ax,'on');
    grid(ax,'off');

    ax.LineWidth = Plot.AxesLineWidth;
    ax.TickDir   = 'in';

    % Note: MATLAB does not provide a separate "tick thickness" property.
    % Tick thickness follows ax.LineWidth visually, so TickLineWidth is kept
    % as a parameter for compatibility.
end

function [t_h, y] = clip_for_auc(t_s, y, tStart_s, tEnd_s, tfinal_s)
% Clip (t,y) to [tStart_s, tEnd_s] (seconds), then convert time to hours.
    if isempty(tEnd_s), tEnd_s = tfinal_s; end
    mask = (t_s >= tStart_s) & (t_s <= tEnd_s);
    if ~any(mask)
        t_h = []; y = [];
        return;
    end
    t_h = t_s(mask) / 3600;
    y   = y(mask);
end

function I = calculate_AUC(t_h, y)
% Right Riemann sum in hour units:
% I = sum_{i=2..N} (t_i - t_{i-1}) * y_i
    I = 0;
    for i = 2:length(t_h)
        I = I + (t_h(i)-t_h(i-1)) * y(i);
    end
end

function plot_auc_experiment_matrix(ax, DurMin, Y, C, lw, ms, labels)
% Plot experimental AUC curves (solid lines) with FILLED markers + legend labels.
    for i = 1:size(Y,1)
        plot(ax, DurMin, Y(i,:), '-o', ...
            'LineWidth', lw, ...
            'MarkerSize', ms, ...
            'Color', C(i,:), ...
            'MarkerFaceColor', C(i,:), ...
            'MarkerEdgeColor', C(i,:), ...
            'DisplayName', labels(i));
    end
end

function plot_auc_sim_matrix(ax, DurMin, Ysim, C, lw, ms)
% Plot simulated AUC curves (dashed) but hide them from legend.
    for i = 1:size(Ysim,1)
        plot(ax, DurMin, Ysim(i,:), '--o', ...
            'LineWidth', lw, ...
            'MarkerSize', ms, ...
            'Color', C(i,:), ...
            'MarkerFaceColor', 'none', ...
            'HandleVisibility', 'off');  % <-- legendには出さない
    end
end

function [chIdx, denomIdx, chLabel] = map_channel(Channel)
% Map channel specification to state index and denominator index.
% In this model:
%   Ch1 -> G_dA(1) / G_tot(1)  (G2C1)
%   Ch2 -> G_dA(3) / G_tot(3)  (G3R1)
%   Ch3 -> G_dA(5) / G_tot(5)  (G1S1)
    if ischar(Channel) || isstring(Channel)
        switch lower(string(Channel))
            case "ch1", chIdx = 1; denomIdx = 1; chLabel = 'Ch1 (G2C1)';
            case "ch2", chIdx = 3; denomIdx = 3; chLabel = 'Ch2 (G3R1)';
            case "ch3", chIdx = 5; denomIdx = 5; chLabel = 'Ch3 (G1S1)';
            otherwise, error('Unknown Channel = %s', Channel);
        end
    else
        switch Channel
            case 1, chIdx = 1; denomIdx = 1; chLabel = 'Ch1 (G2C1)';
            case 2, chIdx = 3; denomIdx = 3; chLabel = 'Ch2 (G3R1)';
            case 3, chIdx = 5; denomIdx = 5; chLabel = 'Ch3 (G1S1)';
            otherwise, error('Channel must be 1,2,3 or ''ch1|ch2|ch3''.');
        end
    end
end

function p = override_total_G3R1_and_dB3(p, G3R1_new_nM)
% Override G3R1 total (nM) and co-increase dB3 total (nM) by the same delta
% so that (dB3_tot - G3R1_tot) stays constant.

    G3R1_old = p.aux.G3R1_tot;
    dB3_old  = p.aux.dB3_tot;

    dG = G3R1_new_nM - G3R1_old;
    dB3_new_nM = dB3_old + dG;

    p.aux.G3R1_tot = G3R1_new_nM;
    p.aux.dB3_tot  = dB3_new_nM;

    % rebuild vectors
    G2C1 = p.aux.G2C1_tot;
    G2C3 = p.aux.G2C3_tot;
    G3R1 = p.aux.G3R1_tot;
    G1R3 = p.aux.G1R3_tot;
    G1S1 = p.aux.G1S1_tot;

    dB2 = p.aux.dB2_tot;
    dB3 = p.aux.dB3_tot;
    dB1 = p.aux.dB1_tot;

    p.G_tot  = [G2C1 G2C3 G3R1 G1R3 G1S1]' * 1e-9;
    p.dB_tot = [dB2  dB3  dB1 ]' * 1e-9;
end

function AUC = build_embedded_AUC_data(Din_vals)
% Build a struct with embedded AUC mean/std arrays for Exp1/2/3.
% All are measured at Din=[20 40 60] (minutes).
    AUC = struct();
    AUC.Din_vals = Din_vals(:).';

    % ======================
    % Exp0: Time-series data
    % ======================
    AUC.exp0.mean_5nM  = [2.74930579412568	2.66620579066230	2.28017948042252	1.39676542435326	1.08574275521072   0.695213494830666];
    AUC.exp0.std_5nM   = [0.100788511315339	0.484170804678733	0.904889362984810	0.457328322955851	0.591191023914664   0.136386991174033];

    % ======================
    % Exp0.5: dose sweep
    % ======================
    AUC.exp0.mean_dose  = [2.72745037326219	2.64401754542818	1.24894291418106];
    AUC.exp0.std_dose   = [0.0775436965437986	0.618621802409835	0.550975215004794];

    % ======================
    % Exp1: G3R1 sweep (5,10,25 nM)
    % ======================
    AUC.exp1.levels_nM = [5 10 25];
    AUC.exp1.G3R1_levels_nM = AUC.exp1.levels_nM;  % backward compatibility

    AUC.exp1.mean_5nM  = [2.56499499837718, 1.96172243338906, 0.864853882077734];
    AUC.exp1.std_5nM   = [0.503811562305211, 0.771135745740686, 0.303357578619457];

    AUC.exp1.mean_10nM = [1.57539880189863, 0.662969895854546, 0.294228900361590];
    AUC.exp1.std_10nM  = [0.335147554254185, 0.576005047016623, 0.129556181103264];

    AUC.exp1.mean_25nM = [0.498801084899379, 0.234218003967067, 0.335176590776887];
    AUC.exp1.std_25nM  = [0.114438607347594, 0.246755060480399, 0.0379034435609612];

    % ======================
    % Exp2: dA1 toehold length = 8,6,4 nt (nominal=8)
    % ======================
    AUC.exp2.toehold_nt = [8 6 4];

    AUC.exp2.mean_8b = [2.62112235459656, 1.75124015886793, 1.38857435443720];
    AUC.exp2.std_8b  = [0.0812330403411200, 0.371381421328752, 0.466800671479998];

    AUC.exp2.mean_6b = [3.03573082836344, 2.90763225699382, 2.45183187679066];
    AUC.exp2.std_6b  = [0.141605278494657, 0.104876195918514, 0.253786613030568];

    AUC.exp2.mean_4b = [3.09994898455949, 2.97298907575560, 2.71276629624907];
    AUC.exp2.std_4b  = [0.102625364660198, 0.125931764802170, 0.182025433716606];

    % ======================
    % Exp3: dB1 toehold length = 8,6,4 nt (nominal=8)
    % ======================
    AUC.exp3.toehold_nt = [8 6 4];

    AUC.exp3.mean_8b = [2.56173667051312, 2.09954690172282, 0.543000279185901];
    AUC.exp3.std_8b  = [0.194379133634788, 0.298704078162605, 0.151383327353484];

    AUC.exp3.mean_6b = [0.523584825432856, 0.487934703430427, -0.0421065906067014];
    AUC.exp3.std_6b  = [0.308899901539074, 0.144229025818429, 0.256954930837777];

    AUC.exp3.mean_4b = [0.184195637493338, -0.232970384315233, -0.0902961430380806];
    AUC.exp3.std_4b  = [0.191360767844190, 0.332762080709247, 0.0852841374464548];

    % --- Package matrices (rows x cols = 3 x 3) ---
    AUC.exp1.Y    = [AUC.exp1.mean_5nM;  AUC.exp1.mean_10nM; AUC.exp1.mean_25nM];
    AUC.exp1.Ystd = [AUC.exp1.std_5nM;   AUC.exp1.std_10nM;  AUC.exp1.std_25nM];

    AUC.exp2.Y    = [AUC.exp2.mean_8b;   AUC.exp2.mean_6b;   AUC.exp2.mean_4b];
    AUC.exp2.Ystd = [AUC.exp2.std_8b;    AUC.exp2.std_6b;    AUC.exp2.std_4b];

    AUC.exp3.Y    = [AUC.exp3.mean_8b;   AUC.exp3.mean_6b;   AUC.exp3.mean_4b];
    AUC.exp3.Ystd = [AUC.exp3.std_8b;    AUC.exp3.std_6b;    AUC.exp3.std_4b];
end

function C = make_hex_ramp(hexEndpoints, n)
% make a color ramp between two hex colors, returns n x 3 RGB
% hexEndpoints: ["#aaaaaa","#bbbbbb"] (2 colors)
    if n <= 1
        C = hex2rgb_mat(hexEndpoints(1));
        return;
    end
    c0 = hex2rgb_mat(hexEndpoints(1));
    c1 = hex2rgb_mat(hexEndpoints(2));
    t = linspace(0,1,n).';
    C = (1-t).*c0 + t.*c1;
end

function rgb = hex2rgb_mat(hex)
% "#RRGGBB" -> 1x3 double
    hex = char(hex);
    if hex(1) == '#', hex = hex(2:end); end
    r = hex2dec(hex(1:2));
    g = hex2dec(hex(3:4));
    b = hex2dec(hex(5:6));
    rgb = [r g b]/255;
end

%% ==================== Model functions (copied) ==========================
function [t, Y, ok] = simulate_IFFL(p, n, m)
% Simulate 4h dynamics with a pulse duration p.T (sec).
    ok = true;
    tfinal = 14400;  % 4 h
    Tpulse = p.T;
    if ~isfinite(Tpulse) || Tpulse < 0, Tpulse = inf; end
    t1_end = min(Tpulse, tfinal);

    Y0 = initial_conditions(p)*1e-9;

    AbsTol_vec = [
        repmat(1e-10,m,1);  % G_dA
        repmat(1e-10,m,1);  % G_dB
        repmat(1e-10,n,1);  % dA
        repmat(1e-10,n,1);  % dB
        repmat(1e-9 ,n,1);  % rC
        repmat(1e-9 ,n,1);  % rR
        repmat(1e-9 ,n,1);  % dR
        repmat(1e-9 ,n,1);  % dA_dR
        repmat(1e-9 ,n,1);  % rI
        1e-9                % S
    ];

    neq = 2*m + 7*n + 1;  % number of ODE states
    odeOpts = odeset('RelTol',1e-7,'AbsTol',AbsTol_vec,'NonNegative',1:neq,'MaxStep',100);

    p1 = p;
    f1 = @(t,y) General_genelet_model_segment(t,y,p1,n,m);
    [t1, Y1] = ode15s(f1, [0 t1_end], Y0, odeOpts);
    if any(~isfinite(Y1(:))) || isempty(t1), t=[]; Y=[]; ok=false; return; end
    if t1_end >= tfinal, t=t1; Y=Y1; return; end

    p2 = p;
    if ~isfield(p2,'added_dR_amount') || isempty(p2.added_dR_amount)
        p2.added_dR_amount = 250e-9;
    end
    Y0_2 = Y1(end,:).';

    % dR starts after [G_dA(m), G_dB(m), dA(n), dB(n), rC(n), rR(n)] => offset = 2m+4n
    dR1_idx = 2*m + 4*n + 1;  % dR(1)

    Y0_2(dR1_idx) = Y0_2(dR1_idx) + p2.added_dR_amount;

    f2 = @(t,y) General_genelet_model_segment(t,y,p2,n,m);
    [t2, Y2] = ode15s(f2, [t1_end tfinal], Y0_2, odeOpts);
    if any(~isfinite(Y2(:))) || isempty(t2), t=[]; Y=[]; ok=false; return; end

    t = [t1; t2(2:end)];
    Y = [Y1; Y2(2:end,:)];
end

function dYdt = General_genelet_model_segment(t,Y,p,n,m) %#ok<INUSD>
    k_gar   = p.k_gar(:);
    k_ar    = p.k_ar(:);
    k_pr    = p.k_pr(:);
    k_ir    = p.k_ir(:);
    k_dr    = p.k_dr(:);
    k_ga    = p.k_ga(:);
    k_gb    = p.k_gb(:);
    k_gab   = p.k_gab(:);
    k_gbc   = p.k_gbc(:);
    k_bc    = p.k_bc(:);
    k_pc    = p.k_pc(:);
    k_dc    = p.k_dc(:);
    k_pi    = p.k_pi(:);
    k_gar_d = p.k_gar_d(:);

    idx_input = [1];
    k_aa_vec = p.k_aa_other * ones(n,1);
    k_ai_vec = p.k_aa_other * ones(n,1);
    k_aa_vec(idx_input) = p.k_aa_input;
    k_ai_vec(idx_input) = p.k_aa_input;

    G_tot  = p.G_tot(:);
    dA_tot = p.dA_tot(:);
    dB_tot = p.dB_tot(:);

    N = 0;
    G_dA   = Y(N+1:N+m); N = N + m;
    G_dB   = Y(N+1:N+m); N = N + m;
    dA     = Y(N+1:N+n); N = N + n;
    dB     = Y(N+1:N+n); N = N + n;
    rC     = Y(N+1:N+n); N = N + n;
    rR     = Y(N+1:N+n); N = N + n;
    dR     = Y(N+1:N+n); N = N + n;
    dA_dR  = Y(N+1:N+n); N = N + n;
    rI     = Y(N+1:N+n); N = N + n;
    S      = Y(N+1);

    M_ac = [1 1 0 0 0; 0 0 1 0 0; 0 0 0 1 1];
    M_rp = [0 0 0 0 0; 0 0 0 0 0; 0 0 1 0 0];
    M_bc = [1 1 0 0 0; 0 0 1 0 0; 0 0 0 1 1];
    M_cp = [0 0 0 0 0; 0 1 0 0 0; 1 0 0 0 0];
    M_ip = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 1 0];

    leak = p.leak;

    dA_rR = dA_tot - dA - M_ac*G_dA - dA_dR;
    dB_rC = dB_tot - dB - M_bc*G_dB;
    G     = G_tot  - G_dA - G_dB;

    dG_dA = (M_ac'*(k_ga.*dA)).*G ...
          - (M_ac'*(k_gar  .*rR)).*G_dA ...
          - (M_ac'*(k_gar_d.*dR)).*G_dA ...
          - (M_bc'*(k_gab.*dB)).*G_dA;

    dG_dB = (M_bc'*(k_gb.*dB)).*G ...
          + (M_bc'*(k_gab.*dB)).*G_dA ...
          - (M_bc'*(k_gbc.*rC)).*G_dB;

    ddA   = k_dr.*dA_rR ...
          - k_ga.*(M_ac*G).*dA ...
          - k_ar.*dA.*rR ...
          - k_ar.*dA.*dR ...
          + k_gab.*(M_bc*G_dA).*dB;

    ddB   = k_dc.*dB_rC ...
          - k_gb.*(M_bc*G).*dB ...
          - k_bc.*dB.*rC ...
          - k_gab.*(M_bc*G_dA).*dB;

    drC   = - k_gbc.*(M_bc*G_dB).*rC ...
          -   k_bc.*dB.*rC ...
          +   M_cp*(k_pc.*G_dA) + M_cp*(leak.*k_pc.*G_dB) ...
          +   k_aa_vec.*dB;

    drR   = - k_gar.*(M_ac*G_dA).*rR ...
          -   k_ar.*dA.*rR ...
          +   M_rp*(k_pr.*G_dA) - k_ir.*rR.*rI + M_rp*(leak.*k_pr.*G_dB) ...
          +   k_ai_vec.*dA;

    ddR   = - k_gar_d.*(M_ac*G_dA).*dR ...
          -   k_ar.*dA.*dR;

    ddA_dR = k_ar.*dA.*dR + k_gar_d.*(M_ac*G_dA).*dR;

    drI   = - k_ir.*rR.*rI + M_ip*(k_pi.*G_dA) + M_ip*(leak.*k_pi.*G_dB);
    dS    = k_pc(m)*G_dA(m) - k_dc(n)*S;

    dYdt = [dG_dA; dG_dB; ddA; ddB; drC; drR; ddR; ddA_dR; drI; dS];
end

function Y0 = initial_conditions(p)
    G2C1_tot = p.aux.G2C1_tot; G2C3_tot = p.aux.G2C3_tot; G3R1_tot = p.aux.G3R1_tot;
    G1R3_tot = p.aux.G1R3_tot; G1S1_tot = p.aux.G1S1_tot;
    dA2_tot  = p.aux.dA2_tot;  dA3_tot  = p.aux.dA3_tot;  dA1_tot  = p.aux.dA1_tot;
    dB2_tot  = p.aux.dB2_tot;  dB3_tot  = p.aux.dB3_tot;  dB1_tot  = p.aux.dB1_tot;

    Y0 = [
        min(dA2_tot*G2C1_tot/(G2C1_tot+G2C3_tot), max(0,G2C1_tot-dB2_tot*G2C1_tot/(G2C1_tot+G2C3_tot)) );
        min(dA2_tot*G2C3_tot/(G2C1_tot+G2C3_tot), max(0,G2C3_tot-dB2_tot*G2C3_tot/(G2C1_tot+G2C3_tot)) );
        min(dA3_tot,  max(0,G3R1_tot-dB3_tot));
        min(dA1_tot*G1R3_tot/(G1R3_tot+G1S1_tot), max(0,G1R3_tot-dB1_tot*G1R3_tot/(G1R3_tot+G1S1_tot)) );
        min(dA1_tot*G1S1_tot/(G1R3_tot+G1S1_tot), max(0,G1S1_tot-dB1_tot*G1S1_tot/(G1R3_tot+G1S1_tot)) );

        min(G2C1_tot,dB2_tot*G2C1_tot/(G2C1_tot+G2C3_tot));
        min(G2C3_tot,dB2_tot*G2C3_tot/(G2C1_tot+G2C3_tot));
        min(G3R1_tot,dB3_tot);
        min(G1R3_tot,dB1_tot*G1R3_tot/(G1R3_tot+G1S1_tot));
        min(G1S1_tot,dB1_tot*G1S1_tot/(G1R3_tot+G1S1_tot));

        max(0,dA2_tot-max(0,(G2C1_tot+G2C3_tot)-dB2_tot));
        max(0,dA3_tot-max(0,G3R1_tot-dB3_tot));
        max(0,dA1_tot-max(0,(G1R3_tot+G1S1_tot)-dB1_tot));

        max(dB2_tot-G2C1_tot-G2C3_tot,0);
        max(dB3_tot-G3R1_tot,0);
        max(dB1_tot-G1R3_tot-G1S1_tot,0);

        0;0;0;
        0;0;0;
        0;0;0;
        0;0;0;
        0;0;0;
        0
    ];
end

function data = read_iffl_dir(dirpath, varargin)
    p = inputParser;
    addParameter(p,'Pattern','Din_%dmin_means.csv');
    addParameter(p,'T_list_min',[]);
    addParameter(p,'TimeCol','time_h');
    addParameter(p,'Ch1Col','ON_mean_ch1');
    addParameter(p,'Ch2Col','ON_mean_ch2');
    addParameter(p,'Ch3Col','ON_mean_ch3');
    addParameter(p,'ClipZero',false);
    parse(p,varargin{:});

    conds = struct('T_sec',{},'time_s',{},'y_ch1',{},'y_ch2',{},'y_ch3',{});

    if isempty(p.Results.T_list_min)
        D = dir(fullfile(dirpath,'Din_*min*.csv'));
        if isempty(D), error('No matching CSVs found in %s',dirpath); end
        Tmins = [];
        for k = 1:numel(D)
            tok = regexp(D(k).name,'Din_(\d+)min','tokens','once');
            if ~isempty(tok), Tmins(end+1) = str2double(tok{1}); end %#ok<AGROW>
        end
        Tmins = unique(Tmins);
    else
        Tmins = p.Results.T_list_min(:)';
    end

    idx = 1;
    for Tm = Tmins
        fname = sprintf(p.Results.Pattern,Tm);
        fpath = fullfile(dirpath,fname);
        if ~exist(fpath,'file')
            warning('File not found for T=%d min: %s (skipped)',Tm,fname);
            continue
        end
        T = readtable(fpath);
        time_s = T.(p.Results.TimeCol)*3600;
        y1 = T.(p.Results.Ch1Col);
        y2 = T.(p.Results.Ch2Col);
        y3 = T.(p.Results.Ch3Col);
        if p.Results.ClipZero
            y1 = max(y1,0); y2 = max(y2,0); y3 = max(y3,0);
        end
        conds(idx).T_sec   = Tm*60; %#ok<AGROW>
        conds(idx).time_s  = time_s; %#ok<AGROW>
        conds(idx).y_ch1   = y1;     %#ok<AGROW>
        conds(idx).y_ch2   = y2;     %#ok<AGROW>
        conds(idx).y_ch3   = y3;     %#ok<AGROW>
        idx = idx+1;
    end

    if isempty(conds)
        error('No conditions loaded. Check directory/pattern/columns.');
    end
    data.conds = conds;
end

function plot_auc_experiment_matrix_with_std(ax, DurMin, Ymean, Ystd, C, lw, ms, alphaFill, labels)
% Plot experimental AUC curves with ±1σ shadow + mean line (FILLED markers).
% Ymean, Ystd: (nRows x nCols)
    for i = 1:size(Ymean,1)
        x = DurMin(:).';
        y = Ymean(i,:);
        s = Ystd(i,:);
        if any(isnan(y)) || any(isnan(s)), continue; end

        ylo = y - s;
        yhi = y + s;

        % shadow band (no legend)
        xx = [x, fliplr(x)];
        yy = [yhi, fliplr(ylo)];
        fill(ax, xx, yy, C(i,:), ...
            'FaceAlpha', alphaFill, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % mean line + filled markers (legend entry)
        plot(ax, x, y, '-o', ...
            'LineWidth', lw, ...
            'MarkerSize', ms, ...
            'Color', C(i,:), ...
            'MarkerFaceColor', C(i,:), ...
            'MarkerEdgeColor', C(i,:), ...
            'DisplayName', labels(i));
    end
end

function plot_auc_single_curve_with_std(ax, x, y, s, color, lw, ms, alphaFill, displayName)
% Single curve version for Tile1 (Exp0-like): shadow + mean line + filled markers.
    x = x(:).';
    y = y(:).';
    s = s(:).';
    if any(isnan(y)) || any(isnan(s)), return; end

    ylo = y - s;
    yhi = y + s;

    xx = [x, fliplr(x)];
    yy = [yhi, fliplr(ylo)];
    fill(ax, xx, yy, color, ...
        'FaceAlpha', alphaFill, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    plot(ax, x, y, '-o', ...
        'LineWidth', lw, ...
        'MarkerSize', ms, ...
        'Color', color, ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color, ...
        'DisplayName', displayName);
end

function p = override_input_G2C1_G2C3(p, G2C1_new_nM, G2C3_new_nM)
% Override input-node genelet totals (assume G8C1/G8C3 correspond to G2C1/G2C3).
% Only updates p.aux.G2C1_tot, p.aux.G2C3_tot and rebuilds p.G_tot.
% (Other totals and dA2/dB2 are kept as fitted values.)

    p.aux.G2C1_tot = G2C1_new_nM;
    p.aux.G2C3_tot = G2C3_new_nM;

    % rebuild vectors (keep others as-is)
    G2C1 = p.aux.G2C1_tot;
    G2C3 = p.aux.G2C3_tot;
    G3R1 = p.aux.G3R1_tot;
    G1R3 = p.aux.G1R3_tot;
    G1S1 = p.aux.G1S1_tot;

    p.G_tot = [G2C1 G2C3 G3R1 G1R3 G1S1]' * 1e-9;
end

function p = apply_alpha_kp_scale_local(p, alpha)
% Apply a uniform scaling factor alpha to production-related rates.
% IMPORTANT:
% - This function must be used ONLY in the Exp0.5 simulation context.
% - It modifies only the local copy "p" passed in (caller must ensure no side effects).
    if ~isfinite(alpha) || alpha <= 0
        error('alpha must be positive and finite. Got alpha=%g', alpha);
    end
    p.k_pr = p.k_pr(:) * alpha;
    p.k_pc = p.k_pc(:) * alpha;
    p.k_pi = p.k_pi(:) * alpha;
end
