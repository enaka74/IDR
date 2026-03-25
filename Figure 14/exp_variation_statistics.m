%% ============================================================
%  debug_IFFL_ODE_reproducibility_with_plots.m
%
%  PURPOSE (minimal / no unnecessary changes):
%   - Extract ONLY the necessary parts that affect ODE simulation + AUC.
%   - Add plotting so you can visually inspect:
%       (A) fixed-p repeatability (AUC run1 vs run2)
%       (B) CV=0 Monte Carlo equivalence (AUC_mc distribution vs AUC_nom)
%       (C) optional time-series overlay (run1 vs run2) for one Din
%
%  REQUIREMENT:
%   - Run AFTER run_iffl_fit_include_autoinhibition() so that
%     IFFL_fit_results exists in the base workspace.
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

% Model size (should match fitting code)
n = 3;
m = 5;

%% ======== 1. Settings for AUC ===========================================
Channel = 'ch3';  % 'ch1' | 'ch2' | 'ch3' or 1|2|3
[chIdx, denomIdx, chLabel] = map_channel(Channel);

Tfinal_s     = 4*3600;
AUC_tStart_s = 0;
AUC_tEnd_s   = [];     % [] means use Tfinal_s

Din_list_min = [10 20 30 40 60 240];  % edit freely

idx_conn = Din_list_min <= 60;
idx_iso  = Din_list_min >= 240;

fprintf('\n=== Debug: ODE reproducibility + CV=0 equivalence (%s) ===\n', chLabel);

%% ======== 2. (A) Fixed-p repeatability check ============================
AUC_run1 = nan(size(Din_list_min));
AUC_run2 = nan(size(Din_list_min));
Nt1      = nan(size(Din_list_min));
Nt2      = nan(size(Din_list_min));

% store an example time-series for optional overlay plot
Din_overlay_min = 40;        % pick one
t_overlay_1 = []; y_overlay_1 = [];
t_overlay_2 = []; y_overlay_2 = [];

fprintf('\n--- (A) Fixed-p repeatability check ---\n');
for j = 1:numel(Din_list_min)
    Din_j = Din_list_min(j);

    p0 = p_fit;
    p0.T = Din_j * 60;

    [t1, Y1, ok1] = simulate_IFFL(p0, n, m);
    [t2, Y2, ok2] = simulate_IFFL(p0, n, m);

    if ~(ok1 && ok2)
        fprintf('Din=%d: solver failed (ok1=%d, ok2=%d)\n', Din_j, ok1, ok2);
        continue;
    end

    frac1 = Y1(:,chIdx) ./ p0.G_tot(denomIdx);
    frac2 = Y2(:,chIdx) ./ p0.G_tot(denomIdx);

    [t1_h, f1] = clip_for_auc(t1, frac1, AUC_tStart_s, AUC_tEnd_s, Tfinal_s);
    [t2_h, f2] = clip_for_auc(t2, frac2, AUC_tStart_s, AUC_tEnd_s, Tfinal_s);

    AUC1 = calculate_AUC(t1_h, f1);
    AUC2 = calculate_AUC(t2_h, f2);

    AUC_run1(j) = AUC1;
    AUC_run2(j) = AUC2;
    Nt1(j) = numel(t1);
    Nt2(j) = numel(t2);

    fprintf('Din=%3d min: AUC1=%.12g, AUC2=%.12g, dAUC=%.3g | Nt1=%d, Nt2=%d\n', ...
        Din_j, AUC1, AUC2, (AUC2-AUC1), numel(t1), numel(t2));

    if Din_j == Din_overlay_min
        t_overlay_1 = t1(:)/3600;
        y_overlay_1 = frac1(:);
        t_overlay_2 = t2(:)/3600;
        y_overlay_2 = frac2(:);
    end
end

%% ======== 3. (B) Monte Carlo =====================
Nmc = 500;     % keep small for debug; raise later
rng(1);       % MC reproducible

CV_kp = 0.1;
CV_kd = 0.1;

kp_fields = {'k_pc','k_pr','k_pi'};
kd_fields = {'k_dc','k_dr'};

AUC_nom_vec   = nan(size(Din_list_min));
AUC_mc_mean   = nan(size(Din_list_min));
AUC_mc_std    = nan(size(Din_list_min));
AUC_mc_maxabs = nan(size(Din_list_min));
maxScaleErr   = nan(size(Din_list_min));

% keep MC samples for plotting (Din x Nmc)
AUC_mc_samples = nan(numel(Din_list_min), Nmc);

fprintf('\n--- (B) CV=0 Monte Carlo equivalence check ---\n');
for j = 1:numel(Din_list_min)
    Din_j = Din_list_min(j);

    p0 = p_fit; p0.T = Din_j * 60;

    AUC_nom = compute_auc_single(p0, n, m, chIdx, denomIdx, AUC_tStart_s, AUC_tEnd_s, Tfinal_s);
    AUC_nom_vec(j) = AUC_nom;

    max_abs_scale_err = 0;

    for k = 1:Nmc
        p = p0;

        [s_kp, s_kd] = draw_lognormal_multipliers(CV_kp, CV_kd);
        max_abs_scale_err = max([max_abs_scale_err, abs(s_kp-1), abs(s_kd-1)]);

        p = scale_fields(p, kp_fields, s_kp);
        p = scale_fields(p, kd_fields, s_kd);

        AUC_mc_samples(j,k) = compute_auc_single(p, n, m, chIdx, denomIdx, AUC_tStart_s, AUC_tEnd_s, Tfinal_s);
    end

    d = AUC_mc_samples(j,:) - AUC_nom;
    AUC_mc_mean(j)   = mean(AUC_mc_samples(j,:), 'omitnan');
    AUC_mc_std(j)    = std(AUC_mc_samples(j,:), 'omitnan');
    AUC_mc_maxabs(j) = max(abs(d), [], 'omitnan');
    maxScaleErr(j)   = max_abs_scale_err;

    fprintf('Din=%3d min: AUC_nom=%.12g | MC mean=%.12g (std=%.3g) | max|d|=%.3g | max|scale-1|=%.3g\n', ...
        Din_j, AUC_nom, AUC_mc_mean(j), AUC_mc_std(j), AUC_mc_maxabs(j), maxScaleErr(j));
end

%% ======== 4. PLOTS (NEW) ================================================

% ======== x-axis "break" mapping for Din (DISPLAY ONLY) ==================
% We compress the gap (60 -> 240) by mapping Din=240 to x=70.
% This keeps earlier Din values unchanged.
[xDin, xTicks, xTickLabels, xBreakMid] = map_Din_with_break(Din_list_min, 60, 240, 70);


% -------------------------
% Fig A: AUC run1 vs run2
% -------------------------
figA = figure('Name','(A) Fixed-p repeatability: AUC run1 vs run2');
ax = axes(figA); hold(ax,'on'); box(ax,'on'); grid(ax,'off');

% 1:1 line
mn = min([AUC_run1(:); AUC_run2(:)], [], 'omitnan');
mx = max([AUC_run1(:); AUC_run2(:)], [], 'omitnan');
pad = 0.05*(mx-mn + eps);
plot(ax, [mn-pad, mx+pad], [mn-pad, mx+pad], '--', 'LineWidth', 2, 'Color', [0 0 0], 'HandleVisibility','off');

% scatter (with labels)
plot(ax, AUC_run1, AUC_run2, 'o', 'LineWidth', 2, 'MarkerSize', 9, 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'DisplayName','AUC pairs');
for j = 1:numel(Din_list_min)
    if isfinite(AUC_run1(j)) && isfinite(AUC_run2(j))
        text(ax, AUC_run1(j), AUC_run2(j), sprintf('  %d', Din_list_min(j)), ...
            'FontSize', 12, 'Color', [0 0 0], 'Interpreter','none');
    end
end

xlabel(ax,'AUC (run 1)');
ylabel(ax,'AUC (run 2)');
title(ax, sprintf('Fixed-p repeatability (%s)', chLabel), 'Interpreter','none');
ax.TickDir = 'in';
ax.LineWidth = 2;
legend(ax,'Location','best','Box','off');

% -------------------------
% Fig A2: dAUC vs Din + solver step counts
% -------------------------
figA2 = figure('Name','(A2) Fixed-p repeatability: dAUC and solver steps');
tl = tiledlayout(figA2, 2, 1, 'TileSpacing','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'off');
dAUC = AUC_run2 - AUC_run1;
plot(ax1, xDin, dAUC, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0]);
set(ax1, 'XTick', xTicks, 'XTickLabel', xTickLabels);
add_xbreak_marker(ax1, xBreakMid);

xlabel(ax1,'D_{in} [min]');
ylabel(ax1,'AUC(run2) - AUC(run1)');
title(ax1,'AUC difference across identical repeats','Interpreter','none');
ax1.TickDir = 'in'; ax1.LineWidth = 2;

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'off');
plot(ax2, xDin, Nt1, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0 0], ...
    'MarkerFaceColor', [0 0 0], 'DisplayName','Nt run1');
plot(ax2, xDin, Nt2, '--o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0 0], ...
    'MarkerFaceColor', 'none', 'DisplayName','Nt run2');
set(ax2, 'XTick', xTicks, 'XTickLabel', xTickLabels);
add_xbreak_marker(ax2, xBreakMid);

xlabel(ax2,'D_{in} [min]');
ylabel(ax2,'# solver time points');
title(ax2,'ode15s adaptive step counts','Interpreter','none');
ax2.TickDir = 'in'; ax2.LineWidth = 2;
legend(ax2,'Location','best','Box','off');

title(tl, 'Fixed-p repeatability diagnostics', 'Interpreter','none');

% -------------------------
% Fig B: MC distribution vs nominal
% -------------------------
figB = figure('Name','(B) CV=0 Monte Carlo: distribution vs nominal');
ax = axes(figB); hold(ax,'on'); box(ax,'on'); grid(ax,'off');

% jittered scatter per Din
for j = 1:numel(Din_list_min)
    x0 = Din_list_min(j);
    y  = AUC_mc_samples(j,:);
    xj = x0 + 0.35*(rand(size(y))-0.5); % small jitter
    %plot(ax, xj, y, 'o', 'MarkerSize', 6, 'LineWidth', 1.5, 'Color', [0 0 0], ...
    %    'MarkerFaceColor', 'none', 'HandleVisibility','off');
end

% nominal line
plot(ax, xDin(idx_conn), AUC_nom_vec(idx_conn), '-o', ...
    'LineWidth', 3, 'MarkerSize', 8, ...
    'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'DisplayName','Nominal');

plot(ax, xDin(idx_iso), AUC_nom_vec(idx_iso), 'o', ...
    'LineWidth', 3, 'MarkerSize', 8, ...
    'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'HandleVisibility','off');

% MC mean ± std as error bars
errorbar(ax, xDin(idx_conn), AUC_mc_mean(idx_conn), AUC_mc_std(idx_conn), 'o', ...
    'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0 0], ...
    'MarkerFaceColor', 'none', 'DisplayName','MC mean ± std');

errorbar(ax, xDin(idx_iso), AUC_mc_mean(idx_iso), AUC_mc_std(idx_iso), 'o', ...
    'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0 0], ...
    'MarkerFaceColor', 'none', 'HandleVisibility','off');

xlabel(ax,'D_{in} [min]');
ylabel(ax, sprintf('AUC (%s)', chLabel));
title(ax, sprintf('CV=0 MC should collapse to nominal (Nmc=%d)', Nmc), 'Interpreter','none');
ax.TickDir = 'in';
ax.LineWidth = 2;
xlim([0,70]);
legend(ax,'Location','best','Box','off');
set(ax, 'XTick', xTicks, 'XTickLabel', xTickLabels);
add_xbreak_marker(ax, xBreakMid);
xlim([min(xDin)-5, max(xDin)+5]);  % or keep [0,70] if you prefer fixed


%% ==================== Local helper functions (minimal) ===================

function AUCval = compute_auc_single(p, n, m, chIdx, denomIdx, tStart_s, tEnd_s, tfinal_s)
    [tSim, YSim, ok] = simulate_IFFL(p, n, m);
    if ~ok || isempty(tSim) || isempty(YSim)
        AUCval = NaN; return;
    end
    frac = YSim(:,chIdx) ./ p.G_tot(denomIdx);
    [t_h, f] = clip_for_auc(tSim, frac, tStart_s, tEnd_s, tfinal_s);
    if isempty(t_h)
        AUCval = NaN; return;
    end
    AUCval = calculate_AUC(t_h, f);
end

function [s_kp, s_kd] = draw_lognormal_multipliers(CV_kp, CV_kd)
% Lognormal multipliers with mean 1 and specified CV.
% sigma = sqrt(log(1+CV^2)), multiplier = exp(sigma * N(0,1)).
    sigma_kp = sqrt(log(1 + CV_kp^2));
    sigma_kd = sqrt(log(1 + CV_kd^2));
    s_kp = exp(sigma_kp * randn());
    s_kd = exp(sigma_kd * randn());
end

function p = scale_fields(p, fieldList, s)
% Multiply each existing field in fieldList by scalar s.
    for ii = 1:numel(fieldList)
        f = fieldList{ii};
        if isfield(p,f) && ~isempty(p.(f))
            p.(f) = p.(f) * s;
        end
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

%% ==================== Model functions (COPIED AS-IS) =====================
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

function [xDin, xTicks, xTickLabels, xBreakMid] = map_Din_with_break(Din_list_min, Din_left, Din_right, x_right_mapped)
% Map Din values to a display x-axis that "breaks" between Din_left and Din_right.
% - Din <= Din_left: x = Din (unchanged)
% - Din >= Din_right: x = x_right_mapped (for Din_right), and shift any larger Din accordingly (if ever used)

    Din_list_min = Din_list_min(:).';
    xDin = Din_list_min;

    % map Din_right -> x_right_mapped
    shift = x_right_mapped - Din_right;
    xDin(Din_list_min >= Din_right) = Din_list_min(Din_list_min >= Din_right) + shift;

    xTicks = xDin;
    xTickLabels = string(Din_list_min);

    % where to draw the "//" marker (midpoint between mapped Din_left and mapped Din_right)
    xLeftMapped  = Din_left;                 % unchanged
    xRightMapped = Din_right + shift;        % equals x_right_mapped
    xBreakMid = 0.5*(xLeftMapped + xRightMapped);
end

function add_xbreak_marker(ax, xmid)
% Draw a small "//" marker on the x-axis area to indicate a broken axis.
% Uses data coordinates (so it scales nicely).

    yl = ax.YLim;
    yr = yl(2) - yl(1);
    if yr <= 0, return; end

    % position a bit above bottom
    y0 = yl(1) + 0.06*yr;

    % marker size
    dx = 0.6;          % adjust if you want thicker/thinner
    dy = 0.03*yr;

    hold(ax, 'on');
    plot(ax, [xmid-dx, xmid-dx/3], [y0-dy, y0+dy], 'k-', 'LineWidth', 2, 'HandleVisibility','off');
    plot(ax, [xmid+dx/3, xmid+dx], [y0-dy, y0+dy], 'k-', 'LineWidth', 2, 'HandleVisibility','off');
end
