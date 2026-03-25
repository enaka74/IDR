function run_iffl_fit()
% RUN_IFFL_FIT — Fit reaction rates (time-series) + optional embedded AUC penalties (Exp1/2/3).
%
% Update (this revision):
% - AUC fraction_ON is FIXED for all Exp: Y(:,5)/G_tot(5) (G1S1)
% - AUC integral uses calculate_AUC() (Right Riemann sum, hour units)
% - Stage1/2/3 GA bounds are clipped by ABSOLUTE bounds defined from Stage0 baseline:
%     absolute bounds = [lb, ub] returned by make_rate_bounds(spreadPow) for relevant parameters.
%   Only tuned parameters are constrained by absolute bounds (others never change after Stage0).

% ===================== USER INPUT SECTION =====================
% --- Common time-series data (existing) ---
T_list_min   = [];
Weights      = [1 1 1];          % [w_ch1, w_ch2, w_ch3] for time-series SSE
spreadPow    = 1;                % +/- spread in log10-space around baselines
filePattern  = 'Din_%dmin_means.csv';
timeCol      = 'time_h';
ch1Col       = 'ON_mean_ch1';    % Ch1 → G_dA(1)/G2C1_tot
ch2Col       = 'ON_mean_ch2';    % Ch2 → G_dA(3)/G3R1_tot
ch3Col       = 'ON_mean_ch3';    % Ch3 → G_dA(5)/G1S1_tot
clipZero     = false;
rngSeed      = 1;

% Fixed totals (nM) and baseline k_p/k_d for initialization.
% Vector layout: [G2C1 G2C3 G3R1 G1R3 G1S1, dA2 dA3 dA1, dB2 dB3 dB1, k_p k_d]
fixed_opt_param = [25 10 5 1 25, 150 50 100, 0 125 37.5, 0.10 0.0015];

% --- Circuit size (fixed) ---
n = 3;   % regulator components
m = 5;   % genelet nodes

% --- Choose optimizer for baseline pipeline (Stage0/0.5 + default) ---
method = 'ga';   % keep: used for building GA options

% --- Choose optimizer for Exp stages (Exp1/2/3) ---
% 'ga' or 'fmincon'
% expOptMethod = 'ga';

% If you want to switch Exp by Exp (expOptMethod is ignored）
expOptMethodByID = containers.Map({1,2,3},{'ga','ga','ga'});


% --- GA options ---
gaUseParallel        = true;
gaPopSize            = 100;
gaEliteCount         = [];
gaCrossoverFraction  = 0.8;
gaMutationScale      = 0.5;
gaMutationShrink     = 0.95;
gaMaxGenerations       = 30;
gaMaxStallGenerations  = 10;
gaFunctionTolerance    = 1e-6;
gaFitnessLimit         = [];
gaMaxTime              = 18000;

% --- AUC integration options ---
SelectedExperiments = [1 2 3 4];   % [] or any subset of [1 2 3]

% AUC weights (independent from time-series channel weights)
wAUC = struct();
wAUC.exp1 = 1.5;
wAUC.exp2 = 1.0;
wAUC.exp3 = 1.0;
wAUC.exp4 = 1.0;

alpha_lb = 0.9;   % Lower bound for alpha (linear scale)
alpha_ub = 1.1;   % Upper bound for alpha (linear scale)

% Embedded AUC Din (minutes): columns correspond to these in AUC matrices
AUC_Din_vals_min = [20 40 60];

% --- Per-Din weights for time-series fitting (condition weights) ---
Din_vals_min     = [10 20 30 40 60 240];
Din_weights      = [2  1  1  1  1  2];   % example

% Toehold levels (nt) used in exp2/exp3 (rows of AUC.exp2.Y / AUC.exp3.Y)
toeholdLengths = [8 6 4];

% Stage 0 includes AUC or not (baseline time-series first is usually safer)
Stage0IncludeAUC = false;

% --- Run control (resume / stop) ---
% StartFromStage:
%   "0"    : run from Stage0 (default)
%   "0.5a" : resume from Stage0.5a (requires IFFL_fit_results in base)
%   "0.5b" : resume from Stage0.5b (requires IFFL_fit_results in base)
%   "1"    : resume from Stage1 (Exp1)
%   "2"    : resume from Stage2 (Exp2)
%   "3"    : resume from Stage3 (Exp3)
StartFromStage = "0";

% StopAfterStage:
%   "none": run through all selected stages
%   "0"   : stop after Stage0
%   "0.5a": stop after Stage0.5a
%   "0.5b": stop after Stage0.5b (common checkpoint)
%   "1"   : stop after Stage1
%   "2"   : stop after Stage2
%   "3"   : stop after Stage3
StopAfterStage = "none";

% ===============================================================

if ~isempty(rngSeed), rng(rngSeed); end

dirpath = fileparts(mfilename('fullpath'));

% Load time-series data
dataTS = read_iffl_dir(dirpath, 'Pattern',filePattern, 'T_list_min',T_list_min, ...
                               'TimeCol',timeCol, 'Ch1Col',ch1Col, 'Ch2Col',ch2Col, 'Ch3Col',ch3Col, ...
                               'ClipZero',clipZero);
dataTS.DinWeightMap = containers.Map(num2cell(Din_vals_min), num2cell(Din_weights));

% Build embedded AUC struct (provided function; embedded below)
AUC = build_embedded_AUC_data(AUC_Din_vals_min);

% Bounds for baseline time-series rates only (log10 space)
[p0,lb,ub,rateNames] = make_rate_bounds(spreadPow);

% ===== MINIMAL CHANGE: widen ONLY k_gar search range (log10 space) =====
%{
kgar_extra_decades = 0.5;    % additional decade
lb(1) = lb(1) - kgar_extra_decades;
ub(1) = ub(1) + kgar_extra_decades;
%}

% Base mapper (time-series baseline): returns p with per-node vectors
base_mapper = make_rates_only_param_mapper(fixed_opt_param, n, m);

% Absolute bounds (from Stage0 baseline) for relevant parameters
absBounds = build_abs_bounds_struct(lb, ub);

% Build base opts
opts = struct;
opts.n = n; opts.m = m;
opts.param_mapper = base_mapper;

opts.dataTS       = dataTS;
opts.DinWeightMap = containers.Map(num2cell(Din_vals_min), num2cell(Din_weights));
opts.Weights      = Weights;
opts.Ch1Index     = 1;

opts.AUC            = AUC;
opts.wAUC           = wAUC;
opts.toeholdLengths = toeholdLengths;

% absolute bounds container for stage builders
opts.absBounds = absBounds;

% Evaluate initial guess (TS-only by default here)
f_init = objective_total(p0, opts, []);
fprintf('\nInitial TS objective (no AUC): %.6g\n', f_init);
for i = 1:numel(rateNames)
    fprintf('init %-12s = %.6g\n', rateNames{i}, 10.^p0(i));
end

% Plot initial overlays (time-series only)
p_init = opts.param_mapper(p0);
figure; tiledlayout(numel(dataTS.conds), 1, 'TileSpacing', 'compact');
for k = 1:numel(dataTS.conds)
    ck = dataTS.conds(k);
    p_init.T = ck.T_sec;
    [tSim0, YSim0] = simulate_IFFL(p_init, n, m);
    frac1_0 = YSim0(:,1) / p_init.G_tot(1);
    frac3_0 = YSim0(:,3) / p_init.G_tot(3);
    frac5_0 = YSim0(:,5) / p_init.G_tot(5);

    nexttile; hold on
    plot(ck.time_s/3600, ck.y_ch1, 'o', 'DisplayName', 'ch1 data');
    plot(ck.time_s/3600, ck.y_ch2, 'o', 'DisplayName', 'ch2 data');
    plot(ck.time_s/3600, ck.y_ch3, 'o', 'DisplayName', 'ch3 data');
    plot(tSim0/3600, frac1_0, '-', 'DisplayName', 'sim ch1 (G2C1)');
    plot(tSim0/3600, frac3_0, '-', 'DisplayName', 'sim ch2 (G3R1)');
    plot(tSim0/3600, frac5_0, '-', 'DisplayName', 'sim ch3 (G1S1)');
    xlabel('time [h]'); ylabel('fraction ON'); legend('Location','best'); grid on
    title(sprintf('Initial guess — T = %d min', round(ck.T_sec/60)));
end

% GA options
gaOpts = [];
if strcmpi(method,'ga')
    if gaUseParallel
        try
            if isempty(gcp('nocreate')), parpool('local'); end
        catch
            warning('Parallel pool unavailable. Falling back to serial GA.');
            gaUseParallel = false;
        end
    end
    gaOpts = optimoptions('ga','Display','iter', ...
                              'UseParallel', gaUseParallel, ...
                              'CrossoverFraction', gaCrossoverFraction, ...
                              'MutationFcn', { @mutationgaussian, gaMutationScale, gaMutationShrink }, ...
                              'MaxGenerations', gaMaxGenerations, ...
                              'MaxStallGenerations', gaMaxStallGenerations, ...
                              'FunctionTolerance', gaFunctionTolerance, ...
                              'MaxTime', gaMaxTime, ...
                              'HybridFcn', []);
    if ~isempty(gaPopSize),    gaOpts = optimoptions(gaOpts,'PopulationSize',gaPopSize); end
    if ~isempty(gaEliteCount), gaOpts = optimoptions(gaOpts,'EliteCount',gaEliteCount); end
    if ~isempty(gaFitnessLimit)
        gaOpts = optimoptions(gaOpts,'FitnessLimit',gaFitnessLimit);
    end
end

% ===================== Resume logic =====================
didResume = false;

if ~strcmpi(StartFromStage, "0")
    try
        R = evalin('base','IFFL_fit_results');
        if ~isfield(R,'stages') || isempty(R.stages)
            error('IFFL_fit_results.stages not found or empty.');
        end

        results = R;  % carry over history

        % --- choose p_prev depending on StartFromStage ---
        S = string(StartFromStage);

        if any(strcmpi(S, ["0.5a","0.5b"]))
            idx0 = find_stage_by_prefix(results.stages, 'Stage0: time-series baseline');
            if isempty(idx0), idx0 = find_stage_by_prefix(results.stages, 'Stage0'); end
            if isempty(idx0)
                error('Could not find Stage0 in IFFL_fit_results.stages. Run Stage0 once (StopAfterStage="0") to create it.');
            end
            p_prev = results.stages(idx0).p_best;
            fprintf('\n[RESUME] StartFromStage=%s => using Stage0 p_best (stage #%d) as p_prev.\n', S, idx0);

        elseif strcmpi(S, "1")
            idx05b = find_stage_by_prefix(results.stages, 'Stage0.5b: TS-only local refine');
            idx05a = find_stage_by_prefix(results.stages, 'Stage0.5a: TS-only');
            idx0   = find_stage_by_prefix(results.stages, 'Stage0: time-series baseline');
            if isempty(idx0), idx0 = find_stage_by_prefix(results.stages, 'Stage0'); end

            if ~isempty(idx05b)
                p_prev = results.stages(idx05b).p_best;
                fprintf('\n[RESUME] StartFromStage=1 => using Stage0.5b p_best (stage #%d) as p_prev.\n', idx05b);
            elseif ~isempty(idx05a)
                p_prev = results.stages(idx05a).p_best;
                fprintf('\n[RESUME] StartFromStage=1 => using Stage0.5a p_best (stage #%d) as p_prev.\n', idx05a);
            elseif ~isempty(idx0)
                p_prev = results.stages(idx0).p_best;
                fprintf('\n[RESUME] StartFromStage=1 => using Stage0 p_best (stage #%d) as p_prev.\n', idx0);
            else
                error('Could not find Stage0/0.5a/0.5b in history. Run from Stage0 once.');
            end

        elseif strcmpi(S, "2")
            idx1 = find_stage_by_prefix(results.stages, 'Stage1:');
            if isempty(idx1)
                error('Could not find Stage1 in history. Run Stage1 once (StopAfterStage="1") before resuming from Stage2.');
            end
            p_prev = results.stages(idx1).p_best;
            fprintf('\n[RESUME] StartFromStage=2 => using Stage1 p_best (stage #%d) as p_prev.\n', idx1);

        elseif strcmpi(S, "3")
            idx2 = find_stage_by_prefix(results.stages, 'Stage2:');
            if isempty(idx2)
                error('Could not find Stage2 in history. Run Stage2 once (StopAfterStage="2") before resuming from Stage3.');
            end
            p_prev = results.stages(idx2).p_best;
            fprintf('\n[RESUME] StartFromStage=3 => using Stage2 p_best (stage #%d) as p_prev.\n', idx2);

        else
            error('Unknown StartFromStage="%s". Use "0","0.5a","0.5b","1","2","3".', S);
        end

        stageIdx = numel(results.stages) + 1;
        didResume = true;

    catch ME
        error('Resume requested (StartFromStage=%s) but could not load IFFL_fit_results: %s', StartFromStage, ME.message);
    end
end


% ===================== Optimization Pipeline =====================
% IMPORTANT: do NOT overwrite results if resuming.
if ~didResume
    results = struct();
    results.stages = struct([]);
end

if strcmpi(StartFromStage, "0")
    % ---- Stage 0 (baseline) ----
    stage = struct();
    stage.name = 'Stage0: time-series baseline';
    stage.mapper = base_mapper;
    stage.p_prev = [];
    stage.selectedExp = [];

    if Stage0IncludeAUC
        stage.selectedExp = SelectedExperiments(:)';
    end

    fprintf('\n=== %s ===\n', stage.name);
    obj0 = @(x) objective_total(x, opts, stage.selectedExp);
    [bestP0,bestF0] = run_ga(obj0, p0, lb, ub, gaOpts);
    p_best0 = stage.mapper(bestP0);

    results.stages(1).name = stage.name;
    results.stages(1).bestP_log10 = bestP0;
    results.stages(1).bestF = bestF0;
    results.stages(1).p_best = p_best0;

    fprintf('Best (Stage0) objective: %.6g\n', bestF0);

    % checkpoint save
    p_prev = p_best0;
    results.final = p_prev;
    assignin('base','IFFL_fit_results',results);

    if strcmpi(StopAfterStage,"0")
        fprintf('\n[STOP] Stopped after Stage0. Inspect IFFL_fit_results and plots, then rerun with StartFromStage="1".\n');
        return;
    end

    stageIdx = 2;
end

% ---- Stage 0.5a: TS + Exp1(AUC), tune node multipliers for k_p/k_d + (k_ar3,k_gar3) ----
run05a = any(strcmpi(StartFromStage, ["0","0.5a"]));
if run05a
    stage = struct();
    stage.name = 'Stage0.5a: TS + Exp1(AUC) (tune kp/kd node multipliers + k_ar3,k_gar3)';
    stage.selectedExp = [1];  % <-- IMPORTANT: include Exp1 AUC in 0.5a
    stage = build_stage05a_kp_kd_plus_ar3_gar3(stage, p_prev, opts.absBounds, spreadPow, opts.n, opts.m);

    fprintf('\n=== %s ===\n', stage.name);

    % objective takes z (log10) -> apply -> p -> evaluate TS + AUCexp1
    obj05 = @(z) objective_total_fromP(stage.apply(z, stage.p_prev), opts, stage.selectedExp);

    % GA on z
    [bestZ05,bestF05] = run_ga(obj05, stage.z0, stage.zlb, stage.zub, gaOpts);
    p_prev = stage.apply(bestZ05, stage.p_prev);

    results.stages(stageIdx).name = stage.name;
    results.stages(stageIdx).bestZ_log10 = bestZ05;
    results.stages(stageIdx).bestF = bestF05;
    results.stages(stageIdx).p_best = p_prev;
    stageIdx = stageIdx + 1;

    fprintf('Best (%s) objective: %.6g\n', stage.name, bestF05);

    % checkpoint save
    results.final = p_prev;
    assignin('base','IFFL_fit_results',results);

    if strcmpi(StopAfterStage,"0.5a")
        fprintf('\n[STOP] Stopped after Stage0.5a. Rerun with StartFromStage="0.5b" or "2" to continue.\n');
        return;
    end
end

% ---- Stage 0.5b: local refine with fmincon (ALL params) ----
run05b = any(strcmpi(StartFromStage, ["0","0.5a","0.5b"]));
if run05b
    stage2 = struct();
    stage2.name = 'Stage0.5b: local refine (fmincon on ALL params; TS + Exp1(AUC))';
    stage2.selectedExp = [1];  % <-- CHANGED (was TS-only)
    stage2.p_prev = p_prev;

    fprintf('\n=== %s ===\n', stage2.name);

    multLo = 0.5;
    multHi = 2.0;
    stage2 = build_stage05b_allparams(stage2, opts.absBounds, multLo, multHi, opts.n, opts.m);

    obj05b = @(x) objective_total_fromP(stage2.apply(x, stage2.p_prev), opts, stage2.selectedExp);

    epsIn = 1e-10;
    x0_fmin = min(max(stage2.x0(:), stage2.xlb(:) + epsIn), stage2.xub(:) - epsIn);

    try
        [bestX05b,bestF05b] = run_fmincon(obj05b, x0_fmin, stage2.xlb, stage2.xub);

        p_prev = stage2.apply(bestX05b, stage2.p_prev);

        results.stages(stageIdx).name = stage2.name;
        results.stages(stageIdx).bestX_log10 = bestX05b;
        results.stages(stageIdx).bestF = bestF05b;
        results.stages(stageIdx).p_best = p_prev;
        stageIdx = stageIdx + 1;

        fprintf('Best (%s) objective: %.6g\n', stage2.name, bestF05b);

    catch ME
        warning('Stage0.5b fmincon failed or unavailable: %s. Skipping Stage0.5b.', ME.message);
    end

    results.final = p_prev;
    assignin('base','IFFL_fit_results',results);

    if strcmpi(StopAfterStage,"0.5b")
        fprintf('\n[STOP] Stopped after Stage0.5b. Rerun with StartFromStage="2" to continue.\n');
        return;
    end
end

% ---- Sequential fine-tuning: Exp2 -> Exp3 ----
startExpID = 2;
if strcmpi(StartFromStage,"2"), startExpID = 2; end
if strcmpi(StartFromStage,"3"), startExpID = 3; end

% Backward-compat: if someone sets StartFromStage="1", treat it as "2"
if strcmpi(StartFromStage,"1")
    warning('StartFromStage="1" is deprecated (Stage1 removed). Starting from Stage2 instead.');
    startExpID = 2;
end

for expID = [2 3]
    if expID < startExpID
        continue;
    end

    if ~any(SelectedExperiments==expID), continue; end

    switch expID
        case 2
            stage.name = 'Stage2: +Exp2 AUC (tune k_ar3,k_gar3 for 8/6/4 nt)';
            stage.selectedExp = [2];
            stage = build_stage_exp2(stage, p_prev, spreadPow, opts.absBounds);
        case 3
            stage.name = 'Stage3: +Exp3 AUC (tune k_bc3,k_gbc3 for 8/6/4 nt)';
            stage.selectedExp = [3];
            stage = build_stage_exp3(stage, p_prev, spreadPow, opts.absBounds);

        otherwise
            error('Unknown expID');
    end

    fprintf('\n=== %s ===\n', stage.name);

    objStage = @(z) objective_total_fromP(stage.apply(z, stage.p_prev), opts, stage.selectedExp);

    % --- decide optimizer for this Exp stage ---
    thisMethod = 'ga';  % default if map not provided
    
    if exist('expOptMethodByID','var') && isa(expOptMethodByID,'containers.Map') && isKey(expOptMethodByID, expID)
        thisMethod = expOptMethodByID(expID);
    end

    % Switch GA / fmincon at each stage
    [thisBestZ,thisBestF] = run_stage_optimizer(objStage, stage.z0, stage.zlb, stage.zub, gaOpts, thisMethod);
    bestZ = thisBestZ;  bestF = thisBestF;

    p_prev = stage.apply(bestZ, stage.p_prev);

    %%%%% DEBUG %%%%%%
    if expID == 3
        % --- Debug: check Exp3 ratios immediately after apply ---
        z_lin = 10.^bestZ(:);
        fprintf('[Stage3 bestZ lin] r_bc6=%.6g r_bc4=%.6g r_gbc6=%.6g r_gbc4=%.6g\n', ...
            z_lin(1), z_lin(2), z_lin(3), z_lin(4));
    
        k_bc8  = p_prev.k_bc(3);
        k_gbc8 = p_prev.k_gbc(3);
    
        fprintf('[Stage3 after apply ratios] r_bc6=%.6g r_bc4=%.6g r_gbc6=%.6g r_gbc4=%.6g\n', ...
            p_prev.tuned.exp3.k_bc3(2)/k_bc8, ...
            p_prev.tuned.exp3.k_bc3(3)/k_bc8, ...
            p_prev.tuned.exp3.k_gbc3(2)/k_gbc8, ...
            p_prev.tuned.exp3.k_gbc3(3)/k_gbc8);
    end


    results.stages(stageIdx).name = stage.name;
    results.stages(stageIdx).bestZ_log10 = bestZ;
    results.stages(stageIdx).bestF = bestF;
    results.stages(stageIdx).p_best = p_prev;
    stageIdx = stageIdx + 1;

    fprintf('Best (%s) objective: %.6g\n', stage.name, bestF);

    % checkpoint save after each Exp stage
    results.final = p_prev;
    assignin('base','IFFL_fit_results',results);

    % StopAfterStage for Exp1/2/3
    if strcmpi(StopAfterStage, string(expID))
        fprintf('\n[STOP] Stopped after Stage%d. Inspect IFFL_fit_results, then rerun to continue.\n', expID);
        return;
    end
end

% ---- Stage4: Exp4 AUC-only (fit alpha; do NOT modify p_prev) ----
if any(SelectedExperiments==4)

    fprintf('\n=== Stage4: Exp4 AUC-only (fit alpha; no side effects) ===\n');

    % Decision variable is z = log10(alpha)
    z0  = log10(1.0);
    zlb = log10(alpha_lb);
    zub = log10(alpha_ub);

    % AUC-only objective for Exp4 (no TS term, no modification to p_prev)
    objStage4 = @(z) opts.wAUC.exp4 * objective_auc_exp4_alpha_only(p_prev, 10.^z(1), opts, 5);

    % Prefer fmincon for 1D unless overridden
    thisMethod = "fmincon";
    if exist('expOptMethodByID','var') && isa(expOptMethodByID,'containers.Map') && isKey(expOptMethodByID, 4)
        thisMethod = lower(string(expOptMethodByID(4)));
    end

    [bestZ4, bestF4] = run_stage_optimizer(objStage4, z0, zlb, zub, gaOpts, thisMethod);

    alpha_best = 10.^bestZ4(1);

    % Save alpha separately (do NOT apply to p_prev)
    stage = struct();
    stage.name = 'Stage4: Exp4 AUC-only (fit alpha; do NOT apply to k_p)';
    stage.selectedExp = [4];

    results.stages(stageIdx).name = stage.name;
    results.stages(stageIdx).bestZ_log10 = bestZ4;
    results.stages(stageIdx).bestF = bestF4;
    results.stages(stageIdx).alpha_best = alpha_best;
    results.stages(stageIdx).p_best = p_prev;   % unchanged
    stageIdx = stageIdx + 1;

    % Also store in a dedicated field for convenience
    results.exp4 = struct();
    results.exp4.alpha_best = alpha_best;
    results.exp4.bestF = bestF4;

    fprintf('Best (Stage4) objective: %.6g\n', bestF4);
    fprintf('Best alpha (stored only) = %.6g\n', alpha_best);

    % Checkpoint save
    results.final = p_prev;  % unchanged final params
    assignin('base','IFFL_fit_results',results);

    if strcmpi(StopAfterStage,"4")
        fprintf('\n[STOP] Stopped after Stage4.\n');
        return;
    end
end


results.final = results.stages(end).p_best;

fprintf('\n==== FINAL PARAMS (after selected stages) ====\n');
print_key_params(results.final);

% Quick overlays final
p_plot = results.final;
figure; tiledlayout(numel(dataTS.conds),1,'TileSpacing','compact');
for k=1:numel(dataTS.conds)
    ck = dataTS.conds(k);
    p_plot.T = ck.T_sec;
    [tSim,YSim] = simulate_IFFL(p_plot, n, m);
    frac1 = YSim(:,1)/p_plot.G_tot(1);
    frac3 = YSim(:,3)/p_plot.G_tot(3);
    frac5 = YSim(:,5)/p_plot.G_tot(5);
    nexttile; hold on
    plot(ck.time_s/3600, ck.y_ch1,'o','DisplayName','ch1 data');
    plot(ck.time_s/3600, ck.y_ch2,'o','DisplayName','ch2 data');
    plot(ck.time_s/3600, ck.y_ch3,'o','DisplayName','ch3 data');
    plot(tSim/3600, frac1, '-', 'DisplayName','sim ch1 (G2C1)');
    plot(tSim/3600, frac3, '-', 'DisplayName','sim ch2 (G3R1)');
    plot(tSim/3600, frac5, '-', 'DisplayName','sim ch3 (G1S1)');
    xlabel('time [h]'); ylabel('fraction ON'); legend('Location','best'); grid on
    title(sprintf('FINAL — T = %d min', round(ck.T_sec/60)));
end

assignin('base','IFFL_fit_results',results);
end

%% ============== Data input (time-series, multi-T, 3ch) ===================
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

%% ============== Bounds & mappers (base rates) ============================
function [p0,lb,ub,names] = make_rate_bounds(spreadPow)
    if nargin<1, spreadPow = 0.5; end
    [r0,names] = baseline_rates();
    logr0 = log10(r0(:));
    p0 = logr0;
    lb = logr0 - spreadPow;
    ub = logr0 + spreadPow;
end

function [r0,names] = baseline_rates()
    p = default_params();
    r0 = [p.k_gar, p.k_ar, p.k_ir, p.k_ga, p.k_gb, ...
          p.k_gab, p.k_gbc, p.k_bc, p.k_p, p.k_d, ...
          p.k_aa_input, p.k_aa_other, p.k_gar_d, p.leak];
    names = {'k_gar','k_ar','k_ir','k_ga','k_gb','k_gab','k_gbc','k_bc','k_p','k_d', ...
             'k_aa_input','k_aa_other','k_gar_d','leak'};
end

function mapper = make_rates_only_param_mapper(fixed_opt_param, n, m)
    if nargin<1
        fixed_opt_param = [25 10 5 1 25 150 50 100, 0 125 37.5, 0.10 0.0015];
    end
    if nargin<2, n=3; end
    if nargin<3, m=5; end

    function p = inner(x_log10)
        p = apply_optimized_params(default_params(), fixed_opt_param);
        x = 10.^x_log10(:);

        auxScalar = struct();
        auxScalar.k_gar       = x(1);
        auxScalar.k_ar        = x(2);
        auxScalar.k_ir        = x(3);
        auxScalar.k_ga        = x(4);
        auxScalar.k_gb        = x(5);
        auxScalar.k_gab       = x(6);
        auxScalar.k_gbc       = x(7);
        auxScalar.k_bc        = x(8);
        auxScalar.k_p         = x(9);
        auxScalar.k_d         = x(10);
        auxScalar.k_aa_input  = x(11);
        auxScalar.k_aa_other  = x(12);
        auxScalar.k_gar_d     = x(13);
        auxScalar.leak        = x(14);

        p.leak        = auxScalar.leak;
        p.k_aa_input  = auxScalar.k_aa_input;
        p.k_aa_other  = auxScalar.k_aa_other;

        % per-component vectors (n)
        p.k_gar   = auxScalar.k_gar   * ones(n,1);
        p.k_ar    = auxScalar.k_ar    * ones(n,1);
        p.k_ir    = auxScalar.k_ir    * ones(n,1);
        p.k_ga    = auxScalar.k_ga    * ones(n,1);
        p.k_gb    = auxScalar.k_gb    * ones(n,1);
        p.k_gab   = auxScalar.k_gab   * ones(n,1);
        p.k_gbc   = auxScalar.k_gbc   * ones(n,1);
        p.k_bc    = auxScalar.k_bc    * ones(n,1);
        p.k_gar_d = auxScalar.k_gar_d * ones(n,1);
        p.k_dr    = auxScalar.k_d     * ones(n,1);
        p.k_dc    = auxScalar.k_d     * ones(n,1);

        % per-node vectors (m)
        p.k_pr = auxScalar.k_p * ones(m,1);
        p.k_pc = auxScalar.k_p * ones(m,1);
        p.k_pi = auxScalar.k_p * ones(m,1);

        p.auxScalar = auxScalar;

        % toehold-specific containers (initialized to baseline)
        p.tuned = struct();
        p.tuned.exp2 = struct('toehold', [8 6 4], ...
                              'k_ar3',  auxScalar.k_ar  * ones(1,3), ...
                              'k_gar3', auxScalar.k_gar * ones(1,3));
        p.tuned.exp3 = struct('toehold', [8 6 4], ...
                              'k_bc3',  auxScalar.k_bc  * ones(1,3), ...
                              'k_gbc3', auxScalar.k_gbc * ones(1,3));
    end
    mapper = @inner;
end

function absB = build_abs_bounds_struct(lb, ub)
% baseline_rates order:
%  1 k_gar, 2 k_ar, 3 k_ir, 4 k_ga, 5 k_gb,
%  6 k_gab, 7 k_gbc, 8 k_bc, 9 k_p, 10 k_d,
%  11 k_aa_input, 12 k_aa_other, 13 k_gar_d, 14 leak

    absB = struct();

    absB.k_gar      = [lb(1);  ub(1)];
    absB.k_ar       = [lb(2);  ub(2)];
    absB.k_ir       = [lb(3);  ub(3)];
    absB.k_ga       = [lb(4);  ub(4)];
    absB.k_gb       = [lb(5);  ub(5)];
    absB.k_gab      = [lb(6);  ub(6)];
    absB.k_gbc      = [lb(7);  ub(7)];
    absB.k_bc       = [lb(8);  ub(8)];
    absB.k_p        = [lb(9);  ub(9)];
    absB.k_d        = [lb(10); ub(10)];
    absB.k_aa_input = [lb(11); ub(11)];
    absB.k_aa_other = [lb(12); ub(12)];
    absB.k_gar_d    = [lb(13); ub(13)];
    absB.leak       = [lb(14); ub(14)];
end

%% ============== Total objectives =========================================
function f = objective_total(x_log10, opts, activeExp)
    p = opts.param_mapper(x_log10);
    f = objective_total_fromP(p, opts, activeExp);
end

function f = objective_total_fromP(p, opts, activeExp)
    if nargin<3 || isempty(activeExp), activeExp = []; end

    fTS = objective_IFFL_multi_fromP(p, opts.dataTS, opts.Weights, opts.Ch1Index, opts.n, opts.m);
    f = fTS;

    % AUC uses fixed fraction_ON = Y(:,5)/G_tot(5) for ALL experiments
    stateIdxAUC = 5;

    if any(activeExp==1)
        f = f + opts.wAUC.exp1 * objective_auc_exp1(p, opts, stateIdxAUC);
    end
    if any(activeExp==2)
        f = f + opts.wAUC.exp2 * objective_auc_exp2(p, opts, stateIdxAUC);
    end
    if any(activeExp==3)
        f = f + opts.wAUC.exp3 * objective_auc_exp3(p, opts, stateIdxAUC);
    end
    
end

function f = objective_IFFL_multi_fromP(p, dataTS, Weights, ch1Idx, n, m)
    W = Weights(:).';
    if numel(W) < 3, W = [1, W]; end
    w1 = W(1); w2 = W(2); w3 = W(3);

    BIG = 1e12;
    err_all = [];

    for k = 1:numel(dataTS.conds)
        ck = dataTS.conds(k);
        
        % ===== NEW: per-Din condition weight =====
        wDin = 1;
        if isfield(dataTS, 'DinWeightMap') && isa(dataTS.DinWeightMap,'containers.Map')
            din_min = round(ck.T_sec/60);
            if isKey(dataTS.DinWeightMap, din_min)
                wDin = dataTS.DinWeightMap(din_min);
            end
        end
        % ========================================
        
        p.T = ck.T_sec;

        try
            [tSim,YSim,ok] = simulate_IFFL(p, n, m);
            if ~ok || any(~isfinite(YSim(:))) || any(~isfinite(tSim))
                f = BIG; return;
            end
        catch
            f = BIG; return;
        end

        frac1 = (YSim(:,ch1Idx) ./ p.G_tot(ch1Idx));
        frac2 = (YSim(:,3)      ./ p.G_tot(3));
        frac3 = (YSim(:,5)      ./ p.G_tot(5));

        if any(~isfinite(frac1)) || any(~isfinite(frac2)) || any(~isfinite(frac3))
            f = BIG; return;
        end

        [tSimU, ia] = unique(tSim, 'stable');
        f1U = frac1(ia); f2U = frac2(ia); f3U = frac3(ia);

        y1 = interp1(tSimU, f1U, ck.time_s, 'linear', 'extrap');
        y2 = interp1(tSimU, f2U, ck.time_s, 'linear', 'extrap');
        y3 = interp1(tSimU, f3U, ck.time_s, 'linear', 'extrap');
        if any(~isfinite(y1)) || any(~isfinite(y2)) || any(~isfinite(y3))
            f = BIG; return;
        end

        % Apply wDin to all channels for this condition
        ww1 = wDin * w1;  ww2 = wDin * w2;  ww3 = wDin * w3;

        if isfield(ck,'y_ch1') && ~isempty(ck.y_ch1)
            err_all = [err_all; ww1*(y1 - ck.y_ch1)]; %#ok<AGROW>
        end
        err_all = [err_all; ww2*(y2 - ck.y_ch2); ww3*(y3 - ck.y_ch3)]; %#ok<AGROW>
    end

    f = sum(err_all.^2);
end

%% ============== AUC objectives (matrix-based) ============================
function f = objective_auc_exp1(p, opts, stateIdx)
% Exp1: rows = G3R1 levels [5 10 25], cols = Din vals [20 40 60]
% Also increase dB3_tot by the same delta as G3R1_tot to keep (dB3_tot - G3R1_tot) constant.
% AUC computed by calculate_AUC(t_h, frac)

    BIG = 1e12;

    Din_min = opts.AUC.Din_vals(:).';

    % --- backward-compatible field name handling ---
    if isfield(opts.AUC.exp1,'G3R1_levels_nM')
        levels = opts.AUC.exp1.G3R1_levels_nM(:);
    elseif isfield(opts.AUC.exp1,'G2C3_levels_nM')
        % fallback if you haven't renamed the field yet
        levels = opts.AUC.exp1.G2C3_levels_nM(:);
    else
        error('Exp1: missing levels field (expected G3R1_levels_nM).');
    end

    Ymeas = opts.AUC.exp1.Y;

    if size(Ymeas,1) ~= numel(levels) || size(Ymeas,2) ~= numel(Din_min)
        error('Exp1 AUC matrix size mismatch.');
    end

    Ypred = nan(size(Ymeas));

    for i = 1:numel(levels)
        for j = 1:numel(Din_min)

            % --- HERE is the key fix: override G3R1_tot and co-adjust dB3_tot ---
            p_ij = override_total_G3R1_and_dB3(p, levels(i));

            p_ij.T = Din_min(j) * 60;

            try
                [tSim,YSim,ok] = simulate_IFFL(p_ij, opts.n, opts.m);
                if ~ok, f = BIG; return; end

                frac = YSim(:,stateIdx) ./ p_ij.G_tot(stateIdx);
                t_h = tSim / 3600;
                Ypred(i,j) = calculate_AUC(t_h, frac);
            catch
                f = BIG; return;
            end
        end
    end

    if any(~isfinite(Ypred(:))), f = BIG; return; end
    f = sum((Ypred(:) - Ymeas(:)).^2);
end


function f = objective_auc_exp2(p, opts, stateIdx)
% Exp2: rows = toehold [8 6 4], cols = Din [20 40 60]
% Apply p.tuned.exp2 (k_ar3,k_gar3) per row.
% AUC computed by calculate_AUC(t_h, frac)

    BIG = 1e12;

    Din_min = opts.AUC.Din_vals(:).';
    thRows  = opts.AUC.exp2.toehold_nt(:);
    Ymeas   = opts.AUC.exp2.Y;

    if size(Ymeas,1) ~= numel(thRows) || size(Ymeas,2) ~= numel(Din_min)
        error('Exp2 AUC matrix size mismatch.');
    end

    Ypred = nan(size(Ymeas));

    for i = 1:numel(thRows)
        idx = find(opts.toeholdLengths==thRows(i), 1);
        if isempty(idx), error('Exp2: unknown toehold %d', thRows(i)); end

        for j = 1:numel(Din_min)
            p_ij = p;
            p_ij.k_ar(3)  = p_ij.tuned.exp2.k_ar3(idx);
            p_ij.k_gar(3) = p_ij.tuned.exp2.k_gar3(idx);
            p_ij.T = Din_min(j) * 60;

            try
                [tSim,YSim,ok] = simulate_IFFL(p_ij, opts.n, opts.m);
                if ~ok, f = BIG; return; end

                frac = YSim(:,stateIdx) ./ p_ij.G_tot(stateIdx);
                t_h = tSim / 3600;
                Ypred(i,j) = calculate_AUC(t_h, frac);
            catch
                f = BIG; return;
            end
        end
    end

    if any(~isfinite(Ypred(:))), f = BIG; return; end
    f = sum((Ypred(:) - Ymeas(:)).^2);
end

function f = objective_auc_exp3(p, opts, stateIdx)
% Exp3: rows = toehold [8 6 4], cols = Din [20 40 60]
% Apply p.tuned.exp3 (k_bc3,k_gbc3) per row.
% AUC computed by calculate_AUC(t_h, frac)

    BIG = 1e12;

    Din_min = opts.AUC.Din_vals(:).';
    thRows  = opts.AUC.exp3.toehold_nt(:);
    Ymeas   = opts.AUC.exp3.Y;

    if size(Ymeas,1) ~= numel(thRows) || size(Ymeas,2) ~= numel(Din_min)
        error('Exp3 AUC matrix size mismatch.');
    end

    Ypred = nan(size(Ymeas));

    for i = 1:numel(thRows)
        idx = find(opts.toeholdLengths==thRows(i), 1);
        if isempty(idx), error('Exp3: unknown toehold %d', thRows(i)); end

        for j = 1:numel(Din_min)
            p_ij = p;
            p_ij.k_bc(3)  = p_ij.tuned.exp3.k_bc3(idx);
            p_ij.k_gbc(3) = p_ij.tuned.exp3.k_gbc3(idx);
            p_ij.T = Din_min(j) * 60;

            try
                [tSim,YSim,ok] = simulate_IFFL(p_ij, opts.n, opts.m);
                if ~ok, f = BIG; return; end

                frac = YSim(:,stateIdx) ./ p_ij.G_tot(stateIdx);
                t_h = tSim / 3600;
                Ypred(i,j) = calculate_AUC(t_h, frac);
            catch
                f = BIG; return;
            end
        end
    end

    if any(~isfinite(Ypred(:))), f = BIG; return; end
    f = sum((Ypred(:) - Ymeas(:)).^2);
end


function f = objective_auc_exp4_alpha_only(p_base, alpha, opts, stateIdx)
% Exp4 (AUC-only):
%   Fit a single scalar alpha that scales k_pr/k_pc/k_pi uniformly,
%   but ONLY within this Exp4 objective evaluation.
%   The base parameter struct p_base is NOT modified (no side effects).

    BIG = 1e12;

    levels  = opts.AUC.exp4.G2C1_G2C3_levels_nM;  % 3x2 (nM)
    Ymeas   = opts.AUC.exp4.Y;                    % 3x1
    Din_min = opts.AUC.exp4.Din_min;

    if size(levels,2) ~= 2
        error('Exp4: expected 3x2 levels array for (G2C1,G2C3).');
    end
    if size(levels,1) ~= numel(Ymeas)
        error('Exp4: levels rows and Y length mismatch.');
    end
    if any(~isfinite(Ymeas))
        error('Exp4: AUC.exp4.Y contains NaN/Inf. Please fill experimental AUC values.');
    end
    if ~isfinite(alpha) || alpha <= 0
        f = BIG; return;
    end

    Ypred = nan(size(Ymeas));

    for i = 1:size(levels,1)
        g2c1 = levels(i,1);
        g2c3 = levels(i,2);

        % Work on a local copy only (no persistent modification)
        p_i = p_base;

        % Override totals for this condition
        p_i = override_total_G2C1_and_G2C3(p_i, g2c1, g2c3);

        % Apply alpha ONLY here (temporary)
        p_i.k_pr = p_i.k_pr(:) * alpha;
        p_i.k_pc = p_i.k_pc(:) * alpha;
        p_i.k_pi = p_i.k_pi(:) * alpha;

        % Set pulse duration
        p_i.T = Din_min * 60;

        try
            [tSim, YSim, ok] = simulate_IFFL(p_i, opts.n, opts.m);
            if ~ok, f = BIG; return; end

            frac = YSim(:,stateIdx) ./ p_i.G_tot(stateIdx);
            t_h  = tSim / 3600;
            Ypred(i) = calculate_AUC(t_h, frac);
        catch
            f = BIG; return;
        end
    end

    if any(~isfinite(Ypred)), f = BIG; return; end
    f = sum((Ypred - Ymeas).^2);
end



%% ============== Stage builder (Stage0.5: node multipliers for k_p, k_d) ===
function stage = build_stage_kp_kd_node_mult(stage, p_prev, n, m)
% Stage0.5a:
%   Tune node-wise multipliers for k_p and k_d using TS-only objective.
%   Multipliers are constrained to [0.5, 2] in linear scale.
%
% Decision variables z (log10):
%   z = [log10(r_kp(1:m)); log10(r_kd(1:n))]

    stage.p_prev = p_prev;

    % Initial guess = all ones (or could use ratios from p_prev if you have them)
    r_kp0 = ones(m,1);
    r_kd0 = ones(n,1);

    stage.z0 = log10([r_kp0; r_kd0]);

    % Bounds: 0.5..2 in linear -> log10 bounds
    lo = log10(0.3);
    hi = log10(3.0);

    stage.zlb = lo * ones(numel(stage.z0),1);
    stage.zub = hi * ones(numel(stage.z0),1);

    stage.apply = @(z,p_in) apply_stage_kp_kd_node_mult(z,p_in,n,m);
end

function stage = build_stage05a_kp_kd_plus_ar3_gar3(stage, p_prev, absB, spreadPow, n, m)
% Stage0.5a (NEW):
%   Tune:
%     - node-wise multipliers for k_p (m vars) in [0.5,2]
%     - component-wise multipliers for k_d (n vars) in [0.5,2]
%     - absolute k_ar3, k_gar3 (2 vars), bounded by local +/-spreadPow clipped by abs bounds
%
% Decision variables z (log10):
%   z = [log10(r_kp(1:m));
%        log10(r_kd(1:n));
%        log10(k_ar3);
%        log10(k_gar3)]

    if nargin<4 || isempty(spreadPow), spreadPow = 1; end
    if nargin<5, n = 3; end
    if nargin<6, m = 5; end

    stage.p_prev = p_prev;

    % --- initial multipliers ---
    r_kp0 = ones(m,1);
    r_kd0 = ones(n,1);

    % --- initial k_ar3/k_gar3 (absolute) ---
    k_ar3_0  = p_prev.k_ar(3);
    k_gar3_0 = p_prev.k_gar(3);

    stage.z0 = log10([r_kp0; r_kd0; k_ar3_0; k_gar3_0]);

    % --- bounds for multipliers ---
    loMult = log10(0.3);
    hiMult = log10(3.0);

    zlb_mult = loMult * ones(m+n,1);
    zub_mult = hiMult * ones(m+n,1);

    % --- bounds for k_ar3/k_gar3: local +/- spreadPow clipped by absolute bounds ---
    z0_ar3  = log10(k_ar3_0);
    z0_gar3 = log10(k_gar3_0);

    zlb_ar3  = max(z0_ar3  - spreadPow, absB.k_ar(1));
    zub_ar3  = min(z0_ar3  + spreadPow, absB.k_ar(2));

    zlb_gar3 = max(z0_gar3 - spreadPow, absB.k_gar(1));
    zub_gar3 = min(z0_gar3 + spreadPow, absB.k_gar(2));

    stage.zlb = [zlb_mult; zlb_ar3; zlb_gar3];
    stage.zub = [zub_mult; zub_ar3; zub_gar3];

    stage.apply = @(z,p_in) apply_stage05a_kp_kd_plus_ar3_gar3(z,p_in,n,m);
end

function p = apply_stage_kp_kd_node_mult(z_log10, p, n, m)
% Apply node-wise multipliers to k_p-related and k_d-related vectors.
%
% NOTE:
%   - k_p in this model is used as per-genelet-node vectors (length m):
%       k_pr, k_pc, k_pi
%   - k_d in this model is used as per-regulator-component vectors (length n):
%       k_dr, k_dc

    z = 10.^z_log10(:);

    r_kp = z(1:m);
    r_kd = z(m+1:m+n);

    % Apply to k_p vectors (m)
    p.k_pr = p.k_pr(:) .* r_kp;
    p.k_pc = p.k_pc(:) .* r_kp;
    p.k_pi = p.k_pi(:) .* r_kp;

    % Apply to k_d vectors (n)
    p.k_dr = p.k_dr(:) .* r_kd;
    p.k_dc = p.k_dc(:) .* r_kd;
end

function p = apply_stage05a_kp_kd_plus_ar3_gar3(z_log10, p, n, m)
% Apply Stage0.5a (NEW) decision vector.
% - multiply existing kp/kd vectors by node/component multipliers
% - set absolute k_ar(3), k_gar(3)
% - (optional but recommended) sync Exp2 8nt baseline in tuned struct to base

    z = 10.^z_log10(:);

    r_kp = z(1:m);
    r_kd = z(m+1:m+n);

    k_ar3_new  = z(m+n+1);
    k_gar3_new = z(m+n+2);

    % --- apply multipliers to kp vectors (m) ---
    p.k_pr = p.k_pr(:) .* r_kp;
    p.k_pc = p.k_pc(:) .* r_kp;
    p.k_pi = p.k_pi(:) .* r_kp;

    % --- apply multipliers to kd vectors (n) ---
    p.k_dr = p.k_dr(:) .* r_kd;
    p.k_dc = p.k_dc(:) .* r_kd;

    % --- set absolute k_ar3/k_gar3 on base vectors ---
    p.k_ar(3)  = k_ar3_new;
    p.k_gar(3) = k_gar3_new;

    % --- safety sync: ensure Exp2 8nt baseline equals updated base ---
    if isfield(p,'tuned') && isfield(p.tuned,'exp2') && isfield(p.tuned.exp2,'toehold')
        idx8 = find(p.tuned.exp2.toehold == 8, 1);
        if isempty(idx8), idx8 = 1; end
        if isfield(p.tuned.exp2,'k_ar3'),  p.tuned.exp2.k_ar3(idx8)  = p.k_ar(3);  end
        if isfield(p.tuned.exp2,'k_gar3'), p.tuned.exp2.k_gar3(idx8) = p.k_gar(3); end
    end
end

%% ============== Stage builders (Exp1/2/3) ================================

function stage = build_stage_exp1(stage, p_prev, spreadPow, absB)
% Stage1 (Exp1): Tune k_ar3, k_gar3 (component 3) with local +/-spreadPow,
% clipped by absolute bounds (Stage0 abs bounds).

    stage.p_prev = p_prev;

    k_ar3_0  = p_prev.k_ar(3);
    k_gar3_0 = p_prev.k_gar(3);

    stage.z0  = log10([k_ar3_0; k_gar3_0]);

    % absB fields are [lb;ub] in log10
    absLB = [absB.k_ar(1);  absB.k_gar(1)];
    absUB = [absB.k_ar(2);  absB.k_gar(2)];

    % local bounds (around current best) clipped by absolute bounds
    stage.zlb = max(stage.z0 - spreadPow, absLB);
    stage.zub = min(stage.z0 + spreadPow, absUB);

    stage.apply = @(z,p_in) apply_stage_exp1(z,p_in);
end


function p = apply_stage_exp1(z_log10, p)
% Apply tuned k_ar3, k_gar3 to base vectors (component 3).
% Also synchronize tuned.exp2 baseline (8nt) to the updated base, for safety.

    z = 10.^z_log10(:);

    p.k_ar(3)  = z(1);
    p.k_gar(3) = z(2);

    % --- Safety sync: ensure exp2 8nt baseline equals base after Stage1 ---
    if isfield(p,'tuned') && isfield(p.tuned,'exp2') && isfield(p.tuned.exp2,'toehold')
        idx8 = find(p.tuned.exp2.toehold == 8, 1);
        if isempty(idx8), idx8 = 1; end  % fallback (your convention is [8 6 4])
        if isfield(p.tuned.exp2,'k_ar3'),  p.tuned.exp2.k_ar3(idx8)  = p.k_ar(3);  end
        if isfield(p.tuned.exp2,'k_gar3'), p.tuned.exp2.k_gar3(idx8) = p.k_gar(3); end
    end
end


function stage = build_stage_exp2(stage, p_prev, spreadPow, absBounds) %#ok<INUSD>
% Exp2: 8nt is the baseline (fixed). Fit ratios for 6nt and 4nt relative to 8nt.
% IMPORTANT: 8nt baseline is taken from the *base* rates optimized up to Stage0.5b:
%   k_ar8  = p_prev.k_ar(3)
%   k_gar8 = p_prev.k_gar(3)
%
% Ratio range: 0.001 .. 1  (log10 range: -3 .. 0)
%
% Decision variables z (log10):
%   z = [log10(r_ar6); log10(r_ar4); log10(r_gar6); log10(r_gar4)]

    stage.p_prev = p_prev;

    % Convention: p_prev.tuned.exp2.toehold = [8 6 4]
    idx8 = 1; idx6 = 2; idx4 = 3;

    % Base (fixed) 8nt values come from base vectors (Stage0.5b-optimized)
    k_ar8  = p_prev.k_ar(3);
    k_gar8 = p_prev.k_gar(3);

    % Current ratios as initial guesses.
    % If tuned values exist, use them to initialize ratios; otherwise fall back to 1.
    if isfield(p_prev,'tuned') && isfield(p_prev.tuned,'exp2') && ...
       isfield(p_prev.tuned.exp2,'k_ar3') && isfield(p_prev.tuned.exp2,'k_gar3')
        r_ar6  = safe_clamp(p_prev.tuned.exp2.k_ar3(idx6)  / k_ar8,  1e-3, 1);
        r_ar4  = safe_clamp(p_prev.tuned.exp2.k_ar3(idx4)  / k_ar8,  1e-3, 1);
        r_gar6 = safe_clamp(p_prev.tuned.exp2.k_gar3(idx6) / k_gar8, 1e-3, 1);
        r_gar4 = safe_clamp(p_prev.tuned.exp2.k_gar3(idx4) / k_gar8, 1e-3, 1);
    else
        r_ar6 = 1; r_ar4 = 1; r_gar6 = 1; r_gar4 = 1;
    end

    stage.z0  = log10([r_ar6; r_ar4; r_gar6; r_gar4]);
    stage.zlb = -3 * ones(4,1);
    stage.zub =  0 * ones(4,1);

    stage.apply = @(z,p_in) apply_stage_exp2(z,p_in);
end


function p = apply_stage_exp2(z_log10, p)
% Apply ratios to 6/4 nt, keep 8 nt baseline fixed.
% IMPORTANT: baseline is taken from base vectors:
%   k_ar8  = p.k_ar(3)
%   k_gar8 = p.k_gar(3)

    z = 10.^z_log10(:);

    idx8 = 1; idx6 = 2; idx4 = 3;

    % Fixed bases (8nt) from base rates (Stage0.5b-optimized)
    k_ar8  = p.k_ar(3);
    k_gar8 = p.k_gar(3);

    % z = [r_ar6; r_ar4; r_gar6; r_gar4]
    p.tuned.exp2.k_ar3(idx6)  = k_ar8  * z(1);
    p.tuned.exp2.k_ar3(idx4)  = k_ar8  * z(2);
    p.tuned.exp2.k_gar3(idx6) = k_gar8 * z(3);
    p.tuned.exp2.k_gar3(idx4) = k_gar8 * z(4);

    % Enforce 8nt tuned baseline to equal base baseline (explicit synchronization)
    p.tuned.exp2.k_ar3(idx8)  = k_ar8;
    p.tuned.exp2.k_gar3(idx8) = k_gar8;
end


function stage = build_stage_exp3(stage, p_prev, spreadPow, absBounds) %#ok<INUSD>
% Exp3: 8nt is the baseline (fixed). Fit ratios for 6nt and 4nt relative to 8nt.
% IMPORTANT: 8nt baseline is taken from the *base* rates optimized up to Stage0.5b:
%   k_bc8  = p_prev.k_bc(3)
%   k_gbc8 = p_prev.k_gbc(3)
%
% Ratio range: 0.001 .. 1  (log10 range: -3 .. 0)
%
% Decision variables z (log10):
%   z = [log10(r_bc6); log10(r_bc4); log10(r_gbc6); log10(r_gbc4)]

    stage.p_prev = p_prev;

    % Convention: p_prev.tuned.exp3.toehold = [8 6 4]
    idx8 = 1; idx6 = 2; idx4 = 3;

    % Base (fixed) 8nt values come from base vectors (Stage0.5b-optimized)
    k_bc8  = p_prev.k_bc(3);
    k_gbc8 = p_prev.k_gbc(3);

    % Current ratios as initial guesses.
    % If tuned values exist, use them to initialize ratios; otherwise fall back to 1.
    if isfield(p_prev,'tuned') && isfield(p_prev.tuned,'exp3') && ...
       isfield(p_prev.tuned.exp3,'k_bc3') && isfield(p_prev.tuned.exp3,'k_gbc3')
        r_bc6  = safe_clamp(p_prev.tuned.exp3.k_bc3(idx6)  / k_bc8,  1e-3, 1);
        r_bc4  = safe_clamp(p_prev.tuned.exp3.k_bc3(idx4)  / k_bc8,  1e-3, 1);
        r_gbc6 = safe_clamp(p_prev.tuned.exp3.k_gbc3(idx6) / k_gbc8, 1e-3, 1);
        r_gbc4 = safe_clamp(p_prev.tuned.exp3.k_gbc3(idx4) / k_gbc8, 1e-3, 1);
    else
        r_bc6 = 1; r_bc4 = 1; r_gbc6 = 1; r_gbc4 = 1;
    end

    stage.z0  = log10([r_bc6; r_bc4; r_gbc6; r_gbc4]);
    stage.zlb = -3 * ones(4,1);
    stage.zub =  0 * ones(4,1);

    stage.apply = @(z,p_in) apply_stage_exp3(z,p_in);
end


function p = apply_stage_exp3(z_log10, p)
% Apply ratios to 6/4 nt, keep 8 nt baseline fixed.
% IMPORTANT: baseline is taken from base vectors:
%   k_bc8  = p.k_bc(3)
%   k_gbc8 = p.k_gbc(3)

    z = 10.^z_log10(:);

    idx8 = 1; idx6 = 2; idx4 = 3;

    % Fixed bases (8nt) from base rates (Stage0.5b-optimized)
    k_bc8  = p.k_bc(3);
    k_gbc8 = p.k_gbc(3);

    % z = [r_bc6; r_bc4; r_gbc6; r_gbc4]
    p.tuned.exp3.k_bc3(idx6)  = k_bc8  * z(1);
    p.tuned.exp3.k_bc3(idx4)  = k_bc8  * z(2);
    p.tuned.exp3.k_gbc3(idx6) = k_gbc8 * z(3);
    p.tuned.exp3.k_gbc3(idx4) = k_gbc8 * z(4);

    % Enforce 8nt tuned baseline to equal base baseline (explicit synchronization)
    p.tuned.exp3.k_bc3(idx8)  = k_bc8;
    p.tuned.exp3.k_gbc3(idx8) = k_gbc8;
end


%% ============== GA runner ===============================================
function [bestX,bestF] = run_ga(obj, x0, lb, ub, gaOpts)
    if isempty(gaOpts)
        gaOpts = optimoptions('ga','Display','iter','UseParallel',true);
    end

    gaOpts = optimoptions(gaOpts, 'MutationFcn', @mutationadaptfeasible);

    nvars = numel(x0);
    [bestX,bestF] = ga(obj,nvars,[],[],[],[],lb,ub,[],gaOpts);
end

%% ============== fmincon runner (local refine) ============================
function [bestX,bestF,exitflag,output] = run_fmincon(obj, x0, lb, ub)
% Local refine with bound constraints.
% Requires Optimization Toolbox (fmincon).

    opts = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','interior-point', ...
        'MaxFunctionEvaluations', 2e4, ...
        'MaxIterations', 500, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-10);

    problem = struct();
    problem.objective = obj;
    problem.x0 = x0(:);
    problem.lb = lb(:);
    problem.ub = ub(:);
    problem.Aineq = [];
    problem.bineq = [];
    problem.Aeq = [];
    problem.beq = [];
    problem.nonlcon = [];
    problem.solver = 'fmincon';
    problem.options = opts;

    [bestX,bestF,exitflag,output] = fmincon(problem);
end

%% ============== Core simulation utilities ===============================
function [t, Y, ok] = simulate_IFFL(p, n, m)
    ok = true;
    tfinal = 14400;            % 4 h
    Tpulse = p.T;
    if ~isfinite(Tpulse) || Tpulse < 0, Tpulse = inf; end
    t1_end = min(Tpulse, tfinal);

    Y0 = initial_conditions(p)*1e-9;

    AbsTol_vec = [
        repmat(1e-10,5,1);
        repmat(1e-10,5,1);
        repmat(1e-10,3,1);
        repmat(1e-10,3,1);
        repmat(1e-9 ,3,1);
        repmat(1e-9 ,3,1);
        repmat(1e-9 ,3,1);
        repmat(1e-9 ,3,1);
        repmat(1e-9 ,3,1);
        1e-9
    ];
    odeOpts = odeset('RelTol',1e-7,'AbsTol',AbsTol_vec,'NonNegative',1:32,'MaxStep',100);

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
    dR1_idx = 23;
    Y0_2(dR1_idx) = Y0_2(dR1_idx) + p2.added_dR_amount;

    f2 = @(t,y) General_genelet_model_segment(t,y,p2,n,m);
    [t2, Y2] = ode15s(f2, [t1_end tfinal], Y0_2, odeOpts);
    if any(~isfinite(Y2(:))) || isempty(t2), t=[]; Y=[]; ok=false; return; end

    t = [t1; t2(2:end)];
    Y = [Y1; Y2(2:end,:)];
end

function dYdt = General_genelet_model_segment(t,Y,p,n,m) %#ok<INUSD>
    k_gar   = p.k_gar;
    k_ar    = p.k_ar;
    k_pr    = p.k_pr;
    k_ir    = p.k_ir;
    k_dr    = p.k_dr;
    k_ga    = p.k_ga;
    k_gb    = p.k_gb;
    k_gab   = p.k_gab;
    k_gbc   = p.k_gbc;
    k_bc    = p.k_bc;
    k_pc    = p.k_pc;
    k_dc    = p.k_dc;
    k_pi    = p.k_pi;
    k_gar_d = p.k_gar_d;

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

%% ============== Initial conditions & params ==============================
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

function p = default_params()
    p.k_gar = 1e5; p.k_ar = 1e4; p.k_ir = 1e5; p.k_ga = 1e4; p.k_gb = 1e5;
    p.k_gab = 1e4; p.k_gbc = 1e3; p.k_bc = 1e4;
    p.k_p = 0.10; p.k_d = 0.0015;

    p.k_aa_input = 1e-3;
    p.k_aa_other = 1e-4;

    p.k_gar_d = 5e5;
    p.leak = 0.02;

    p.G_tot  = zeros(5,1);
    p.dA_tot = zeros(3,1);
    p.dB_tot = zeros(3,1);
    p.aux = struct();

    p.T = NaN;
end

function p = apply_optimized_params(p,opt)
    G2C1 = opt(1);  G2C3 = opt(2);  G3R1 = opt(3);  G1R3 = opt(4);  G1S1 = opt(5);
    dA2  = opt(6);  dA3  = opt(7);  dA1  = opt(8);
    dB2  = opt(9);  dB3  = opt(10); dB1  = opt(11);

    p.G_tot  = [G2C1 G2C3 G3R1 G1R3 G1S1]'*1e-9;
    p.dA_tot = [dA2 dA3 dA1]'*1e-9;
    p.dB_tot = [dB2 dB3 dB1]'*1e-9;

    p.aux = struct('G2C1_tot',G2C1,'G2C3_tot',G2C3,'G3R1_tot',G3R1,'G1R3_tot',G1R3,'G1S1_tot',G1S1, ...
                   'dA2_tot',dA2,'dA3_tot',dA3,'dA1_tot',dA1, ...
                   'dB2_tot',dB2,'dB3_tot',dB3,'dB1_tot',dB1);
end

function p = override_total_G3R1_and_dB3(p, G3R1_new_nM)
% Override G3R1_tot and also increase dB3_tot by the same delta so that
% free-floating dB3 (dB3_tot - G3R1_tot) stays constant.

    % current totals (nM) from aux
    G3R1_old = p.aux.G3R1_tot;
    dB3_old  = p.aux.dB3_tot;

    % delta and updated dB3
    dG = G3R1_new_nM - G3R1_old;
    dB3_new_nM = dB3_old + dG;

    % write back
    p.aux.G3R1_tot = G3R1_new_nM;
    p.aux.dB3_tot  = dB3_new_nM;

    % rebuild vectors (nM -> M)
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

function p = override_total_G2C1_and_G2C3(p, G2C1_new_nM, G2C3_new_nM)
% Override G2C1_tot and G2C3_tot (in nM) and rebuild p.G_tot (in M).
% NOTE: This function does NOT co-adjust dA/dB totals. If your experiment
% requires coupled changes, modify here accordingly.

    p.aux.G2C1_tot = G2C1_new_nM;
    p.aux.G2C3_tot = G2C3_new_nM;

    % Rebuild G_tot vector (nM -> M)
    G2C1 = p.aux.G2C1_tot;
    G2C3 = p.aux.G2C3_tot;
    G3R1 = p.aux.G3R1_tot;
    G1R3 = p.aux.G1R3_tot;
    G1S1 = p.aux.G1S1_tot;

    p.G_tot = [G2C1 G2C3 G3R1 G1R3 G1S1]' * 1e-9;
end



%% ============== Embedded AUC data (user-provided) ========================
function AUC = build_embedded_AUC_data(Din_vals)
% Build a struct with embedded AUC mean arrays for Exp1/2/3.
% All are measured at Din=[20 40 60] (minutes).

AUC = struct();
AUC.Din_vals = Din_vals(:).';

% Exp1: (FIX) G3R1_tot = 5,10,25 nM
AUC.exp1.G3R1_levels_nM = [5 10 25];
AUC.exp1.mean_5nM  = [2.56499499837718, 1.96172243338906, 0.864853882077734];
AUC.exp1.mean_10nM = [1.57539880189863, 0.662969895854546, 0.294228900361590];
AUC.exp1.mean_25nM = [0.498801084899379, 0.234218003967067, 0.335176590776887];

% Exp2: dA1 toehold length = 8,6,4 nt (nominal=8)
AUC.exp2.toehold_nt = [8 6 4];
AUC.exp2.mean_8b = [2.62112235459656, 1.75124015886793, 1.38857435443720];
AUC.exp2.mean_6b = [3.03573082836344, 2.90763225699382, 2.45183187679066];
AUC.exp2.mean_4b = [3.09994898455949, 2.97298907575560, 2.71276629624907];

% Exp3: dB1 toehold length = 8,6,4 nt (nominal=8)
AUC.exp3.toehold_nt = [8 6 4];
AUC.exp3.mean_8b = [2.56173667051312, 2.09954690172282, 0.543000279185901];
AUC.exp3.mean_6b = [0.523584825432856, 0.487934703430427, -0.0421065906067014];
AUC.exp3.mean_4b = [0.184195637493338, -0.232970384315233, -0.0902961430380806];

% --- Exp4: AUC vs (G2C1_tot, G2C3_tot) conditions (3 conditions) ---
% In this implementation, we map [G8C1_tot, G8C3_tot] onto the model totals
% [G2C1_tot, G2C3_tot] by overriding p.aux.G2C1_tot and p.aux.G2C3_tot.
AUC.exp4.G2C1_G2C3_levels_nM = [10  2;
                               25  5;
                               50 10];   % 3x2 (nM)

% Experimental AUC values (3x1). Replace NaNs with your measured values.
AUC.exp4.Y = [2.72745037326219	2.64401754542818	1.24894291418106];

% Pulse duration (Din) for Exp4 AUC evaluation (minutes). Adjust to match your experiment.
AUC.exp4.Din_min = 30;


% Package in consistent matrices (3 levels x 3 Din)
AUC.exp1.Y = [AUC.exp1.mean_5nM;  AUC.exp1.mean_10nM; AUC.exp1.mean_25nM];
AUC.exp2.Y = [AUC.exp2.mean_8b;   AUC.exp2.mean_6b;  AUC.exp2.mean_4b];
AUC.exp3.Y = [AUC.exp3.mean_8b;   AUC.exp3.mean_6b;  AUC.exp3.mean_4b];
end

%% ============== Reporting ===============================================
function print_key_params(p)
    fprintf('k_bc (n=3):   [%.3g %.3g %.3g]\n', p.k_bc(1), p.k_bc(2), p.k_bc(3));
    fprintf('k_gbc (n=3):  [%.3g %.3g %.3g]\n', p.k_gbc(1), p.k_gbc(2), p.k_gbc(3));
    fprintf('k_ar  (n=3):  [%.3g %.3g %.3g]\n', p.k_ar(1), p.k_ar(2), p.k_ar(3));
    fprintf('k_gar (n=3):  [%.3g %.3g %.3g]\n', p.k_gar(1), p.k_gar(2), p.k_gar(3));
    fprintf('leak: %.3g\n', p.leak);

    if isfield(p,'tuned') && isfield(p.tuned,'exp2')
        fprintf('Exp2 tuned (k_ar3 by toehold [8 6 4]):  [%.3g %.3g %.3g]\n', p.tuned.exp2.k_ar3);
        fprintf('Exp2 tuned (k_gar3 by toehold [8 6 4]): [%.3g %.3g %.3g]\n', p.tuned.exp2.k_gar3);
    end
    if isfield(p,'tuned') && isfield(p.tuned,'exp3')
        fprintf('Exp3 tuned (k_bc3 by toehold [8 6 4]):  [%.3g %.3g %.3g]\n', p.tuned.exp3.k_bc3);
        fprintf('Exp3 tuned (k_gbc3 by toehold [8 6 4]): [%.3g %.3g %.3g]\n', p.tuned.exp3.k_gbc3);
    end
end

function I = calculate_AUC(t_h, y)
% Right Riemann sum in hour units: sum_{i=2..N} (t_i - t_{i-1}) * y_i
    I = 0;
    for i = 2:length(t_h)
        I = I + (t_h(i)-t_h(i-1)) * y(i);
    end
end

function y = safe_clamp(x, lo, hi)
% Clamp x into [lo, hi], guarding NaN/Inf.
    if ~isfinite(x) || x <= 0
        x = 1; % fall back to neutral ratio
    end
    y = min(max(x, lo), hi);
end

function stage = build_stage05b_allparams(stage, absB, multLo, multHi, n, m)
% Stage0.5b: optimize ALL baseline rates + node-wise kp + component-wise kd (log10 variables).
%
% Decision vector x (log10):
%   x = [ log10(k_gar_common);          % common (replicated to all n)
%         log10(k_ar_common);           % common (replicated to all n)
%         log10(k_ar3);                 % %% NEW: component-3 specific
%         log10(k_gar3);                % %% NEW: component-3 specific
%         log10(k_ir);
%         log10(k_ga);
%         log10(k_gb);
%         log10(k_gab);
%         log10(k_gbc);
%         log10(k_bc);
%         log10(k_aa_input);
%         log10(k_aa_other);
%         log10(k_gar_d);
%         log10(leak);
%         log10(kp_nodes(1:m));
%         log10(kd_comps(1:n)) ]
%
% Bounds:
%   - "other" rates: Stage0 abs bounds
%   - kp,kd: Stage0 abs bounds expanded by multLo/multHi (absolute)

    p0 = stage.p_prev;

    % --- extract current values from p_prev (Stage0.5a result) ---
    % common scalar-ish rates (take element 1 of vectors)
    k_gar_common = p0.k_gar(1);
    k_ar_common  = p0.k_ar(1);

    % %% NEW: keep component-3 specific values as separate decision vars
    k_ar3  = p0.k_ar(3);
    k_gar3 = p0.k_gar(3);

    k_ir  = p0.k_ir(1);
    k_ga  = p0.k_ga(1);
    k_gb  = p0.k_gb(1);
    k_gab = p0.k_gab(1);
    k_gbc = p0.k_gbc(1);
    k_bc  = p0.k_bc(1);

    k_aa_input = p0.k_aa_input;
    k_aa_other = p0.k_aa_other;
    k_gar_d    = p0.k_gar_d(1);
    leak       = p0.leak;

    % node-wise kp: use k_pr (length m)
    kp_nodes = p0.k_pr(:);
    % component-wise kd: use k_dr (length n)
    kd_comps = p0.k_dr(:);

    % --- build x0 (log10) ---
    x0 = log10([k_gar_common;
                k_ar_common;
                k_ar3;        % %% NEW
                k_gar3;       % %% NEW
                k_ir;
                k_ga;
                k_gb;
                k_gab;
                k_gbc;
                k_bc;
                k_aa_input;
                k_aa_other;
                k_gar_d;
                leak;
                kp_nodes;
                kd_comps]);

    % --- bounds for the 14 "other" rates (log10) from absB ---
    % absB fields: [lb;ub] in log10 for each rate

    xlb_other = [absB.k_gar(1);
                 absB.k_ar(1);
                 absB.k_ar(1);     % %% NEW: k_ar3 uses same absolute bounds as k_ar
                 absB.k_gar(1);    % %% NEW: k_gar3 uses same absolute bounds as k_gar
                 absB.k_ir(1);
                 absB.k_ga(1);
                 absB.k_gb(1);
                 absB.k_gab(1);
                 absB.k_gbc(1);
                 absB.k_bc(1);
                 absB.k_aa_input(1);
                 absB.k_aa_other(1);
                 absB.k_gar_d(1);
                 absB.leak(1)];

    xub_other = [absB.k_gar(2);
                 absB.k_ar(2);
                 absB.k_ar(2);     % %% NEW
                 absB.k_gar(2);    % %% NEW
                 absB.k_ir(2);
                 absB.k_ga(2);
                 absB.k_gb(2);
                 absB.k_gab(2);
                 absB.k_gbc(2);
                 absB.k_bc(2);
                 absB.k_aa_input(2);
                 absB.k_aa_other(2);
                 absB.k_gar_d(2);
                 absB.leak(2)];

    % --- bounds for kp/kd in absolute scale (log10 expanded by multLo/multHi) ---
    kp_lb = absB.k_p(1) + log10(multLo);
    kp_ub = absB.k_p(2) + log10(multHi);

    kd_lb = absB.k_d(1) + log10(multLo);
    kd_ub = absB.k_d(2) + log10(multHi);

    xlb = [xlb_other;
           kp_lb * ones(m,1);
           kd_lb * ones(n,1)];

    xub = [xub_other;
           kp_ub * ones(m,1);
           kd_ub * ones(n,1)];

    stage.x0  = x0;
    stage.xlb = xlb;
    stage.xub = xub;

    stage.apply = @(x,p_in) apply_stage05b_allparams(x,p_in,n,m);
end

function p = apply_stage05b_allparams(x_log10, p, n, m)
% Apply the full decision vector to parameter struct p.
% Must match build_stage05b_allparams ordering.

    x = 10.^x_log10(:);

    % unpack (must match build_stage05b_allparams)
    idx = 1;

    k_gar_common = x(idx); idx=idx+1;
    k_ar_common  = x(idx); idx=idx+1;

    % %% NEW: component-3 specific values
    k_ar3  = x(idx); idx=idx+1;
    k_gar3 = x(idx); idx=idx+1;

    k_ir  = x(idx); idx=idx+1;
    k_ga  = x(idx); idx=idx+1;
    k_gb  = x(idx); idx=idx+1;
    k_gab = x(idx); idx=idx+1;
    k_gbc = x(idx); idx=idx+1;
    k_bc  = x(idx); idx=idx+1;

    k_aa_input = x(idx); idx=idx+1;
    k_aa_other = x(idx); idx=idx+1;
    k_gar_d    = x(idx); idx=idx+1;
    leak       = x(idx); idx=idx+1;

    kp_nodes = x(idx:idx+m-1); idx = idx + m;
    kd_comps = x(idx:idx+n-1); idx = idx + n;

    % apply scalars (vectors are replicated)
    p.k_gar = k_gar_common * ones(n,1);
    p.k_ar  = k_ar_common  * ones(n,1);

    % %% NEW: preserve component-3 specificity (do NOT get overwritten)
    p.k_ar(3)  = k_ar3;
    p.k_gar(3) = k_gar3;

    p.k_ir  = k_ir  * ones(n,1);
    p.k_ga  = k_ga  * ones(n,1);
    p.k_gb  = k_gb  * ones(n,1);
    p.k_gab = k_gab * ones(n,1);
    p.k_gbc = k_gbc * ones(n,1);
    p.k_bc  = k_bc  * ones(n,1);

    p.k_aa_input = k_aa_input;
    p.k_aa_other = k_aa_other;

    p.k_gar_d = k_gar_d * ones(n,1);
    p.leak    = leak;

    % node-wise kp (m): apply to all three
    p.k_pr = kp_nodes(:);
    p.k_pc = kp_nodes(:);
    p.k_pi = kp_nodes(:);

    % component-wise kd (n): apply to both
    p.k_dr = kd_comps(:);
    p.k_dc = kd_comps(:);
end


function idx = find_stage_by_prefix(stages, prefix)
% Return first index i such that stages(i).name starts with prefix. If not found, [].

    idx = [];
    if isempty(stages), return; end

    for ii = 1:numel(stages)
        if ~isfield(stages(ii),'name'), continue; end

        nm = stages(ii).name;

        if ischar(nm)
            if strncmp(nm, prefix, length(prefix))
                idx = ii;
                return;
            end
        elseif isstring(nm)
            if startsWith(nm, prefix)
                idx = ii;
                return;
            end
        end
    end
end

function [bestZ,bestF] = run_stage_optimizer(obj, z0, zlb, zub, gaOpts, method)
% Run stage optimizer using 'ga' or 'fmincon' on decision vector z (log10 vars).

    if nargin < 6 || isempty(method), method = 'ga'; end
    method = lower(string(method));

    switch method
        case "ga"
            [bestZ,bestF] = run_ga(obj, z0, zlb, zub, gaOpts);

        case "fmincon"
            % interior-point likes a strictly feasible start
            epsIn = 1e-10;
            z0_fmin = min(max(z0(:), zlb(:) + epsIn), zub(:) - epsIn);

            % run_fmincon is already defined in your file
            [bestZ,bestF] = run_fmincon(obj, z0_fmin, zlb, zub);

        otherwise
            error('Unknown expOptMethod="%s". Use "ga" or "fmincon".', method);
    end
end
