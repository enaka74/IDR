clear all
close all

% ============================
% Experiment-specific settings
% ============================

% Control data file (used for baseline/min normalization)
ctr_filename   = 'PD_R_0x_20251218_control_ave.xlsx';
ctr_sheetname  = 'PD_R_0x_20251218_control_ave';
ctr_min        = 2;
ctr_baseline   = 1;

% Sample (main) data file
filename = 'PD_R_05_20251203_60min_3rep.xlsx';
sheetname = 'PD_R_05_20251203_60min_3rep';

% Experiment layout
N_s = 9;  % number of samples (tiles)
N_n = 3;  % number of channels
N   = 9;  % number of samples excluding controls

Sample_label = ["dA1\_8b, Rep#1" ...
                "dA1\_8b, Rep#2" ...
                "dA1\_8b, Rep#3" ...
                "dA1\_6b, Rep#1" ...
                "dA1\_6b, Rep#2" ...
                "dA1\_6b, Rep#3" ...
                "dA1\_4b, Rep#1" ...
                "dA1\_4b, Rep#2" ...
                "dA1\_4b, Rep#3"];

cond_labels = ["dA1\_8b", "dA1\_6b", "dA1\_4b"];

Channels = ["G_{in}C_{out} (TYE563)", "G_{r}R_{out} (TEX615)", "G_{out}S (6-FAM)"];
interval = 2; % minutes per data point

% Common time axis aligned at T7+1: 0..240 min (121 points at 2-min interval)
t_T7_axis = (0:120) * interval / 60;  % [hours]

% Well/row layout for sample data
% well_list = ["A6" "B6" "C6" "D6" "E6" "F6" "G6" "H6" "I6"];
well_list = ["A14" "D14" "G14" "B14" "E14" "H14" "C14" "F14" "I14"];
row_list  = ["333" "522" "139" "328" "527" "716"];  % 2 rows per channel

% Well/row layout for control data
% ctr_well_list = ["M1" "N1"];
% ctr_row_list  = ["382" "621" "138" "377" "626" "865"]; % 2 rows per channel
%ctr_well_list = ["B11" "F11"];
%ctr_row_list  = ["400" "678" "117" "395" "683" "961"]; % 2 rows per channel
ctr_well_list = ["A8" "B8"];
ctr_row_list  = ["340" "558" "117" "335" "563" "781"]; % 12/11/2025 or 12/18/2025 ave, 2 rows per channel

% For 'PD_R_0x_20251218_control_ave.xlsx'. Samples were loaded into A3:H5, 

%ctr_well_list = ["A4" "B4"];
%ctr_row_list  = ["340" "558" "117" "335" "563" "781"]; % 12/18/2025, 2 rows per channel


% -------------------------------
% Reagent addition (index points)
% -------------------------------
t_dA   = 16;
t_T7   = 38;
t_dR   = 69;             % samples
ctr_t_dA = 16;
ctr_t_T7 = 38;
ctr_t_dR = 69;           % controls

% For samples: [dA, T7, dR]; later rows may have NaN for no further additions
time_pts = [t_dA,t_T7,t_dR;
            t_dA,t_T7,t_dR;
            t_dA,t_T7,t_dR;
            t_dA,t_T7,t_dR;
            t_dA,t_T7,t_dR;
            t_dA,t_T7,t_dR;
            t_dA,t_T7,t_dR;
            t_dA,t_T7,t_dR;
            t_dA,t_T7,t_dR];

% Template for segment-wise smoothing (not used directly below)
time_pts_norm = repmat([t_dA,t_T7], N_s, 1);

% Reference points used in normalization (index into vectors)
% [ref_for_ch1, ref_for_ch2, ref_for_ch3]
% Here: use dA for ch1, T7+1 for ch2 and ch3
time_pts_ref = [t_dA, t_T7+1, t_T7+1];

% Control: define segment change points (NaN to stop smoothing safely)
ctr_time_pts = [ctr_t_dA, ctr_t_T7, NaN;   % baseline row
                ctr_t_dA, NaN,      NaN];  % min row

% Optional: if you need a 'two-point' normalization for controls
ctr_time_pts_norm = [ctr_t_T7, NaN;
                     ctr_t_dA, NaN];

% -----------------------
% Smoothing configuration
% -----------------------
% 'window' and 'method' are per-channel arrays (3 channels)
window = [1,1,1;
          1,1,1;
          0.2,0.2,0.2];

m1="movmean"; m2="movmedian"; m3="gaussian"; m4="lowess"; m5="loess"; m6="rlowess"; m7="rloess"; m8="sgolay";
method = [m4,m4,m4];

% Baseline-only smoothing window override (per channel).
% If you want to use a different window only for the control baseline row:
%  - No override: leave as []
%  - One value for all channels: 0.5 (scalar or vector)
%  - Per-channel cell: { [5 5], 3, [0.5 0.5 0.5] }
%    (inside the cell, [] or NaN falls back to window(j,:))
ctr_baseline_window = [1,1,1;
                       1,1,1;
                       1,1,1];

% ============================
% Run common analysis pipeline
% ============================
plate_analysis_core_ratio_NormExample;
