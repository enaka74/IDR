% import IDD score data
sheetname = 'Sheet1';
edges = logspace(-6,2,45);

data_IFFL = readmatrix('IFFL_GRN_parameter_screening_RandomParameter_rescaled_yFloor_u=0.1.xlsx', ...
    "Sheet",sheetname,"range","A1:A10000");
data_NFL  = readmatrix('NFL_GRN_parameter_screening_RandomParameter_rescaled_yFloor_u=0.1.xlsx', ...
    "Sheet",sheetname,"range","A1:A10000");

figure; hold on

% --- IFFL histogram ---
h1 = histogram(data_IFFL, edges, 'Normalization','probability');
h1.FaceColor = [240, 178, 122]/256;
h1.FaceAlpha = 0.75;
h1.EdgeColor = 'none';
h1.LineWidth = 2;

% --- NFL histogram ---
h2 = histogram(data_NFL, edges, 'Normalization','probability');
h2.FaceColor = [119, 202, 188]/256;
h2.FaceAlpha = 0.75;
h2.EdgeColor = 'none';
h2.LineWidth = 2;

% axes styling
set(gca, 'Xscale', 'log');
xlabel("IDD Score")
ylabel("Fraction of parameter sets")
xlim([1e-6,1e2])
ylim([0,0.5e-2])
xticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2])
set(gca, 'FontSize', 30, 'LineWidth', 3);
box on

ax = gca;
ax.Units = 'normalized';
ax.Position = [0.15 0.25 0.85 0.65];
ax.LooseInset = ax.TightInset;

legend({'IFFL','NFL'},'Location','best','Box','off')

% ============================================================
% Change proposal #1: tail fraction annotation at T
% ============================================================
T = 0.1;

pTail_IFFL = mean(data_IFFL >= T);
pTail_NFL  = mean(data_NFL  >= T);

xline(T,'k--','LineWidth',2);

yl = ylim;


% ============================================================
% Optional: inset zoom panel for [0.1, 100]
% (Comment out this whole block if you don't use it.)
% ============================================================
%{
ax_main = gca;

% Inset position: [left bottom width height] in normalized figure units
ax_in = axes('Position',[0.62 0.58 0.30 0.30]);  % <- tweak if needed
hold(ax_in,'on')

histogram(ax_in, data_IFFL, edges, 'Normalization','probability', ...
    'FaceColor',[240, 178, 122]/256, 'FaceAlpha',0.35, 'EdgeColor','none', 'LineWidth',2);
histogram(ax_in, data_NFL,  edges, 'Normalization','probability', ...
    'FaceColor',[119, 202, 188]/256, 'FaceAlpha',0.35, 'EdgeColor','none', 'LineWidth',2);

set(ax_in,'XScale','log');
xlim(ax_in,[0.1, 100])
xticks(ax_in,[0.1, 1, 10, 100])
box(ax_in,'on')
set(ax_in,'FontSize',18,'LineWidth',2)

% Draw the same threshold line in inset
xline(ax_in, T, 'k--', 'LineWidth', 1.5);

% Optional: make inset y-limits tight
ylim(ax_in,'auto')

% Return focus to main axes
axes(ax_main);
%}

% export
set(gcf, 'Renderer', 'painters');
print('filename.eps', '-depsc', '-r300');

hold off
