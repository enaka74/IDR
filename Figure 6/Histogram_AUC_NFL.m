% import IDD score data
filename = 'NFL_signaling_parameter_screening_RandomParameter_minscale.xlsx';
sheetname = 'Sheet1';
data = readmatrix(filename,"Sheet",sheetname,"range","B1:B5000"); 

% Create histogram
edges = logspace(-6,2,45);
%edges = [0:0.25:20];
h = histogram(data,edges);
h.FaceColor = [93, 173, 226]/256;
h.LineWidth = 2;
set(gca, 'Xscale', 'log');
xlabel("IDA Score")
ylabel("# of parameter sets")
xlim([1e-6,1e2])
ylim([0,30])
xticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2])
set(gca, 'FontSize', 30, 'LineWidth', 3);
box on

set(gcf, 'Renderer', 'painters');
print('filename.eps', '-depsc', '-r300');