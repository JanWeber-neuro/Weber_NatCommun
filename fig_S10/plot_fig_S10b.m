function plot_fig_S10b(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S10b');

time = table2array(tbl(1,2:end));

%% Plot PFC

ga_pfc = table2array(tbl(3:16,2:end));

figure;
a = shadedErrorBar(time, mean(ga_pfc), ...
    std(ga_pfc) / sqrt(size(ga_pfc,1)));
a.mainLine.Color = 'r'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'r'; a.patch.EdgeColor =' r'; a.patch.FaceAlpha = 0.5;

xlabel('Time to BR [s]');
ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 3, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off
set(gcf, 'Position', [2604 438 287 222]);
    
%% Plot motor cortex

ga_pfc = table2array(tbl(20:29,2:end));

figure;
a = shadedErrorBar(time, mean(ga_pfc), ...
    std(ga_pfc) / sqrt(size(ga_pfc,1)));
a.mainLine.Color = 'b'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'b'; a.patch.EdgeColor = 'b'; a.patch.FaceAlpha = 0.5;

xlabel('Time to BR [s]');
ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 3, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off
set(gcf, 'Position', [2604 438 287 222]);
    
