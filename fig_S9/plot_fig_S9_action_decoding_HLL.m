function plot_fig_S9_action_decoding_HLL(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S9');

time = table2array(tbl(1,2:end));

%% Plot for PFC

ga_pfc = table2array(tbl(2:11,2:end));

figure;
a = shadedErrorBar(time, mean(ga_pfc), ...
    std(ga_pfc) / sqrt(size(ga_pfc,1)));
a.mainLine.Color = 'r'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'r'; a.patch.EdgeColor = 'r'; a.patch.FaceAlpha = 0.5;

xlabel('Time to HLL [s]');
ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 1, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off

set(gcf, 'Position', [357 341 433 275]);

%% Plot for motor cortex

ga_motor = table2array(tbl(15:21,2:end));

figure;
a = shadedErrorBar(time, mean(ga_motor), ...
    std(ga_motor) / sqrt(size(ga_motor,1)));
a.mainLine.Color = 'b'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'b'; a.patch.EdgeColor = 'b'; a.patch.FaceAlpha = 0.5;

xlabel('Time to HLL [s]');
ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 1, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off

set(gcf, 'Position', [357 341 433 275]);

