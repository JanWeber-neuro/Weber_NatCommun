function plot_fig_S6ab_MDD(path)

%% Plot MDD PFC

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6a left');

time = table2array(tbl(1,2:end));
ga = table2array(tbl(3:19,2:end));

figure;
    
a = shadedErrorBar(time, mean(ga), std(ga) / sqrt(size(ga,1)));
a.mainLine.Color = 'r'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'r'; a.patch.EdgeColor = 'r'; a.patch.FaceAlpha = 0.5;

xlim([time(1) time(end)])
ylim([-1 6]);
ylabel('MDD', 'FontSize', 13);
set(gca, 'fontsize', 13, 'linewidth', 3);
box off
set(gcf, 'Position', [2602 472 205 151]);
    
%% Plot MDD motor cortex

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6b left');

time = table2array(tbl(1,2:end));
ga = table2array(tbl(3:16,2:end));

figure;
    
a = shadedErrorBar(time, mean(ga), std(ga) / sqrt(size(ga,1)));
a.mainLine.Color = 'b'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'b'; a.patch.EdgeColor = 'b'; a.patch.FaceAlpha = 0.5;

xlim([time(1) time(end)])
ylim([-4 7]);
ylabel('MDD', 'FontSize', 13);
set(gca, 'fontsize', 13, 'linewidth', 3);
box off
set(gcf, 'Position', [2602 472 205 151]);
