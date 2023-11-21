function plot_fig_1d(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_1d');

%% Define colors for plotting

color.m1    = [0 0 1];
color.pfc   = [1 0 0];

%% Extract data from file

time = table2array(tbl(1,2:end));

GA_active   = [];
GA_inactive = [];

GA_active.frontal   = table2array(tbl(3:20,2:end));
GA_inactive.frontal = table2array(tbl(25:42,2:end));
  
GA_active.motor     = table2array(tbl(47:64,2:end));
GA_inactive.motor   = table2array(tbl(69:86,2:end));

% remove patients without data
GA_active.frontal(find(all(isnan(GA_active.frontal),2)),:) = [];
GA_inactive.frontal(find(all(isnan(GA_active.frontal),2)),:) = [];

GA_active.motor(find(all(isnan(GA_active.motor),2)),:) = [];
GA_inactive.motor(find(all(isnan(GA_active.motor),2)),:) = [];

% plot data
figure;
% active
a = shadedErrorBar(time, nanmean(GA_active.frontal)*100, nanstd(GA_active.frontal)*100 ...
    / sqrt(size(GA_active.frontal,1)));
a.mainLine.Color = color.pfc; a.mainLine.LineWidth = 2;
a.patch.FaceColor = color.pfc; a.patch.EdgeColor = color.pfc; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = color.pfc; a.edge(2).Color = color.pfc;
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

hold on

% inactive
plot(time, nanmean(GA_inactive.frontal)*100, 'Color', color.pfc, 'LineStyle', '--', 'LineWidth', 2);

% plot motor

% active
a = shadedErrorBar(time, nanmean(GA_active.motor)*100, nanstd(GA_active.motor)*100 ...
    / sqrt(size(GA_active.motor,1)));
a.mainLine.Color = color.m1; a.mainLine.LineWidth = 2;
a.patch.FaceColor = color.m1; a.patch.EdgeColor = color.m1; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = color.m1; a.edge(2).Color = color.m1;
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

hold on

% % inactive
plot(time, nanmean(GA_inactive.motor)*100, 'Color', color.m1, 'LineStyle', '--', 'LineWidth', 2);

% make plot nice
xline(0, 'k--', 'linewidth', 0.5); yline(0, 'k--', 'linewidth', 0.5);
xlim([time(1) time(end)]); ylabel('%EV'); 
xlabel('Time to HLL [s]')

set(gca, 'fontsize', 13, 'linewidth', 3); box off;
xticks([-0.4 -0.2 0 0.2]);

set(gcf, 'Position', [500 380 252 168]);
