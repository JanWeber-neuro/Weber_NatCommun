function plot_fig_S6cd_PC1(path)

% colors for plotting
color = {'#008450', '#EFB700', '#B81D13'};

%% Plot population activity PFC

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6c left');

time = table2array(tbl(1,2:end));

ga_pfc_0  = table2array(tbl(3:19,2:end));
ga_pfc_25 = table2array(tbl(23:39,2:end));
ga_pfc_75 = table2array(tbl(43:59,2:end));

figure;

% plot first PC frontal for different conditions
a = shadedErrorBar(time, mean(ga_pfc_0), ...
    std(ga_pfc_0) / sqrt(size(ga_pfc_0,1)));
a.mainLine.Color = color{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{1}; a.patch.EdgeColor = color{1}; a.patch.FaceAlpha = 0.5;
hold on
a = shadedErrorBar(time, mean(ga_pfc_25), ...
    std(ga_pfc_25) / sqrt(size(ga_pfc_25,1)));
a.mainLine.Color = color{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{2}; a.patch.EdgeColor = color{2}; a.patch.FaceAlpha = 0.5;
hold on
a = shadedErrorBar(time, mean(ga_pfc_75), ...
    std(ga_pfc_75) / sqrt(size(ga_pfc_75,1)));
a.mainLine.Color = color{3}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{3}; a.patch.EdgeColor = color{3}; a.patch.FaceAlpha = 0.5;

ylim([-3.5 7]);

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)]); box off; hold on;

set(gcf, 'Position', [530 333 348 277]);

%% Plot population activity motor cortex

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6d left');

time = table2array(tbl(1,2:end));

ga_motor_0  = table2array(tbl(3:16,2:end));
ga_motor_25 = table2array(tbl(20:33,2:end));
ga_motor_75 = table2array(tbl(37:50,2:end));

figure;

% plot first PC frontal for different conditions
a = shadedErrorBar(time, mean(ga_motor_0), ...
    std(ga_motor_0) / sqrt(size(ga_motor_0,1)));
a.mainLine.Color = color{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{1}; a.patch.EdgeColor = color{1}; a.patch.FaceAlpha = 0.5;
hold on
a = shadedErrorBar(time, mean(ga_motor_25), ...
    std(ga_motor_25) / sqrt(size(ga_motor_25,1)));
a.mainLine.Color = color{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{2}; a.patch.EdgeColor = color{2}; a.patch.FaceAlpha = 0.5;
hold on
a = shadedErrorBar(time, mean(ga_motor_75), ...
    std(ga_motor_75) / sqrt(size(ga_motor_75,1)));
a.mainLine.Color = color{3}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{3}; a.patch.EdgeColor = color{3}; a.patch.FaceAlpha = 0.5;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)]); box off; hold on;

set(gcf, 'Position', [530 333 348 277]);
