function plot_fig_S11_decoding(path)

% colors
cmap    = cbrewer('div', 'Spectral', 200, 'PCHIP');

cb       = cell(2,1);
cb{1}    = cmap(60,:);
cb{2}    = cmap(190,:);

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S11');

time = table2array(tbl(1,2:end));

%% Plot decoding for PFC

ga_pfc_context = table2array(tbl(3:19,2:end));
ga_pfc_action = table2array(tbl(24:40,2:end));

figure;
a = shadedErrorBar(time, mean(ga_pfc_context), ...
    std(ga_pfc_context) / sqrt(size(ga_pfc_context,1)));
a.mainLine.Color = cb{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{1}; a.patch.EdgeColor = cb{1}; a.patch.FaceAlpha = 0.5;

hold on

a = shadedErrorBar(time, mean(ga_pfc_action), ...
    std(ga_pfc_action) / sqrt(size(ga_pfc_action,1)));
a.mainLine.Color = cb{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{2}; a.patch.EdgeColor = cb{2}; a.patch.FaceAlpha = 0.5;


xlabel('Time to BR [s]');
ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 1, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off

set(gcf, 'Position', [411 372 422 264]);

%% Plot decoding for motor cortex

ga_motor_context = table2array(tbl(45:58,2:end));
ga_motor_action = table2array(tbl(62:75,2:end));

figure;
a = shadedErrorBar(time, mean(ga_motor_context), ...
    std(ga_motor_context) / sqrt(size(ga_motor_context,1)));
a.mainLine.Color = cb{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{1}; a.patch.EdgeColor = cb{1}; a.patch.FaceAlpha = 0.5;

hold on

a = shadedErrorBar(time, mean(ga_motor_action), ...
    std(ga_motor_action) / sqrt(size(ga_motor_action,1)));
a.mainLine.Color = cb{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{2}; a.patch.EdgeColor = cb{2}; a.patch.FaceAlpha = 0.5;


xlabel('Time to BR [s]');
ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 1, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off

set(gcf, 'Position', [411 372 422 264]);

