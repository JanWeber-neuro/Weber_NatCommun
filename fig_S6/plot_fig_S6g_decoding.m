function plot_fig_S6g_decoding(path)

% colors
cmap    = cbrewer('div', 'Spectral', 200, 'PCHIP');

cb       = cell(2,1);
cb{1}    = cmap(60,:);
cb{2}    = cmap(190,:);

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6g');

time = table2array(tbl(1,2:end));

%% Plot decoding for PFC

ga_pfc_context = table2array(tbl(3:20,2:end));
ga_pfc_context(all(isnan(ga_pfc_context),2),:) = [];

ga_pfc_action = table2array(tbl(24:41,2:end));
ga_pfc_action(all(isnan(ga_pfc_action),2),:) = [];

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
yline(0.33, 'LineWidth', 3, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off

set(gcf, 'Position', [596 417 229 198]);

%% Plot decoding for motor cortex

ga_motor_context = table2array(tbl(45:62,2:end));
ga_motor_context(all(isnan(ga_motor_context),2),:) = [];

ga_motor_action = table2array(tbl(66:83,2:end));
ga_motor_action(all(isnan(ga_motor_action),2),:) = [];

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
yline(0.33, 'LineWidth', 3, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off

set(gcf, 'Position', [596 417 229 198]);
