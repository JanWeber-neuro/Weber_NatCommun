function plot_fig_5c_decoding(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_5c');

% colors
cmap    = cbrewer('div', 'Spectral', 200, 'PCHIP');

cb       = cell(2,1);
cb{1}    = cmap(60,:);
cb{2}    = cmap(190,:);

time = table2array(tbl(1,2:end));

%% Plot PFC

pfc_context = table2array(tbl(3:20,2:end)); 
pfc_context(all(isnan(pfc_context),2),:) = [];
pfc_action  = table2array(tbl(24:41,2:end)); 
pfc_action(all(isnan(pfc_action),2),:) = [];

% plot data
figure;
a = shadedErrorBar(time, mean(pfc_context), ...
    std(pfc_context) / sqrt(size(pfc_context,1)));
a.mainLine.Color = cb{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{1}; a.patch.EdgeColor = cb{1}; a.patch.FaceAlpha = 0.5;

hold on

a = shadedErrorBar(time, mean(pfc_action), ...
    std(pfc_action) / sqrt(size(pfc_action,1)));
a.mainLine.Color = cb{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{2}; a.patch.EdgeColor = cb{2}; a.patch.FaceAlpha = 0.5;

xlabel('Time to BR [s]');

ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 3, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.22 0.47]); 
box off

set(gcf, 'Position', [596 417 229 198]);

%% Plot motor cortex

motor_context = table2array(tbl(45:61,2:end)); 
motor_context(all(isnan(motor_context),2),:) = [];
motor_action  = table2array(tbl(67:82,2:end)); 
motor_action(all(isnan(motor_action),2),:) = [];

% plot data
figure;
a = shadedErrorBar(time, mean(motor_context), ...
    std(motor_context) / sqrt(size(motor_context,1)));
a.mainLine.Color = cb{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{1}; a.patch.EdgeColor = cb{1}; a.patch.FaceAlpha = 0.5;

hold on

a = shadedErrorBar(time, mean(motor_action), ...
    std(motor_action) / sqrt(size(motor_action,1)));
a.mainLine.Color = cb{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{2}; a.patch.EdgeColor = cb{2}; a.patch.FaceAlpha = 0.5;

xlabel('Time to BR [s]');

ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 3, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.22 0.47]); 
box off

set(gcf, 'Position', [596 417 229 198]);
