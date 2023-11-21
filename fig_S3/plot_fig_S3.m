function plot_fig_S3(path)

tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S3');

time = table2array(tbl(1,2:end));

%% Plot PFC

ga_pfc_1 = table2array(tbl(3:18,2:end));
ga_pfc_2 = table2array(tbl(22:37,2:end));

figure;
% old
a = shadedErrorBar(time, mean(ga_pfc_1)*100, std(ga_pfc_1)*100 ...
    / sqrt(size(ga_pfc_1,1)));
a.mainLine.Color = [0.5 0.5 0.5]; a.mainLine.LineWidth = 2;
a.patch.FaceColor = [0.5 0.5 0.5]; a.patch.EdgeColor = [0.5 0.5 0.5]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0.5 0.5 0.5]; a.edge(2).Color = [0.5 0.5 0.5];
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

hold on

% new
a = shadedErrorBar(time, mean(ga_pfc_2)*100, std(ga_pfc_2)*100 ...
    / sqrt(size(ga_pfc_2,1)));
a.mainLine.Color = [0 0 0]; a.mainLine.LineWidth = 2;
a.patch.FaceColor = [0 0 0]; a.patch.EdgeColor = [0 0 0]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0 0 0]; a.edge(2).Color = [0 0 0];
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

% make plot nice
xline(0, 'k--', 'linewidth', 0.5); yline(0, 'k--', 'linewidth', 0.5);
xlim([time(1) time(end)]); ylabel('PEV'); xlabel('Time to HLL [s]');
set(gca, 'fontsize', 13, 'linewidth', 3); box off;
xticks([-0.4 -0.2 0 0.2]);

set(gcf, 'Position',[718 461 290 197]);

%% Plot motor cortex

ga_motor_1 = table2array(tbl(41:51,2:end));
ga_motor_2 = table2array(tbl(55:70,2:end));
ga_motor_2(all(isnan(ga_motor_2),2),:) = [];

figure;
% old
a = shadedErrorBar(time, mean(ga_motor_1)*100, std(ga_motor_1)*100 ...
    / sqrt(size(ga_motor_1,1)));
a.mainLine.Color = [0.5 0.5 0.5]; a.mainLine.LineWidth = 2;
a.patch.FaceColor = [0.5 0.5 0.5]; a.patch.EdgeColor = [0.5 0.5 0.5]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0.5 0.5 0.5]; a.edge(2).Color = [0.5 0.5 0.5];
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

hold on

% new
a = shadedErrorBar(time, mean(ga_motor_2)*100, std(ga_motor_2)*100 ...
    / sqrt(size(ga_motor_2,1)));
a.mainLine.Color = [0 0 0]; a.mainLine.LineWidth = 2;
a.patch.FaceColor = [0 0 0]; a.patch.EdgeColor = [0 0 0]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0 0 0]; a.edge(2).Color = [0 0 0];
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

% make plot nice
xline(0, 'k--', 'linewidth', 0.5); yline(0, 'k--', 'linewidth', 0.5);
xlim([time(1) time(end)]); ylabel('PEV'); xlabel('Time to HLL [s]');
set(gca, 'fontsize', 13, 'linewidth', 3); box off;
xticks([-0.4 -0.2 0 0.2]);

set(gcf, 'Position',[718 461 290 197]);
