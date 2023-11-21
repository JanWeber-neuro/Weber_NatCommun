function plot_fig_2cd_singleTrials_HFA(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2cd');
time = linspace(-0.5, 0.3, 411);

%% Plot single trials PFC

figure;

a = plot(time, tbl.Var1(2:end), 'LineWidth', 3);
a.Color = '#008450'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

figure;

a = plot(time, tbl.Var2(2:end), 'LineWidth', 3);
a.Color = '#EFB700'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

figure;

a = plot(time, tbl.Var3(2:end), 'LineWidth', 3);
a.Color = '#B81D13'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

%% Plot single trials PFC

figure;

a = plot(time, tbl.Var5(2:end), 'LineWidth', 3);
a.Color = '#008450'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

figure;

a = plot(time, tbl.Var6(2:end), 'LineWidth', 3);
a.Color = '#EFB700'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

figure;

a = plot(time, tbl.Var7(2:end), 'LineWidth', 3);
a.Color = '#B81D13'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

