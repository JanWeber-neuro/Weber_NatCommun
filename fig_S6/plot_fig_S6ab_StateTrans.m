function plot_fig_S6ab_StateTrans(path)

%% Plot state transition PFC

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6a right');
fsample = 512;

time = table2array(tbl(1,2:end));
ga_pfc_0  = table2array(tbl(3:19,2:end));
ga_pfc_25 = table2array(tbl(23:39,2:end));
ga_pfc_75 = table2array(tbl(43:59,2:end));

figure;
plot(time, smoothdata(mean(ga_pfc_0), 'gaussian', round(fsample*0.01)), 'Color', '#008450', 'LineWidth', 4);
hold on
plot(time, smoothdata(mean(ga_pfc_25), 'gaussian', round(fsample*0.01)), 'Color', '#EFB700', 'LineWidth', 4);
hold on
plot(time, smoothdata(mean(ga_pfc_75), 'gaussian', round(fsample*0.01)), 'Color', '#B81D13', 'LineWidth', 4);

ylabel('State Change');

a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)], ...
    'Ylim', [a(1) a(2)+0.001]); box off; hold on;
set(gcf, 'Position', [552 342 205 151]);

%% Plot state transition motor cortex

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6b right');
fsample = 512;

time = table2array(tbl(1,2:end));
ga_motor_0  = table2array(tbl(3:16,2:end));
ga_motor_25 = table2array(tbl(21:34,2:end));
ga_motor_75 = table2array(tbl(38:51,2:end));

figure;
plot(time, smoothdata(mean(ga_motor_0), 'gaussian', round(fsample*0.01)), 'Color', '#008450', 'LineWidth', 4);
hold on
plot(time, smoothdata(mean(ga_motor_25), 'gaussian', round(fsample*0.01)), 'Color', '#EFB700', 'LineWidth', 4);
hold on
plot(time, smoothdata(mean(ga_motor_75), 'gaussian', round(fsample*0.01)), 'Color', '#B81D13', 'LineWidth', 4);

ylabel('State Change');

a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)], ...
    'Ylim', [a(1) a(2)+0.001]); box off; hold on;
set(gcf, 'Position', [552 342 205 151]);
