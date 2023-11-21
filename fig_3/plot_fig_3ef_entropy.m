function plot_fig_3ef_entropy(path)

%% Plot sample entropy PFC

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_3e right');

time = table2array(tbl(1,2:end));
fsample = 512;

ga_pfc_0  = table2array(tbl(3:18,2:end));
ga_pfc_25 = table2array(tbl(22:37,2:end));
ga_pfc_75 = table2array(tbl(41:56,2:end));

% plot data
figure;
smoothwin_plot = 0.002;
plot(time, smoothdata(nanmean(ga_pfc_0), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#008450', 'LineWidth', 4);
hold on
plot(time, smoothdata(nanmean(ga_pfc_25), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#EFB700', 'LineWidth', 4);
hold on
plot(time, smoothdata(nanmean(ga_pfc_75), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#B81D13', 'LineWidth', 4);

% xlabel('Time to HLL [s]');
ylabel('Sample entropy');

% plot significant timepoints
a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)], ...
    'Ylim', [a(1)-0.001 a(2)+0.001]); box off; hold on;

set(gcf, 'Position', [728 436 304 210]);

clear tbl

%% Plot sample entropy motor cortex

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_3f right');

time = table2array(tbl(1,2:end));
fsample = 512;

ga_motor_0  = table2array(tbl(2:12,2:end));
ga_motor_25 = table2array(tbl(16:26,2:end));
ga_motor_75 = table2array(tbl(30:40,2:end));

% plot data
figure;
smoothwin_plot = 0.002;
plot(time, smoothdata(nanmean(ga_motor_0), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#008450', 'LineWidth', 4);
hold on
plot(time, smoothdata(nanmean(ga_motor_25), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#EFB700', 'LineWidth', 4);
hold on
plot(time, smoothdata(nanmean(ga_motor_75), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#B81D13', 'LineWidth', 4);

% xlabel('Time to HLL [s]');
ylabel('Sample entropy');

% plot significant timepoints
a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)], ...
    'Ylim', [a(1)-0.001 a(2)+0.001]); box off; hold on;

set(gcf, 'Position', [728 436 304 210]);
