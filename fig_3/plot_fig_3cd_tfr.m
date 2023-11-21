function plot_fig_3cd_tfr(path)

% color for conditons
color = {'#008450', '#EFB700', '#B81D13'};

%% Plot PFC

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_3c_right');

freq  = table2array(tbl(1,:));
ga_0  = table2array(tbl(2,:));
ga_25 = table2array(tbl(3,:));
ga_75 = table2array(tbl(4,:));

figure;

plot(freq, ga_0, 'Color', color{1}, 'LineWidth', 4);
hold on
plot(freq, ga_25, 'Color', color{2}, 'LineWidth', 4);
hold on
plot(freq, ga_75, 'Color', color{3}, 'LineWidth', 4);

ylabel('Power [z]');
xlabel('Frequency [Hz]');
set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [freq(1) freq(end)]); 
box off; 

set(gca, 'xscale', 'log', 'xtick', [4 8 16 32 64 128], 'ylim', [-3.2 2.2]);

set(gcf, 'Position', [479 351 199 184]);

%% Plot motor cortex

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_3d_right');

freq  = table2array(tbl(1,:));
ga_0  = table2array(tbl(2,:));
ga_25 = table2array(tbl(3,:));
ga_75 = table2array(tbl(4,:));

figure;

plot(freq, ga_0, 'Color', color{1}, 'LineWidth', 4);
hold on
plot(freq, ga_25, 'Color', color{2}, 'LineWidth', 4);
hold on
plot(freq, ga_75, 'Color', color{3}, 'LineWidth', 4);

ylabel('Power [z]');
xlabel('Frequency [Hz]');
set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [freq(1) freq(end)]); 
box off; 

set(gca, 'xscale', 'log', 'xtick', [4 8 16 32 64 128], 'ylim', [-3 4.2]);

set(gcf, 'Position', [479 351 199 184]);
