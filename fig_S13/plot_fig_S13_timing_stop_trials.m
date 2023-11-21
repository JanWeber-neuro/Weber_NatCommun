function plot_fig_S13_timing_stop_trials(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S13');

figure;
a = histogram(tbl.SSD, 'Normalization', 'probability');
a.FaceColor = 'k';
a.FaceAlpha = 1;

hold on

[f,xi] = ksdensity(tbl.SSD);

plot(xi, f*0.012, 'r', 'LineWidth', 2)

ylabel('probability');
xlabel('Time [s]');
box off

set(gca, 'fontsize', 13);
