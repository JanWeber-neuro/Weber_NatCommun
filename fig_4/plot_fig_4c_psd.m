function plot_fig_4c_psd(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_4c');

freq  = table2array(tbl(1,2:end));

ga_0  = table2array(tbl(3:18,2:end));
ga_25 = table2array(tbl(24:39,2:end));
ga_75 = table2array(tbl(43:58,2:end));

figure;

color = {'#008450', '#EFB700', '#B81D13'};

plot(freq, mean(ga_0), 'Color', color{1}, 'LineWidth', 4);
hold on
plot(freq, mean(ga_25), 'Color', color{2}, 'LineWidth', 4);
hold on
plot(freq, mean(ga_75), 'Color', color{3}, 'LineWidth', 4);

% make plot look nice
xlabel('Frequency [Hz]');
ylabel('Power')
set(gca, 'linewidth', 3, 'fontsize', 13, 'xlim', [1 15]);
box off
set(gcf, 'Position', [489 361 270 201]);
