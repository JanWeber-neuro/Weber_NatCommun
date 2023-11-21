function plot_fig_S1_RTDist(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S1');

% color 2 use for plotting
cb = {'#008450', '#EFB700', '#B81D13'};

figure;
histogram(tbl.Var1(2:find(isnan(tbl.Var1),1,'first')), 'Normalization', 'pdf', 'FaceColor', cb{1}, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);
[f,xi] = ksdensity(tbl.Var1(2:find(isnan(tbl.Var1),1,'first')));
hold on
plot(xi,f,'Color', cb{1}, 'LineWidth', 1.5)
hold on
histogram(tbl.Var2(2:end), 'Normalization', 'pdf', 'FaceColor', cb{2}, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);
[f,xi] = ksdensity(tbl.Var2);
plot(xi,f,'Color', cb{2}, 'LineWidth', 1.5)
hold on 
histogram(tbl.Var3(2:find(isnan(tbl.Var3),1,'first')), 'Normalization', 'pdf', 'FaceColor', cb{3}, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);
[f,xi] = ksdensity(tbl.Var3(2:find(isnan(tbl.Var3),1,'first')));
plot(xi,f,'Color', cb{3}, 'LineWidth', 1.5)

xlim([0 0.2]);

xlabel('Movement Onset');
xticks([0 0.1 0.2]);
xticklabels({'HLL', '0.1', '0.2'});
ylabel('PDF');

box off
set(gca, 'fontsize', 13);

set(gcf, 'Position', [505 452 265 247]);
