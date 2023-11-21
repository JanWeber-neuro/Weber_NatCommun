function plot_fig_S14_pca_shuffling(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S14');

figure;
h = histogram(tbl.PermDist);
h.FaceColor = 'k';
h.FaceAlpha = 1;

ylabel('Count');
xlabel('% explained variance');

box off
set(gca, 'fontsize', 13);

xline(tbl.ChanceLevel(1), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');

xline(tbl.ExpVar(tbl.Significant_PCs(9)), 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);

xline(tbl.ExpVar(tbl.Significant_PCs(9)+1), 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);

xlim([0.8 1.5]);

figure;
plot(cumsum(tbl.ExpVar(1:9)), 'k.', 'MarkerSize', 20);
hold on
plot(cumsum(tbl.ExpVar(1:9)), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

xlabel('Principal Component');
ylabel('Cumulative %EV');

box off
set(gca, 'fontsize', 13);
