function plot_fig_5d_connectivity(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_5d_upper_panel');

figure;
a = histfit(tbl.Var1);
a(1).FaceColor = 'k';
a(1).EdgeColor = a(1).FaceColor;
a(1).FaceAlpha = 1;
a(1).EdgeAlpha = 1;
a(1).LineWidth = 3;

xlabel('rho [zval]');
ylabel('Count');
xline(0, 'LineStyle', '--', 'LineWidth', 1.5);

set(gca, 'linew', 2, 'fontsize', 13); box off;
set(gcf, 'Position', [547 457 227 197]);

clear tbl

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_5d_lower_panel');

figure;
a = jw_scatterplot({tbl.Connectivity}, 0, 0, {'k'});
a.MarkerFaceColor = 'k';
a.MarkerEdgeColor = 'k';
a.MarkerEdgeAlpha = 0.5;
a.SizeData        = 80;

xlim([0.6 1.5]); ylim([-2 2]); yline(0, 'LineWidth', 2);
xlabel('Pow Corr.'); xticks('');
ylabel('rho [z]');
set(gca, 'linew', 2, 'fontsize', 13); box off;
set(gcf, 'Position', [547 457 227 197]);
