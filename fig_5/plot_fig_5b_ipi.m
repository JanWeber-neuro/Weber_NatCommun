function plot_fig_5b_ipi(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_5b');

figure;
a = histfit(tbl.GA_frontal);
a(1).FaceColor = 'r';
a(1).EdgeColor = a(1).FaceColor;
a(1).FaceAlpha = 1;
a(1).EdgeAlpha = 1;
a(1).LineWidth = 1;
a(2).Color     = 'k';
a(2).Color(4) = 0;

xlabel('IPI [Hz]');
ylabel('Count');

set(gca, 'linew', 0.5, 'fontsize', 13); box off;
set(gcf, 'Position', [524 394 230 193]);
            
figure;
a = histfit(tbl.GA_motor);
a(1).FaceColor = 'b';
a(1).EdgeColor = a(1).FaceColor;
a(1).FaceAlpha = 1;
a(1).EdgeAlpha = 1;
a(1).LineWidth = 1;
a(2).Color     = 'k';
a(2).Color(4) = 0;

xlabel('IPI [Hz]');
ylabel('Count');

set(gca, 'linew', 0.5, 'fontsize', 13); box off;
set(gcf, 'Position', [524 394 230 193]);
                 