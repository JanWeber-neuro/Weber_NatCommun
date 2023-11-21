function plot_fig_4e_ipi(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_4e');

figure;
[f,xi] = ksdensity(tbl.LH0);
a = plot(xi, f,'LineWidth', 3); a.Color(4) = 0.8; a.Color = '#008450';
[~, idx] = max(f); hold on
a = xline(xi(idx), 'LineWidth', 1, 'LineStyle', '--');
a.Color = '#008450';
hold on
[f,xi] = ksdensity(tbl.LH25);
a = plot(xi, f,'LineWidth', 3); a.Color(4) = 0.8; a.Color = '#EFB700';
[~, idx] = max(f); hold on
a = xline(xi(idx), 'LineWidth', 1, 'LineStyle', '--');
a.Color = '#EFB700';
hold on
[f,xi] = ksdensity(tbl.LH75);
a = plot(xi, f,'LineWidth', 3); a.Color(4) = 0.8; a.Color = '#B81D13';
[~, idx] = max(f); hold on
a = xline(xi(idx), 'LineWidth', 1, 'LineStyle', '--');
a.Color = '#B81D13';

set(gca, 'fontsize', 13, 'linewidth', 3);
xlabel('IPI [Hz]');
ylabel('Probability');
box off
set(gcf, 'Position', [510 356 291 230]);
