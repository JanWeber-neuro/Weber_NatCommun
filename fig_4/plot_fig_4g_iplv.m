function plot_fig_4g_iplv(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_4g');

freq = table2array(tbl(1,2:end));

ga_0  = table2array(tbl(3:16,2:end));
ga_25 = table2array(tbl(20:33,2:end));
ga_75 = table2array(tbl(37:50,2:end));

figure;
h = shadedErrorBar(freq, nanmean(ga_0), ...
    nanstd(ga_0)/sqrt(size(ga_0,1)));
h.mainLine.Color = '#008450'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#008450'; h.patch.EdgeColor = '#008450'; h.patch.FaceAlpha = 0.5;

hold on

h = shadedErrorBar(freq, nanmean(ga_25), ...
    nanstd(ga_25)/sqrt(size(ga_25,1)));
h.mainLine.Color = '#EFB700'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#EFB700'; h.patch.EdgeColor = '#EFB700'; h.patch.FaceAlpha = 0.5;

hold on

h = shadedErrorBar(freq, nanmean(ga_75), ...
    nanstd(ga_75)/sqrt(size(ga_75,1)));
h.mainLine.Color = '#B81D13'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#B81D13'; h.patch.EdgeColor = '#B81D13'; h.patch.FaceAlpha = 0.5;

ylabel('normalized iPLV [z]'); xlabel('Frequency [Hz]'); xlim([freq(1) freq(end)]);
set(gca, 'linewidth', 3, 'fontsize', 13, 'xtick', [4 8 16 32], 'xscale', 'log');

set(gcf, 'Position', [760 459 261 229]);
