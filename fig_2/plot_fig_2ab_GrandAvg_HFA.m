function plot_fig_2ab_GrandAvg_HFA(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2a lower left');

%% Plot data for frontal cortex

time = table2array(tbl(1,2:end));

ga_0_pfc  = table2array(tbl(3:20,2:end));
ga_25_pfc = table2array(tbl(25:42,2:end));
ga_75_pfc = table2array(tbl(47:64,2:end));

figure;

h = shadedErrorBar(time, nanmean(ga_0_pfc), nanstd(ga_0_pfc) / ...
    sqrt(size(ga_0_pfc,1)));
h.mainLine.Color = '#008450'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#008450'; h.patch.EdgeColor = '#008450'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#008450'; h.edge(2).Color = '#008450';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

hold on

h = shadedErrorBar(time, nanmean(ga_25_pfc), nanstd(ga_25_pfc) / ...
    sqrt(size(ga_25_pfc,1)));
h.mainLine.Color = '#EFB700'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#EFB700'; h.patch.EdgeColor = '#EFB700'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#EFB700'; h.edge(2).Color = '#EFB700';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

hold on

h = shadedErrorBar(time, nanmean(ga_75_pfc), nanstd(ga_75_pfc) / ...
    sqrt(size(ga_75_pfc,1)));
h.mainLine.Color = '#B81D13'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#B81D13'; h.patch.EdgeColor = '#B81D13'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#B81D13'; h.edge(2).Color = '#B81D13';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

ylabel('HFA [z]');
xlabel('Time to HLL');

a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)], ...
    'Ylim', [a(1)-0.1 a(2)+0.1]); box off; hold on;
set(gcf, 'Position', [491 320 201 175]);

clear tbl

%% Plot data for motor cortex

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2b lower left');

time = table2array(tbl(1,2:end));

ga_0_motor  = table2array(tbl(3:20,2:end));
ga_25_motor = table2array(tbl(25:42,2:end));
ga_75_motor = table2array(tbl(47:64,2:end));

figure;

h = shadedErrorBar(time, nanmean(ga_0_motor), nanstd(ga_0_motor) / ...
    sqrt(size(ga_0_motor,1)));
h.mainLine.Color = '#008450'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#008450'; h.patch.EdgeColor = '#008450'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#008450'; h.edge(2).Color = '#008450';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

hold on

h = shadedErrorBar(time, nanmean(ga_25_motor), nanstd(ga_25_motor) / ...
    sqrt(size(ga_25_motor,1)));
h.mainLine.Color = '#EFB700'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#EFB700'; h.patch.EdgeColor = '#EFB700'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#EFB700'; h.edge(2).Color = '#EFB700';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

hold on

h = shadedErrorBar(time, nanmean(ga_75_motor), nanstd(ga_75_motor) / ...
    sqrt(size(ga_75_motor,1)));
h.mainLine.Color = '#B81D13'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#B81D13'; h.patch.EdgeColor = '#B81D13'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#B81D13'; h.edge(2).Color = '#B81D13';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

ylabel('HFA [z]');
xlabel('Time to HLL');

a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)], ...
    'Ylim', [a(1)-0.1 a(2)+0.1]); box off; hold on;
set(gcf, 'Position', [491 320 201 175]);

