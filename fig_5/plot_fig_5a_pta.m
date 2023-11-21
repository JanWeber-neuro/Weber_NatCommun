function plot_fig_5a_pta(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_5a');

time = table2array(tbl(1,2:end));

%% Plot peak-triggered average PFC

pta_pfc = table2array(tbl(3:20,2:end));
pta_pfc(all(isnan(pta_pfc),2),:) = [];

% plot
figure;
a = shadedErrorBar(time, mean(pta_pfc), ...
    std(pta_pfc) / sqrt(size(pta_pfc,1)));
a.mainLine.Color = 'r'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'r'; a.patch.EdgeColor = 'r'; a.patch.FaceAlpha = 0.5;

ylabel('Activity');
xlabel('Time from Peak');

set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)]); 
box off

set(gcf, 'Position', [488 381 331 252]);

hold on

% sine fit
sineparams = sineFit(time,mean(pta_pfc), 0);
plot(time,sineparams(1)+sineparams(2)*sin(2*pi*sineparams(3)*time+sineparams(4)), 'k', 'LineWidth', 1);

xlim([-0.3 0.3]);

%% Plot peak-triggered average motor cortex

pta_motor = table2array(tbl(25:42,2:end));
pta_motor(all(isnan(pta_motor),2),:) = [];

% plot
figure;
a = shadedErrorBar(time, mean(pta_motor), ...
    std(pta_motor) / sqrt(size(pta_motor,1)));
a.mainLine.Color = 'b'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'b'; a.patch.EdgeColor =' b'; a.patch.FaceAlpha = 0.5;

ylabel('Activity');
xlabel('Time from Peak');

set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)]); 
box off

set(gcf, 'Position', [488 381 331 252]);

hold on

% sine fit
sineparams = sineFit(time,mean(pta_motor), 0);
plot(time,sineparams(1)+sineparams(2)*sin(2*pi*sineparams(3)*time+sineparams(4)), 'k', 'LineWidth', 1);

xlim([-0.3 0.3]);
