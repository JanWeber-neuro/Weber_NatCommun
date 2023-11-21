function plot_fig_S6f_ExplVar(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6f');
nPCs = 5;
GA_frontal = table2array(tbl(1:17,:));
GA_motor = table2array(tbl(21:33,:));

figure;
a = shadedErrorBar(1:nPCs, mean(GA_frontal), std(GA_frontal) / sqrt(size(GA_frontal,1) - sum(isnan(GA_frontal(:,1)))));
a.mainLine.Color = 'r'; a.mainLine.LineWidth = 2; 
a.patch.FaceColor = 'r'; a.patch.EdgeColor = 'r'; a.patch.FaceAlpha = 0.5;
a.mainLine.Marker = 'o'; a.mainLine.MarkerFaceColor = 'r'; a.mainLine.MarkerEdgeColor = 'r';

hold on

a = shadedErrorBar(1:nPCs, mean(GA_motor), std(GA_motor) / sqrt(size(GA_motor,1) - sum(isnan(GA_motor(:,1)))));
a.mainLine.Color = 'b'; a.mainLine.LineWidth = 2;
a.patch.FaceColor = 'b'; a.patch.EdgeColor = 'b'; a.patch.FaceAlpha = 0.5;
a.mainLine.Marker = 'o'; a.mainLine.MarkerFaceColor = 'b'; a.mainLine.MarkerEdgeColor = 'b';

set(gca, 'linewidth', 3, 'fontsize', 13);
ylabel('Cumulative Sum');
xlim([1 nPCs]);
xlabel('PCs');
xticks(1:nPCs);
box off

set(gcf, 'Position', [515 450 243 173]);
