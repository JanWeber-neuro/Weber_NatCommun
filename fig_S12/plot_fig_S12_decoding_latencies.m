function plot_fig_S12_decoding_latencies(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S12');

cmap = viridis(100);

data2plot = {tbl.PFCContext, tbl.PFCAction', ...
    tbl.MotorAction'}';

% color 2 use for plotting
cb = {cmap(1,:), cmap(20,:), cmap(80,:)};

figure; 
h = rm_raincloud(data2plot, [1 0 0], 0, 'ks', 0.15);

% admin
for Icond = 1:numel(data2plot)
    
    % change density plot color
    h.p{Icond,1}.FaceColor        = cb{Icond};
    h.p{Icond,1}.EdgeColor        = cb{Icond};
    h.p{Icond,1}.FaceAlpha        = 0.8;

    % change color for indiv. dots
    h.s{Icond,1}.MarkerFaceColor  = cb{Icond};
    h.s{Icond,1}.MarkerFaceAlpha  = 0.8;
    h.s{Icond,1}.SizeData         = 40;

    % change mean dots
    h.m(Icond).MarkerFaceColor    = cb{Icond};
    h.m(Icond).SizeData           = 80;
    
end

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];
    
% general figure settings
yticks('');
xline(0, 'LineStyle', '--');
set(gca, 'fontsize', 13);
box off
view([0 90]);

set(gcf, 'Position', [2666 403 299 182]);
