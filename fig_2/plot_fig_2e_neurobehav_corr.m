function plot_fig_2e_neurobehav_corr(path)

%% Correlation RT & Peak Amplitude

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2e upper right');

cmap = cbrewer('seq', 'YlOrRd', length(tbl.RT), 'PCHIP');

figure;
h = scatter(tbl.RT, tbl.Peak_Amplitude, 300, cmap, 'Marker', '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
hold on
a = lsline;
a.LineWidth = 2;

ylabel('HFA_m_a_x');
xlabel('RT');
set(gca, 'linewidth', 1, 'fontsize', 13);
box off

set(gcf, 'Position', [484 351 228 156]);

clear tbl

%% Correlation RT & Peak Latency

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2e lower right');

figure;
h = scatter(tbl.RT, tbl.Peak_Latency, 300, cmap, 'Marker', '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
hold on
a = lsline;
a.LineWidth = 2;

ylabel('HFA_p_e_a_k_t_i_m_e');
xlabel('RT');
set(gca, 'linewidth', 1, 'fontsize', 13);
box off

set(gcf, 'Position', [484 351 228 156]);

