function plot_fig_S8_PC1_to_PC5(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S8');

time = table2array(tbl(1,2:end));

ga_pc1 = table2array(tbl(2:18,2:end));
ga_pc2 = table2array(tbl(22:38,2:end));
ga_pc3 = table2array(tbl(42:58,2:end));
ga_pc4 = table2array(tbl(62:78,2:end));
ga_pc5 = table2array(tbl(82:98,2:end));

ga_all = NaN(size(ga_pc1,1),5,size(ga_pc1,2));
ga_all(:,1,:) = ga_pc1;
ga_all(:,2,:) = ga_pc2;
ga_all(:,3,:) = ga_pc3;
ga_all(:,4,:) = ga_pc4;
ga_all(:,5,:) = ga_pc5;

cmap = viridis(size(ga_all,2));

figure;
for jj = 1:size(ga_all,2)

    subplot(2,3,jj)

    jw_shadedErrorBar(time, squeeze(mean(ga_all(:,jj,:),1)), ...
        std(squeeze(ga_all(:,jj,:))), size(ga_all,1), cmap(jj,:));

    xlim([time(1) time(end)]);
    set(gca, 'fontsize', 13, 'linewidth', 1, 'box', 'off');

    xlabel('Time to HLL [s]');    
    ylabel(sprintf('Population Act. [PC%d]', jj));
    xline(0, '--', 'LineWidth', 1.5)
    xlim([time(1) time(end)]);
    set(gca, 'fontsize', 13, 'linewidth', 1, 'box', 'off');
    
end

set(gcf, 'Position', [381 356 587 307]);

%% - PLOT THE CHANGE IN ACTIVATION OVER TIME - %%

figure;
for jj = 1:size(ga_all,2)

    subplot(2,3,jj)

    tmp = diff(smoothdata(squeeze(ga_all(:,jj,:)),2,'gaussian'),1,2);

    jw_shadedErrorBar(time(2:end), squeeze(mean(tmp,1)), ...
        std(squeeze(tmp)), size(tmp,1), cmap(jj,:));

    xlabel('Time to HLL [s]');    
    ylabel(sprintf('change in activity [PC%d]', jj));
    xlim([time(2) time(end)]);
    set(gca, 'fontsize', 13, 'linewidth', 1, 'box', 'off');
    
end

set(gcf, 'Position', [381 356 587 307]);
