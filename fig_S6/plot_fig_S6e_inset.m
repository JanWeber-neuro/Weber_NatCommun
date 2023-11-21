function plot_fig_S6e_inset(path)

c2use = {'#008450', '#EFB700', '#B81D13'};

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S6e inset');

avg_zval_0  = tbl.Var1(2:end);
avg_zval_25 = tbl.Var2(2:end);
avg_zval_75 = tbl.Var3(2:end);

figure;

data2use = {avg_zval_0, avg_zval_25, avg_zval_75};

for ii = 1:3
        
    vertical_errorbars(ii, 0.1, mean(data2use{ii}), std(data2use{ii}) / sqrt(length(data2use{ii})), c2use{ii}, 2)    

    hold on
    
    scatter(ones(size(data2use{ii}))*ii, data2use{ii}, 'filled', 'MarkerFaceColor', ...
        c2use{ii}, 'jitter','on', 'jitterAmount',0.04);

end
    
xticks(1:3);
xticklabels({'0%', '25%', '75%'}); xtickangle(45);
ylabel('rho');

set(gca, 'fontsize', 13, 'linewidth', 2);

box off

set(gcf, 'Position', [511 383 258 184]);
