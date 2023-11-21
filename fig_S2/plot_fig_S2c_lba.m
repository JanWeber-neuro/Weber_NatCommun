function plot_fig_S2c_lba(path)

hexcolors = {'#EFB700', '#B81D13'};

colors2use = NaN(2,3);
colorBins  = [2:3; 4:5; 6:7];

for iHex = 1:numel(hexcolors)
    for iBin = 1:size(colorBins)
        colors2use(iHex,iBin) = hex2dec(hexcolors{iHex}(colorBins(iBin,1):colorBins(iBin,end)));
    end
end

% scale by dividing through 256 (Matlab conversion)
colors2use = colors2use/ 256;

%% Plot shift in starting point

tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_S2c');

figure;

% plot data
jw_raincloud_pairplot([tbl.Diff25_Vs0_, tbl.Diff75_Vs0_], {colors2use(1,:), colors2use(2,:)});

set(gca, 'ytick', []);
xlabel('starting point bias');

h = gca; h.YAxis.Visible = 'off';

set(gca, 'fontsize', 13);

box off
