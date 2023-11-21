%% Figure S8

%%

clear
ft_defaults

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_data = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Round_3/Data';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_S6/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

% add colormaps
addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/colormaps/;

%% - LOAD DATA - %%

% load HFA data locked to the LL
load(fullfile(path_data, 'hfa_hll_locked.mat'));

%% - SUBSELECT TO PFC AND PERFORM PCA - %%

% load ROI info
load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/Electrode_Information_ROI/elec_info_ROI.mat;

% preallocate memory
pcadata = cell(numel(allData),1);
pcainfo = cell(numel(allData),1);

for Isub = 1:numel(pcadata)

    fprintf('%d/%d\n', Isub, numel(pcadata));

    if isfield(ROIs(Isub).frontal, 'elecidx')

        % only use PFC subset
        cfg = [];
        cfg.channel   = ROIs(Isub).frontal.elecidx;
        pcadata{Isub} = ft_selectdata(cfg, allData{Isub});

        % perform PCA on data
        [pcadata{Isub}.trial, pcainfo{Isub}] = rh_pca(pcadata{Isub}.trial);

    end

end%Isub

% remove empty cells
pcadata = pcadata(~cellfun(@isempty, pcadata));
pcainfo = pcainfo(~cellfun(@isempty, pcainfo));

%% - PERFORM PCA ON DATA - %%

% numer of PCs to consider
nPCs = 5;

dat2plot = cell(numel(pcadata),1);

for Isub = 1:numel(dat2plot)

    if numel(pcadata{Isub}.label) >= nPCs

        % get only relevant PCs
        cfg = [];
        cfg.channel    = 1:nPCs;
        dat2plot{Isub} = ft_timelockanalysis(cfg, pcadata{Isub});
        dat2plot{Isub}.label = {'PC1', 'PC2', 'PC3', 'PC4', 'PC5'};
    
        dat2plot{Isub}.cfg = [];

    end

end

dat2plot = dat2plot(~cellfun(@isempty, dat2plot));

%% - COMPUTE THE VARIANCE EXPLAINED BY THE FIRST 5 PCs - %%

ev = NaN(numel(pcadata),1);

for Isub = 1:size(ev,1)

    if numel(pcadata{Isub}.label) >= nPCs

        tmp      = cumsum(pcainfo{Isub}.explained(1:nPCs));
        ev(Isub) = tmp(end); clear tmp

    end

end

fprintf('Explained Variance [%%] for PC1 to PC%d = %.2f%% ± %.2f%% [mean ± SD]\n', nPCs, nanmean(ev), nanstd(ev));

%% - PERFORM GRAND AVERAGE - %%

close all;

cfg = [];
cfg.keepindividual = 'yes';

ga = ft_timelockgrandaverage(cfg, dat2plot{:});

cmap = viridis(size(ga.individual,2));

figure;
for jj = 1:size(ga.individual,2)

    subplot(2,3,jj)
%     subplot(2,1,1);

    jw_shadedErrorBar(ga.time, squeeze(mean(ga.individual(:,jj,:),1)), ...
        std(squeeze(ga.individual(:,jj,:))), size(ga.individual,1), cmap(jj,:));

    xlim([ga.time(1) ga.time(end)]);
    set(gca, 'fontsize', 13, 'linewidth', 1, 'box', 'off');
    set(gcf, 'Position', [494 299 368 282]);

%     subplot(2,1,2);
%     plot(ga.time, squeeze(ga.individual(:,jj,:))', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
%     
    xlabel('Time to HLL [s]');    
    ylabel(sprintf('Population Act. [PC%d]', jj));
    xline(0, '--', 'LineWidth', 1.5)
    xlim([ga.time(1) ga.time(end)]);
    set(gca, 'fontsize', 13, 'linewidth', 1, 'box', 'off');
    
end

set(gcf, 'Position', [381 356 587 307]);

cd(path_out);
saveas(gcf, 'pca_1to5_pfc_locked_HLL.pdf')

%% - SAVE SOURCE DATA FOR PCA 1-5 - %%

% go into results folder
mkdir(path_source, 'Fig_S6');
cd(fullfile(path_source, 'Fig_S6'));

for iPC = 1:5

    sourcedata = [];
    
    time = ga.time;
    
    % - save source data for PC1 - %
    sourcedata = squeeze(ga.individual(:,iPC,:));
    
    tmp = NaN(size(sourcedata,1)+1, size(sourcedata,2)+1);
    tmp(1,2:end)      = time;
    tmp(2:end,1)      = (1:17)';
    tmp(2:end,2:end)  = sourcedata;
    
    table_data = array2table(tmp); clear tmp
    
    % write data into excelsheet
    filename = ['fig_S6_frontal_PC' num2str(iPC) '.xlsx'];
    writetable(table_data,filename, 'Sheet', 1);

    clear table_data

end

%% - PLOT THE CHANGE IN ACTIVATION OVER TIME - %%

cfg = [];
cfg.keepindividual = 'yes';

ga = ft_timelockgrandaverage(cfg, dat2plot{:});

cmap = viridis(size(ga.individual,2));

figure;
for jj = 1:size(ga.individual,2)

    subplot(2,3,jj)

%     subplot(2,1,1);
    % get the derivative
    tmp = diff(smoothdata(squeeze(ga.individual(:,jj,:)),2,'gaussian'),1,2);

    jw_shadedErrorBar(ga.time(2:end), squeeze(mean(tmp,1)), ...
        std(squeeze(tmp)), size(tmp,1), cmap(jj,:));

    xlim([ga.time(2) ga.time(end)]);
    set(gca, 'fontsize', 13, 'linewidth', 1, 'box', 'off');
    set(gcf, 'Position', [494 299 368 282]);

%     subplot(2,1,2);
%     plot(ga.time(2:end), tmp', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
% 
    xlabel('Time to HLL [s]');    
    ylabel(sprintf('change in activity [PC%d]', jj));
    xlim([ga.time(2) ga.time(end)]);
    set(gca, 'fontsize', 13, 'linewidth', 1, 'box', 'off');
    
end

set(gcf, 'Position', [381 356 587 307]);

cd(path_out);
saveas(gcf, 'derivative_pca_1to5_pfc_locked_HLL.pdf')

%% - PLOT SINGLE SUBJECT PC1 - %%

figure;
% plot single subjects 
cmap = viridis(size(ga.individual,1));

for Isub = 1:size(ga.individual,1)

    subplot(3,6,Isub)

    plot(ga.time, squeeze(ga.individual(Isub,1,:)), 'Color', cmap(Isub,:), 'LineWidth', 2);

    if Isub == 1
        xlabel('Time to HLL [s]');    
        ylabel('PC1');
    end
    xlim([ga.time(1) ga.time(end)]);
    set(gca, 'fontsize', 13, 'linewidth', 1, 'box', 'off');
    
end

sgtitle('PFC Population Responses Single Subject', 'fontsize', 13);

cd(path_out);
saveas(gcf, 'single_subjects_pc1_pfc_locked_HLL.pdf')
