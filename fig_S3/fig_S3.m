%% Figure S3
% this script rungs additional analysis required for the plots of figure
% S3.

%%

clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_S2/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

path_elec   = '/Volumes/IMPECOG/elecs_and_info';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% Load Data

% preallocate memory
allData_old = cell(numel(subjects),1); % old data without accounting for accuracy
allData_new = cell(numel(subjects),1); % new data after accounting for accuracy

for Isub = 1:numel(subjects)
        
    % get indicated subject ID

    subID = subjects{Isub};

    % - [LOAD "OLD" DATA] -%
    load(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data', [subID '_activedata_ROI_pca.mat']));
                        
    data = roidata; clear roidata;
        
    allData_old{Isub} = data; clear data
        
    % - [LOAD "NEW" DATA] -%
    load(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data', [subID '_activedata_ROI_pca_w_main_effect_accuracy.mat']));
                        
    data = roidata; clear roidata;
        
    allData_new{Isub} = data; clear data

end

%% Extract signals within ROIs and average for each subject to avoid bias due to different n electrodes

% name for regions of interest within structure 
rois = {'motor', 'motor'};

% load example data to store time 
load /Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/OSL13/Data/HFA/PCA/OSL13_TaskActive_PCAexpvar_w_main_effect_accuracy.mat

npts_hll  = pcadata.hll.pca.time;
npts_br   = pcadata.br.pca.time;

% initialize storage for grand average data
ga_old = struct();
ga_new = struct();

for Isub = 1:numel(subjects)

    %% HLL PFC "OLD"
    
    % get data for motor active/inactive electrodes
    if ~isempty(allData_old{Isub}.hll.frontal.active)

        % get the mean over active electrodes
        ga_old.hll.frontal.active(Isub,:) = nanmean(allData_old{Isub}.hll.frontal.active.pca.trial,1);

    else

        ga_old.hll.frontal.active(Isub,:) = NaN(1, length(npts_hll));

    end
            
    if ~isempty(allData_old{Isub}.hll.frontal.inactive)

        % get the mean over active electrodes
        ga_old.hll.frontal.inactive(Isub,:) = nanmean(allData_old{Isub}.hll.frontal.inactive.pca.trial,1);

    else

        ga_old.hll.frontal.inactive(Isub,:) = NaN(1, length(npts_hll));

    end

    %% BR PFC "OLD"
    
    % get data for motor active/inactive electrodes
    if ~isempty(allData_old{Isub}.br.frontal.active)

        % get the mean over active electrodes
        ga_old.br.frontal.active(Isub,:) = nanmean(allData_old{Isub}.br.frontal.active.pca.trial,1);

    else

        ga_old.br.frontal.active(Isub,:) = NaN(1, length(npts_br));

    end
            
    if ~isempty(allData_old{Isub}.br.frontal.inactive)

        % get the mean over active electrodes
        ga_old.br.frontal.inactive(Isub,:) = nanmean(allData_old{Isub}.br.frontal.inactive.pca.trial,1);

    else

        ga_old.br.frontal.inactive(Isub,:) = NaN(1, length(npts_br));

    end

    %% HLL PFC "NEW"
    
    % get data for frontal active/inactive electrodes
    if ~isempty(allData_new{Isub}.hll.frontal.active)

        % get the mean over active electrodes
        ga_new.hll.frontal.active(Isub,:) = nanmean(allData_new{Isub}.hll.frontal.active.pca.trial,1);

    else

        ga_new.hll.frontal.active(Isub,:) = NaN(1, length(npts_hll));

    end
            
    if ~isempty(allData_new{Isub}.hll.frontal.inactive)

        % get the mean over active electrodes
        ga_new.hll.frontal.inactive(Isub,:) = nanmean(allData_new{Isub}.hll.frontal.inactive.pca.trial,1);

    else

        ga_new.hll.frontal.inactive(Isub,:) = NaN(1, length(npts_hll));

    end

    %% BR PFC "NEW"
    
    % get data for frontal active/inactive electrodes
    if ~isempty(allData_new{Isub}.br.frontal.active)

        % get the mean over active electrodes
        ga_new.br.frontal.active(Isub,:) = nanmean(allData_new{Isub}.br.frontal.active.pca.trial,1);

    else

        ga_new.br.frontal.active(Isub,:) = NaN(1, length(npts_br));

    end
            
    if ~isempty(allData_new{Isub}.br.frontal.inactive)

        % get the mean over active electrodes
        ga_new.br.frontal.inactive(Isub,:) = nanmean(allData_new{Isub}.br.frontal.inactive.pca.trial,1);

    else

        ga_new.br.frontal.inactive(Isub,:) = NaN(1, length(npts_br));

    end

    %% HLL M1 "OLD"
    
    % get data for motor active/inactive electrodes
    
    if ~isempty(allData_old{Isub}.hll.motor.active)

        % get the mean over active electrodes
        ga_old.hll.motor.active(Isub,:) = nanmean(allData_old{Isub}.hll.motor.active.pca.trial,1);

    else

        ga_old.hll.motor.active(Isub,:) = NaN(1, length(npts_hll));

    end
            
            
    if ~isempty(allData_old{Isub}.hll.motor.inactive)

        % get the mean over active electrodes
        ga_old.hll.motor.inactive(Isub,:) = nanmean(allData_old{Isub}.hll.motor.inactive.pca.trial,1);

    else

        ga_old.hll.motor.inactive(Isub,:) = NaN(1, length(npts_hll));

    end

    %% BR M1 "OLD"
    
    % get data for motor active/inactive electrodes
    
    if ~isempty(allData_old{Isub}.br.motor.active)

        % get the mean over active electrodes
        ga_old.br.motor.active(Isub,:) = nanmean(allData_old{Isub}.br.motor.active.pca.trial,1);

    else

        ga_old.br.motor.active(Isub,:) = NaN(1, length(npts_br));

    end
            
            
    if ~isempty(allData_old{Isub}.br.motor.inactive)

        % get the mean over active electrodes
        ga_old.br.motor.inactive(Isub,:) = nanmean(allData_old{Isub}.br.motor.inactive.pca.trial,1);

    else

        ga_old.br.motor.inactive(Isub,:) = NaN(1, length(npts_br));

    end

    %% HLL M1 "NEW"
    
    % get data for motor active/inactive electrodes
    
    if ~isempty(allData_new{Isub}.hll.motor.active)

        % get the mean over active electrodes
        ga_new.hll.motor.active(Isub,:) = nanmean(allData_new{Isub}.hll.motor.active.pca.trial,1);

    else

        ga_new.hll.motor.active(Isub,:) = NaN(1, length(npts_hll));

    end
                        
    if ~isempty(allData_new{Isub}.hll.motor.inactive)

        % get the mean over active electrodes
        ga_new.hll.motor.inactive(Isub,:) = nanmean(allData_new{Isub}.hll.motor.inactive.pca.trial,1);

    else

        ga_new.hll.motor.inactive(Isub,:) = NaN(1, length(npts_hll));

    end

    %% BR M1 "NEW"
    
    % get data for motor active/inactive electrodes
    
    if ~isempty(allData_new{Isub}.br.motor.active)

        % get the mean over active electrodes
        ga_new.br.motor.active(Isub,:) = nanmean(allData_new{Isub}.br.motor.active.pca.trial,1);

    else

        ga_new.br.motor.active(Isub,:) = NaN(1, length(npts_br));

    end
                        
    if ~isempty(allData_new{Isub}.br.motor.inactive)

        % get the mean over active electrodes
        ga_new.br.motor.inactive(Isub,:) = nanmean(allData_new{Isub}.br.motor.inactive.pca.trial,1);

    else

        ga_new.br.motor.inactive(Isub,:) = NaN(1, length(npts_br));

    end  

end
 
%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_S2');
cd(fullfile(path_source, 'Fig_S2'));

time = npts_hll;

% save data before orthogonalization
tmp = NaN(size(ga_old.hll.frontal.active,1)+1, size(ga_old.hll.frontal.active,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(ga_old.hll.frontal.active,1))';
tmp(2:end,2:end)  = ga_old.hll.frontal.active;

tmpTable = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S2_frontal_without_ortho.xlsx';
writetable(tmpTable,filename, 'Sheet', 1); clear tmpTable

% save data after orthogonalization
tmp = NaN(size(ga_new.hll.frontal.active,1)+1, size(ga_new.hll.frontal.active,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(ga_new.hll.frontal.active,1))';
tmp(2:end,2:end)  = ga_new.hll.frontal.active;

tmpTable = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S2_frontal_with_ortho.xlsx';
writetable(tmpTable,filename, 'Sheet', 1); clear tmpTable

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_S2');
cd(fullfile(path_source, 'Fig_S2'));

time = npts_hll;

% save data before orthogonalization
tmp = NaN(size(ga_old.hll.motor.active,1)+1, size(ga_old.hll.motor.active,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(ga_old.hll.motor.active,1))';
tmp(2:end,2:end)  = ga_old.hll.motor.active;

tmpTable = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S2_motor_without_ortho.xlsx';
writetable(tmpTable,filename, 'Sheet', 1); clear tmpTable

% save data after orthogonalization
tmp = NaN(size(ga_new.hll.motor.active,1)+1, size(ga_new.hll.motor.active,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(ga_new.hll.motor.active,1))';
tmp(2:end,2:end)  = ga_new.hll.motor.active;

tmpTable = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S2_motor_with_ortho.xlsx';
writetable(tmpTable,filename, 'Sheet', 1); clear tmpTable

%% RUN CLUSTER STATS TO TEST FOR DIFFERENCES IN EXPLAINED VARIANCE WITH AND WITHOUT ACCOUNTING FOR ACCURACY

% select the epoch 

selcond = input('Select Epoch: (1) HLL, (2) BR: ');

if selcond == 1
    % get data
    time = npts_hll;
    
    data2use_old.frontal  = ga_old.hll.frontal.active;
    data2use_old.motor    = ga_old.hll.motor.active;
        
    data2use_new.frontal  = ga_new.hll.frontal.active;
    data2use_new.motor    = ga_new.hll.motor.active;
else
    % get data
    time = npts_br;
    
    data2use_old.frontal  = ga_old.br.frontal.active;
    data2use_old.motor    = ga_old.br.motor.active;
        
    data2use_new.frontal  = ga_new.br.frontal.active;
    data2use_new.motor    = ga_new.br.motor.active;
end

% remove NaNs

idx = find(all(isnan(data2use_old.frontal),2) | all(isnan(data2use_new.frontal),2));
data2use_old.frontal(idx,:)   = [];
data2use_new.frontal(idx,:)   = [];

clear idx

idx = find(all(isnan(data2use_old.motor),2) | all(isnan(data2use_new.motor),2));
data2use_old.motor(idx,:)   = [];
data2use_new.motor(idx,:)   = [];

% Create FT structure
statsdata = struct();

statsdata.frontal.old  = cell(size(data2use_old.frontal,1),1);
statsdata.frontal.new  = cell(size(data2use_new.frontal,1),1);

statsdata.motor.old    = cell(size(data2use_old.motor,1),1);
statsdata.motor.new    = cell(size(data2use_new.motor,1),1);

% Get data for frontal cortex

for Isub = 1:numel(statsdata.frontal.old)
    
    % get "old" data
    statsdata.frontal.old{Isub}.avg      = data2use_old.frontal(Isub,:);
    statsdata.frontal.old{Isub}.dimord   = 'chan_time';
    statsdata.frontal.old{Isub}.time     = time;
    statsdata.frontal.old{Isub}.label    = {'X'};
    statsdata.frontal.old{Isub}.cfg      = [];

    % get "new" data
    statsdata.frontal.new{Isub}.avg      = data2use_new.frontal(Isub,:);
    statsdata.frontal.new{Isub}.dimord   = 'chan_time';
    statsdata.frontal.new{Isub}.time     = time;
    statsdata.frontal.new{Isub}.label    = {'X'};
    statsdata.frontal.new{Isub}.cfg      = [];
        
end

% Get data for motor cortex

for Isub = 1:numel(statsdata.motor.old)
    
    % get "old" data
    statsdata.motor.old{Isub}.avg      = data2use_old.motor(Isub,:);
    statsdata.motor.old{Isub}.dimord   = 'chan_time';
    statsdata.motor.old{Isub}.time     = time;
    statsdata.motor.old{Isub}.label    = {'X'};
    statsdata.motor.old{Isub}.cfg      = [];

    % get "new" data
    statsdata.motor.new{Isub}.avg      = data2use_new.motor(Isub,:);
    statsdata.motor.new{Isub}.dimord   = 'chan_time';
    statsdata.motor.new{Isub}.time     = time;
    statsdata.motor.new{Isub}.label    = {'X'};
    statsdata.motor.new{Isub}.cfg      = [];
        
end

stat      = []; % for frontal and motor cortex
cohensd   = []; % "

% stats frontal cortex

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = 'all';
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.alpha               = 0.05;
cfg.correctm            = 'cluster';
cfg.correcttail         = 'prob';
cfg.computestat         = 'yes'; 
cfg.computecritval      = 'yes'; 
cfg.computeprob         = 'yes'; 
cfg.numrandomization    = 10000;  
cfg.spmversion          = 'spm12';
    
% design matrix for frontal cortex

subj = numel(statsdata.frontal.old);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

stat.frontal = ft_timelockstatistics(cfg, statsdata.frontal.old{:}, statsdata.frontal.new{:});    

% compute cohen's d for frontal cortex
cfg            = [];
cfg.method     = 'analytic';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

cohensd.frontal = ft_timelockstatistics(cfg, statsdata.frontal.old{:}, statsdata.frontal.new{:});

% stats motor cortex

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = 'all';
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.alpha               = 0.05;
cfg.correctm            = 'cluster';
cfg.correcttail         = 'prob';
cfg.computestat         = 'yes'; 
cfg.computecritval      = 'yes'; 
cfg.computeprob         = 'yes'; 
cfg.numrandomization    = 10000;  
cfg.spmversion          = 'spm12';

subj = numel(statsdata.motor.old);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

stat.motor = ft_timelockstatistics(cfg, statsdata.motor.old{:}, statsdata.motor.new{:});    

% compute cohen's d for frontal cortex
cfg            = [];
cfg.method     = 'analytic';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

cohensd.motor = ft_timelockstatistics(cfg, statsdata.motor.old{:}, statsdata.motor.new{:});

% get grand average

GA_old = [];
GA_new = [];

cfg = [];
cfg.keepindividual = 'yes';

GA_old.frontal  = ft_timelockgrandaverage(cfg, statsdata.frontal.old{:});
GA_new.frontal  = ft_timelockgrandaverage(cfg, statsdata.frontal.new{:});

GA_old.motor    = ft_timelockgrandaverage(cfg, statsdata.motor.old{:});
GA_new.motor    = ft_timelockgrandaverage(cfg, statsdata.motor.new{:});
  
%% PLOT STATS FOR PFC

figure;
% old
a = shadedErrorBar(time, nanmean(squeeze(GA_old.frontal.individual))*100, std(squeeze(GA_old.frontal.individual))*100 ...
    / sqrt(size(GA_old.frontal.individual,1)));
a.mainLine.Color = [0.5 0.5 0.5]; a.mainLine.LineWidth = 2;
a.patch.FaceColor = [0.5 0.5 0.5]; a.patch.EdgeColor = [0.5 0.5 0.5]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0.5 0.5 0.5]; a.edge(2).Color = [0.5 0.5 0.5];
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

hold on

% new
a = shadedErrorBar(time, nanmean(squeeze(GA_new.frontal.individual))*100, std(squeeze(GA_new.frontal.individual))*100 ...
    / sqrt(size(GA_new.frontal.individual,1)));
a.mainLine.Color = [0 0 0]; a.mainLine.LineWidth = 2;
a.patch.FaceColor = [0 0 0]; a.patch.EdgeColor = [0 0 0]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0 0 0]; a.edge(2).Color = [0 0 0];
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

% plot significant timepoints
a = ylim;

island = bwconncomp(stat.frontal.mask);
for jj = 1:island.NumObjects
    h = line(time(island.PixelIdxList{jj}), ...
        repmat(a(1)-0.01, 1, numel(island.PixelIdxList{jj})), 'linewidth', 3);
    b = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 3);
    b.Color = [0 0 0];
    hold on
    b = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 3);
    b.Color = [0 0 0];
end
clear island

hold on

% make plot nice
xline(0, 'k--', 'linewidth', 0.5); yline(0, 'k--', 'linewidth', 0.5);
xlim([npts_hll(1) npts_hll(end)]); ylabel('PEV'); xlabel('Time to BR [s]');
set(gca, 'fontsize', 13, 'linewidth', 3); box off;
xticks([-0.4 -0.2 0 0.2]);

set(gcf, 'Position', [500 272 379 276]);

allposthoc = {stat.frontal};
allcohens  = {cohensd.frontal};
posthocidx = {'Frontal'};

for Ipair = 1:numel(allposthoc)

    % extract start point of cluster
    if isfield(allposthoc{Ipair}, 'posclusters')

        for Icluster = 1:numel(allposthoc{Ipair}.posclusters)

            fprintf('\n');

            fprintf('p-value cluster %d = %.3f\n', Icluster, allposthoc{Ipair}.posclusters(Icluster).prob);
            fprintf('clustermass cluster %d = %.3f\n', Icluster, allposthoc{Ipair}.posclusters(Icluster).clusterstat);
            
            begeffect = allposthoc{Ipair}.time...
                (find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'first'));
            endeffect = allposthoc{Ipair}.time...
                (find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'last'));

            fprintf('Effect (%s): cluster %d from %.3f to %.3f ms\n', posthocidx{Ipair}, Icluster, begeffect, endeffect);

            fprintf('Cohens d cluster %d = %.3f\n', Icluster, mean(allcohens{Ipair}.cohensd(...
                find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'first'):...
                find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'last'))));
            
            fprintf('\n');

        end

    end
    
end
   
% go to output folder to save figures
cd(path_out)

if selcond == 1
    print(gcf,'w2_trace_contrast_PFC_HLL.pdf','-dpdf','-r400');   
else
    print(gcf,'w2_trace_contrast_PFC_BR.pdf','-dpdf','-r400');   
end

%% PLOT STATS FOR MOTOR CORTEX

figure;

% old
a = shadedErrorBar(time, nanmean(squeeze(GA_old.motor.individual))*100, std(squeeze(GA_old.motor.individual))*100 ...
    / sqrt(size(GA_old.motor.individual,1)));
a.mainLine.Color = [0.5 0.5 0.5]; a.mainLine.LineWidth = 2;
a.patch.FaceColor = [0.5 0.5 0.5]; a.patch.EdgeColor = [0.5 0.5 0.5]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0.5 0.5 0.5]; a.edge(2).Color = [0.5 0.5 0.5];
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

hold on

% new
a = shadedErrorBar(time, nanmean(squeeze(GA_new.motor.individual))*100, std(squeeze(GA_new.motor.individual))*100 ...
    / sqrt(size(GA_new.motor.individual,1)));
a.mainLine.Color = [0 0 0]; a.mainLine.LineWidth = 2;
a.patch.FaceColor = [0 0 0]; a.patch.EdgeColor = [0 0 0]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0 0 0]; a.edge(2).Color = [0 0 0];
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

% plot significant timepoints
a = ylim;

island = bwconncomp(stat.motor.mask);
for jj = 1:island.NumObjects
    h = line(time(island.PixelIdxList{jj}), ...
        repmat(a(1)-0.01, 1, numel(island.PixelIdxList{jj})), 'linewidth', 3);
    b = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 3);
    b.Color = [0 0 0];
    hold on
    b = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 3);
    b.Color = [0 0 0];
end
clear island

hold on

% make plot nice
xline(0, 'k--', 'linewidth', 0.5); yline(0, 'k--', 'linewidth', 0.5);
xlim([npts_hll(1) npts_hll(end)]); ylabel('PEV'); xlabel('Time to BR [s]');
set(gca, 'fontsize', 13, 'linewidth', 3); box off;
xticks([-0.4 -0.2 0 0.2]);

set(gcf, 'Position', [500 272 379 276]);

allposthoc = {stat.motor};
allcohens  = {cohensd.motor};
posthocidx = {'Motor'};

for Ipair = 1:numel(allposthoc)

    % extract start point of cluster
    if isfield(allposthoc{Ipair}, 'posclusters')

        for Icluster = 1:numel(allposthoc{Ipair}.posclusters)

            fprintf('\n');

            fprintf('p-value cluster %d = %.3f\n', Icluster, allposthoc{Ipair}.posclusters(Icluster).prob);
            fprintf('clustermass cluster %d = %.3f\n', Icluster, allposthoc{Ipair}.posclusters(Icluster).clusterstat);
            
            begeffect = allposthoc{Ipair}.time...
                (find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'first'));
            endeffect = allposthoc{Ipair}.time...
                (find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'last'));

            fprintf('Effect (%s): cluster %d from %.3f to %.3f ms\n', posthocidx{Ipair}, Icluster, begeffect, endeffect);

            fprintf('Cohens d cluster %d = %.3f\n', Icluster, mean(allcohens{Ipair}.cohensd(...
                find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'first'):...
                find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'last'))));
            
            fprintf('\n');

        end

    end
    
end

% go to output folder to save figures
cd(path_out)

if selcond == 1
    print(gcf,'w2_trace_contrast_M1_HLL.pdf','-dpdf','-r400');   
else
    print(gcf,'w2_trace_contrast_M1_BR.pdf','-dpdf','-r400');   
end

