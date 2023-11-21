%% Figure 1d
% this script rungs additional analysis required for the plots of figure 1d

%%

clear
ft_defaults

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_1/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/;

allData     = cell(numel(subjects),1);

%% - LOAD DATA - %%

for Isub = 1:numel(subjects)
        
    % get indicated subject ID

    subID = subjects{Isub};

    % load data
    load(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data', [subID '_activedata_ROI_pca.mat']));
                        
    data = roidata; clear roidata;
        
    allData{Isub} = data; clear data
        
end

%% - EXTRACT SIGNAL PER ROI - %%

% name for regions of interest within structure 
rois = {'frontal', 'motor'};

% load example data to store time 
load /Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/OSL13/Data/HFA/PCA/OSL13_TaskActive_PCAexpvar.mat;

npts_hll    = pcadata.hll.pca.time;
npts_br     = pcadata.br.pca.time;

% initialize storage for grand average data
ga = struct();

for Isub = 1:numel(subjects)
   
    %% - HLL FRONTAL - %%
    
    % get data for frontal active/inactive electrodes
    
    if ~isempty(allData{Isub}.hll.frontal.active)

        % get the mean over active electrodes
        ga.hll.frontal.active(Isub,:) = nanmean(allData{Isub}.hll.frontal.active.pca.trial,1);

    else

        ga.hll.frontal.active(Isub,:) = NaN(1, length(npts_hll));

    end
            
            
    if ~isempty(allData{Isub}.hll.frontal.inactive)

        % get the mean over active electrodes
        ga.hll.frontal.inactive(Isub,:) = nanmean(allData{Isub}.hll.frontal.inactive.pca.trial,1);

    else

        ga.hll.frontal.inactive(Isub,:) = NaN(1, length(npts_hll));

    end

    
    %% - HLL MOTOR - %%
    
    % get data for motor active/inactive electrodes
    
    if ~isempty(allData{Isub}.hll.motor.active)

        % get the mean over active electrodes
        ga.hll.motor.active(Isub,:) = nanmean(allData{Isub}.hll.motor.active.pca.trial,1);

    else

        ga.hll.motor.active(Isub,:) = NaN(1, length(npts_hll));

    end
            
            
    if ~isempty(allData{Isub}.hll.motor.inactive)

        % get the mean over active electrodes
        ga.hll.motor.inactive(Isub,:) = nanmean(allData{Isub}.hll.motor.inactive.pca.trial,1);

    else

        ga.hll.motor.inactive(Isub,:) = NaN(1, length(npts_hll));

    end
            
    %% - BR FRONTAL - %%
    
    % get data for frontal active/inactive electrodes
    
    if ~isempty(allData{Isub}.br.frontal.active)

        % get the mean over active electrodes
        ga.br.frontal.active(Isub,:) = nanmean(allData{Isub}.br.frontal.active.pca.trial,1);

    else

        ga.br.frontal.active(Isub,:) = NaN(1, length(npts_br));

    end
            
            
    if ~isempty(allData{Isub}.br.frontal.inactive)

        % get the mean over active electrodes
        ga.br.frontal.inactive(Isub,:) = nanmean(allData{Isub}.br.frontal.inactive.pca.trial,1);

    else

        ga.br.frontal.inactive(Isub,:) = NaN(1, length(npts_br));

    end

    
    %% - BR MOTOR - %%
    
    % get data for motor active/inactive electrodes
    
    if ~isempty(allData{Isub}.br.motor.active)

        % get the mean over active electrodes
        ga.br.motor.active(Isub,:) = nanmean(allData{Isub}.br.motor.active.pca.trial,1);

    else

        ga.br.motor.active(Isub,:) = NaN(1, length(npts_br));

    end
            
            
    if ~isempty(allData{Isub}.br.motor.inactive)

        % get the mean over active electrodes
        ga.br.motor.inactive(Isub,:) = nanmean(allData{Isub}.br.motor.inactive.pca.trial,1);

    else

        ga.br.motor.inactive(Isub,:) = NaN(1, length(npts_br));

    end
       
end
    
%% Save Source Data

% go into results folder
mkdir(path_source, 'Fig_1d');
cd(fullfile(path_source, 'Fig_1d'));

time = npts_hll;

% -- FRONTAL CORTEX - %%
tmp = NaN(size(ga.hll.frontal.active,1)*2+1, size(ga.hll.frontal.active,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = repmat(1:18,1,2)';
tmp(2:end,2:end)  = [ga.hll.frontal.active; ga.hll.frontal.inactive] * 100;

table_frontal = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'pev_frontal.xlsx';
writetable(table_frontal,filename, 'Sheet', 1);

% -- MOTOR CORTEX - %%
tmp = NaN(size(ga.hll.motor.active,1)*2, size(ga.hll.motor.active,2)+1);
tmp(:,1)     = repmat(1:18,1,2)';
tmp(:,2:end) = [ga.hll.motor.active; ga.hll.motor.inactive] * 100;

table_motor  = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'pev_motor.xlsx';
writetable(table_motor,filename, 'Sheet', 1);

%% Define colors for plotting

color.m1    = [0 0 1];
color.pfc   = [1 0 0];

%% Compute Stats Between Regions

% select the epoch 

selcond = input('Select Epoch: (1) HLL, (2) BR: ');

if selcond == 1
    time = npts_hll;
    data2use_frontal = ga.hll.frontal.active;
    data2use_motor   = ga.hll.motor.active;
else
    time = npts_br;
    data2use_frontal = ga.br.frontal.active;
    data2use_motor   = ga.br.motor.active;
end

statsdata.frontal = cell(numel(subjects),1);
statsdata.motor   = cell(numel(subjects),1);

for Isub = 1:numel(subjects)
    
    % get data frontal cortex
    if all(isnan(data2use_frontal(Isub,:)),2)
        continue;
    else
        statsdata.frontal{Isub}.avg    = data2use_frontal(Isub,:);
        statsdata.frontal{Isub}.dimord = 'chan_time';
        statsdata.frontal{Isub}.time   = time;
        statsdata.frontal{Isub}.label  = {'X'};
        statsdata.frontal{Isub}.cfg    = [];
    end
    
    % get data motor cortex
    if all(isnan(data2use_motor(Isub,:)),2)
        continue;
    else
        statsdata.motor{Isub}.avg      = data2use_motor(Isub,:);
        statsdata.motor{Isub}.dimord   = 'chan_time';
        statsdata.motor{Isub}.time     = time;
        statsdata.motor{Isub}.label    = {'X'};
        statsdata.motor{Isub}.cfg      = [];
    end
    
    
end

% remove empty cells
statsdata.frontal = statsdata.frontal(~cellfun(@isempty, statsdata.frontal));
statsdata.motor   = statsdata.motor(~cellfun(@isempty, statsdata.motor));

% do the stats

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = 'all';
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.alpha               = 0.05;
cfg.correctm            = 'cluster';
cfg.correcttail         = 'prob';
cfg.computestat         = 'yes'; 
cfg.computecritval      = 'yes'; 
cfg.computeprob         = 'yes'; 
cfg.numrandomization    = 10000;  
cfg.spmversion          = 'spm12';
    
subj  = numel(statsdata.frontal);
design_frontal = zeros(1,1*subj);
for i = 1:subj
design_frontal(1,i) = 1;
end
clear subj
subj  = numel(statsdata.motor);
design_motor = zeros(1,1*subj);
for i = 1:subj
design_motor(1,i) = 2;
end
design = horzcat(design_frontal, design_motor);

cfg.design = design;
cfg.ivar  = 1;

stat = ft_timelockstatistics(cfg, statsdata.frontal{:}, statsdata.motor{:});    

% plot the stats

plot(stat.time, stat.stat, 'k', 'LineWidth', 3);
hold on
plot(stat.time, stat.mask);

if selcond == 1
    xlabel('Time to HLL [s]');
else
    xlabel('Time to BR [s]');
end

ylabel('stats [t]');

set(gca, 'fontsize', 13, 'linewidth', 3);
box off

% check cluster stats and print
if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    fprintf('First positive cluster: p = %.3f\n', stat.posclusters(1).prob);
end

% check cluster stats and print
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
    fprintf('First negative cluster: p = %.3f\n', stat.negclusters(1).prob);
end

%% Compute Stats Within a Region (Active vs. Inactive Electrodes)

% select the epoch 

selcond = input('Select Epoch: (1) HLL, (2) BR: ');

data2use_active   = struct();
data2use_inactive = struct();

if selcond == 1
    
    time = npts_hll;
    data2use_active.frontal     = ga.hll.frontal.active;
    data2use_inactive.frontal   = ga.hll.frontal.inactive;
    data2use_active.motor       = ga.hll.motor.active;
    data2use_inactive.motor     = ga.hll.motor.inactive;
    
else
    
    time = npts_br;
    data2use_active.frontal     = ga.br.frontal.active;
    data2use_inactive.frontal   = ga.br.frontal.inactive;
    data2use_active.motor       = ga.br.motor.active;
    data2use_inactive.motor     = ga.br.motor.inactive;

end

% remove NaNs

idx = find(all(isnan(data2use_active.frontal),2) | all(isnan(data2use_inactive.frontal),2));
data2use_active.frontal(idx,:)   = [];
data2use_inactive.frontal(idx,:) = [];

clear idx

idx = find(all(isnan(data2use_active.motor),2) | all(isnan(data2use_inactive.motor),2));
data2use_active.motor(idx,:)     = [];
data2use_inactive.motor(idx,:)   = [];


% Create FT structure
statsdata = struct();

statsdata.frontal.active     = cell(size(data2use_active.frontal,1),1);
statsdata.frontal.inactive   = cell(size(data2use_active.frontal,1),1);

statsdata.motor.active       = cell(size(data2use_active.motor,1),1);
statsdata.motor.inactive     = cell(size(data2use_active.motor,1),1);

% Get data for frontal cortex

for Isub = 1:numel(statsdata.frontal.active)
    
        
    % get data for active elecs
    statsdata.frontal.active{Isub}.avg        = data2use_active.frontal(Isub,:);
    statsdata.frontal.active{Isub}.dimord     = 'chan_time';
    statsdata.frontal.active{Isub}.time       = time;
    statsdata.frontal.active{Isub}.label      = {'X'};
    statsdata.frontal.active{Isub}.cfg        = [];

    % get data for inactive elecs
    statsdata.frontal.inactive{Isub}.avg      = data2use_inactive.frontal(Isub,:);
    statsdata.frontal.inactive{Isub}.dimord   = 'chan_time';
    statsdata.frontal.inactive{Isub}.time     = time;
    statsdata.frontal.inactive{Isub}.label    = {'X'};
    statsdata.frontal.inactive{Isub}.cfg      = [];
        
end


% Same for motor cortex
for Isub = 1:numel(statsdata.motor.active)
    
        
    % get data for active elecs
    statsdata.motor.active{Isub}.avg        = data2use_active.motor(Isub,:);
    statsdata.motor.active{Isub}.dimord     = 'chan_time';
    statsdata.motor.active{Isub}.time       = time;
    statsdata.motor.active{Isub}.label      = {'X'};
    statsdata.motor.active{Isub}.cfg        = [];

    % get data for inactive elecs
    statsdata.motor.inactive{Isub}.avg      = data2use_inactive.motor(Isub,:);
    statsdata.motor.inactive{Isub}.dimord   = 'chan_time';
    statsdata.motor.inactive{Isub}.time     = time;
    statsdata.motor.inactive{Isub}.label    = {'X'};
    statsdata.motor.inactive{Isub}.cfg      = [];
        
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

subj = numel(statsdata.frontal.active);
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

stat.frontal = ft_timelockstatistics(cfg, statsdata.frontal.active{:}, statsdata.frontal.inactive{:});    

% compute cohen's d for frontal cortex
cfg            = [];
cfg.method     = 'analytic';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

cohensd.frontal = ft_timelockstatistics(cfg, statsdata.frontal.active{:}, statsdata.frontal.inactive{:});

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

subj = numel(statsdata.motor.active);
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

stat.motor = ft_timelockstatistics(cfg, statsdata.motor.active{:}, statsdata.motor.inactive{:});    

% compute cohen's d for frontal cortex
cfg            = [];
cfg.method     = 'analytic';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

cohensd.motor = ft_timelockstatistics(cfg, statsdata.motor.active{:}, statsdata.motor.inactive{:});

% get grand average

GA_active   = [];
GA_inactive = [];

cfg = [];
cfg.keepindividual = 'yes';

GA_active.frontal    = ft_timelockgrandaverage(cfg, statsdata.frontal.active{:});
GA_inactive.frontal  = ft_timelockgrandaverage(cfg, statsdata.frontal.inactive{:});

GA_active.motor      = ft_timelockgrandaverage(cfg, statsdata.motor.active{:});
GA_inactive.motor    = ft_timelockgrandaverage(cfg, statsdata.motor.inactive{:});
  
% plot data with stats

% active
a = shadedErrorBar(time, nanmean(squeeze(GA_active.frontal.individual))*100, std(squeeze(GA_active.frontal.individual))*100 ...
    / sqrt(size(GA_active.frontal.individual,1)));
a.mainLine.Color = color.pfc; a.mainLine.LineWidth = 2;
a.patch.FaceColor = color.pfc; a.patch.EdgeColor = color.pfc; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = color.pfc; a.edge(2).Color = color.pfc;
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

hold on

% % inactive
plot(time, nanmean(squeeze(GA_inactive.frontal.individual))*100, 'Color', color.pfc, 'LineStyle', '--', 'LineWidth', 2);

% plot motor

% active
a = shadedErrorBar(time, nanmean(squeeze(GA_active.motor.individual))*100, std(squeeze(GA_active.motor.individual))*100 ...
    / sqrt(size(GA_active.motor.individual,1)));
a.mainLine.Color = color.m1; a.mainLine.LineWidth = 2;
a.patch.FaceColor = color.m1; a.patch.EdgeColor = color.m1; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = color.m1; a.edge(2).Color = color.m1;
a.edge(1).LineWidth = 0.1; a.edge(2).LineWidth = 0.1;

hold on

% % inactive
plot(time, nanmean(squeeze(GA_inactive.motor.individual))*100, 'Color', color.m1, 'LineStyle', '--', 'LineWidth', 2);

% plot significant timepoints
a = ylim;

island = bwconncomp(stat.frontal.mask);
for jj = 1:island.NumObjects
    h = line(time(island.PixelIdxList{jj}), ...
        repmat(a(1)-0.01, 1, numel(island.PixelIdxList{jj})), 'linewidth', 3);
    b = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 3);
    b.Color = color.pfc;
    hold on
    b = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 3);
    b.Color = color.pfc;
end
clear island

hold on

island = bwconncomp(stat.motor.mask);
for jj = 1:island.NumObjects
    h = line(time(island.PixelIdxList{jj}), ...
        repmat(a(1)-0.2, 1, numel(island.PixelIdxList{jj})), 'linewidth', 3);
    b = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 3);
    b.Color = color.m1;
    hold on
    b = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 3);
    b.Color = color.m1;
end
clear island

% make plot nice
xline(0, 'k--', 'linewidth', 0.5); yline(0, 'k--', 'linewidth', 0.5);
xlim([npts_hll(1) npts_hll(end)]); ylabel('%EV'); 

if selcond == 1
    xlabel('Time to HLL [s]')
else
    xlabel('Time to BR [s]')
end

set(gca, 'fontsize', 13, 'linewidth', 3); box off;
xticks([-0.4 -0.2 0 0.2]);

set(gcf, 'Position', [500 380 252 168]);
