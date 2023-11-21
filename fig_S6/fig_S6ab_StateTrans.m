%% Figure S6a & b (right panel)
% this script rungs additional analysis required for the plots of figure
% S6a & b (right panel)

%%

clear
ft_defaults

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_data = '/Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/HFA/ExpVariance_Approach/HFA_state_space_dynamics';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_5/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% - LOAD DATA - %%

load(fullfile(path_data, 'state_space_data_br_all_elecs.mat'));

%% - DEFINE COLORS - %%

color    = {'#008450', '#EFB700', '#B81D13'}; % color for condition plotting
roicolor = {'r', 'b'}; % color per ROI

%% - COMPUTE EUCLIDEAN DISTANCE - %%

% load single subject to get time
load /Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/OSL24/Data/HFA/PCA/PCA_ROI_Data/OSL24_activedata_ROI_pca.mat

fsample       = 512;
timewin       = roidata.hll.frontal.active.time; % full time window
stepsize      = round(fsample*0.02);             % stepsize for sliding window
winsize       = round(fsample*0.05);             % window size over which data will be averaged to smooth

% create Nx2 sample vector (with overlap of stepsize)
samples2use = [1:stepsize:length(timewin)-winsize; ...
    1+winsize:stepsize:length(timewin)]';            

time = roidata.hll.frontal.active.time(samples2use(:,1));

%% - COMPUTE STATE CHANGES - %%

% initialize structure for data
data_frontal = struct(); data_motor = struct();

for Isub = 1:numel(subjects)
    
    % -----------------------
    % Frontal Cortex
    % -----------------------

    if isfield(sp_data{Isub}, 'frontal')
        
        % activity change for 0%
        data_frontal.ac.lh_0{Isub}.avg       = mean(sp_data{Isub}.frontal.ActChange.lh_0);
        data_frontal.ac.lh_0{Isub}.basenorm  = mean(sp_data{Isub}.frontal.ActChange.lh_0) ...
            - mean(sp_data{Isub}.frontal.ActChange.lh_0(:,1));
        data_frontal.ac.lh_0{Isub}.time      = time;
        data_frontal.ac.lh_0{Isub}.dimord    = 'chan_time';
        data_frontal.ac.lh_0{Isub}.label     = {'frontal'};
        data_frontal.ac.lh_0{Isub}.cfg       = [];
           
        % activity change for 25%
        data_frontal.ac.lh_25{Isub}.avg      = mean(sp_data{Isub}.frontal.ActChange.lh_25);
        data_frontal.ac.lh_25{Isub}.basenorm = mean(sp_data{Isub}.frontal.ActChange.lh_25) ...
            - mean(sp_data{Isub}.frontal.ActChange.lh_25(:,1));
        data_frontal.ac.lh_25{Isub}.time     = time;
        data_frontal.ac.lh_25{Isub}.dimord   = 'chan_time';
        data_frontal.ac.lh_25{Isub}.label    = {'frontal'};
        data_frontal.ac.lh_25{Isub}.cfg      = [];
        
        % activity change for 75%
        data_frontal.ac.lh_75{Isub}.avg      = mean(sp_data{Isub}.frontal.ActChange.lh_75);
        data_frontal.ac.lh_75{Isub}.basenorm = mean(sp_data{Isub}.frontal.ActChange.lh_75) ...
            - mean(sp_data{Isub}.frontal.ActChange.lh_75(:,1));
        data_frontal.ac.lh_75{Isub}.time     = time;
        data_frontal.ac.lh_75{Isub}.dimord   = 'chan_time';
        data_frontal.ac.lh_75{Isub}.label    = {'frontal'};
        data_frontal.ac.lh_75{Isub}.cfg      = [];
                
    end
    
    % -----------------------
    % Motor Cortex
    % -----------------------

    if isfield(sp_data{Isub}, 'motor')
        
        % activity change for 0%
        data_motor.ac.lh_0{Isub}.avg       = mean(sp_data{Isub}.motor.ActChange.lh_0);
        data_motor.ac.lh_0{Isub}.basenorm  = mean(sp_data{Isub}.motor.ActChange.lh_0) ...
            - mean(sp_data{Isub}.motor.ActChange.lh_0(:,1));
        data_motor.ac.lh_0{Isub}.time      = time;
        data_motor.ac.lh_0{Isub}.dimord    = 'chan_time';
        data_motor.ac.lh_0{Isub}.label     = {'motor'};
        data_motor.ac.lh_0{Isub}.cfg       = [];
           
        % activity change for 25%
        data_motor.ac.lh_25{Isub}.avg      = mean(sp_data{Isub}.motor.ActChange.lh_25);
        data_motor.ac.lh_25{Isub}.basenorm = mean(sp_data{Isub}.motor.ActChange.lh_25) ...
            - mean(sp_data{Isub}.motor.ActChange.lh_25(:,1));
        data_motor.ac.lh_25{Isub}.time     = time;
        data_motor.ac.lh_25{Isub}.dimord   = 'chan_time';
        data_motor.ac.lh_25{Isub}.label    = {'motor'};
        data_motor.ac.lh_25{Isub}.cfg      = [];
        
        % activity change for 75%
        data_motor.ac.lh_75{Isub}.avg      = mean(sp_data{Isub}.motor.ActChange.lh_75);
        data_motor.ac.lh_75{Isub}.basenorm = mean(sp_data{Isub}.motor.ActChange.lh_75) ...
            - mean(sp_data{Isub}.motor.ActChange.lh_75(:,1));
        data_motor.ac.lh_75{Isub}.time     = time;
        data_motor.ac.lh_75{Isub}.dimord   = 'chan_time';
        data_motor.ac.lh_75{Isub}.label    = {'motor'};
        data_motor.ac.lh_75{Isub}.cfg      = [];
                
    end
        
end %Isub
    
%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_5a_right');
cd(fullfile(path_source, 'Fig_5a_right'));

sourcedata = struct;

time = data_frontal.ac.lh_0{1}.time;

for Isub = 1:numel(data_frontal.ac.lh_0)
    
    if isempty(data_frontal.ac.lh_0{Isub})
        
        sourcedata.lh_0(Isub,:)  = NaN(1, length(time));
        sourcedata.lh_25(Isub,:) = NaN(1, length(time));
        sourcedata.lh_75(Isub,:) = NaN(1, length(time));

    else
        
        sourcedata.lh_0(Isub,:)  = data_frontal.ac.lh_0{Isub}.avg;
        sourcedata.lh_25(Isub,:) = data_frontal.ac.lh_25{Isub}.avg;
        sourcedata.lh_75(Isub,:) = data_frontal.ac.lh_75{Isub}.avg;

    end        
        
end%Isub

% - STORE 0% LH - %
tmp = NaN(size(sourcedata.lh_0,1)+1, size(sourcedata.lh_0,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_frontal.lh_0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5a_frontal_StateTransition_lh0.xlsx';
writetable(table_frontal.lh_0,filename, 'Sheet', 1);

% - STORE 25% LH - %
tmp = NaN(size(sourcedata.lh_25,1)+1, size(sourcedata.lh_25,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_frontal.lh_25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5a_frontal_StateTransition_lh25.xlsx';
writetable(table_frontal.lh_25,filename, 'Sheet', 1);

% - STORE 75% LH - %
tmp = NaN(size(sourcedata.lh_75,1)+1, size(sourcedata.lh_75,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.lh_75;

table_frontal.lh_75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5a_frontal_StateTransition_lh75.xlsx';
writetable(table_frontal.lh_75,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_5b_right');
cd(fullfile(path_source, 'Fig_5b_right'));

sourcedata = struct;

time = data_motor.ac.lh_0{1}.time;

for Isub = 1:numel(data_motor.ac.lh_0)
    
    if isempty(data_motor.ac.lh_0{Isub})
        
        sourcedata.lh_0(Isub,:)  = NaN(1, length(time));
        sourcedata.lh_25(Isub,:) = NaN(1, length(time));
        sourcedata.lh_75(Isub,:) = NaN(1, length(time));

    else
        
        sourcedata.lh_0(Isub,:)  = data_motor.ac.lh_0{Isub}.avg;
        sourcedata.lh_25(Isub,:) = data_motor.ac.lh_25{Isub}.avg;
        sourcedata.lh_75(Isub,:) = data_motor.ac.lh_75{Isub}.avg;

    end        
        
end%Isub

% - STORE 0% LH - %
tmp = NaN(size(sourcedata.lh_0,1)+1, size(sourcedata.lh_0,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_motor.lh_0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5b_motor_StateTransition_lh0.xlsx';
writetable(table_motor.lh_0,filename, 'Sheet', 1);

% - STORE 25% LH - %
tmp = NaN(size(sourcedata.lh_25,1)+1, size(sourcedata.lh_25,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_motor.lh_25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5b_motor_StateTransition_lh25.xlsx';
writetable(table_motor.lh_25,filename, 'Sheet', 1);

% - STORE 75% LH - %
tmp = NaN(size(sourcedata.lh_75,1)+1, size(sourcedata.lh_75,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_75;

table_motor.lh_75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5b_motor_StateTransition_lh75.xlsx';
writetable(table_motor.lh_75,filename, 'Sheet', 1);

%%

% remove empty cells
data_frontal.ac.lh_0     = data_frontal.ac.lh_0(~cellfun(@isempty, data_frontal.ac.lh_0));
data_frontal.ac.lh_25    = data_frontal.ac.lh_25(~cellfun(@isempty, data_frontal.ac.lh_25));
data_frontal.ac.lh_75    = data_frontal.ac.lh_75(~cellfun(@isempty, data_frontal.ac.lh_75));

data_motor.ac.lh_0       = data_motor.ac.lh_0(~cellfun(@isempty, data_motor.ac.lh_0));
data_motor.ac.lh_25      = data_motor.ac.lh_25(~cellfun(@isempty, data_motor.ac.lh_25));
data_motor.ac.lh_75      = data_motor.ac.lh_75(~cellfun(@isempty, data_motor.ac.lh_75));

%% - COMPUTE STATS - %%

% Get data of interest
q       = 'Statistic: Frontal = 1; Motor = 2: ';
selroi = input(q);

statsdata = struct();

if selroi == 1
    
    statsdata.lh_0   = data_frontal.ac.lh_0;
    statsdata.lh_25  = data_frontal.ac.lh_25;
    statsdata.lh_75  = data_frontal.ac.lh_75;

elseif selroi == 2
    
    statsdata.lh_0   = data_motor.ac.lh_0;
    statsdata.lh_25  = data_motor.ac.lh_25;
    statsdata.lh_75  = data_motor.ac.lh_75;

end
    
% Run stats

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = 'all';
cfg.tail                = 1;
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesFunivariate';
cfg.alpha               = 0.05;
cfg.correctm            = 'cluster';
cfg.numrandomization    = 10000;  
cfg.spmversion          = 'spm12';
cfg.computestat         = 'yes'; 
cfg.computecritval      = 'yes'; 
cfg.computeprob         = 'yes'; 

% make design matrix
subj = numel(statsdata.lh_0);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
for i = 1:subj
design(1,subj+subj+i) = i;
end
design(2,1:subj)             = 1;
design(2,subj+1:2*subj)      = 2;
design(2,subj+subj+1:3*subj) = 3;

% clarify independent and dependent variable
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

fstat = ft_timelockstatistics(cfg, statsdata.lh_0{:}, statsdata.lh_25{:}, statsdata.lh_75{:});    

%% - PAIRWISE COMPARISONS - %%

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = 'all';
cfg.tail                = 0;
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.alpha               = 0.05;
cfg.correctm            = 'cluster';
cfg.correcttail         = 'prob';
cfg.numrandomization    = 10000;  
cfg.spmversion          = 'spm12';
cfg.computestat         = 'yes'; 
cfg.computecritval      = 'yes'; 
cfg.computeprob         = 'yes'; 

subj = numel(statsdata.lh_0);
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

stat = cell(3,1);
stat{1}  = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_0{:});
stat{2}  = ft_timelockstatistics(cfg, statsdata.lh_25{:}, statsdata.lh_0{:});
stat{3}  = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_25{:});

%% - PLOT STATS - %%

% get current grand averages

cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter      = 'avg';

tmpGA_lh0     = ft_timelockgrandaverage(cfg, statsdata.lh_0{:});
tmpGA_lh25    = ft_timelockgrandaverage(cfg, statsdata.lh_25{:});
tmpGA_lh75    = ft_timelockgrandaverage(cfg, statsdata.lh_75{:});

% plot data

plot(tmpGA_lh0.time, smoothdata(nanmean(squeeze(tmpGA_lh0.individual)), 'gaussian', round(fsample*0.01)), 'Color', '#008450', 'LineWidth', 4);
hold on
plot(tmpGA_lh25.time, smoothdata(nanmean(squeeze(tmpGA_lh25.individual)), 'gaussian', round(fsample*0.01)), 'Color', '#EFB700', 'LineWidth', 4);
hold on
plot(tmpGA_lh75.time, smoothdata(nanmean(squeeze(tmpGA_lh75.individual)), 'gaussian', round(fsample*0.01)), 'Color', '#B81D13', 'LineWidth', 4);

ylabel('State Change');
% xlabel('Time to BR');

% plot significant timepoints
a = ylim;

hold on

% mask f-stats
island = bwconncomp(fstat.mask);
for jj = 1:island.NumObjects
    h = line(time(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.001, 1, numel(island.PixelIdxList{jj})), 'linewidth', 5);
    plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', 'k', 'LineWidth', 5);
    hold on
    plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', 'k', 'LineWidth', 5);
end 
clear island

if isfield(fstat, 'posclusters') && fstat.posclusters(1).prob < 0.05

    % mask 0 vs. 75%
    island = bwconncomp(stat{1}.mask);
    for jj = 1:island.NumObjects
        h = line(time(island.PixelIdxList{jj}), ...
            repmat(a(2)+0.002, 1, numel(island.PixelIdxList{jj})), 'linewidth', 5);
        plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#008450', 'LineWidth', 5);
        hold on
        plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#B81D13', 'LineWidth', 5);
    end
    clear island

    % mask 0 vs. 25%
    island = bwconncomp(stat{2}.mask);
    for jj = 1:island.NumObjects
        h = line(time(island.PixelIdxList{jj}), ...
            repmat(a(2)+0.003, 1, numel(island.PixelIdxList{jj})), 'linewidth', 5);
        plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#008450', 'LineWidth', 5);
        hold on
        plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#EFB700', 'LineWidth', 5);
    end
    clear island

    % mask 25 vs. 75%
    island = bwconncomp(stat{3}.mask);
    for jj = 1:island.NumObjects
        h = line(time(island.PixelIdxList{jj}), ...
            repmat(a(2)+0.004, 1, numel(island.PixelIdxList{jj})), 'linewidth', 5);
        plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#EFB700', 'LineWidth', 5);
        hold on
        plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#B81D13', 'LineWidth', 5);
    end
    clear island
    
end

a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time(1) time(end)], ...
    'Ylim', [a(1) a(2)+0.001]); box off; hold on;
set(gcf, 'Position', [552 342 205 151]);


%% - SAVE FIGURES - %%

if selroi == 1
    print(gcf,'fig_5_frontal_state_change.pdf','-dpdf','-r400');    
else
    print(gcf,'fig_5_motor_state_change.pdf','-dpdf','-r400');    
end
