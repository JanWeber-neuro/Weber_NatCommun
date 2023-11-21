%% Figure S6a & b (left panel)
% this script rungs additional analysis required for the plots of figure
% S6a & b (left panel)

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

% initialize structure for data
data_frontal = struct(); data_motor = struct();

% loop over subjects to get data and make FT structure for stats
for Isub = 1:numel(subjects)
    
    % -----------------------
    % Frontal Cortex
    % -----------------------

    if isfield(sp_data{Isub}, 'frontal')
        
        % euclidean distance between 0 and 25%
        data_frontal.eucl.lh_0_25{Isub}.avg      = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh0_25)))';
        data_frontal.eucl.lh_0_25{Isub}.avgnorm  = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh0_25)))' - ...
            mean(sp_data{Isub}.frontal.euclDist.lh0_25, 'all');
        data_frontal.eucl.lh_0_25{Isub}.basenorm = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh0_25)))' - ...
            mean(squeeze(sp_data{Isub}.frontal.euclDist.lh0_25(:,:,1)), 'all');     
        data_frontal.eucl.lh_0_25{Isub}.time     = time;
        data_frontal.eucl.lh_0_25{Isub}.dimord   = 'chan_time';
        data_frontal.eucl.lh_0_25{Isub}.label    = {'frontal'};
        data_frontal.eucl.lh_0_25{Isub}.cfg      = [];
           
        % euclidean distance between 0 and 75%
        data_frontal.eucl.lh_0_75{Isub}.avg      = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh0_75)))';
        data_frontal.eucl.lh_0_75{Isub}.avgnorm  = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh0_75)))' - ...
            mean(sp_data{Isub}.frontal.euclDist.lh0_75, 'all');
        data_frontal.eucl.lh_0_75{Isub}.basenorm = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh0_75)))' - ...
            mean(squeeze(sp_data{Isub}.frontal.euclDist.lh0_75(:,:,1)), 'all');     
        data_frontal.eucl.lh_0_75{Isub}.time     = time;
        data_frontal.eucl.lh_0_75{Isub}.dimord   = 'chan_time';
        data_frontal.eucl.lh_0_75{Isub}.label    = {'frontal'};
        data_frontal.eucl.lh_0_75{Isub}.cfg      = [];
        
        % euclidean distance between 25 and 75%
        data_frontal.eucl.lh_25_75{Isub}.avg      = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh25_75)))';
        data_frontal.eucl.lh_25_75{Isub}.avgnorm  = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh25_75)))' - ...
            mean(sp_data{Isub}.frontal.euclDist.lh25_75, 'all');
        data_frontal.eucl.lh_25_75{Isub}.basenorm = squeeze(mean(mean(...
            sp_data{Isub}.frontal.euclDist.lh25_75)))' - ...
            mean(squeeze(sp_data{Isub}.frontal.euclDist.lh25_75(:,:,1)), 'all');     
        data_frontal.eucl.lh_25_75{Isub}.time     = time;
        data_frontal.eucl.lh_25_75{Isub}.dimord   = 'chan_time';
        data_frontal.eucl.lh_25_75{Isub}.label    = {'frontal'};
        data_frontal.eucl.lh_25_75{Isub}.cfg      = [];
        
        % euclidean distance between all conditions
        data_frontal.eucl.all{Isub}.avg           = squeeze(mean(mean(mean(...
            sp_data{Isub}.frontal.euclDist.all))))';
        data_frontal.eucl.all{Isub}.avgnorm       = squeeze(mean(mean(mean(...
            sp_data{Isub}.frontal.euclDist.all))))' - ...
            mean(sp_data{Isub}.frontal.euclDist.all, 'all');
        data_frontal.eucl.all{Isub}.basenorm      = squeeze(mean(mean(mean(...
            sp_data{Isub}.frontal.euclDist.all))))' - ...
            mean(squeeze(sp_data{Isub}.frontal.euclDist.all(:,:,1)), 'all');     
        data_frontal.eucl.all{Isub}.time          = time;
        data_frontal.eucl.all{Isub}.dimord        = 'chan_time';
        data_frontal.eucl.all{Isub}.label         = {'frontal'};
        data_frontal.eucl.all{Isub}.cfg           = [];        
            
                        
    end
    
    % -----------------------
    % Motor Cortex
    % -----------------------
    
    if isfield(sp_data{Isub}, 'motor')
        
        % euclidean distance between 0 and 25%
        data_motor.eucl.lh_0_25{Isub}.avg      = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh0_25)))'; % get the mean between pairwise trial differences 
        data_motor.eucl.lh_0_25{Isub}.avgnorm  = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh0_25)))' - ...
            mean(sp_data{Isub}.motor.euclDist.lh0_25, 'all');
        data_motor.eucl.lh_0_25{Isub}.basenorm = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh0_25)))' - ...
            mean(squeeze(sp_data{Isub}.motor.euclDist.lh0_25(:,:,1)), 'all');     
        data_motor.eucl.lh_0_25{Isub}.time     = time;
        data_motor.eucl.lh_0_25{Isub}.dimord   = 'chan_time';
        data_motor.eucl.lh_0_25{Isub}.label    = {'motor'};
        data_motor.eucl.lh_0_25{Isub}.cfg      = [];
           
        % euclidean distance between 0 and 75%
        data_motor.eucl.lh_0_75{Isub}.avg      = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh0_75)))';
        data_motor.eucl.lh_0_75{Isub}.avgnorm  = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh0_75)))' - ...
            mean(sp_data{Isub}.motor.euclDist.lh0_75, 'all');
        data_motor.eucl.lh_0_75{Isub}.basenorm = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh0_75)))' - ...
            mean(squeeze(sp_data{Isub}.motor.euclDist.lh0_75(:,:,1)), 'all');     
        data_motor.eucl.lh_0_75{Isub}.time     = time;
        data_motor.eucl.lh_0_75{Isub}.dimord   = 'chan_time';
        data_motor.eucl.lh_0_75{Isub}.label    = {'motor'};
        data_motor.eucl.lh_0_75{Isub}.cfg      = [];
        
        % euclidean distance between 25 and 75%
        data_motor.eucl.lh_25_75{Isub}.avg      = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh25_75)))';
        data_motor.eucl.lh_25_75{Isub}.avgnorm  = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh25_75)))' - ...
            mean(sp_data{Isub}.motor.euclDist.lh25_75, 'all');
        data_motor.eucl.lh_25_75{Isub}.basenorm = squeeze(mean(mean(...
            sp_data{Isub}.motor.euclDist.lh25_75)))' - ...
            mean(squeeze(sp_data{Isub}.motor.euclDist.lh25_75(:,:,1)), 'all');     
        data_motor.eucl.lh_25_75{Isub}.time     = time;
        data_motor.eucl.lh_25_75{Isub}.dimord   = 'chan_time';
        data_motor.eucl.lh_25_75{Isub}.label    = {'motor'};
        data_motor.eucl.lh_25_75{Isub}.cfg      = [];
        
        % euclidean distance between all conditions
        data_motor.eucl.all{Isub}.avg           = squeeze(mean(mean(mean(...
            sp_data{Isub}.motor.euclDist.all))))';
        data_motor.eucl.all{Isub}.avgnorm       = squeeze(mean(mean(mean(...
            sp_data{Isub}.motor.euclDist.all))))' - ...
            mean(sp_data{Isub}.motor.euclDist.all, 'all');
        data_motor.eucl.all{Isub}.basenorm      = squeeze(mean(mean(mean(...
            sp_data{Isub}.motor.euclDist.all))))' - ...
            mean(squeeze(sp_data{Isub}.motor.euclDist.all(:,:,1)), 'all');     
        data_motor.eucl.all{Isub}.time          = time;
        data_motor.eucl.all{Isub}.dimord        = 'chan_time';
        data_motor.eucl.all{Isub}.label         = {'motor'};
        data_motor.eucl.all{Isub}.cfg           = [];        
        
                
    end
        
end %Isub
    
%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_5a_left');
cd(fullfile(path_source, 'Fig_5a_left'));

sourcedata = struct;

time = data_frontal.eucl.all{1}.time;

for Isub = 1:numel(data_frontal.eucl.all)
    
    if isempty(data_frontal.eucl.all{Isub})
        
        sourcedata.all(Isub,:) = NaN(1, length(time));

    else
        
        sourcedata.all(Isub,:) = data_frontal.eucl.all{Isub}.basenorm;

    end        
        
end%Isub

% STORE IN TABLE %%
tmp = NaN(size(sourcedata.all,1)+1, size(sourcedata.all,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.all;

table_frontal.all = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5a_frontal_MDD.xlsx';
writetable(table_frontal.all,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_5b_left');
cd(fullfile(path_source, 'Fig_5b_left'));

sourcedata = struct;

time = data_motor.eucl.all{1}.time;

for Isub = 1:numel(data_motor.eucl.all)
    
    if isempty(data_motor.eucl.all{Isub})
        
        sourcedata.all(Isub,:) = NaN(1, length(time));

    else
        
        sourcedata.all(Isub,:) = data_motor.eucl.all{Isub}.basenorm;

    end        
        
end%Isub

% STORE IN TABLE %%
tmp = NaN(size(sourcedata.all,1)+1, size(sourcedata.all,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.all;

table_motor.all = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5b_motor_MDD.xlsx';
writetable(table_motor.all,filename, 'Sheet', 1);

%%

% remove empty cells
data_frontal.eucl.lh_0_25  = data_frontal.eucl.lh_0_25(~cellfun(@isempty, data_frontal.eucl.lh_0_25));
data_frontal.eucl.lh_0_75  = data_frontal.eucl.lh_0_75(~cellfun(@isempty, data_frontal.eucl.lh_0_75));
data_frontal.eucl.lh_25_75 = data_frontal.eucl.lh_25_75(~cellfun(@isempty, data_frontal.eucl.lh_25_75));
data_frontal.eucl.all      = data_frontal.eucl.all(~cellfun(@isempty, data_frontal.eucl.all));

data_motor.eucl.lh_0_25    = data_motor.eucl.lh_0_25(~cellfun(@isempty, data_motor.eucl.lh_0_25));
data_motor.eucl.lh_0_75    = data_motor.eucl.lh_0_75(~cellfun(@isempty, data_motor.eucl.lh_0_75));
data_motor.eucl.lh_25_75   = data_motor.eucl.lh_25_75(~cellfun(@isempty, data_motor.eucl.lh_25_75));
data_motor.eucl.all        = data_motor.eucl.all(~cellfun(@isempty, data_motor.eucl.all));

%% - COMPUTE STATISTICS - %%

% go to output folder
cd(path_out);

statsdata = struct();
 
% ---------------------
% Frontal Cortex
% ---------------------

statsdata.frontal   = data_frontal.eucl.all;

% create a "null-distribution"
statsdata.frontal_surro = statsdata.frontal;

for Isub = 1:numel(statsdata.frontal)
    
    statsdata.frontal_surro{Isub}.avg(:)      = 0;
    statsdata.frontal_surro{Isub}.avgnorm(:)  = 0;
    statsdata.frontal_surro{Isub}.basenorm(:) = 0;
    
end

% ---------------------
% Motor Cortex
% ---------------------

statsdata.motor   = data_motor.eucl.all;

% create a "null-distribution"
statsdata.motor_surro = statsdata.motor;

for Isub = 1:numel(statsdata.motor)
    
    statsdata.motor_surro{Isub}.avg(:)      = 0;
    statsdata.motor_surro{Isub}.avgnorm(:)  = 0;
    statsdata.motor_surro{Isub}.basenorm(:) = 0;
    
end


% ------------------------------
% Merge Data 
% ------------------------------

mergeData = {statsdata.frontal, statsdata.frontal_surro; ...
             statsdata.motor, statsdata.motor_surro};
         
stat = cell(size(mergeData,1),1);

for Iroi = 1:size(mergeData,1)

    cfg                     = [];
    cfg.channel             = 'all';
    cfg.latency             = 'all';
    cfg.tail                = 0;
    cfg.parameter           = 'basenorm';
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

    subj = numel(mergeData{Iroi,1});
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

    stat{Iroi} = ft_timelockstatistics(cfg, mergeData{Iroi,1}{:}, mergeData{Iroi,2}{:});

    % compute grand average
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter      = 'basenorm';

    tmpGA_all = ft_timelockgrandaverage(cfg, mergeData{Iroi,1}{:});

    % ---------------
    % Plot Data
    % ---------------
    
    figure(Iroi);
    
    a = shadedErrorBar(tmpGA_all.time, mean(squeeze(tmpGA_all.individual)), ...
        std(squeeze(tmpGA_all.individual)) / sqrt(size(tmpGA_all.individual,1)));
    a.mainLine.Color = roicolor{Iroi}; a.mainLine.LineWidth = 3;
    a.patch.FaceColor = roicolor{Iroi}; a.patch.EdgeColor = roicolor{Iroi}; a.patch.FaceAlpha = 0.5;

    hold on

    % find significant timepoints
    sigpnts = time(logical(stat{Iroi}.posclusterslabelmat == 1));
    
    if ~isempty(sigpnts)
        
        if Iroi == 1
            scatter(sigpnts, ones(1, length(sigpnts))*-0.5, 70, 'filled', 'MarkerFaceColor', roicolor{Iroi}, 'MarkerEdgeColor', roicolor{Iroi}, ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2);
        else
            scatter(sigpnts, ones(1, length(sigpnts))*-3.5, 70, 'filled', 'MarkerFaceColor', roicolor{Iroi}, 'MarkerEdgeColor', roicolor{Iroi}, ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2);
        end
        
    end
    
    xlim([tmpGA_all.time(1) tmpGA_all.time(end)])
    if Iroi == 1
        ylim([-1 6]);
    else
        ylim([-4 7]);
    end
            
    ylabel('MDD', 'FontSize', 13);
    set(gca, 'fontsize', 13, 'linewidth', 3);
    box off
    set(gcf, 'Position', [552 342 205 151]);
        
end



