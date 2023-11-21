%% Figure S6g
% this script rungs additional analysis required for the plots of figure
% S6g

%%

clear
ft_defaults

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_data = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_5/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

%%

fsample   = 512;

smoothwin = 0.025; % window for smoothing

%%

stimdata       = struct();
actdata        = struct();
surrodata_stim = struct();
surrodata_act  = struct();

for Isub = 1:numel(subjects)
    
    subID = subjects{Isub};
    
    if isfile(fullfile(path_data, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']))
        
        load(fullfile(path_data, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']));
        
        % --------------------------------
        % Frontal Cortex
        % --------------------------------

        if isfield(all_data, 'frontal') && ~isempty(all_data.frontal.stimdim)
            
            stimdata.frontal{Isub}.avg          = smoothdata(all_data.frontal.DecAcc_stim(all_data.frontal.stimdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
            stimdata.frontal{Isub}.time         = all_data.frontal.time;
            stimdata.frontal{Isub}.dimord       = 'chan_time';
            stimdata.frontal{Isub}.label        = {'Frontal'};
            stimdata.frontal{Isub}.cfg          = [];
            
            surrodata_stim.frontal{Isub}.avg    = ones(1, length(all_data.frontal.time))*0.33;
            surrodata_stim.frontal{Isub}.time   = all_data.frontal.time;
            surrodata_stim.frontal{Isub}.dimord = 'chan_time';
            surrodata_stim.frontal{Isub}.label  = {'Frontal'};
            surrodata_stim.frontal{Isub}.cfg    = [];
            
        end
        
        if isfield(all_data, 'frontal') && ~isempty(all_data.frontal.actdim)
            
            actdata.frontal{Isub}.avg          = smoothdata(all_data.frontal.DecAcc_act(all_data.frontal.actdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
            actdata.frontal{Isub}.time         = all_data.frontal.time;
            actdata.frontal{Isub}.dimord       = 'chan_time';
            actdata.frontal{Isub}.label        = {'Frontal'};
            actdata.frontal{Isub}.cfg          = [];
            
            surrodata_act.frontal{Isub}.avg    = ones(1, length(all_data.frontal.time))*0.33;
            surrodata_act.frontal{Isub}.time   = all_data.frontal.time;
            surrodata_act.frontal{Isub}.dimord = 'chan_time';
            surrodata_act.frontal{Isub}.label  = {'Frontal'};
            surrodata_act.frontal{Isub}.cfg    = [];
        
        end
        
        % --------------------------------
        % Motor Cortex
        % --------------------------------

        if isfield(all_data, 'motor') && ~isempty(all_data.motor.stimdim)
            
            stimdata.motor{Isub}.avg          = smoothdata(all_data.motor.DecAcc_stim(all_data.motor.stimdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
            stimdata.motor{Isub}.time         = all_data.motor.time;
            stimdata.motor{Isub}.dimord       = 'chan_time';
            stimdata.motor{Isub}.label        = {'motor'};
            stimdata.motor{Isub}.cfg          = [];
            
            surrodata_stim.motor{Isub}.avg    = ones(1, length(all_data.motor.time))*0.33;
            surrodata_stim.motor{Isub}.time   = all_data.motor.time;
            surrodata_stim.motor{Isub}.dimord = 'chan_time';
            surrodata_stim.motor{Isub}.label  = {'motor'};
            surrodata_stim.motor{Isub}.cfg    = [];
            
        end
        
        if isfield(all_data, 'motor') && ~isempty(all_data.motor.actdim)
            
            actdata.motor{Isub}.avg          = smoothdata(all_data.motor.DecAcc_act(all_data.motor.actdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
            actdata.motor{Isub}.time         = all_data.motor.time;
            actdata.motor{Isub}.dimord       = 'chan_time';
            actdata.motor{Isub}.label        = {'motor'};
            actdata.motor{Isub}.cfg          = [];
            
            surrodata_act.motor{Isub}.avg    = ones(1, length(all_data.motor.time))*0.33;
            surrodata_act.motor{Isub}.time   = all_data.motor.time;
            surrodata_act.motor{Isub}.dimord = 'chan_time';
            surrodata_act.motor{Isub}.label  = {'motor'};
            surrodata_act.motor{Isub}.cfg    = [];
        
        end
                        
    end
    
end %Isub
      
%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_5g');
cd(fullfile(path_source, 'Fig_5g'));

sourcedata = struct;

time = stimdata.frontal{1}.time;

% - GET CONTEXT CLASSIFICATION ACCURACY - %
for Isub = 1:numel(stimdata.frontal)
    
    if isempty(stimdata.frontal{Isub})
        
        sourcedata.stim(Isub,:)  = NaN(1, length(time));
        
    else
        
        sourcedata.stim(Isub,:)  = stimdata.frontal{Isub}.avg;

    end        
        
end%Isub

% - GET ACTION CLASSIFICATION ACCURACY - %
for Isub = 1:numel(actdata.frontal)
    
    if isempty(actdata.frontal{Isub})
        
        sourcedata.act(Isub,:)  = NaN(1, length(time));
        
    else
        
        sourcedata.act(Isub,:)  = actdata.frontal{Isub}.avg;

    end        
        
end%Isub

% - STORE CONTEXT DECODING - %
tmp = NaN(size(sourcedata.stim,1)+1, size(sourcedata.stim,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.stim;

table_frontal.stim = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5g_frontal_context_decoding.xlsx';
writetable(table_frontal.stim,filename, 'Sheet', 1);

% - STORE ACTION DECODING - %
tmp = NaN(size(sourcedata.act,1)+1, size(sourcedata.act,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.act;

table_frontal.act = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5g_frontal_action_decoding.xlsx';
writetable(table_frontal.act,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_5g');
cd(fullfile(path_source, 'Fig_5g'));

sourcedata = struct;

time = stimdata.motor{1}.time;

% - GET CONTEXT CLASSIFICATION ACCURACY - %
for Isub = 1:numel(stimdata.motor)
    
    if isempty(stimdata.motor{Isub})
        
        sourcedata.stim(Isub,:)  = NaN(1, length(time));
        
    else
        
        sourcedata.stim(Isub,:)  = stimdata.motor{Isub}.avg;

    end        
        
end%Isub

% - GET ACTION CLASSIFICATION ACCURACY - %
for Isub = 1:numel(actdata.motor)
    
    if isempty(actdata.motor{Isub})
        
        sourcedata.act(Isub,:)  = NaN(1, length(time));
        
    else
        
        sourcedata.act(Isub,:)  = actdata.motor{Isub}.avg;

    end        
        
end%Isub

% - STORE CONTEXT DECODING - %
tmp = NaN(size(sourcedata.stim,1)+1, size(sourcedata.stim,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.stim;

table_motor.stim = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5h_motor_context_decoding.xlsx';
writetable(table_motor.stim,filename, 'Sheet', 1);

% - STORE ACTION DECODING - %
tmp = NaN(size(sourcedata.act,1)+1, size(sourcedata.act,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.act;

table_motor.act = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5h_motor_action_decoding.xlsx';
writetable(table_motor.act,filename, 'Sheet', 1);

%% Remove Empty Cells

stimdata.frontal       = stimdata.frontal(~cellfun(@isempty, stimdata.frontal));
actdata.frontal        = actdata.frontal(~cellfun(@isempty, actdata.frontal));
surrodata_stim.frontal = surrodata_stim.frontal(~cellfun(@isempty, surrodata_stim.frontal));
surrodata_act.frontal  = surrodata_act.frontal(~cellfun(@isempty, surrodata_act.frontal));

stimdata.motor         = stimdata.motor(~cellfun(@isempty, stimdata.motor));
actdata.motor          = actdata.motor(~cellfun(@isempty, actdata.motor));
surrodata_stim.motor   = surrodata_stim.motor(~cellfun(@isempty, surrodata_stim.motor));
surrodata_act.motor    = surrodata_act.motor(~cellfun(@isempty, surrodata_act.motor));

%% Merge Data 

% select ROI
selroi = input('Select ROI: (1) Frontal, (2) Motor: ');

% select ROI for stimulus dimension

if selroi == 1 
    
    statsdata_stim  = stimdata.frontal;
    statssurro_stim = surrodata_stim.frontal;
    
    statsdata_act   = actdata.frontal;
    statssurro_act  = surrodata_act.frontal;

elseif selroi == 2
    
    statsdata_stim  = stimdata.motor;
    statssurro_stim = surrodata_stim.motor;
    
    statsdata_act   = actdata.motor;
    statssurro_act  = surrodata_act.motor;
            
end

    
%% Statistics Prediction Axis

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

subj = numel(statsdata_stim);
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

stat_stim = ft_timelockstatistics(cfg, statsdata_stim{:}, statssurro_stim{:});

%% Statistics Action Axis

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

subj = numel(statsdata_act);
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

stat_act = ft_timelockstatistics(cfg, statsdata_act{:}, statssurro_act{:});

%% Get Grand Average and Plot

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/colormaps/cbrewer/

% colors
cmap    = cbrewer('div', 'Spectral', 200, 'PCHIP');

cb       = cell(2,1);
cb{1}    = cmap(60,:);
cb{2}    = cmap(190,:);

% get grand average
cfg = [];
cfg.keepindividual = 'yes';

GA_Stim  = ft_timelockgrandaverage(cfg, statsdata_stim{:});
GA_Act   = ft_timelockgrandaverage(cfg, statsdata_act{:});

% plot data
time = GA_Stim.time;

a = shadedErrorBar(time, squeeze(mean(GA_Stim.individual)), ...
    std(squeeze(GA_Stim.individual)) / sqrt(size(GA_Stim.individual,1)));
a.mainLine.Color = cb{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{1}; a.patch.EdgeColor = cb{1}; a.patch.FaceAlpha = 0.5;

hold on

a = shadedErrorBar(time, squeeze(mean(GA_Act.individual)), ...
    std(squeeze(GA_Act.individual)) / sqrt(size(GA_Act.individual,1)));
a.mainLine.Color = cb{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{2}; a.patch.EdgeColor = cb{2}; a.patch.FaceAlpha = 0.5;

xlabel('Time to BR [s]');
ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 3, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off

% mask stats
if isfield(stat_stim, 'posclusters')

    ncluster = length(stat_stim.posclusters);

    for Icluster = 1:ncluster    

        if stat_stim.posclusters(Icluster).prob < 0.025

            idx = find(stat_stim.posclusterslabelmat == Icluster);

            % shade grey where significant
            sigpnts = stat_stim.time(idx);

            hold on

            scatter(sigpnts, ones(1, length(sigpnts))*0.29, 30, 'filled',...
                'MarkerFaceColor', cb{1}, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', cb{1}, 'MarkerEdgeAlpha', 0.2);
            clear sigpnts

        end

    end

end

if isfield(stat_act, 'posclusters')
    
    ncluster = length(stat_act.posclusters);

    for Icluster = 1:ncluster    

        if stat_act.posclusters(Icluster).prob < 0.025

            idx = find(stat_act.posclusterslabelmat == Icluster);

            % shade grey where significant
            sigpnts = stat_act.time(idx);

            hold on

            scatter(sigpnts, ones(1, length(sigpnts))*0.31, 30, 'filled',...
                'MarkerFaceColor', cb{2}, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', cb{2}, 'MarkerEdgeAlpha', 0.2);
            clear sigpnts

        end

    end

end

set(gcf, 'Position', [596 417 229 198]);

if selroi == 1
    print(gcf,'fig_5_frontal_DecAccuracy.pdf','-dpdf','-r400');    
else
    print(gcf,'fig_5_motor_DecAccuracy.pdf','-dpdf','-r400');    
end


