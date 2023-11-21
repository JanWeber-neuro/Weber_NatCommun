%% Figure S9

%%

clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_data = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_S7/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

%% Determine whether decoding accuracy based on all or only task-informative electrodes should be plotted

selcond = input('Locked to (1) HLL, (2) BR: ');

if selcond == 1
    data2use = '_stimcoding_all_elecs_HLL_w_surro.mat';
elseif selcond == 2
    data2use = '_stimcoding_all_elecs_BR_w_surro.mat';
end

%%

fsample   = 512;

smoothwin = 0.035; % window for smoothing

%%

stimdata       = struct();
actdata        = struct();
surrodata_stim = struct();
surrodata_act  = struct();

for Isub = 1:numel(subjects)
    
    subID = subjects{Isub};
    
    if isfile(fullfile(path_data, subID, 'Data', 'State_Space', [subID data2use]))
        
        load(fullfile(path_data, subID, 'Data', 'State_Space', [subID data2use]));
        
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
        
%% Remove Empty Cells

stimdata.frontal       = stimdata.frontal(~cellfun(@isempty, stimdata.frontal));
actdata.frontal        = actdata.frontal(~cellfun(@isempty, actdata.frontal));
surrodata_stim.frontal = surrodata_stim.frontal(~cellfun(@isempty, surrodata_stim.frontal));
surrodata_act.frontal  = surrodata_act.frontal(~cellfun(@isempty, surrodata_act.frontal));

stimdata.motor         = stimdata.motor(~cellfun(@isempty, stimdata.motor));
actdata.motor          = actdata.motor(~cellfun(@isempty, actdata.motor));
surrodata_stim.motor   = surrodata_stim.motor(~cellfun(@isempty, surrodata_stim.motor));
surrodata_act.motor    = surrodata_act.motor(~cellfun(@isempty, surrodata_act.motor));

%% Save source data PFC action decoding

% go into results folder
mkdir(path_source, 'Fig_S7');
cd(fullfile(path_source, 'Fig_S7'));

sourcedata = struct;

time = actdata.frontal{1}.time;

for Isub = 1:numel(actdata.frontal)
    
    if isempty(actdata.frontal{Isub})
        
        sourcedata.pfc(Isub,:)  = NaN(1, length(time));

    else
        
        sourcedata.pfc(Isub,:)  = actdata.frontal{Isub}.avg(1,:);

    end        
        
end%Isub

tmp = NaN(size(sourcedata.pfc,1)+1, size(sourcedata.pfc,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata.pfc,1))';
tmp(2:end,2:end)  = sourcedata.pfc;

table_frontal = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S7_frontal_action.xlsx';
writetable(table_frontal,filename, 'Sheet', 1);

%% Save source data Motor action decoding

% go into results folder
mkdir(path_source, 'Fig_S7');
cd(fullfile(path_source, 'Fig_S7'));

sourcedata = struct;

time = actdata.motor{1}.time;

for Isub = 1:numel(actdata.motor)
    
    if isempty(actdata.motor{Isub})
        
        sourcedata.m1(Isub,:)  = NaN(1, length(time));

    else
        
        sourcedata.m1(Isub,:)  = actdata.motor{Isub}.avg(1,:);

    end        
        
end%Isub

tmp = NaN(size(sourcedata.m1,1)+1, size(sourcedata.m1,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata.m1,1))';
tmp(2:end,2:end)  = sourcedata.m1;

table_motor = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S7_motor_action.xlsx';
writetable(table_motor,filename, 'Sheet', 1);

%% Merge Data 

% select ROI
sel_stim = input('Statistic for Stimulus Dimension: (1) Frontal, (2) Motor: ');
sel_act  = input('Statistic for Action Dimension: (1) Frontal, (2) Motor: ');

% select ROI for stimulus dimension

if sel_stim == 1 
    
    statsdata_stim  = stimdata.frontal;
    statssurro_stim = surrodata_stim.frontal;

elseif sel_stim == 2
    
    statsdata_stim  = stimdata.motor;
    statssurro_stim = surrodata_stim.motor;
        
end

% select ROI for action dimension

if sel_act == 1 
    
    statsdata_act   = actdata.frontal;
    statssurro_act  = surrodata_act.frontal;
    
elseif sel_act == 2
    
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

% compute cohen's d as effect size
cfg            = [];
cfg.parameter  = 'avg';
cfg.method     = 'analytic';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

d_stim         = ft_timelockstatistics(cfg, statsdata_stim{:}, statssurro_stim{:});


% check positive clusters

if isfield(stat_stim, 'posclusters')

    for Icluster = 1:numel(stat_stim.posclusters)

        fprintf('%%---------------------------------------------------\n');
        fprintf('%% POSITIVE CLUSTERS---------------------------------\n');
        fprintf('%%---------------------------------------------------\n');

        fprintf('\n');

        fprintf('p-value cluster %d = %.3f\n', Icluster, stat_stim.posclusters(Icluster).prob);
        fprintf('clustermass cluster %d = %.3f\n', Icluster, stat_stim.posclusters(Icluster).clusterstat);

        begeffect = stat_stim.time...
            (find(stat_stim.posclusterslabelmat == Icluster, 1, 'first'));
        endeffect = stat_stim.time...
            (find(stat_stim.posclusterslabelmat == Icluster, 1, 'last'));

        fprintf('Stimulus Decoding Effect: cluster %d from %.3f to %.3f ms\n', Icluster, begeffect, endeffect);

        fprintf('Cohens d cluster %d = %.3f\n', Icluster, mean(d_stim.cohensd(...
            find(stat_stim.posclusterslabelmat == Icluster, 1, 'first'):...
            find(stat_stim.posclusterslabelmat == Icluster, 1, 'last'))));

        fprintf('\n');

    end

end


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

% compute cohen's d as effect size
cfg            = [];
cfg.parameter  = 'avg';
cfg.method     = 'analytic';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

d_act         = ft_timelockstatistics(cfg, statsdata_act{:}, statssurro_act{:});

% check positive clusters

if isfield(stat_act, 'posclusters')

    for Icluster = 1:numel(stat_act.posclusters)

        fprintf('%%---------------------------------------------------\n');
        fprintf('%% POSITIVE CLUSTERS---------------------------------\n');
        fprintf('%%---------------------------------------------------\n');

        fprintf('\n');

        fprintf('p-value cluster %d = %.3f\n', Icluster, stat_act.posclusters(Icluster).prob);
        fprintf('clustermass cluster %d = %.3f\n', Icluster, stat_act.posclusters(Icluster).clusterstat);

        begeffect = stat_act.time...
            (find(stat_act.posclusterslabelmat == Icluster, 1, 'first'));
        endeffect = stat_act.time...
            (find(stat_act.posclusterslabelmat == Icluster, 1, 'last'));

        fprintf('Action Decoding Effect: cluster %d from %.3f to %.3f ms\n', Icluster, begeffect, endeffect);

        fprintf('Cohens d cluster %d = %.3f\n', Icluster, nanmean(d_act.cohensd(...
            find(stat_act.posclusterslabelmat == Icluster, 1, 'first'):...
            find(stat_act.posclusterslabelmat == Icluster, 1, 'last'))));

        fprintf('\n');

    end

end

%% Get Grand Average and Plot

% colors
addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/colormaps/cbrewer/

cmap    = cbrewer('div', 'Spectral', 200, 'PCHIP');

cb       = cell(2,1);
if sel_act == 1
    cb{1} = 'r';
else
    cb{1} = 'b';
end

% cb{2}    = cmap(190,:);

% get grand average
cfg = [];
cfg.keepindividual = 'yes';

GA_Stim  = ft_timelockgrandaverage(cfg, statsdata_stim{:});
GA_Act   = ft_timelockgrandaverage(cfg, statsdata_act{:});

% plot data
time = GA_Stim.time;

% a = shadedErrorBar(time, squeeze(mean(GA_Stim.individual)), ...
%     std(squeeze(GA_Stim.individual)) / sqrt(size(GA_Stim.individual,1)));
% a.mainLine.Color = cb{1}; a.mainLine.LineWidth = 3;
% a.patch.FaceColor = cb{1}; a.patch.EdgeColor = cb{1}; a.patch.FaceAlpha = 0.5;
% 
% hold on

a = shadedErrorBar(time, squeeze(mean(GA_Act.individual)), ...
    std(squeeze(GA_Act.individual)) / sqrt(size(GA_Act.individual,1)));
a.mainLine.Color = cb{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = cb{1}; a.patch.EdgeColor = cb{1}; a.patch.FaceAlpha = 0.5;

if selcond == 1
    xlabel('Time to HLL [s]');
else
    xlabel('Time to BR [s]');

end
ylabel('Classifier Accuracy');
yline(0.33, 'LineWidth', 1, 'LineStyle', '--');
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
box off

% % mask stats
% if isfield(stat_stim, 'posclusters')
% 
%     ncluster = length(stat_stim.posclusters);
% 
%     for Icluster = 1:ncluster    
% 
%         if stat_stim.posclusters(Icluster).prob < 0.025
% 
%             idx = find(stat_stim.posclusterslabelmat == Icluster);
% 
%             % shade grey where significant
%             sigpnts = stat_stim.time(idx);
% 
%             hold on
% 
%             scatter(sigpnts, ones(1, length(sigpnts))*0.29, 30, 'filled',...
%                 'MarkerFaceColor', cb{1}, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', cb{1}, 'MarkerEdgeAlpha', 0.2);
%             clear sigpnts
% 
%         end
% 
%     end
% 
% end

if isfield(stat_act, 'posclusters')
    
    ncluster = length(stat_act.posclusters);

    for Icluster = 1:ncluster    

        if stat_act.posclusters(Icluster).prob < 0.025

            idx = find(stat_act.posclusterslabelmat == Icluster);

            % shade grey where significant
            sigpnts = stat_act.time(idx);

            hold on

            scatter(sigpnts, ones(1, length(sigpnts))*0.31, 30, 'filled',...
                'MarkerFaceColor', cb{1}, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', cb{1}, 'MarkerEdgeAlpha', 0.2);
            clear sigpnts

        end

    end

end

set(gcf, 'Position', [357 341 433 275]);

cd(path_out);

if sel_act == 1
    print(gcf,'Decoding_Action_PFC_Locked_HLL.pdf','-dpdf','-r400');  
else
    print(gcf,'Decoding_Action_M1_Locked_HLL.pdf','-dpdf','-r400');  
end


