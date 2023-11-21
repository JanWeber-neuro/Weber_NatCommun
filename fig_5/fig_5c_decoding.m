%% - Figure 5c
% this script rungs additional analysis required for the plots of figure
% 5c.

%%

clear
close all;
ft_defaults;

%% -SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_5/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% - GENERAL SETTINGS - %%

nPad      = 2;     % zero paddings in seconds for PSD analysis

fsample   = 512;   % sampling frequency

smoothwin = 0.025; % window for smoothing decoding timeseries

%% - LOAD DATA - %%

% preallocate data for data with theta dimension
ThetaDim  = struct();

% preallocate data for data with action dimension
ActionDim = struct();

% preallocate memory for 1/f corrected PSDs
PSD       = struct();

for Isub = 1:numel(subjects)
    
    % subject ID
    subID = subjects{Isub};
    
    fprintf('%% ----------------------------------------------------------- %%\n');
    
    fprintf('%% Computation running for %s\n', subID);

    fprintf('%% ----------------------------------------------------------- %%\n');

    
    if isfile(fullfile(path_in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']))
        
        load(fullfile(path_in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']));

        % #########################
        % FRONTAL CORTEX
        % #########################
        
        if isfield(all_data, 'frontal')
            
                                    
            % ------------------------- 
            % Bring data into FT format
            % -------------------------

            ftdata = struct();
            ftdata.trial  = all_data.frontal.trial;
            ftdata.label  = sprintfc('PC%d',1:size(ftdata.trial,2));
            ftdata.dimord = all_data.frontal.dimord; % dimensions
            ftdata.time   = all_data.frontal.time;
            ftdata.cfg    = [];
            
            % ------------------------- 
            % Compute IRASA
            % ------------------------- 
            
            PadData = struct();              

            for chani = 1:numel(ftdata.label)

                % zero padding
                tmptrial  = squeeze(ftdata.trial(:,chani,:));

                % make FT structure
                PadData.trial(:,chani,:)   = [zeros(size(tmptrial,1),fsample*nPad) tmptrial zeros(size(tmptrial,1),fsample*nPad)];
                PadData.time               = 0:1/fsample:1/fsample*(size(PadData.trial,3)-1);
                PadData.dimord             = 'rpt_chan_time';
                PadData.label(chani)       = ftdata.label(chani);

            end

            % actual analysis

            cfg               = [];
            cfg.foilim        = [2 13];
            cfg.channel       = PadData.label;
            cfg.taper         = 'hanning';
            cfg.pad           = 'nextpow2';
            cfg.keeptrials    = 'yes';
            cfg.method        = 'irasa';
            powspec.frac      = ft_freqanalysis(cfg, PadData); % fractal part
            cfg.method        = 'mtmfft';
            powspec.orig      = ft_freqanalysis(cfg, PadData); % original spectrum

            clear PadData

            % correct for 1/f by subtracting fractal from original spectrum
            cfg               = [];
            cfg.parameter     = 'powspctrm';
            cfg.operation     = 'x2-x1';
            powspec.osci      = ft_math(cfg, powspec.frac, powspec.orig);
            powspec.osci. ...
                powspctrm(powspec.osci.powspctrm < 0) = 0;

            % ------------------------- 
            % Identify Peak 
            % ------------------------- 
            
            % loop over PCs
            
            peakamp = NaN(size(powspec.osci.powspctrm,2),1);
            
            for Icomp = 1:size(powspec.osci.powspctrm,2)
            
                % find peaks
                [pks,locs] = findpeaks(squeeze(mean(powspec.osci.powspctrm(:,Icomp,:))));
                
                peakamp(Icomp) = max(pks);
                
            end
            
            % find PC with max peak
            [~, peakPC] = max(peakamp); clear peakamp
            
            
            % admin
            if ~isempty(peakPC)
                
                all_data.frontal.peakdim = peakPC;
                
            else
                
                all_data.frontal.peakdim = [];
                
            end
            
            % pass data
            ThetaDim.frontal{Isub} = all_data.frontal;
            
            PSD.frontal{Isub}      = powspec.osci; clear powspec
            
        end
        
        % #########################
        % MOTOR CORTEX
        % #########################
        
        if isfield(all_data, 'motor')
                        
            % ------------------------- 
            % Bring data into FT format
            % -------------------------

            ftdata = struct();
            ftdata.trial  = all_data.motor.trial;
            ftdata.label  = sprintfc('PC%d',1:size(ftdata.trial,2));
            ftdata.dimord = all_data.motor.dimord; % dimensions
            ftdata.time   = all_data.motor.time;
            ftdata.cfg    = [];
            
            % ------------------------- 
            % Compute IRASA
            % ------------------------- 
            
            PadData = struct();              

            for chani = 1:numel(ftdata.label)

                % zero padding
                tmptrial  = squeeze(ftdata.trial(:,chani,:));

                % make FT structure
                PadData.trial(:,chani,:)   = [zeros(size(tmptrial,1),fsample*nPad) tmptrial zeros(size(tmptrial,1),fsample*nPad)];
                PadData.time               = 0:1/fsample:1/fsample*(size(PadData.trial,3)-1);
                PadData.dimord             = 'rpt_chan_time';
                PadData.label(chani)       = ftdata.label(chani);

            end

            % actual analysis

            cfg               = [];
            cfg.foilim        = [2 13];
            cfg.channel       = PadData.label;
            cfg.taper         = 'hanning';
            cfg.pad           = 'nextpow2';
            cfg.keeptrials    = 'yes';
            cfg.method        = 'irasa';
            powspec.frac      = ft_freqanalysis(cfg, PadData); % fractal part
            cfg.method        = 'mtmfft';
            powspec.orig      = ft_freqanalysis(cfg, PadData); % original spectrum

            clear PadData

            % correct for 1/f by subtracting fractal from original spectrum
            cfg               = [];
            cfg.parameter     = 'powspctrm';
            cfg.operation     = 'x2-x1';
            powspec.osci      = ft_math(cfg, powspec.frac, powspec.orig);
            powspec.osci. ...
                powspctrm(powspec.osci.powspctrm < 0) = 0;

            % ------------------------- 
            % Identify Peak 
            % ------------------------- 
            
            % loop over PCs
            
            peakamp = NaN(size(powspec.osci.powspctrm,2),1);
            
            for Icomp = 1:size(powspec.osci.powspctrm,2)
            
                % find peaks
                [pks,locs] = findpeaks(squeeze(mean(powspec.osci.powspctrm(:,Icomp,:))));
                
                peakamp(Icomp) = max(pks);
                
            end
            
            % find PC with max peak
            [~, peakPC] = max(peakamp); clear peakamp
            
            
            % admin
            if ~isempty(peakPC)
                
                all_data.motor.peakdim = peakPC;
                
            else
                
                all_data.motor.peakdim = [];
                
            end
            
            % pass data
            ThetaDim.motor{Isub} = all_data.motor;
            
            PSD.motor{Isub}      = powspec.osci; clear powspec

            clear all_data
            
        end
            
    end
    
end %Isub
           
% %% - REMOVE EMPTY CELLS - %%
% 
% ThetaDim.frontal  = ThetaDim.frontal(~cellfun(@isempty, ThetaDim.frontal));
% ThetaDim.motor    = ThetaDim.motor(~cellfun(@isempty, ThetaDim.motor));
% 
% PSD.frontal       = PSD.frontal(~cellfun(@isempty, PSD.frontal));
% PSD.motor         = PSD.motor(~cellfun(@isempty, PSD.motor));

%% - COMPUTE DECODING ACCURACY - %%

stimdata       = struct(); % stimulus data
actdata        = struct(); % action data
surrodata_stim = struct(); % respective surrogates
surrodata_act  = struct();

% #########################
% FRONTAL CORTEX
% #########################

for Isub = 1:numel(ThetaDim.frontal)

    % ------------------------- 
    % Get Stimulus Data
    % ------------------------- 
    
    if ~isempty(ThetaDim.frontal{Isub}) && ~isempty(ThetaDim.frontal{Isub}.stimdim)
        
        stimdata.frontal{Isub}.avg          = smoothdata(ThetaDim.frontal{Isub}.DecAcc_stim...
            (ThetaDim.frontal{Isub}.peakdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
        stimdata.frontal{Isub}.time         = ThetaDim.frontal{Isub}.time;
        stimdata.frontal{Isub}.dimord       = 'chan_time';
        stimdata.frontal{Isub}.label        = {'Frontal'};
        stimdata.frontal{Isub}.cfg          = [];

        surrodata_stim.frontal{Isub}.avg    = ones(1, length(ThetaDim.frontal{Isub}.time))*0.33;
        surrodata_stim.frontal{Isub}.time   = ThetaDim.frontal{Isub}.time;
        surrodata_stim.frontal{Isub}.dimord = 'chan_time';
        surrodata_stim.frontal{Isub}.label  = {'Frontal'};
        surrodata_stim.frontal{Isub}.cfg    = [];
        
    end

    % ------------------------- 
    % Get Action Data
    % ------------------------- 
    
    if ~isempty(ThetaDim.frontal{Isub}) && ~isempty(ThetaDim.frontal{Isub}.actdim)
        
        actdata.frontal{Isub}.avg          = smoothdata(ThetaDim.frontal{Isub}.DecAcc_act...
            (ThetaDim.frontal{Isub}.peakdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
        actdata.frontal{Isub}.time         = ThetaDim.frontal{Isub}.time;
        actdata.frontal{Isub}.dimord       = 'chan_time';
        actdata.frontal{Isub}.label        = {'Frontal'};
        actdata.frontal{Isub}.cfg          = [];

        surrodata_act.frontal{Isub}.avg    = ones(1, length(ThetaDim.frontal{Isub}.time))*0.33;
        surrodata_act.frontal{Isub}.time   = ThetaDim.frontal{Isub}.time;
        surrodata_act.frontal{Isub}.dimord = 'chan_time';
        surrodata_act.frontal{Isub}.label  = {'Frontal'};
        surrodata_act.frontal{Isub}.cfg    = [];

    end
       
end
     

% #########################
% MOTOR CORTEX
% #########################

for Isub = 1:numel(ThetaDim.motor)

    % ------------------------- 
    % Get Stimulus Data
    % ------------------------- 
    
    if ~isempty(ThetaDim.motor{Isub}) && ~isempty(ThetaDim.motor{Isub}.stimdim)
        
        stimdata.motor{Isub}.avg          = smoothdata(ThetaDim.motor{Isub}.DecAcc_stim...
            (ThetaDim.motor{Isub}.peakdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
        stimdata.motor{Isub}.time         = ThetaDim.motor{Isub}.time;
        stimdata.motor{Isub}.dimord       = 'chan_time';
        stimdata.motor{Isub}.label        = {'motor'};
        stimdata.motor{Isub}.cfg          = [];

        surrodata_stim.motor{Isub}.avg    = ones(1, length(ThetaDim.motor{Isub}.time))*0.33;
        surrodata_stim.motor{Isub}.time   = ThetaDim.motor{Isub}.time;
        surrodata_stim.motor{Isub}.dimord = 'chan_time';
        surrodata_stim.motor{Isub}.label  = {'motor'};
        surrodata_stim.motor{Isub}.cfg    = [];
        
    end

    % ------------------------- 
    % Get Action Data
    % ------------------------- 
    
    if ~isempty(ThetaDim.motor{Isub}) && ~isempty(ThetaDim.motor{Isub}.actdim)
        
        actdata.motor{Isub}.avg          = smoothdata(ThetaDim.motor{Isub}.DecAcc_act...
            (ThetaDim.motor{Isub}.peakdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
        actdata.motor{Isub}.time         = ThetaDim.motor{Isub}.time;
        actdata.motor{Isub}.dimord       = 'chan_time';
        actdata.motor{Isub}.label        = {'motor'};
        actdata.motor{Isub}.cfg          = [];

        surrodata_act.motor{Isub}.avg    = ones(1, length(ThetaDim.motor{Isub}.time))*0.33;
        surrodata_act.motor{Isub}.time   = ThetaDim.motor{Isub}.time;
        surrodata_act.motor{Isub}.dimord = 'chan_time';
        surrodata_act.motor{Isub}.label  = {'motor'};
        surrodata_act.motor{Isub}.cfg    = [];

    end
       
end %Isub

%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_6c');
cd(fullfile(path_source, 'Fig_6c'));

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
filename = 'fig_6c_frontal_context_decoding.xlsx';
writetable(table_frontal.stim,filename, 'Sheet', 1);

% - STORE ACTION DECODING - %
tmp = NaN(size(sourcedata.act,1)+1, size(sourcedata.act,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.act;

table_frontal.act = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_6c_frontal_action_decoding.xlsx';
writetable(table_frontal.act,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_6c');
cd(fullfile(path_source, 'Fig_6c'));

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
filename = 'fig_6c_motor_context_decoding.xlsx';
writetable(table_motor.stim,filename, 'Sheet', 1);

% - STORE ACTION DECODING - %
tmp = NaN(size(sourcedata.act,1)+1, size(sourcedata.act,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.act;

table_motor.act = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_6c_motor_action_decoding.xlsx';
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
set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [time(1) time(end)], 'ylim', [0.22 0.47]); 
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

            scatter(sigpnts, ones(1, length(sigpnts))*0.24, 30, 'filled',...
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

            scatter(sigpnts, ones(1, length(sigpnts))*0.24, 30, 'filled',...
                'MarkerFaceColor', cb{2}, 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', cb{2}, 'MarkerEdgeAlpha', 0.2);
            clear sigpnts

        end

    end

end

set(gcf, 'Position', [596 417 229 198]);

cd(path_out);
if sel_stim == 1 && sel_act == 1
    print(gcf,'fig_6_frontal_ThetaDim_DecAccuracy.pdf','-dpdf','-r400');    
elseif sel_stim == 2 && sel_act == 2
    print(gcf,'fig_6_motor_ThetaDim_DecAccuracy.pdf','-dpdf','-r400');    
end

            
