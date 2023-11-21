%% Figure 5a
% this script rungs additional analysis required for the plots of figure
% 5a.

%% - HOUSE KEEPING - %%

clear
close all;

ft_defaults;

%% -SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path.in  = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';
path.fun = '/Users/janweber/Documents/MATLAB/general_functions/';
path.out = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_5/revision_natcomm_final/';
path.source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath(path.fun); % add functions

%% - GENERAL SETTINGS - %%

nPad      = 2;     % zero paddings in seconds for PSD analysis

fsample   = 512;   % sampling frequency

%% - LOAD ELECTRODE LOCATION - %%

load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/Electrode_Information_ROI/elec_info_ROI.mat

%% - COMPUTE PC DIMENSION WITH STRONGEST THETA POWER - %%

% initialize peak components

peakcomp.frontal = NaN(numel(subjects),1);
peakcomp.motor   = NaN(numel(subjects),1);

for Isub = 1:numel(subjects)
    
    % subject ID
    subID = subjects{Isub};
    
    fprintf('%% ----------------------------------------------------------- %%\n');
    
    fprintf('%% Computation running for %s\n', subID);

    fprintf('%% ----------------------------------------------------------- %%\n');

    
    if isfile(fullfile(path.in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']))
        
        load(fullfile(path.in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']));

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
                
                peakcomp.frontal(Isub) = peakPC;
                                
            end
                        
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
                
                peakcomp.motor(Isub) = peakPC;
                                
            end
                        
        end
            
    end
    
end %Isub

%% - COMPUTE PEAK-TRIGGERED AVERAGE - %%

% preallocate memory for frontal and motor peak-locked data
ThetaLockedData.frontal = cell(numel(subjects),1);
ThetaLockedData.motor   = cell(numel(subjects),1);

CtxLockedData.frontal   = cell(numel(subjects),1);
CtxLockedData.motor     = cell(numel(subjects),1);

% store peak components for frontal and motor cortex for looping purpose
roinames = fieldnames(peakcomp);

for Isub = 1:numel(subjects)
    
    % get ID
    subID = subjects{Isub};
        
    % ---------------------------
    % Load Data
    % ---------------------------

    load(fullfile(path.in, subID, 'Data', 'HFA', ...
        [subID '_HFAavgoverbins.mat']));

    % load trial information table
    load(fullfile(path.in, subID, 'Data', [subID '_Trial_Info_clean.mat']));
    
    % ---------------------------
    % Extract ROI Electrodes
    % ---------------------------

    frontal_elec = [];
    if isfield(ROIs(Isub).frontal, 'elecidx') && length(ROIs(Isub).frontal.elecidx) > 1
        frontal_elec = ROIs(Isub).frontal.elecidx;
    end
    
    motor_elec = [];
    if isfield(ROIs(Isub).motor, 'elecidx') && length(ROIs(Isub).motor.elecidx) > 1
        motor_elec = ROIs(Isub).motor.elecidx;
    end
    
    mtl_elec = [];
    if isfield(ROIs(Isub).mtl, 'elecidx') && length(ROIs(Isub).mtl.elecidx) > 1
        mtl_elec = ROIs(Isub).mtl.elecidx;
    end
    
    % merge electrodes
    all_elecs = {frontal_elec, motor_elec};
     
    % ----------------------------
    % Loop over ROIs
    % ----------------------------
    
    % for Iroi = 1:numel(all_elecs)
    for Iroi = 1:numel(all_elecs)
        
        if ~isnan(peakcomp.(roinames{Iroi})(Isub))

            data = struct();

            if ~isempty(all_elecs{Iroi})

                fprintf('Running ID %d/%d, ROI %d/%d\n', Isub, numel(subjects), Iroi, numel(all_elecs));

                % ---------------------------------------
                % select electrodes for ROI
                % ---------------------------------------

                cfg = [];
                cfg.channel = hfa.label(all_elecs{Iroi});

                data        = ft_selectdata(cfg, hfa);

                % ---------------------------------------
                % Redefine trial
                % ---------------------------------------

                % get index for stop trials to remove
                IdxStopTrl = find(trl_mod_table.TrlType == 2 | isnan(trl_mod_table.TrlType)); % trialtype 2 = stop trial

                % remove stop trials from table
                tmptable               = trl_mod_table;
                tmptable(IdxStopTrl,:) = [];

                % get trials

                trl = NaN(size(data.sampleinfo,1),3);

                for triali = 1:size(data.trialinfo,1)

                    curr_onset = data.sampleinfo(triali,1); % current onset = baseline
                    trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonRelease(triali) - fsample*2 ...
                        curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonRelease(triali) + fsample*2 -fsample*2];

                end    

                trl(IdxStopTrl,:,:)    = [];

                cfg = [];
                cfg.trl = trl;

                data    = ft_redefinetrial(cfg, data);

                % ---------------------------------------
                % Perform PCA on Data
                % ---------------------------------------

                [pcadata, pcainfo] = rh_pca(data.trial);

                % time window for peak detection
                twin = nearest(data.time, -0.5):nearest(data.time, 0.3);

                % ---------------------------------------
                % Select PC with strongest theta power
                % ---------------------------------------

                cfg = [];
                cfg.channel = data.label{peakcomp.(roinames{Iroi})};

                data2use    = ft_selectdata(cfg, data);

                % remove trial information
                data2use = rmfield(data2use, 'trialinfo');

                % initialize structure for new trials according to PC peaks
                newtrl = [];

                for Itrial = 1:size(data2use.trial,1)

                    % find peaks per trial
                    [~,timeidx] = findpeaks(squeeze(data2use.trial(Itrial,1,twin)));

                    % add the samples of the time window
                    timeidx = timeidx + twin(1);

                    % get first sample of trial X
                    firstsample = data.sampleinfo(Itrial,1);

                    % concatenate trials (locked to the peak)
                    newtrl = vertcat(newtrl, ...
                        [firstsample+timeidx - round(fsample/2), ...
                         firstsample+timeidx + round(fsample/2), ...
                         repmat(-round(fsample/2), length(timeidx), 1)]);
                    
                end %Itrial

                % ---------------------------------------
                % Redefine trials according to peak
                % ---------------------------------------

                cfg     = [];
                cfg.trl = newtrl;
                 
                ThetaLockedData.(roinames{Iroi}){Isub} = ft_redefinetrial(cfg, data2use);
                
                ThetaLockedData.(roinames{Iroi}){Isub} = ft_timelockanalysis([], ThetaLockedData.(roinames{Iroi}){Isub});
                
                ThetaLockedData.(roinames{Iroi}){Isub}.label = roinames(Iroi);
                
                clear data2use
                
            end
            
        end
        
    end %Iroi
    
end %Isub
     
%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path.source, 'Fig_6a');
cd(fullfile(path.source, 'Fig_6a'));

time = ThetaLockedData.frontal{1}.time;

sourcedata = NaN(numel(subjects), length(time));

for Isub = 1:numel(ThetaLockedData.frontal)
    
    if isempty(ThetaLockedData.frontal{Isub})
        
        sourcedata(Isub,:)  = NaN(1, length(time));

    else
        
        sourcedata(Isub,:)  = ThetaLockedData.frontal{Isub}.avg;

    end        
        
end%Isub

% - STORE DATA - %
tmp = NaN(size(sourcedata,1)+1, size(sourcedata,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata;

table_frontal = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_6a_frontal_PTA.xlsx';
writetable(table_frontal,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path.source, 'Fig_6a');
cd(fullfile(path.source, 'Fig_6a'));

time = ThetaLockedData.motor{1}.time;

sourcedata = NaN(numel(subjects), length(time));

for Isub = 1:numel(ThetaLockedData.motor)
    
    if isempty(ThetaLockedData.motor{Isub})
        
        sourcedata(Isub,:)  = NaN(1, length(time));

    else
        
        sourcedata(Isub,:)  = ThetaLockedData.motor{Isub}.avg;

    end        
        
end%Isub

% - STORE DATA - %
tmp = NaN(size(sourcedata,1)+1, size(sourcedata,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata;

table_motor = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_6a_motor_PTA.xlsx';
writetable(table_motor,filename, 'Sheet', 1);

%% - REMOVE EMPTY CELLS - %%

ThetaLockedData.(roinames{1})  = ThetaLockedData.(roinames{1})(~cellfun(@isempty, ThetaLockedData.(roinames{1})));
ThetaLockedData.(roinames{2})  = ThetaLockedData.(roinames{2})(~cellfun(@isempty, ThetaLockedData.(roinames{2})));

%% - GRAND AVERAGE AND PLOTTING - %%

cfg = [];
cfg.keepindividual = 'yes';

GA_frontal.theta   = ft_timelockgrandaverage(cfg, ThetaLockedData.frontal{:});
GA_motor.theta     = ft_timelockgrandaverage(cfg, ThetaLockedData.motor{:});

% de-meaning
GA_frontal.theta.individual(:,1,:) = squeeze(GA_frontal.theta.individual(:,1,:)) - ...
    squeeze(mean(GA_frontal.theta.individual(:,1,:), 'all'));
GA_motor.theta.individual(:,1,:)   = squeeze(GA_motor.theta.individual(:,1,:)) - ...
    squeeze(mean(GA_motor.theta.individual(:,1,:), 'all'));

% plot
figure;
a = shadedErrorBar(GA_frontal.theta.time, squeeze(mean(GA_frontal.theta.individual)), ...
    std(squeeze(GA_frontal.theta.individual)) / sqrt(size(GA_frontal.theta.individual,1)));
a.mainLine.Color = 'r'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'r'; a.patch.EdgeColor = 'r'; a.patch.FaceAlpha = 0.5;

ylabel('Activity');
xlabel('Time from Peak');

set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [GA_motor.theta.time(1) GA_motor.theta.time(end)]); 
box off

set(gcf, 'Position', [488 381 331 252]);

hold on

% sine fit
time = GA_frontal.theta.time;
sineparams = sineFit(time,squeeze(nanmean(GA_frontal.theta.individual))', 0);
plot(time,sineparams(1)+sineparams(2)*sin(2*pi*sineparams(3)*time+sineparams(4)), 'k', 'LineWidth', 1);

xlim([-0.3 0.3]);

cd(path.out);
print(gcf,'fig_6_thetadim_frontal_pta.pdf','-dpdf','-r400');

figure;
a = shadedErrorBar(GA_motor.theta.time, squeeze(mean(GA_motor.theta.individual)), ...
    std(squeeze(GA_motor.theta.individual)) / sqrt(size(GA_motor.theta.individual,1)));
a.mainLine.Color = 'b'; a.mainLine.LineWidth = 3;
a.patch.FaceColor = 'b'; a.patch.EdgeColor =' b'; a.patch.FaceAlpha = 0.5;

ylabel('Activity');
xlabel('Time from Peak');

set(gca, 'fontsize', 13, 'linewidth', 3, 'xlim', [GA_motor.theta.time(1) GA_motor.theta.time(end)]); 
box off

set(gcf, 'Position', [488 381 331 252]);

hold on

% sine fit
sineparams = sineFit(time,squeeze(nanmean(GA_frontal.theta.individual))', 0);
plot(time,sineparams(1)+sineparams(2)*sin(2*pi*sineparams(3)*time+sineparams(4)), 'k', 'LineWidth', 1);

xlim([-0.3 0.3]);

cd(path.out);
print(gcf,'fig_6_thetadim_motor_pta.pdf','-dpdf','-r400');

            
            
            
            
            
            