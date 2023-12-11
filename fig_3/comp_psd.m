%% PSD on Raw Data during HLL and BR

% ########################################################################

% Steps:
% 1. Load Raw Data
% 2. Perform PSD using IRASA
% 3. Save PSD for each subject

% ########################################################################

%%
clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_elec   = '/Volumes/IMPECOG/elecs_and_info';

addpath('/Volumes/IMPECOG/iEEG_Data_JW/functions/');

%% Input for PSD

fsample = 512; % sampling freq
nPad    = 1;   % zero padding to increase freq res (in sec.)

%% 

for Isub = 1:numel(subjects)
    
    % get subject ID
    subID = subjects{Isub};
    
    fprintf('PSD for %s\n', subID);
    
    % go into folder and get clean data
    load(fullfile(path_in, subID, 'Data', [subID '_clean.mat']));
        
    % load trial information
    load(fullfile(path_in, subID, 'Data', [subID '_Trial_Info_clean.mat']));
    
    %% Get Index for Stop Trials
   
    % get index for stop trials to remove
    IdxStopTrl = find(trl_mod_table.TrlType == 2 | isnan(trl_mod_table.TrlType)); % trialtype 2 = stop trial

    % remove stop trials from table
    tmptable               = trl_mod_table;
    tmptable(IdxStopTrl,:) = [];
   
    %% Segment Data from Start -0.5 to +0.3 sec. relative to HLL (plus some zero padding to increase freq resolution)

    trl = NaN(size(data_clean.sampleinfo,1),3);

    % segmentation: from -0.5 pre HLL to 0.3 sec post-HLL
    for triali = 1:size(data_clean.trialinfo,1)

        curr_onset = data_clean.sampleinfo(triali,1); % current onset = baseline
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.HitLowLim(triali) - fsample/2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.HitLowLim(triali) + round(fsample * 0.3) -fsample/2];

    end
        
    trl(IdxStopTrl,:,:)    = [];

    % redefine trial to probe onset
    cfg          = [];
    cfg.trl      = trl;
    taskdata_hll = ft_redefinetrial(cfg, data_clean);
    
    clear trl
    
    %% Segment Data from Start -2 to +2 sec. relative to BR

    trl = NaN(size(data_clean.sampleinfo,1),3);

    % segmentation: from -0.5 pre BL to 0.3 sec post-BL
    for triali = 1:size(data_clean.trialinfo,1)

        curr_onset = data_clean.sampleinfo(triali,1); % current onset = baseline  
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonRelease(triali) - fsample/2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonRelease(triali) + round(fsample * 0.3) -fsample/2];

    end 
        
    trl(IdxStopTrl,:,:)    = [];

    % redefine trial to probe onset
    cfg          = [];
    cfg.trl      = trl;
    taskdata_br  = ft_redefinetrial(cfg, data_clean);
    
    clear trl
    
    %% Perform PSD on HLL
    
    % Compute PSD (1/f corrected using IRASA)

    PadData = struct();              

    for chani = 1:numel(taskdata_hll.label)

        % zero padding
        tmptrial  = squeeze(taskdata_hll.trial(:,chani,:));

        % make FT structure
        PadData.trial(:,chani,:)   = [zeros(size(tmptrial,1),fsample*nPad) tmptrial zeros(size(tmptrial,1),fsample*nPad)];
        PadData.time               = 0:1/fsample:1/fsample*(size(PadData.trial,3)-1);
        PadData.dimord             = 'rpt_chan_time';
        PadData.label(chani)       = taskdata_hll.label(chani);

    end

    % actual analysis

    cfg               = [];
    cfg.foilim        = [1 80];
    cfg.channel       = PadData.label;
    cfg.taper         = 'hanning';
    cfg.pad           = 'nextpow2';
    cfg.keeptrials    = 'yes';
    cfg.method        = 'irasa';
    powspec.hll.frac  = ft_freqanalysis(cfg, PadData); % fractal part
    cfg.method        = 'mtmfft';
    powspec.hll.orig  = ft_freqanalysis(cfg, PadData); % original spectrum

    clear PadData
    
    % correct for 1/f by subtracting fractal from original spectrum
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2-x1';
    powspec.hll.osci  = ft_math(cfg, powspec.hll.frac, powspec.hll.orig);
    powspec.hll.osci.powspctrm(powspec.hll.osci.powspctrm < 0) = 0;
        
    %% Perform PSD on BR
    
    % Compute PSD (1/f corrected using IRASA)

    PadData = struct();              

    for chani = 1:numel(taskdata_br.label)

        % zero padding
        tmptrial  = squeeze(taskdata_br.trial(:,chani,:));

        % make FT structure
        PadData.trial(:,chani,:)   = [zeros(size(tmptrial,1),fsample*nPad) tmptrial zeros(size(tmptrial,1),fsample*nPad)];
        PadData.time               = 0:1/fsample:1/fsample*(size(PadData.trial,3)-1);
        PadData.dimord             = 'rpt_chan_time';
        PadData.label(chani)       = taskdata_br.label(chani);

    end

    % actual analysis

    cfg               = [];
    cfg.foilim        = [1 80];
    cfg.channel       = PadData.label;
    cfg.taper         = 'hanning';
    cfg.pad           = 'nextpow2';
    cfg.keeptrials    = 'yes';
    cfg.method        = 'irasa';
    powspec.br.frac  = ft_freqanalysis(cfg, PadData); % fractal part
    cfg.method        = 'mtmfft';
    powspec.br.orig  = ft_freqanalysis(cfg, PadData); % original spectrum

    clear PadData
    
    % correct for 1/f by subtracting fractal from original spectrum
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2-x1';
    powspec.br.osci  = ft_math(cfg, powspec.br.frac, powspec.br.orig);
    powspec.br.osci.powspctrm(powspec.br.osci.powspctrm < 0) = 0;
   
    
    %% Save Data
    
    % save PSD for HLL & BR
    savingname = [path_in filesep subID filesep 'Data' filesep 'PSD' filesep subID '_PSD.mat'];
    save(savingname, 'powspec', '-v7.3');

end
    
    
    
    
    
    
    