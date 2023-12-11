%% Time-Frequency on Raw Data

% #######################################################################

% Steps:
% 1. Load Raw Data
% 2. Perform TFR
% 3. Save TFR for each subject

% #######################################################################


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

%% Load ROI Information 

% electrode information containing electrode index for ROIs
load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/Electrode_Information_ROI/elec_info_ROI.mat;

% sampling frequency
fsample = 512;

%% 

for Isub = 1:numel(subjects)
    
    % get subject ID
    subID = subjects{Isub};
    
    fprintf('TFR for %s\n', subID);
    
    % go into folder and get clean data
    cd([path_in filesep subID filesep 'Data'])
    
    load([subID '_clean.mat']);
    
    % load trial information
    load([subID '_Trial_Info_clean.mat']);
    
    %% Get Index for Stop Trials
   
    % get index for stop trials to remove
    IdxStopTrl = find(trl_mod_table.TrlType == 2 | isnan(trl_mod_table.TrlType)); % trialtype 2 = stop trial

    % remove stop trials from table
    tmptable               = trl_mod_table;
    tmptable(IdxStopTrl,:) = [];
    
    %% Get Baseline Data ± 2 sec.
      
    trl = NaN(size(data_clean.sampleinfo,1),3);

    for triali = 1:size(data_clean.trialinfo,1)

        curr_onset = data_clean.sampleinfo(triali,1); % current onset = baseline
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) - fsample*2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + fsample*2 -fsample*2];

    end   
    
    trl(IdxStopTrl,:,:)    = [];

    % redefine trial to probe onset
    cfg          = [];
    cfg.trl      = trl;
    basedata     = ft_redefinetrial(cfg, data_clean);
    
    clear trl
    
    %% Segment Data from Start -2 to +2 sec. relative to HLL

    trl = NaN(size(data_clean.sampleinfo,1),3);

    for triali = 1:size(data_clean.trialinfo,1)

        curr_onset = data_clean.sampleinfo(triali,1); % current onset = baseline
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.HitLowLim(triali) - fsample*2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.HitLowLim(triali) + fsample*2 -fsample*2];

    end    
        
    trl(IdxStopTrl,:,:)    = [];

    % redefine trial to probe onset
    cfg          = [];
    cfg.trl      = trl;
    taskdata_hll = ft_redefinetrial(cfg, data_clean);
    
    clear trl
    
    %% Segment Data from Start -2 to +2 sec. relative to BR

    trl = NaN(size(data_clean.sampleinfo,1),3);

    for triali = 1:size(data_clean.trialinfo,1)

        curr_onset = data_clean.sampleinfo(triali,1); % current onset = baseline
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonRelease(triali) - fsample*2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonRelease(triali) + fsample*2 -fsample*2];

    end    
        
    trl(IdxStopTrl,:,:)    = [];

    % redefine trial to probe onset
    cfg          = [];
    cfg.trl      = trl;
    taskdata_br  = ft_redefinetrial(cfg, data_clean);
    
    clear trl
    
    %% Perform TFR
    
    % prep log spacing for freq analysis on the data
    toi             = taskdata_br.time(1):0.01:taskdata_br.time(end);
    foi_center      = 2.^(-1:1/4:7);  % Center frequencies
    octave          = 1/2;            % Frequency resolution
    foi_min         = 2^(-octave/2)*foi_center;
    foi_max         = 2^(octave/2)*foi_center;
    foi             = (foi_min+foi_max)/2;
    delta_freq      = foi_max-foi_min;
    delta_time      = 0.2;
    n_taper_all     = max(1,round(delta_freq.*delta_time-1));
    foi_center      = round(foi_center*10)/10;
    delta_freq_true = (n_taper_all+1)./delta_time;

    cfg             = [];
    cfg.method      = 'mtmconvol';
    cfg.output      = 'pow';
    cfg.toi         = toi;
    cfg.foi         = foi_center;
    cfg.tapsmofrq   = delta_freq_true./2; 
    cfg.keeptrials  = 'yes';
    cfg.keeptapers  = 'no';
    cfg.t_ftimwin   = ones(1,length(cfg.tapsmofrq))*delta_time;
    cfg.taper       = 'dpss';
    cfg.pad         = 'nextpow2';    
    
    % TFR on baseline activity
    freq_base      = ft_freqanalysis(cfg, basedata);
    
    % TFR on task-related activity HLL
    freq_task_hll  = ft_freqanalysis(cfg, taskdata_hll);
   
    % TFR on task-related activity HLL
    freq_task_br   = ft_freqanalysis(cfg, taskdata_br);
    
    % z-normalize task-related data using bootstrapping
    
    % HLL
    tmpfreq = jw_zbaseline_2(freq_base, freq_task_hll, 0, 0.4, 500);
    freq_task_hll.powspctrm = tmpfreq.zspctrm;
    clear tmpfreq
    
    % BR
    tmpfreq = jw_zbaseline_2(freq_base, freq_task_br, 0, 0.4, 500);
    freq_task_br.powspctrm = tmpfreq.zspctrm;
    clear tmpfreq
    
    
    %% Save Data
    
    % save TFR for HLL
    savingname_hll = [path_in filesep subID filesep 'Data' filesep 'TFR' filesep subID '_HLL_TFR.mat'];
    save(savingname_hll, 'freq_task_hll', '-v7.3');
    
    % save TFR for BR
    savingname_br = [path_in filesep subID filesep 'Data' filesep 'TFR' filesep subID '_BR_TFR.mat'];
    save(savingname_br, 'freq_task_br', '-v7.3');
    
end
    
    
    
    
    
    
    