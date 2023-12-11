%% PLV Analysis

% Compute the PLV between a seed motor electrode that explains most of the
% behavioral variance and electrodes within frontal cortex

% ##########################################################################
% Input:
% - cleaned and segmented data and HFA signal after ANOVA in order to
% quanitfy which motor electrode explains most behavioral variance by means
% of the HFA signal

% Output:
% - PLV/Imaginary PLV between channel pairs

% ##########################################################################

%%

clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_data  = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out   = '/Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/PLV';

addpath /Users/janweber/Documents/MATLAB/general_functions/

addpath /Volumes/IMPECOG/iEEG_Data_JW/Code/NeuralAnalysis/Raw_Data/PLV/

%% Load ROI Information 

% Load all Electrodes
load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/Electrode_Information_ROI/elec_info_ROI.mat

%% Some Infos

% sampling frequency
fsample = 512;

% number of iterations for random shuffling
Niterations = 500;

%% Frequencies for PLV

minfreq = 3;
maxfreq = 32;

centerfreq = logspace(log10(minfreq), log10(maxfreq), (maxfreq - minfreq)*2^(1/8));
foi = [centerfreq - (centerfreq/4); centerfreq + (centerfreq/4)]';

%% 

% preallocate memory
PLV     = cell(numel(subjects),1);
imagPLV = cell(numel(subjects),1);

for Isub = 1:numel(subjects)
    
    % get subject ID
    subID = subjects{Isub};
    
    fprintf('Computing PLV for Subject %s\n', subID);
        
    % load data after ANOVA
    load(fullfile(path_data, subID, 'Data', 'HFA',  [subID '_TaskActive_ExpVariance.mat']))
    
    % load raw, but segmented data
    load(fullfile(path_data, subID, 'Data', [subID '_clean.mat']));
    
    % load trial information
    load(fullfile(path_data, subID, 'Data', [subID '_Trial_Info_clean.mat']));
    
    % ------------------
    % Segment Data
    % ------------------

    % Get index for stop trials to remove
   
    IdxStopTrl = find(trl_mod_table.TrlType == 2 | isnan(trl_mod_table.TrlType)); % trialtype 2 = stop trial

    % remove stop trials from table
    tmptable               = trl_mod_table;
    tmptable(IdxStopTrl,:) = [];
    
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
    data.hll     = ft_redefinetrial(cfg, data_clean);
    
    time2use     = [dsearchn(data.hll.time', -0.5') ...
                     dsearchn(data.hll.time', 0.3')];
        
    clear trl data_clean
        
    % -------------------------------------------------------------------
    % Determine Motor Electrode That Explains Behavioral Variance
    % -------------------------------------------------------------------
    
    if isfield(ROIs(Isub).motor, 'elecidx')
        
        % get channels corresponding to motor electrodes
        motor_elec = ROIs(Isub).motor.elecidx;

        % get data
        seedData = squeeze(anovadata.hll.trial(:,motor_elec,:));

        % get time
        time     = anovadata.hll.time;

        % get behavior
        behav    = anovadata.hll.trialinfo(:,11);

        % get index for channel that explains most variance of the data
        % --- function explVar_motorseed --- %
        [maxPredictor, Rsquared] = explVar_motorseed(seedData, time, behav); 

        chan2use = motor_elec(maxPredictor);

        % get data for motor elecrode
        data_motor   = data.hll.trial(:,chan2use,:);

        % number of motor seeds (just for good coding practice)
        nseed = numel(chan2use);

        % get condition indices
        PLV{Isub}.idx_lh0      = find(data.hll.trialinfo(:,13) == 0);
        PLV{Isub}.idx_lh25     = find(data.hll.trialinfo(:,13) == 25);
        PLV{Isub}.idx_lh75     = find(data.hll.trialinfo(:,13) == 75);

        % get condition indices
        imagPLV{Isub}.idx_lh0  = find(data.hll.trialinfo(:,13) == 0);
        imagPLV{Isub}.idx_lh25 = find(data.hll.trialinfo(:,13) == 25);
        imagPLV{Isub}.idx_lh75 = find(data.hll.trialinfo(:,13) == 75);

    else % if no motor channels present
        
        continue;
        
    end
    
    % check if subject has electrodes within PFC to
    % compute inter-regional PLV
    
    % -------------------------------------------------------------------
    % PLV Frontal & Motor
    % -------------------------------------------------------------------

    if isfield(ROIs(Isub).frontal, 'elecidx')

        fprintf('%%------------------------------------------------\n');
        fprintf('\n');
        fprintf('Compting PLV Frontal - Motor\n');
        fprintf('\n');
        fprintf('%%------------------------------------------------\n');
        
        % get channels corresponding to frontal electrodes
        frontal_elec = ROIs(Isub).frontal.elecidx;

        % get data for frontal electrodes
        data_frontal = data.hll.trial(:,frontal_elec,:);

        % -----------------------------------------
        % Compute Task-PLV
        % -----------------------------------------
        method = 'trialshuffling'; % method to use for surrogate distribution 
        
        [all_plv, all_iplv] = channelwise_plv(data_motor, data_frontal, foi, nseed, Niterations, ...
            fsample, 1, time2use, method);        

        % store data in cell structure
        PLV{Isub}.f_motor.true           = all_plv.tmp_plv;
        PLV{Isub}.f_motor.z_true         = all_plv.plv_z; 
        PLV{Isub}.f_motor.surro          = all_plv.tmp_plv_surro;
        clear all_plv

        % store data in cell structure
        imagPLV{Isub}.f_motor.true       = all_iplv.tmp_imagplv;
        imagPLV{Isub}.f_motor.z_true     = all_iplv.iplv_z; 
        imagPLV{Isub}.f_motor.surro      = all_iplv.tmp_imagplv_surro;
        clear all_iplv

    end
    
end %Isub

%% Save Data
   
save(fullfile(path_out, 'plv_trialshuffling'), 'PLV')
save(fullfile(path_out, 'imag_plv_trialshuffling'), 'imagPLV')


