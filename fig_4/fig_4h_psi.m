%% Phase-Slope Index using the Peak-Frequencies from the PLV Spectrum

% Compute the PSI between a seed motor electrode that explains most of the
% behavioral variance and electrodes within frontal cortex

% ########################################################################

% Input:
% - cleaned and segmented data and HFA signal after ANOVA in order to
% quanitfy which motor electrode explains most behavioral variance by means
% of the HFA signal

% Steps:
% - Get seed electrode in motor cortex
% - Get the peak frequency for the PLV per electrode and only use
%   electrodes that have a peak value (no NaNs)
% - Compute the PSI based on the peak frequency (center frequency) ± 2Hz
% - Create a bootstrap distribution to account for the different # trials

% Output:
% - Z-scored PSI between channel pairs

% ########################################################################

%%

clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_data  = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_plv   = '/Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/PLV/';

path_save  = '/Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/PSI';

addpath('/Volumes/IMPECOG/iEEG_Data_JW/Code/NeuralAnalysis/Raw_Data/PSI/');

addpath('/Users/Janweber/Documents/MATLAB/general_functions/');

%% Load ROI Information 

% electrode information containing electrode index for ROIs
load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/Electrode_Information_ROI/elec_info_ROI.mat

%% Load Peak Frequency Files

datapeak = 'peakfreq_plv_trialshuffling.mat';

load(fullfile(path_plv, datapeak));

%% Frequencies for PSI

% ± 3 Hz from peak frequency
freqedges = 3;

%% Zero-Padding

fsample = 512; % sampling frequency
nPad    = 2; % zero padding in sec per sites

%% 

% initialize structure for phase-slope index
PSI = cell(numel(subjects),1);

for Isub = 1:numel(subjects)

    % get subject ID
    subID = subjects{Isub};

    fprintf('Computing PSI for Subject %s\n', subID);

    % load HFA data that was saved after ANOVA (already segmented)
    load(fullfile(path_data, subID, 'Data', 'HFA', [subID '_TaskActive_ExpVariance.mat']));

    % load raw, but segmented data
    load(fullfile(path_data, subID, 'Data', [subID '_clean_segmented.mat']));

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
        [maxPredictor, Rsquared, p] = explVar_motorseed(seedData, time, behav); 

        chan2use = motor_elec(maxPredictor);

        % get data for motor elecrode
        tmp_motor   = data_clean.hll.trial(:,chan2use,:);

        % zero-pad data to increase frequeency resolution (trials x channel x time)
        data_motor = NaN(size(tmp_motor,1), size(tmp_motor,2), fsample*nPad*2+size(tmp_motor,3));
        for jj = 1:size(tmp_motor,2)
            data_motor(:,jj,:) = [zeros(size(tmp_motor,1),fsample*nPad) squeeze(tmp_motor(:,jj,:)) zeros(size(tmp_motor,1),fsample*nPad)];
        end        

        % get condition indices
        idx_lh0      = find(data_clean.hll.trialinfo(:,13) == 0);
        idx_lh25     = find(data_clean.hll.trialinfo(:,13) == 25);
        idx_lh75     = find(data_clean.hll.trialinfo(:,13) == 75);

    else % if no motor channels present

        continue;

    end

    alltrlidx = {idx_lh0, idx_lh25, idx_lh75, 1:size(data_clean.hll.trialinfo,1)};

    mintrl = min([length(idx_lh0) length(idx_lh25) length(idx_lh75)]);

    % check if subject has either active electrodes within PFC to
    % compute inter-regional PSI

    for Icond = 1:numel(alltrlidx)

        if Icond == 3 
            nboot = 1;
        else
            nboot = 1000;
        end

        % -------------------------------------------------------------------
        % PSI Between Frontal & Motor Cortex
        % -------------------------------------------------------------------

        if isfield(ROIs(Isub).frontal, 'elecidx')

            % -----------------------------------------------------
            % Get Electrodes that showed a distinct peak in the PLV
            % -----------------------------------------------------
            elecs2use = find(~isnan(peak_plv_f_motor{Isub}));

            % get frontal electrodes
            frontal_elec = ROIs(Isub).frontal.elecidx(elecs2use);

            if ~isempty(frontal_elec)

                % -----------------------------------------------------------------
                % Compute PSI using Bootstrapping to Account for TrialDifferences
                % -----------------------------------------------------------------

                % permutation x channel sender x channel receiver x center frequencies
                PSI_Dist = NaN(nboot, length(frontal_elec));

                for Iboot = 1:nboot

                    fprintf('Bootstrapping Frontal-Motor, Iteration: %d/%d, Condition: %d/%d, ID: %d/%d\n', Iboot, nboot, Icond, numel(alltrlidx), Isub, numel(subjects));

                    % get random trial indices
                    rndidx = randsample(alltrlidx{Icond}, mintrl);

                    for Ichan = 1:length(frontal_elec)

                        % get the peak frequency for that channel
                        tmp_PeakFreq = peak_plv_f_motor{Isub}(elecs2use(Ichan));

                        % get the frequency vector
                        freqbins = [round(tmp_PeakFreq - freqedges) round(tmp_PeakFreq + freqedges)];

                        % get data for frontal electrodes
                        tmp_frontal = data_clean.hll.trial(rndidx,frontal_elec(Ichan),:);

                        % zero-pad data to increase frequeency resolution (trials x channel x time)
                        data_frontal = NaN(size(tmp_frontal,1), size(tmp_frontal,2), size(data_motor,3));
                        for jj = 1:size(tmp_frontal,2)
                            data_frontal(:,jj,:) = [zeros(size(tmp_frontal,1),fsample*nPad) squeeze(tmp_frontal(:,jj,:)) zeros(size(tmp_frontal,1),fsample*nPad)];
                        end

                        % concatenate data for PSI calculation (channels x time x trials)
                        data2use  = permute(cat(2, data_motor(rndidx,:,:), data_frontal), [2 3 1]);

                        % compute PSI using function "data2psiX" (needs input channels x time x trials)
                        tmpPSI    = data2psiX(data2use,fsample,freqbins,1);

                        % get important cells within the matrix (first row = motor
                        % electrode, 2:end = frontal electrodes)
                        tmpPSI    = tmpPSI(2:end,1,:);

                        PSI_Dist(Iboot,Ichan) = tmpPSI; clear tmpPSI

                    end

                end

                tmpPSI = squeeze(mean(PSI_Dist, 'all')); clear PSI_Dist

                if Icond == 1
                    PSI{Isub}.f_motor.lh_0.z_task   = tmpPSI; clear tmpPSI
                elseif Icond == 2
                    PSI{Isub}.f_motor.lh_25.z_task  = tmpPSI; clear tmpPSI
                elseif Icond == 3
                    PSI{Isub}.f_motor.lh_75.z_task  = tmpPSI; clear tmpPSI
                elseif Icond == 4
                    PSI{Isub}.f_motor.all.z_task    = tmpPSI; clear tmpPSI
                end

            end

        end

    end

end

PSI = PSI(~cellfun(@isempty, PSI));

%% Save Data

save(fullfile(path_save, 'psi_all_elec_individual_peakfreq_3Hz'), 'PSI')
    

