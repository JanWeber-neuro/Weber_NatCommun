%% Figure 4f
% this script rungs additional analysis required for the plots of figure
% 4f.

%%

clear
ft_defaults

rng('default');

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_4/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

%% - GENERAL SETTINGS - %%

nboot = 1000;

fsample = 512;

%% - LOAD DATA - %%

for Isub = 1:numel(subjects)
        
    % get indicated subject ID
    subID = subjects{Isub};
        
    % load HFA data
    load(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data', [subID '_activedata_ROI_pca.mat']))
    
    tmpData    = {roidata.hll, roidata.br};
    
    % Get index for conditions
    
    if ~isempty(tmpData{1}.frontal.inactive)
        
        idx_l0  = find(tmpData{1}.frontal.inactive.trialinfo(:,13) == 0);
        idx_l25 = find(tmpData{1}.frontal.inactive.trialinfo(:,13) == 25);
        idx_l75 = find(tmpData{1}.frontal.inactive.trialinfo(:,13) == 75);
        
    elseif ~isempty(tmpData{1}.motor.inactive)
        
        idx_l0  = find(tmpData{1}.motor.inactive.trialinfo(:,13) == 0);
        idx_l25 = find(tmpData{1}.motor.inactive.trialinfo(:,13) == 25);
        idx_l75 = find(tmpData{1}.motor.inactive.trialinfo(:,13) == 75);
                
    end
    
    alltrialidx = {idx_l0, idx_l25, idx_l75}; % all trials concatenated
    
    % get minimum number of trials within all conditions
    min_numtrl  = min(cellfun('size',alltrialidx,1));
            
    peakdistance  = struct(); 
    peakamplitude = struct();        
    
    % loop over trials with different likelihood of stop
    for Icond = 1:numel(alltrialidx)
    
        fprintf('ID: %d/%d\n', Isub, numel(subjects));
        
        % ------------------------------------------------------------
        % Frontal Cortex
        % ------------------------------------------------------------
        
        if ~isempty(tmpData{1}.frontal.active)

            % storage variable for peak interval
            peakinterval = NaN(nboot, size(tmpData{1}.frontal.active.trial,2), min_numtrl);
            % storage variable for peak amplitudes
            p_amplitude  = NaN(size(peakinterval));

            % create bootstrap distribution
            for Iboot = 1:nboot

                % get time for epoch
                time  = tmpData{1}.frontal.active.time;

                % get a random index to account for trial number
                rndidx = randsample(numel(alltrialidx{Icond}), min_numtrl);
                % get data for current likelihood of stop
                dat   = tmpData{1}.frontal.active.trial...
                    (alltrialidx{Icond}(rndidx),:,:);

                % make sure dimensions are correct
                if size(dat,1) == length(time)
                    dat = permute(dat, [2 1]);
                end

                % loop over channel
                for Ichan = 1:size(dat,2)

                    % loop over trials (restrict trial number to minimum
                    % trial number across trials to not bias kurtosis)
                    for Itrial = 1:min_numtrl

                        % current data for the trial
                        trialdata = squeeze(dat(Itrial,Ichan,:));

                        [val, idx] = findpeaks(trialdata);

                        if ~isempty(idx)

                            tmp = fsample ./ diff(idx); % get frequency between peaks
                            tmp(isoutlier(tmp)) = [];   % remove outliers
                            peakinterval(Iboot,Ichan,Itrial) = mean(tmp);
                            val(isoutlier(val)) = [];
                            p_amplitude(Iboot,Ichan,Itrial)  = mean(val);

                        else

                            peakinterval(Iboot,Ichan,Itrial) = NaN;
                            p_amplitude(Iboot,Ichan,Itrial)  = NaN;

                        end

                    end % Itrial

                    idxoutlier = isoutlier(squeeze(peakinterval(Iboot,Ichan,:)), 'median');
                    peakinterval(Iboot,Ichan,idxoutlier) = NaN;

                end % Ichan
                
            end %Iboot

            peakdistance.frontal{Icond}  = mean(peakinterval, 'all', 'omitnan');
            peakamplitude.frontal{Icond} = mean(p_amplitude, 'all', 'omitnan');

        else % if no data for ROI

            peakdistance.frontal{Icond} = NaN;
            peakamplitude.frontal{Icond} = NaN;
            
        end
        
        % ------------------------------------------------------------
        % Motor Cortex
        % ------------------------------------------------------------
        
        if ~isempty(tmpData{1}.motor.active)

            % storage variable for peak interval
            peakinterval = NaN(nboot, size(tmpData{1}.motor.active.trial,2), min_numtrl);
            % storage variable for peak amplitudes
            p_amplitude  = NaN(size(peakinterval));

            % create bootstrap distribution
            for Iboot = 1:nboot

                % get time for epoch
                time  = tmpData{1}.motor.active.time;
                
                % get a random index to account for trial number
                rndidx = randsample(numel(alltrialidx{Icond}), min_numtrl);
                % get data for current likelihood of stop
                dat   = tmpData{1}.motor.active.trial...
                    (alltrialidx{Icond}(rndidx),:,:);

                % make sure dimensions are correct
                if size(dat,1) == length(time)
                    dat = permute(dat, [2 1]);
                end

                % loop over channel
                for Ichan = 1:size(dat,2)

                    % loop over trials (restrict trial number to minimum
                    % trial number across trials to not bias kurtosis)
                    for Itrial = 1:min_numtrl

                        % current data for the trial
                        trialdata = squeeze(dat(Itrial,Ichan,:));

                        [val, idx] = findpeaks(trialdata);

                        if ~isempty(idx)

                            tmp = fsample ./ diff(idx); % get frequency between peaks
                            tmp(isoutlier(tmp)) = [];   % remove outliers
                            peakinterval(Iboot,Ichan,Itrial) = mean(tmp);
                            val(isoutlier(val)) = [];
                            p_amplitude(Iboot,Ichan,Itrial)  = mean(val);

                        else

                            peakinterval(Iboot,Ichan,Itrial) = NaN;
                            p_amplitude(Iboot,Ichan,Itrial)  = NaN;

                        end

                    end % Itrial

                    idxoutlier = isoutlier(squeeze(peakinterval(Iboot,Ichan,:)), 'median');
                    peakinterval(Iboot,Ichan,idxoutlier) = NaN;

                end % Ichan
                
            end %Iboot

            peakdistance.motor{Icond}  = mean(peakinterval, 'all', 'omitnan');
            peakamplitude.motor{Icond} = mean(p_amplitude, 'all', 'omitnan');

        else % if no data for ROI

            peakdistance.motor{Icond} = NaN;
            peakamplitude.motor{Icond} = NaN;
            
        end 
                    
    end %Icond
        
    %% Bring Data to structure

    % extract mean of distribution
                
    hfa_rhythm.hll.frontal.lh0.mu(Isub)  = peakdistance.frontal{1};
    hfa_rhythm.hll.frontal.lh25.mu(Isub) = peakdistance.frontal{2};
    hfa_rhythm.hll.frontal.lh75.mu(Isub) = peakdistance.frontal{3};

    hfa_rhythm.hll.motor.lh0.mu(Isub)    = peakdistance.motor{1};
    hfa_rhythm.hll.motor.lh25.mu(Isub)   = peakdistance.motor{2};
    hfa_rhythm.hll.motor.lh75.mu(Isub)   = peakdistance.motor{3};
   

end %Isub

%% - MERGE DATA - %%

MergedData = {hfa_rhythm.hll};

ga_frontal_lh0  = cell(numel(MergedData),1);
ga_frontal_lh25 = cell(numel(MergedData),1);
ga_frontal_lh75 = cell(numel(MergedData),1);

ga_motor_lh0    = cell(numel(MergedData),1);
ga_motor_lh25   = cell(numel(MergedData),1);
ga_motor_lh75   = cell(numel(MergedData),1);
    
for Iepoch = 1:numel(MergedData)

    %% Frontal Cortex
        
    ga_frontal_lh0{1}.active.mu     = MergedData{1}.frontal.lh0.mu;
    ga_frontal_lh25{1}.active.mu    = MergedData{1}.frontal.lh25.mu;
    ga_frontal_lh75{1}.active.mu    = MergedData{1}.frontal.lh75.mu;

    %% Motor Cortex
        
    ga_motor_lh0{1}.active.mu       = MergedData{1}.motor.lh0.mu;
    ga_motor_lh25{1}.active.mu      = MergedData{1}.motor.lh25.mu;
    ga_motor_lh75{1}.active.mu      = MergedData{1}.motor.lh75.mu;
    
end

%% - PLOT DATA - %%

% go to output folder
cd(path_out);

cb =  {'#008450', '#EFB700', '#B81D13'};

% select epoch for plotting
selcond = input('Select Epoch: (1) HLL, (2) BR: ');

% merge data for looping
mergeData = {ga_frontal_lh0{selcond}.active.mu, ga_frontal_lh25{selcond}.active.mu, ga_frontal_lh75{selcond}.active.mu; 
    ga_motor_lh0{selcond}.active.mu, ga_motor_lh25{selcond}.active.mu, ga_motor_lh75{selcond}.active.mu};

for Iroi = 1:size(mergeData,1)
    
    figure(Iroi); 

    cb =  {'#008450', '#EFB700', '#B81D13'};

    h = rm_raincloud({mergeData{Iroi,1}, mergeData{Iroi,2}, mergeData{Iroi,3}}', [1 0 0], 0, 'ks', 0.2);
    
    % admin
    for Icond = 1:3

        % change density plot color
        h.p{Icond,1}.FaceColor        = cb{Icond};
        h.p{Icond,1}.EdgeColor        = cb{Icond};
        h.p{Icond,1}.FaceAlpha        = 0.8;

        % change color for indiv. dots
        h.s{Icond,1}.MarkerFaceColor  = cb{Icond};
        h.s{Icond,1}.MarkerFaceAlpha  = 0.8;
        h.s{Icond,1}.SizeData         = 40;

        % change mean dots
        h.m(Icond).MarkerFaceColor    = cb{Icond};
        h.m(Icond).SizeData           = 80;

    end

    h.l(1).Color    = [0.7 0.7 0.7];
    h.l(2).Color    = [0.7 0.7 0.7];

    % general figure settings
    tmpxticks = yticks; % due to rotation is different
    yticks(tmpxticks);
    yticklabels(fliplr({'0%', '25%', '75%'}));
    ylabel('Likelihood of Stop')
    xlabel('IPI [Hz]')
    set(gca, 'linewidth', 3, 'fontsize', 13);
    box off

    set(gcf, 'Position', [510 356 291 230]);
    
    % save figure 
    if Iroi == 1
        print(gcf,'fig_4_ipi_group_level_frontal_hll.pdf','-dpdf','-r400');        
    else
        print(gcf,'fig_4_ipi_group_level_motor_hll.pdf','-dpdf','-r400');        
    end
    
end %Iroi
        
            
%% Extract Excel Sheet

% get only subject number without "OSL"
xx  = regexp(subjects,'\d+(\.)?(\d+)?','match');
subjectnumber =str2double([xx{:}])';      

statsdata = struct();

%% HLL Mean of Distribution

% include frontal data
statsdata.hll(:,4) = horzcat(hfa_rhythm.hll.frontal.lh0.mu, ...
                             hfa_rhythm.hll.frontal.lh25.mu, ...
                             hfa_rhythm.hll.frontal.lh75.mu, ...
                             hfa_rhythm.hll.motor.lh0.mu, ...
                             hfa_rhythm.hll.motor.lh25.mu, ...
                             hfa_rhythm.hll.motor.lh75.mu)';

statsdata.hll(1:end,1) = repmat(subjectnumber, 6, 1);

% 1 = Frontal, 2 = Motor
statsdata.hll(1:end,2) = vertcat(ones(size(statsdata.hll, 1) / 2, 1), ...
                                   repmat(3, size(statsdata.hll, 1) / 2, 1));

% document factor conditions (0, 25, 75% likelihood of Stop)                               
statsdata.hll(1:end,3) = repmat(vertcat(zeros(size(statsdata.hll, 1) / 6, 1), ...
                                   repmat(25, size(statsdata.hll, 1) / 6, 1), ... 
                                   repmat(75, size(statsdata.hll, 1) / 6, 1)), 2 , 1);

%% Write Data into Table

% go into results folder
mkdir(path_source, 'Fig_4f');
cd(fullfile(path_source, 'Fig_4f'));

statstable_hll    = array2table(statsdata.hll, 'VariableNames', {'ID', 'ROI', 'Prediction', 'Mean'});

% write data into excelsheet
filename = 'fig_4f.xlsx';

writetable(statstable_hll,filename, 'Sheet', 'HLL', 'WriteRowNames', true);

    