%% Figure 2a & b (peak metrics)
% this script rungs additional analysis required for the plots of figure
% 2a & b.

%%

clear
ft_defaults

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_2/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

%% - LOAD DATA - %%

% number of bootstrap distributions
nboot = 1000; 

% preallocate memory
peakdata_active     = struct();

for Isub = 1:numel(subjects)
        
    % get indicated subject ID
    subID = subjects{Isub};
        
    % load data
    load(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data', [subID '_activedata_ROI_pca.mat']));
  
    % merge data
    tmpData  = {roidata.hll, roidata.br};
    
    % --------------------------
    % Get Index for Conditions
    % --------------------------

    if ~isempty(tmpData{1}.frontal.active)
        
        idx_l0  = find(tmpData{1}.frontal.active.trialinfo(:,13) == 0);
        idx_l25 = find(tmpData{1}.frontal.active.trialinfo(:,13) == 25);
        idx_l75 = find(tmpData{1}.frontal.active.trialinfo(:,13) == 75);
        
    elseif ~isempty(tmpData{1}.motor.active)
        
        idx_l0  = find(tmpData{1}.motor.active.trialinfo(:,13) == 0);
        idx_l25 = find(tmpData{1}.motor.active.trialinfo(:,13) == 25);
        idx_l75 = find(tmpData{1}.motor.active.trialinfo(:,13) == 75);
        
    elseif ~isempty(tmpData{1}.mtl.active)
        
        idx_l0  = find(tmpData{1}.mtl.active.trialinfo(:,13) == 0);
        idx_l25 = find(tmpData{1}.mtl.active.trialinfo(:,13) == 25);
        idx_l75 = find(tmpData{1}.mtl.active.trialinfo(:,13) == 75);
        
    end
    
    % merge trials for looping
    alltrialidx = {idx_l0, idx_l25, idx_l75};
    
    % extract the minimum number of trials
    mintrl = min([numel(idx_l0), numel(idx_l25), numel(idx_l75)]);
    
    % ------------------------------------
    % Extract Peak Latency and Amplitude
    % ------------------------------------
                
    fprintf('Compute Peaks %d/%d\n', Isub, numel(subjects));
    
    % loop over trials with different likelihood of stop
    for jj = 1:numel(alltrialidx)
        
        % ################
        % Frontal Cortex
        % ################

        if ~isempty(tmpData{1}.frontal.active)
            
            % get time for epoch
            time  = tmpData{1}.frontal.active.time;
            
            % preallocate memory for peak latency and amplitude
            peaklat = NaN(numel(tmpData{1}.frontal.active.label), nboot, mintrl);
            ampl    = NaN(numel(tmpData{1}.frontal.active.label), nboot, mintrl);

            for Ichan = 1:numel(tmpData{1}.frontal.active.label)
                
                % -----------------------------
                % Create Bootstrap Distribution
                % -----------------------------
                
                for Iboot = 1:nboot
                
                    % get random index for trials
                    rndidx = randsample(numel(alltrialidx{jj}), mintrl);
                    
                    % loop over trials                    
                    for Itrial = 1:size(rndidx,1)

                        % get data for current likelihood of stop
                        dat   = squeeze(tmpData{1}.frontal.active.trial...
                            (alltrialidx{jj}(rndidx(Itrial)),Ichan,:));

                        % make sure dimensions are correct
                        if size(dat,1) == length(time)
                            dat = permute(dat, [2 1]);
                        end

                        [val, idx] = findpeaks(dat);

                        if ~isempty(idx)

                            peaklat(Ichan,Iboot,Itrial) = time(idx(val == max(val)));
                            ampl(Ichan,Iboot,Itrial)    = max(val);

                        end

                    end %Itrial 

                    % remove outliers
                    idxout_lat = isoutlier(squeeze(peaklat(Ichan,Iboot,:)), 'percentiles', [5 95]);
                    idxout_amp = isoutlier(squeeze(ampl(Ichan,Iboot,:)), 'percentiles', [5 95]);
                    
                end %Iboot 

            end %Ichan 

            peaklatency.frontal{jj} = reshape(mean(peaklat,2), 1, [])';
            peakamp.frontal{jj}     = reshape(mean(ampl,2), 1, [])';

            % extract also number of channels that were active
            nchan.frontal = Ichan;
                            
            clear ampl peaklat
            
        else
            
            peaklatency.frontal{jj} = NaN;
            peakamp.frontal{jj}     = NaN;
            nchan.frontal           = NaN;

        end
        
        % ################
        % Motor Cortex
        % ################

        if ~isempty(tmpData{1}.motor.active)
            
            % get time for epoch
            time  = tmpData{1}.motor.active.time;
            
            % preallocate memory for peak latency and amplitude
            peaklat = NaN(numel(tmpData{1}.motor.active.label), nboot, mintrl);
            ampl    = NaN(numel(tmpData{1}.motor.active.label), nboot, mintrl);

            for Ichan = 1:numel(tmpData{1}.motor.active.label)
                
                % -----------------------------
                % Create Bootstrap Distribution
                % -----------------------------
                
                for Iboot = 1:nboot
                
                    % get random index for trials
                    rndidx = randsample(numel(alltrialidx{jj}), mintrl);
                    
                    % loop over trials                    
                    for Itrial = 1:size(rndidx,1)

                        % get data for current likelihood of stop
                        dat   = squeeze(tmpData{1}.motor.active.trial...
                            (alltrialidx{jj}(rndidx(Itrial)),Ichan,:));

                        % make sure dimensions are correct
                        if size(dat,1) == length(time)
                            dat = permute(dat, [2 1]);
                        end

                        [val, idx] = findpeaks(dat);

                        if ~isempty(idx)

                            peaklat(Ichan,Iboot,Itrial) = time(idx(val == max(val)));
                            ampl(Ichan,Iboot,Itrial)    = max(val);

                        end

                    end %Itrial 

                    % remove outliers
                    idxout_lat = isoutlier(squeeze(peaklat(Ichan,Iboot,:)), 'percentiles', [5 95]);
                    idxout_amp = isoutlier(squeeze(ampl(Ichan,Iboot,:)), 'percentiles', [5 95]);
                    
                end %Iboot 

            end %Ichan 

            peaklatency.motor{jj} = reshape(mean(peaklat,2), 1, [])';
            peakamp.motor{jj}     = reshape(mean(ampl,2), 1, [])';

            % extract also number of channels that were active
            nchan.motor = Ichan;
                            
            clear ampl peaklat
            
        else
            
            peaklatency.motor{jj} = NaN;
            peakamp.motor{jj}     = NaN;
            nchan.motor           = NaN;

        end

    end %jj
    
    % --------------------------
    % Bring Data to structure
    % --------------------------
        
    % ------------------
    % Latency
    % ------------------

    % frontal
    peakdata_active.hll.frontal.latency_lh0{Isub}  = peaklatency.frontal{1};
    peakdata_active.hll.frontal.latency_lh25{Isub} = peaklatency.frontal{2};
    peakdata_active.hll.frontal.latency_lh75{Isub} = peaklatency.frontal{3};

    % motor
    peakdata_active.hll.motor.latency_lh0{Isub}    = peaklatency.motor{1};
    peakdata_active.hll.motor.latency_lh25{Isub}   = peaklatency.motor{2};
    peakdata_active.hll.motor.latency_lh75{Isub}   = peaklatency.motor{3};

    % ------------------
    % Amplitude
    % ------------------
    
    % frontal
    peakdata_active.hll.frontal.amp_lh0{Isub}  = peakamp.frontal{1};
    peakdata_active.hll.frontal.amp_lh25{Isub} = peakamp.frontal{2};
    peakdata_active.hll.frontal.amp_lh75{Isub} = peakamp.frontal{3};
    
    % motor
    peakdata_active.hll.motor.amp_lh0{Isub}    = peakamp.motor{1};
    peakdata_active.hll.motor.amp_lh25{Isub}   = peakamp.motor{2};
    peakdata_active.hll.motor.amp_lh75{Isub}   = peakamp.motor{3};
            
    clear peakamp peaklatency
        
end %Isub
    
%% - EXTRACT AVERAGES - %%

MergedData = {peakdata_active.hll};

ga_frontal_lh0  = cell(numel(MergedData),1);
ga_frontal_lh25 = cell(numel(MergedData),1);
ga_frontal_lh75 = cell(numel(MergedData),1);

ga_motor_lh0    = cell(numel(MergedData),1);
ga_motor_lh25   = cell(numel(MergedData),1);
ga_motor_lh75   = cell(numel(MergedData),1);

%% Frontal
    
% Extract all latencies and amplitudes for grand average
ga_frontal_lh0{1}.active.latency       = cellfun(@mean, MergedData{1}.frontal.latency_lh0);
ga_frontal_lh0{1}.active.amplitude     = cellfun(@mean, MergedData{1}.frontal.amp_lh0);

% Extract all latencies and amplitudes for grand average
ga_frontal_lh25{1}.active.latency       = cellfun(@mean, MergedData{1}.frontal.latency_lh25);
ga_frontal_lh25{1}.active.amplitude     = cellfun(@mean, MergedData{1}.frontal.amp_lh25);
    
% Extract all latencies and amplitudes for grand average
ga_frontal_lh75{1}.active.latency       = cellfun(@mean, MergedData{1}.frontal.latency_lh75);
ga_frontal_lh75{1}.active.amplitude     = cellfun(@mean, MergedData{1}.frontal.amp_lh75);

%% Motor
    
% Extract all latencies and amplitudes for grand average
ga_motor_lh0{1}.active.latency       = cellfun(@mean, MergedData{1}.motor.latency_lh0);
ga_motor_lh0{1}.active.amplitude     = cellfun(@mean, MergedData{1}.motor.amp_lh0);

% Extract all latencies and amplitudes for grand average
ga_motor_lh25{1}.active.latency       = cellfun(@mean, MergedData{1}.motor.latency_lh25);
ga_motor_lh25{1}.active.amplitude     = cellfun(@mean, MergedData{1}.motor.amp_lh25);

% Extract all latencies and amplitudes for grand average
ga_motor_lh75{1}.active.latency       = cellfun(@mean, MergedData{1}.motor.latency_lh75);
ga_motor_lh75{1}.active.amplitude     = cellfun(@mean, MergedData{1}.motor.amp_lh75);
  
%% - COLORS - %%

cb = {'#008450', '#EFB700', '#B81D13'};

%% - SELECT EPOCH - %%

selcond = input('Select Epoch: (1) HLL: ');

%% - PLOT LATENCY - %%

cd(path_out);

mergeData                 = struct();
mergeData.frontal.latency = {ga_frontal_lh0{selcond}.active.latency, ga_frontal_lh25{selcond}.active.latency, ...
                             ga_frontal_lh75{selcond}.active.latency};
mergeData.motor.latency   = {ga_motor_lh0{selcond}.active.latency, ga_motor_lh25{selcond}.active.latency, ...
                             ga_motor_lh75{selcond}.active.latency};
       
                         
% ------------------------
% Plot Frontal Cortex
% ------------------------

figure; 
h = rm_raincloud(mergeData.frontal.latency', [1 0 0], 0, 'ks', 0.03);

% admin
for Icond = 1:numel(mergeData.frontal.latency)
    
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
ylabel('LL Stop'); % x and y axis are flipped in the function
if selcond == 1
    xlabel('Time to HLL [s]');
else
    xlabel('Time to BR [s]');
end
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [441 334 266 126]);
xlim([-0.5 0.3]);

view([0 90]);

% save figure
print(gcf,'fig_2_peaklatency_frontal_HLL.pdf','-dpdf','-r400');

% ------------------------
% Plot Motor Cortex
% ------------------------

figure; 
h = rm_raincloud(mergeData.motor.latency', [1 0 0], 0, 'ks', 0.03);

% admin
for Icond = 1:numel(mergeData.motor.latency)
    
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
ylabel('LL Stop'); % x and y axis are flipped in the function
if selcond == 1
    xlabel('Time to HLL [s]');
else
    xlabel('Time to BR [s]');
end
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [441 334 266 126]);
xlim([-0.5 0.3]);

view([0 90]);
  
% save figure
print(gcf,'fig_2_peaklatency_motor_HLL.pdf','-dpdf','-r400');

%% - PLOT AMPLITUDE - %%

cd(path_out);

mergeData             = struct();
mergeData.frontal.amp = {ga_frontal_lh0{selcond}.active.amplitude, ga_frontal_lh25{selcond}.active.amplitude, ...
                             ga_frontal_lh75{selcond}.active.amplitude};
mergeData.motor.amp   = {ga_motor_lh0{selcond}.active.amplitude, ga_motor_lh25{selcond}.active.amplitude, ...
                             ga_motor_lh75{selcond}.active.amplitude};
                            
% ------------------------
% Plot Frontal Cortex
% ------------------------

figure; 
h = rm_raincloud(mergeData.frontal.amp', [1 0 0], 0, 'ks', 0.9);

% admin
for Icond = 1:numel(mergeData.frontal.amp)
    
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
ylabel('LL Stop'); % x and y axis are flipped in the function
xlabel('HFA Peak [z]');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [544 272 186 189]);

% save figure
print(gcf,'fig_2_peakamp_frontal_HLL.pdf','-dpdf','-r400');

% ------------------------
% Plot Motor Cortex
% ------------------------

figure; 
h = rm_raincloud(mergeData.motor.amp', [1 0 0], 0, 'ks', 1);

% admin
for Icond = 1:numel(mergeData.motor.amp)
    
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
ylabel('LL Stop'); % x and y axis are flipped in the function
xlabel('HFA Peak [z]');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [544 272 186 189]);

% save figure
print(gcf,'fig_2_peakamp_motor_HLL.pdf','-dpdf','-r400');

%% Bring Data into Excelsheet

% get only subject number without "OSL"
xx  = regexp(subjects,'\d+(\.)?(\d+)?','match');
subjectnumber =str2double([xx{:}])';      

statsdata = struct();
    
%% HLL

% include frontal data
statsdata.hll(:,4) = vertcat(ga_frontal_lh0{1}.active.latency', ...
                         ga_frontal_lh25{1}.active.latency', ...
                         ga_frontal_lh75{1}.active.latency', ...
                         ...
                         ga_motor_lh0{1}.active.latency', ...
                         ga_motor_lh25{1}.active.latency', ...
                         ga_motor_lh75{1}.active.latency');
                     
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
mkdir(path_source, 'fig_2a_b_upper_left_lower_right');
cd(fullfile(path_source, 'fig_2a_b_upper_left_lower_right'));

statstable_hll    = array2table(statsdata.hll, 'VariableNames', {'ID', 'ROI', 'Prediction', 'PeakLatency'});

% write data into excelsheet
filename = 'fig_2_a_b_upper_left.xlsx';

writetable(statstable_hll,filename, 'Sheet', 'HLL', 'WriteRowNames', true);

%% Extract Data for Stats in Python (Amplitude)

% get only subject number without "OSL"
xx  = regexp(subjects,'\d+(\.)?(\d+)?','match');
subjectnumber =str2double([xx{:}])';      

statsdata = struct();
                      
%% HLL

% include frontal data
statsdata.hll(:,4) = vertcat(ga_frontal_lh0{1}.active.amplitude', ...
                         ga_frontal_lh25{1}.active.amplitude', ...
                         ga_frontal_lh75{1}.active.amplitude', ...
                         ...
                         ga_motor_lh0{1}.active.amplitude', ...
                         ga_motor_lh25{1}.active.amplitude', ...
                         ga_motor_lh75{1}.active.amplitude');
                     
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
mkdir(path_source, 'fig_2a_b_upper_left_lower_right');
cd(fullfile(path_source, 'fig_2a_b_upper_left_lower_right'));

statstable_hll    = array2table(statsdata.hll, 'VariableNames', {'ID', 'ROI', 'Prediction', 'PeakAmplitude'});

% write data into excelsheet
filename = 'fig_2_a_b_lower_right.xlsx';

writetable(statstable_hll,filename, 'Sheet', 'HLL', 'WriteRowNames', true);

    