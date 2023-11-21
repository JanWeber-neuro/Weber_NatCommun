%% Stacked Single Trials Motor Cortex

%%

clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_in = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/ % add functions

%% Parameters for Thresholding

min_t   = 100;             % ms to consider significance
cutoff  = norminv(0.01);   % z-score at p = 0.05 two-tailed
fsample = 512;             % sampling frequency

%% Initialize Data Storage

allData = cell(numel(subjects),1);

%% Load Data

for Isub = 1:numel(subjects)
        
    % get indicated subject ID

    subID = subjects{Isub};
    
    fprintf('Load Data %s\n', subID);
    
    % load HFA data in ROIs
    load(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data', [subID '_activedata_ROI_pca.mat']));
        
    % load trial information table     
    load(fullfile(path_in, subID, 'Data', [subID '_Trial_Info_clean.mat']));
    
    % modify trial table (remove stop trials or NaN trials)
    trl_mod_table(trl_mod_table.TrlType == 2 | isnan(trl_mod_table.TrlType),:) = [];   
                
    % get relevant data (data for HLL and BR)
    tmpData = roidata.hll;
    
    % extract RT Data
    RTs        = trl_mod_table.RT; % get RTs

    [~, idx]   = sort(RTs);        % sort RT
    RTs_sorted = RTs(idx,:);       % sorted RTs

    %% Sort HFA Data According to Reaction Time

    if ~isempty(tmpData.motor.active)

        time = tmpData.motor.active.time;

        % sort HFA data and average across channels
        sortdata = squeeze(nanmean(tmpData.motor.active.trial(idx,:,:),2));

        % Store Data
        allData{Isub}.hll.RTs_sorted = RTs_sorted;
        allData{Isub}.hll.HFA_sorted = sortdata; 

    end
        
    clear RTs_sorted sortdata
                        
end % Isub
        
%% Concatenate Trials Over Subjects and Plot

% colormap brewer
addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/colormaps/cbrewer/
cmap = cbrewer('seq', 'YlOrRd', 2000, 'PCHIP');

% Initialize variables for HFA, RT and LHs
Stacked_HFA_HLL  = [];
Stacked_RT_HLL   = [];

% concatenate all subject data
for Isub = 1:numel(allData)
    
   if isfield(allData{Isub}, 'hll')
        
        Stacked_HFA_HLL = cat(1, Stacked_HFA_HLL, allData{Isub}.hll.HFA_sorted); % HFA
        Stacked_RT_HLL  = cat(1, Stacked_RT_HLL, allData{Isub}.hll.RTs_sorted);  % reaction time        
   
   end
   
end
    
% sort stacked RT and the rest accordingly
[~, idx]                = sort(Stacked_RT_HLL); % sort RT
Stacked_RT_HLL_sort     = Stacked_RT_HLL(idx);
Stacked_HFA_HLL_sort    = Stacked_HFA_HLL(idx,:); clear idx

% remove outliers
[row, ~] = find(Stacked_HFA_HLL_sort > 15);
Stacked_HFA_HLL_sort(unique(row),:) = [];
Stacked_RT_HLL_sort(unique(row),:) = [];

% smooth a bit
Stacked_HFA_HLL_sort  = movmean(Stacked_HFA_HLL_sort, [20 20], 1);

figure(1);
% plot all conditions
imagesc(time, 1:length(Stacked_RT_HLL_sort), smoothdata(Stacked_HFA_HLL_sort, 'gaussian')); axis xy;
hold on
plot(Stacked_RT_HLL_sort, 1:length(Stacked_RT_HLL_sort), 'k', 'LineWidth', 3);
a = xline(0, 'k--', 'LineWidth', 1);    
colormap(cmap)
yticks([]);
ylabel('Single Trials [N]');
xlabel('Time to HLL [s]');
y = colorbar;
ylabel(y, 'HFA [z]', 'FontSize', 13);
set(gca, 'Clim', [2 5], 'linewidth', 3, 'fontsize', 13); box off

set(gcf, 'Position', [442 335 287 239]);
  
%% Linear Regression Peak Amplitude/Latency and RT

% Initialize variables for HFA, RT and LHs
Stacked_HFA_HLL  = [];
Stacked_RT_HLL   = [];

% concatenate all subject data
for Isub = 1:numel(allData)
    
   if isfield(allData{Isub}, 'hll')
        
        Stacked_HFA_HLL = cat(1, Stacked_HFA_HLL, allData{Isub}.hll.HFA_sorted); % HFA
        Stacked_RT_HLL  = cat(1, Stacked_RT_HLL, allData{Isub}.hll.RTs_sorted);  % reaction time        
   
   end
   
end
    
% sort stacked RT and the rest accordingly
[~, idx]                = sort(Stacked_RT_HLL); % sort RT
Stacked_RT_HLL_sort     = Stacked_RT_HLL(idx);
Stacked_HFA_HLL_sort    = Stacked_HFA_HLL(idx,:); clear idx

Stacked_HFA_HLL_sort  = movmean(Stacked_HFA_HLL_sort, [0 1], 1);

% Extract Peak Power and Latency for Multiple Regression on RT

% initialize peak time and amplitude
peakamp  = NaN(size(Stacked_HFA_HLL_sort,1),1);
peaktime = NaN(size(Stacked_HFA_HLL_sort,1),1);

for ii = 1:size(Stacked_HFA_HLL_sort,1)
    
    [val,idx] = findpeaks(Stacked_HFA_HLL_sort(ii,:));
    if ~isempty(val)
        peaktime(ii) = idx(val == max(val));
        peakamp(ii)  = max(val);
    end
    
end
    
%% Linear Regression Full Model

X   = [peakamp, peaktime];
mdl = fitlm(X, Stacked_RT_HLL_sort);

anova(mdl, 'summary')
fprintf('\nAdjusted Rsquared = %.4f\n', mdl.Rsquared.Adjusted);

%% Linear Regression Partial Model (Peak Time)

X   = peaktime;
mdl = fitlm(X, Stacked_RT_HLL_sort);

anova(mdl, 'summary')
fprintf('\nAdjusted Rsquared = %.4f\n', mdl.Rsquared.Adjusted);

% save source data
mkdir(path_source, 'Table_S1');
cd(fullfile(path_source, 'Table_S1'));

table_data = array2table([time(peaktime)' Stacked_RT_HLL_sort]);

% write data into excelsheet
filename = 'table_S1_motor_peaktiming.xlsx';
writetable(table_data,filename, 'Sheet', 1);

clear table_data

%% Linear Regression Partial Model (Peak Amplitude)

X   = peakamp;
mdl = fitlm(X, Stacked_RT_HLL_sort);

anova(mdl, 'summary')
fprintf('\nAdjusted Rsquared = %.4f\n', mdl.Rsquared.Adjusted);

% save source data
mkdir(path_source, 'Table_S1');
cd(fullfile(path_source, 'Table_S1'));

table_data = array2table([peakamp Stacked_RT_HLL_sort]);

% write data into excelsheet
filename = 'table_S1_motor_peakamp.xlsx';
writetable(table_data,filename, 'Sheet', 1);


