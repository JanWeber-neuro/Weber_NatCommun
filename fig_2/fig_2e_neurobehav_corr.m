%% Figure 2e
% this script rungs additional analysis required for the plots of figure
% 2e.

%%

clear
ft_defaults

%% - SUBJECTS - %%

subjects = {'OSL24'};

%% - PATH SETTINGS - %%

path_in   = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_2/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/ % add functions

%% - LOAD DATA - %%

fsample = 512; % sampling frequency

allData = cell(numel(subjects),1);

for Isub = 1:numel(subjects)
        
    % get indicated subject ID

    subID = subjects{Isub};
    
    fprintf('Load Data %s\n', subID);
    
    % load HFA data in ROIs
    load(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data', [subID '_activedata_ROI_pca.mat']));
        
    % load trial information table     
    load(fullfile(path_in, subID, 'Data', [subID '_Trial_Info_clean.mat']));
    
    % modify trial table (remove stop trials or NaN trials)
    trl_mod_table(trl_mod_table.TrlType == 2 | isnan(trl_mod_table.TrlType),:)  = [];   
            
    % get likelihoods of stopping
    idx_LHs = trl_mod_table.LikehoodStop;
    
    % get relevant data (data for HLL and BR)
    tmpData = {roidata.hll};
    
    for Iepoch = 1:numel(tmpData)
        
        % extract RT Data
        RTs        = trl_mod_table.GoRT - 0.56; % get RTs
        
        [~, idx]   = sort(RTs);        % sort RT
        RTs_sorted = RTs(idx,:);       % sorted RTs
        LHs_sorted = idx_LHs(idx);     % sort LHs

        % Sort HFA Data According to Reaction Time
        
        if ~isempty(tmpData{Iepoch}.frontal.active)
            
            time = tmpData{Iepoch}.frontal.active.time;
            
            % sort HFA data and average across channels
            sortdata = squeeze(nanmean(tmpData{Iepoch}.frontal.active.trial(idx,:,:),2));
                                                        
            % Store Data
            allData{Isub}.hll.RTs_sorted = RTs_sorted;
            allData{Isub}.hll.LHs_sorted = LHs_sorted;
            allData{Isub}.hll.HFA_sorted = sortdata; 

        end
        
        clear RTs_sorted sortdata
                
    end % Iepoch
        
end % Isub
        
%% - CONCATENATE TRIALS ACROSS SUBJECTS - %%

% go into results folder
mkdir(path_source, 'Fig_2e');
cd(fullfile(path_source, 'Fig_2e'));

% colormap brewer
addpath /Users/Janweber/Documents/MATLAB/general_functions/colormaps/cbrewer/
cmap = cbrewer('seq', 'YlOrRd', 1000, 'PCHIP');

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

Stacked_HFA_HLL_sort  = smoothdata(Stacked_HFA_HLL_sort, 1, 'gaussian', [2 2]);

figure;
% plot all conditions
imagesc(time, 1:length(Stacked_RT_HLL_sort), Stacked_HFA_HLL_sort); axis xy;
hold on
plot(Stacked_RT_HLL_sort, 1:length(Stacked_RT_HLL_sort), 'k', 'LineWidth', 3);
a = xline(0, 'k--', 'LineWidth', 1);    
colormap(cmap)
yticks([]);
ylabel('Single Trials [N]');
xlabel('Time to HLL [s]');
y = colorbar;
ylabel(y, 'HFA [z]', 'FontSize', 13);
set(gca, 'Clim', [2 9.5], 'linewidth', 3, 'fontsize', 13); box off
  
set(gcf, 'Position', [558 321 253 253]);
 
%% - SAVE FIGURE - %%

print(gcf,'fig_2_stacked_trials_frontal_single_subject.pdf','-dpdf','-r400');

%% - COMPUTE THE CORRELATION 

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
    
%% Linear Regression Partial Model (Peak Time)

X   = peaktime;
mdl = fitlm(X, Stacked_RT_HLL_sort);

anova(mdl, 'summary')
fprintf('\nAdjusted Rsquared = %.4f\n', mdl.Rsquared.Adjusted);

% compute correlation as well
[r,p] = corr(peaktime, Stacked_RT_HLL_sort, 'type', 'Spearman');

fprintf('Correlation peaktime ~ RT: r = %.3f, p = %.3f\n', r, p);

% scatter peaktime ~ RT

cmap = cbrewer('seq', 'YlOrRd', length(X), 'PCHIP');

figure;
h = scatter(Stacked_RT_HLL_sort, time(peaktime), 300, cmap, 'Marker', '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
hold on
a = lsline;
a.LineWidth = 2;

ylabel('HFA_p_e_a_k_t_i_m_e');
xlabel('RT');
set(gca, 'linewidth', 1, 'fontsize', 13);
box off

set(gcf, 'Position', [484 351 228 156]);

% bring into excel file

tmpTable = array2table([Stacked_RT_HLL_sort time(peaktime)'], 'VariableNames', {'RT', 'PeakTime'});

% write data into excelsheet
filename = 'fig_2e_lower_right.xlsx';
writetable(tmpTable,filename, 'Sheet', 1); clear tmpTable

%% - SAVE FIGURE - %%

print(gcf,'fig_2_corr_frontal_peaktime_RT_single_subject.pdf','-dpdf','-r400');

%% Linear Regression Partial Model (Peak Amplitude)

X   = peakamp;
mdl = fitlm(X, Stacked_RT_HLL_sort);

anova(mdl, 'summary')
fprintf('\nAdjusted Rsquared = %.4f\n', mdl.Rsquared.Adjusted);

% compute correlation as well
[r,p] = corr(peakamp, Stacked_RT_HLL_sort, 'type', 'Spearman');

fprintf('Correlation peakamp ~ RT: r = %.3f, p = %.3f\n', r, p);

% scatter data

cmap = cbrewer('seq', 'YlOrRd', length(X), 'PCHIP');

figure;
h = scatter(Stacked_RT_HLL_sort, peakamp, 300, cmap, 'Marker', '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
hold on
a = lsline;
a.LineWidth = 2;

ylabel('HFA_m_a_x');
xlabel('RT');
set(gca, 'linewidth', 1, 'fontsize', 13);
box off

set(gcf, 'Position', [484 351 228 156]);

% bring into excel file

tmpTable = array2table([Stacked_RT_HLL_sort peakamp], 'VariableNames', {'RT', 'PeakAmplitude'});

% write data into excelsheet
filename = 'fig_2e_upper_right.xlsx';
writetable(tmpTable,filename, 'Sheet', 1); clear tmpTable

%% - SAVE FIGURE - %%

print(gcf,'fig_2_corr_frontal_peakamp_RT_single_subject.pdf','-dpdf','-r400');

