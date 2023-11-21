%% Figure S1 a & b
% this script rungs additional analysis required for the plots of figure
% S1a & b.

%%

clear 
rng('default')

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL18', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in   = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_S1/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

%%

% preallocate memory for d'prime
dprime           = struct();

dprime.lh_25     = NaN(numel(subjects),1); 
dprime.lh_75     = NaN(numel(subjects),1); 

% preallocate memory for criterion
criterion        = struct();

criterion.lh_25  = NaN(numel(subjects),1); 
criterion.lh_75  = NaN(numel(subjects),1); 

%% - COMPUTE PARAMETERS - %%

for Isub = 1:numel(subjects)
    
    % subject ID
    subID = subjects{Isub};
    
    % patient OSL18 does only have the log file
    if strcmp(subID, 'OSL18')

        % go into subject directory
        filename = fullfile(path_in, subID, 'Logsheet', [subID '.xlsx']);

        logsheet = xlsread(filename);

        % columns containing relevant information
        origtrlnum = logsheet(:,20); % 20 - original trial number
        Acc        = logsheet(:,22); % 22 - ACC 
        ReleaseRT  = logsheet(:,27); % 27 - Collect Release RT
        GoRT       = logsheet(:,32); % 32 - Go RT (NaNs indicate that the trial was a stop trial where subjects had to withdraw a response)
        hstoppoint = logsheet(:,34) - logsheet(:,35); % 34-35 - is h stop point in px before hitting lower limit
        probstop   = logsheet(:,41); % 41 - probability
        ssd        = logsheet(:,45); % 45 - SSD (Stop Signal Delay)

        logsheet_mod            = logsheet(:,[20 22 27 32 34 35 41 45]); % get relevant columns
        logsheet_mod            = [logsheet_mod(:,1:4), logsheet_mod(:,5)-logsheet_mod(:,6), logsheet_mod(:,7:end)]; % perform operation column 34-35 
        logsheet_mod(:,3:4)     = logsheet_mod(:,3:4) / 1000; % convert ms to sec
        logsheet_mod(:,6)       = logsheet_mod(:,6) * 100; % probability in %
        logsheet_mod(:,7)       = logsheet_mod(:,7) / 1000;

        logsheet_mod_table = array2table(logsheet_mod, 'VariableNames', {'TrlNum', 'Acc', 'ReleaseRT', 'GoRT', 'HStopPoint', 'LikehoodStop', 'SSD'});

        % ------------------------------------
        % 25% Likelihood of Stop
        % ------------------------------------

        % compute parameters
        Hit  = sum(isnan(logsheet_mod_table.SSD) & logsheet_mod_table.Acc == 1 & logsheet_mod_table.LikehoodStop == 25);
        Miss = sum(isnan(logsheet_mod_table.SSD) & logsheet_mod_table.Acc == 0 & logsheet_mod_table.LikehoodStop == 25);
        FA   = sum(~isnan(logsheet_mod_table.SSD) & logsheet_mod_table.Acc == 0 & logsheet_mod_table.LikehoodStop == 25);
        CR   = sum(~isnan(logsheet_mod_table.SSD) & logsheet_mod_table.Acc == 1 & logsheet_mod_table.LikehoodStop == 25);

        % compute dprime and criterion
        [dp,lambda] = SDT(Hit,FA,Miss,CR);

        if ~isinf(dp)
            dprime.lh_25(Isub) = dp;
        else
            dprime.lh_25(Isub) = NaN;
        end

        criterion.lh_25(Isub)  = lambda; clear dp lambda

        % ------------------------------------
        % 75% Likelihood of Stop
        % ------------------------------------

        % compute parameters
        Hit  = sum(isnan(logsheet_mod_table.SSD) & logsheet_mod_table.Acc == 1 & logsheet_mod_table.LikehoodStop == 75);
        Miss = sum(isnan(logsheet_mod_table.SSD) & logsheet_mod_table.Acc == 0 & logsheet_mod_table.LikehoodStop == 75);
        FA   = sum(~isnan(logsheet_mod_table.SSD) & logsheet_mod_table.Acc == 0 & logsheet_mod_table.LikehoodStop == 75);
        CR   = sum(~isnan(logsheet_mod_table.SSD) & logsheet_mod_table.Acc == 1 & logsheet_mod_table.LikehoodStop == 75);

    else
        
        % load trial information 
        load(fullfile(path_in, subID, 'Data', [subID '_Trial_Info.mat']));

        % ------------------------------------
        % 25% Likelihood of Stop
        % ------------------------------------

        % compute parameters
        Hit  = sum(trl_mod_table.TrlType == 4 & trl_mod_table.Acc == 1 & trl_mod_table.LikehoodStop == 25);
        Miss = sum(trl_mod_table.TrlType == 4 & trl_mod_table.Acc == 0 & trl_mod_table.LikehoodStop == 25);
        FA   = sum(trl_mod_table.TrlType == 2 & trl_mod_table.Acc == 0 & trl_mod_table.LikehoodStop == 25);
        CR   = sum(trl_mod_table.TrlType == 2 & trl_mod_table.Acc == 1 & trl_mod_table.LikehoodStop == 25);

        % compute dprime and criterion
        [dp,lambda] = SDT(Hit,FA,Miss,CR);

        if ~isinf(dp)
            dprime.lh_25(Isub) = dp;
        else
            dprime.lh_25(Isub) = NaN;
        end

        criterion.lh_25(Isub)  = lambda; clear dp lambda

        % ------------------------------------
        % 75% Likelihood of Stop
        % ------------------------------------

        % compute parameters
        Hit  = sum(trl_mod_table.TrlType == 4 & trl_mod_table.Acc == 1 & trl_mod_table.LikehoodStop == 75);
        Miss = sum(trl_mod_table.TrlType == 4 & trl_mod_table.Acc == 0 & trl_mod_table.LikehoodStop == 75);
        FA   = sum(trl_mod_table.TrlType == 2 & trl_mod_table.Acc == 0 & trl_mod_table.LikehoodStop == 75);
        CR   = sum(trl_mod_table.TrlType == 2 & trl_mod_table.Acc == 1 & trl_mod_table.LikehoodStop == 75);

    end
    
    % compute dprime and criterion
    [dp,lambda] = SDT(Hit,FA,Miss,CR);
    
    if ~isinf(dp)
        dprime.lh_75(Isub) = dp;
    else
        dprime.lh_75(Isub) = NaN;
    end
    
    criterion.lh_75(Isub)  = lambda; clear dp lambda
      
end %Isub    
    
%% - SAVE SOURCE DATA - %%

% go into results folder
mkdir(path_source, 'Fig_S1a_b');
cd(fullfile(path_source, 'Fig_S1a_b'));

% - STORE CRITERION FOR 25 & 75% LH - %
tmp          = NaN(numel(subjects),3);
tmp(:,2:end) = horzcat(criterion.lh_25, criterion.lh_75);
tmp(:,1)     = (1:numel(subjects))';

table.c = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'Fig_S1_criterion.xlsx';
writetable(table.c,filename, 'Sheet', 1);
  
% - STORE DPRIME FOR 25 & 75% LH - %
tmp          = NaN(numel(subjects),3);
tmp(:,2:end) = horzcat(dprime.lh_25, dprime.lh_75);
tmp(:,1)     = (1:numel(subjects))';

table.dprime = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'Fig_S1_dprime.xlsx';
writetable(table.dprime,filename, 'Sheet', 1);

%% -- PLOT D'PRIME -- %%

hexcolors = {'#EFB700', '#B81D13'};

colors2use = NaN(2,3);
colorBins = [2:3; 4:5; 6:7];

for iHex = 1:numel(hexcolors)
    for iBin = 1:size(colorBins)
        colors2use(iHex,iBin) = hex2dec(hexcolors{iHex}(colorBins(iBin,1):colorBins(iBin,end)));
    end
end

% scale by dividing through 256 (Matlab conversion)
colors2use = colors2use/ 256;

figure;

% data to plot
plot_dprime = [dprime.lh_25, dprime.lh_75];

jw_raincloud_pairplot(plot_dprime, {colors2use(1,:), colors2use(2,:)});

set(gca, 'ytick', []);
xlabel('d-prime');

h = gca; h.YAxis.Visible = 'off';

set(gca, 'fontsize', 13);

box off

%% - SAVE FIGURE - %%

cd(path_out);
print(gcf,'supplementary_fig_1_dprime.pdf','-dpdf','-r400');    

%% - PLOT CRITERION - %%

figure;

% data to plot
plot_lambda = [criterion.lh_25, criterion.lh_75];

jw_raincloud_pairplot(plot_lambda, {colors2use(1,:), colors2use(2,:)});

set(gca, 'ytick', []);
xlabel('criterion');

h = gca; h.YAxis.Visible = 'off';

set(gca, 'fontsize', 13);

box off

%% - SAVE FIGURE - %%

cd(path_out);
print(gcf,'supplementary_fig_1_criterion.pdf','-dpdf','-r400');    

    