%% Figure 1c
% this script rungs additional analysis required for the plots of figure 1c

%%

clear
close all;

%% Set path

path_main = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';
cd(path_main);

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_1/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% Subjects

subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL18', 'OSL21', 'OSL24', 'OSL26', ...
    'OSL27', 'OSL28', 'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', ...
    'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Loop over subjects and extract reaction time and accuracy

% initilize outcome variables reaction time (RT) and accuracy (ACC)
RT_0 = cell(numel(subjects),1); RT_25 = cell(numel(subjects),1); RT_75 = cell(numel(subjects),1);
ACC_0 = cell(numel(subjects),1); ACC_25 = cell(numel(subjects),1); ACC_75 = cell(numel(subjects),1);

for subji = 1:numel(subjects)
    
    subID = subjects{subji};
        
    % go into subject directory
    cd([path_main subID filesep 'Logsheet']);
    
    % read logsheet
    if exist([subID '.xlsx'], 'file') && ~strcmp(subID, 'OSL39')
        logsheet = xlsread([subID '.xlsx']);
    elseif exist([subID '_modified.xlsx'], 'file')
        logsheet = xlsread([subID '_modified.xlsx']);
    else 
        continue;
    end
    
    fprintf('Extract Data %s\n', subID)
    
    if strcmp(subID, 'OSL39') % logsheet for OSL 39 diverges a bit 
        
        % columns containing relevant information
        origtrlnum = logsheet(:,1); % 1 - original trial number
        Acc        = logsheet(:,2); % 2 - ACC 
        ReleaseRT  = logsheet(:,3); % 3 - Collect Release RT
        GoRT       = logsheet(:,4); % 4 - Go RT (NaNs indicate that the trial was a stop trial where subjects had to withdraw a response)
        probstop   = logsheet(:,6); % 6 - probability

        logsheet_mod            = logsheet(:,[1 2 3 4 6]); % get relevant columns
        logsheet_mod(:,3:4)     = logsheet_mod(:,3:4) / 1000;   % convert ms to sec
        logsheet_mod(:,5)       = logsheet_mod(:,5) * 100;      % probability in %

        % create a table to better index
        logsheet_mod_table = array2table(logsheet_mod, 'VariableNames', {'TrlNum', 'ACC', 'ReleaseRT', 'GoRT', 'LikehoodStop'});
        
    else
        
        % columns containing relevant information
        origtrlnum = logsheet(:,20); % 20 - original trial number
        Acc        = logsheet(:,22); % 22 - ACC 
        ReleaseRT  = logsheet(:,27); % 27 - Collect Release RT
        GoRT       = logsheet(:,32); % 32 - Go RT (NaNs indicate that the trial was a stop trial where subjects had to withdraw a response)
        probstop   = logsheet(:,41); % 41 - probability

        logsheet_mod            = logsheet(:,[20 22 27 32 41]); % get relevant columns
        logsheet_mod(:,3:4)     = logsheet_mod(:,3:4) / 1000;   % convert ms to sec
        logsheet_mod(:,5)       = logsheet_mod(:,5) * 100;      % probability in %

        % create a table to better index
        logsheet_mod_table = array2table(logsheet_mod, 'VariableNames', {'TrlNum', 'ACC', 'ReleaseRT', 'GoRT', 'LikehoodStop'});
        
    end

    % extract RT for the different likelihoods of stopping
    
    RT_0{subji}    = logsheet_mod_table.GoRT(logical(logsheet_mod_table.LikehoodStop == 0 ...
        & ~isnan(logsheet_mod_table.GoRT))) - 0.57; % RT for 0% stopping
    RT_0{subji,2}  = subID;
    
    RT_25{subji}   = logsheet_mod_table.GoRT(logical(logsheet_mod_table.LikehoodStop == 25 ...
        & ~isnan(logsheet_mod_table.GoRT)))- 0.57; % RT for 25% stopping
    RT_25{subji,2} = subID;
    
    RT_75{subji}   = logsheet_mod_table.GoRT(logical(logsheet_mod_table.LikehoodStop == 75 ...
        & ~isnan(logsheet_mod_table.GoRT)))- 0.57; % RT for 75% stopping
    RT_75{subji,2} = subID;
    
    % extract ACC for the different likelihoods of stopping
    
    ACC_0{subji}    = logsheet_mod_table.ACC(logical(logsheet_mod_table.LikehoodStop == 0));   % ACC for 0% stopping
    ACC_0{subji,2}  = subID;
    
    ACC_25{subji}   = logsheet_mod_table.ACC(logical(logsheet_mod_table.LikehoodStop == 25));  % ACC for 25% stopping
    ACC_25{subji,2} = subID;
    
    ACC_75{subji}   = logsheet_mod_table.ACC(logical(logsheet_mod_table.LikehoodStop == 75));  % ACC for 75% stopping
    ACC_75{subji,2} = subID;
    
end

%% Remove potential empty cells

RT_0   = reshape(RT_0(~cellfun(@isempty, RT_0)), size(RT_0(~cellfun(@isempty, RT_0)),1)/2, 2);
RT_25  = reshape(RT_25(~cellfun(@isempty, RT_25)), size(RT_25(~cellfun(@isempty, RT_25)),1)/2, 2);
RT_75  = reshape(RT_75(~cellfun(@isempty, RT_75)), size(RT_75(~cellfun(@isempty, RT_75)),1)/2, 2);

ACC_0  = reshape(ACC_0(~cellfun(@isempty, ACC_0)), size(ACC_0(~cellfun(@isempty, ACC_0)),1)/2, 2);
ACC_25 = reshape(ACC_25(~cellfun(@isempty, ACC_25)), size(ACC_25(~cellfun(@isempty, ACC_25)),1)/2, 2);
ACC_75 = reshape(ACC_75(~cellfun(@isempty, ACC_75)), size(ACC_75(~cellfun(@isempty, ACC_75)),1)/2, 2);

%% Group Level Results

% initilize outcome variables reaction time (RT) and accuracy (ACC)
avg_RT_0 = NaN(size(RT_0,1),1); avg_RT_25 = NaN(size(RT_25,1),1); avg_RT_75 = NaN(size(RT_75,1),1);
avg_ACC_0 = NaN(size(ACC_0,1),1); avg_ACC_25 = NaN(size(ACC_25,1),1); avg_ACC_75 = NaN(size(ACC_75,1),1);

nboot = 1000;

for subji = 1:size(avg_RT_0,1)
    
    fprintf('Computing Bootstrap Distribution %d/%d\n', subji, size(avg_RT_0,1));
    
    % get minimum # trials 
    mintrl = numel(RT_75{subji});
    
    % ----------------------
    % extract mean RT
    % ----------------------

    % perform bootstrapping
    boot_0  = NaN(nboot,1);
    boot_25 = NaN(nboot,1);
    for Iboot = 1:nboot
        
        rndidx         = randsample(numel(RT_0{subji}), mintrl); % randomly sample as many trials as in the lowest condition
        boot_0(Iboot)  = mean(RT_0{subji}(rndidx)); clear rndidx
        
        rndidx         = randsample(numel(RT_25{subji}), mintrl);
        boot_25(Iboot) = mean(RT_25{subji}(rndidx)); clear rndidx
        
    end
            
    avg_RT_0(subji)   = mean(boot_0);
    avg_RT_25(subji)  = mean(boot_25);
    avg_RT_75(subji)  = mean(RT_75{subji});

    % ----------------------
    % extract mean ACC
    % ----------------------

    % perform bootstrapping
    boot_0  = NaN(nboot,1);
    boot_25 = NaN(nboot,1);
    for Iboot = 1:nboot
        
        rndidx         = randsample(numel(ACC_0{subji}), mintrl);
        boot_0(Iboot)  = mean(ACC_0{subji}(rndidx)); clear rndidx
        
        rndidx         = randsample(numel(ACC_25{subji}), mintrl);
        boot_25(Iboot) = mean(ACC_25{subji}(rndidx)); clear rndidx
        
    end
    
    avg_ACC_0(subji)  = mean(boot_0);
    avg_ACC_25(subji) = mean(boot_25);
    avg_ACC_75(subji) = mean(ACC_75{subji});
    
end

% group data into matrix
all_avg_RTs  = horzcat(avg_RT_0, avg_RT_25, avg_RT_75);     % matrix for RTs
all_avg_ACCs = horzcat(avg_ACC_0, avg_ACC_25, avg_ACC_75);  % matrix for ACCs

fprintf('RTs:\n\n0%%  = %.3f ± %.3fms\n25%% = %.3f ± %.3fms\n75%% = %.3f ± %.3fms\n', ...
    nanmean(all_avg_RTs(:,1))*1000, nanstd(all_avg_RTs(:,1))*1000, ...
    nanmean(all_avg_RTs(:,2))*1000, nanstd(all_avg_RTs(:,2))*1000, ...
    nanmean(all_avg_RTs(:,3))*1000, nanstd(all_avg_RTs(:,3))*1000)

fprintf('ACCs:\n\n0%%   = %.3f ± %.3f%%\n25%% = %.3f ± %.3f%%\n75%% = %.3f ± %.3f%%\n', ...
    nanmean(all_avg_ACCs(:,1))*100, nanstd(all_avg_ACCs(:,1))*100, ...
    nanmean(all_avg_ACCs(:,2))*100, nanstd(all_avg_ACCs(:,2))*100, ...
    nanmean(all_avg_ACCs(:,3))*100, nanstd(all_avg_ACCs(:,3))*100)

%% Plot Reaction Time

% create cell array
plotRT = {avg_RT_0, avg_RT_25, avg_RT_75}';

% color 2 use for plotting
cb = {'#008450', '#EFB700', '#B81D13'};

figure(1); 
h = rm_raincloud(plotRT, [1 0 0], 0, 'ks');

%% Make Plot Look Nice

% change density plot color
h.p{1,1}.FaceColor = cb{1};
h.p{2,1}.FaceColor = cb{2};
h.p{3,1}.FaceColor = cb{3};
h.p{1,1}.EdgeColor = cb{1};
h.p{2,1}.EdgeColor = cb{2};
h.p{3,1}.EdgeColor = cb{3};
h.p{1,1}.FaceAlpha = 0.8;
h.p{2,1}.FaceAlpha = 0.8;
h.p{3,1}.FaceAlpha = 0.8;

% change color for indiv. dots
h.s{1,1}.MarkerFaceColor = cb{1};
h.s{2,1}.MarkerFaceColor = cb{2};
h.s{3,1}.MarkerFaceColor = cb{3};
h.s{1,1}.MarkerFaceAlpha = 0.8;
h.s{2,1}.MarkerFaceAlpha = 0.8;
h.s{3,1}.MarkerFaceAlpha = 0.8;
h.s{1,1}.SizeData = 60;
h.s{2,1}.SizeData = 60;
h.s{3,1}.SizeData = 60;

% change mean dots
h.m(1).MarkerFaceColor = cb{1};
h.m(2).MarkerFaceColor = cb{2};
h.m(3).MarkerFaceColor = cb{3};
h.m(1).SizeData = 100;
h.m(2).SizeData = 100;
h.m(3).SizeData = 100;

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
xticks([0 0.15 0.25]);
ylabel('Likelihood of Stop'); % x and y axis are flipped in the function
xlabel('RT [sec.]');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [425 276 290 202]);

cd(path_out);
saveas(gcf, 'fig_1_RT.pdf')

% compute confidence interval
CI_0 = jm_computeConfidenceInterval(avg_RT_0);
fprintf('RT 0%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_0(1), CI_0(2));

CI_25 = jm_computeConfidenceInterval(avg_RT_25);
fprintf('RT 25%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_25(1), CI_25(2));

CI_75 = jm_computeConfidenceInterval(avg_RT_75);
fprintf('RT 75%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_75(1), CI_75(2));

%% Same for Accuracy

% create cell array
plotAcc = {avg_ACC_0, avg_ACC_25, avg_ACC_75}';

figure(2); 
h = rm_raincloud(plotAcc, [1 0 0], 0, 'ks');

%% Make Plot Look Nice

% change density plot color
h.p{1,1}.FaceColor = cb{1};
h.p{2,1}.FaceColor = cb{2};
h.p{3,1}.FaceColor = cb{3};
h.p{1,1}.EdgeColor = cb{1};
h.p{2,1}.EdgeColor = cb{2};
h.p{3,1}.EdgeColor = cb{3};
h.p{1,1}.FaceAlpha = 0.8;
h.p{2,1}.FaceAlpha = 0.8;
h.p{3,1}.FaceAlpha = 0.8;

% change color for indiv. dots
h.s{1,1}.MarkerFaceColor = cb{1};
h.s{2,1}.MarkerFaceColor = cb{2};
h.s{3,1}.MarkerFaceColor = cb{3};
h.s{1,1}.MarkerFaceAlpha = 0.8;
h.s{2,1}.MarkerFaceAlpha = 0.8;
h.s{3,1}.MarkerFaceAlpha = 0.8;
h.s{1,1}.SizeData = 60;
h.s{2,1}.SizeData = 60;
h.s{3,1}.SizeData = 60;

% change mean dots
h.m(1).MarkerFaceColor = cb{1};
h.m(2).MarkerFaceColor = cb{2};
h.m(3).MarkerFaceColor = cb{3};
h.m(1).SizeData = 100;
h.m(2).SizeData = 100;
h.m(3).SizeData = 100;

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('Likelihood of Stop'); % x and y axis are flipped in the function
xlabel('Accuracy');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [425 276 290 202]);

% compute confidence interval
CI_0 = jm_computeConfidenceInterval(avg_ACC_0);
fprintf('ACC 0%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_0(1), CI_0(2));

CI_25 = jm_computeConfidenceInterval(avg_ACC_25);
fprintf('ACC 25%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_25(1), CI_25(2));

CI_75 = jm_computeConfidenceInterval(avg_ACC_75);
fprintf('ACC 75%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_75(1), CI_75(2));

cd(path_out);
saveas(gcf, 'fig_1_Acc.pdf')

%% Compute Inter-Quartile Range

% initilize outcome variables for IQR
avg_iqr_0 = NaN(size(RT_0,1),1); avg_iqr_25 = NaN(size(RT_25,1),1); avg_iqr_75 = NaN(size(RT_75,1),1);

for subji = 1:size(avg_iqr_0,1)
    
    fprintf('Bootstrapping Subject %d/%d\n', subji, size(avg_iqr_0,1));
    
    % get minimum number of trials across all conditions
    mintrial = min([numel(RT_0{subji}), numel(RT_25{subji}), numel(RT_75{subji})]);
    
    % create bootstrap distribution
    
    bootdata = {RT_0{subji}, RT_25{subji}}; % data that needs to be bootstrapped
    nperm    = 10000; % number of bootstrap permutations
    
    for Idata = 1:numel(bootdata)
        
        % assign memory for temporary IQR storage
        tmpiqr = NaN(nperm,1); 

        for Iboot = 1:nperm
            rndidx = randsample(numel(bootdata{Idata}),mintrial);
            tmpiqr(Iboot) = iqr(bootdata{Idata}(rndidx));
        end
        
        if Idata == 1
            avg_iqr_0(subji)  = mean(tmpiqr);
        else
            avg_iqr_25(subji) = mean(tmpiqr);
        end
        
    end
    
    % 75% condition does not need to be bootstrapped since it is the
    % condition with the lowest # trials
    avg_iqr_75(subji)  = iqr(RT_75{subji});
    
end

% group data into matrix
all_avg_IQRs  = horzcat(avg_iqr_0, avg_iqr_25, avg_iqr_75); % matrix for RTs

%% Plot IQR

% create cell array
plotIQR = {avg_iqr_0, avg_iqr_25, avg_iqr_75}';

% color to use for plotting
cb = {'#008450', '#EFB700', '#B81D13'};

figure(3); 
h = rm_raincloud(plotIQR, [1 0 0], 0, 'ks', 0.009);

%% Make Plot Look Nice

% change density plot color
h.p{1,1}.FaceColor = cb{1};
h.p{2,1}.FaceColor = cb{2};
h.p{3,1}.FaceColor = cb{3};
h.p{1,1}.EdgeColor = cb{1};
h.p{2,1}.EdgeColor = cb{2};
h.p{3,1}.EdgeColor = cb{3};
h.p{1,1}.FaceAlpha = 0.8;
h.p{2,1}.FaceAlpha = 0.8;
h.p{3,1}.FaceAlpha = 0.8;

% change color for indiv. dots
h.s{1,1}.MarkerFaceColor = cb{1};
h.s{2,1}.MarkerFaceColor = cb{2};
h.s{3,1}.MarkerFaceColor = cb{3};
h.s{1,1}.MarkerFaceAlpha = 0.8;
h.s{2,1}.MarkerFaceAlpha = 0.8;
h.s{3,1}.MarkerFaceAlpha = 0.8;
h.s{1,1}.SizeData = 60;
h.s{2,1}.SizeData = 60;
h.s{3,1}.SizeData = 60;

% change mean dots
h.m(1).MarkerFaceColor = cb{1};
h.m(2).MarkerFaceColor = cb{2};
h.m(3).MarkerFaceColor = cb{3};
h.m(1).SizeData = 100;
h.m(2).SizeData = 100;
h.m(3).SizeData = 100;

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('Likelihood of Stop'); % x and y axis are flipped in the function
xlabel('IQR (RT)');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [425 276 290 202]);

% compute confidence interval
CI_0 = jm_computeConfidenceInterval(avg_iqr_0);
fprintf('IQR 0%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_0(1), CI_0(2));

CI_25 = jm_computeConfidenceInterval(avg_iqr_25);
fprintf('IQR 25%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_25(1), CI_25(2));

CI_75 = jm_computeConfidenceInterval(avg_iqr_75);
fprintf('IQR 75%% LH: 95%% Confidence Interval: [%.3f %.3f]\n', CI_75(1), CI_75(2));

cd(path_out);
saveas(gcf, 'fig_1_IQR.pdf')

%% Bring Data into Excelsheet

% go into results folder
mkdir(path_source, 'Fig_1c');
cd(fullfile(path_source, 'Fig_1c'));

% RT data into table
RTs        = [repmat(1:size(all_avg_RTs,1),1,3)' ...
    [ones(size(all_avg_RTs,1),1);ones(size(all_avg_RTs,1),1)*2; ones(size(all_avg_RTs,1),1)*3] ...
    reshape(all_avg_RTs, [], 1)];

table_RTs  = array2table(RTs, 'VariableNames', {'ID', 'Prediction', 'RT'});

% write data into excelsheet
filename = 'Behavior_RT.xlsx';
writetable(table_RTs,filename, 'Sheet', 1);

% ACC data into table
ACCs       = [repmat(1:size(all_avg_ACCs,1),1,3)' ...
    [ones(size(all_avg_ACCs,1),1);ones(size(all_avg_ACCs,1),1)*2; ones(size(all_avg_ACCs,1),1)*3] ...
    reshape(all_avg_ACCs, [], 1)];

table_ACCs = array2table(ACCs, 'VariableNames', {'ID', 'Prediction', 'ACC'});

% write data into excelsheet
filename = 'Behavior_ACC.xlsx';
writetable(table_ACCs,filename, 'Sheet', 1);

% IQR data into table
IQRs       = [repmat(1:size(all_avg_IQRs,1),1,3)' ...
    [ones(size(all_avg_IQRs,1),1);ones(size(all_avg_IQRs,1),1)*2; ones(size(all_avg_IQRs,1),1)*3] ...
    reshape(all_avg_IQRs, [], 1)];

table_IQRs = array2table(IQRs, 'VariableNames', {'ID', 'Prediction', 'IQR'});

% write data into excelsheet
filename = 'Behavior_IQR.xlsx';
writetable(table_IQRs,filename, 'Sheet', 1);

