%% Figure 2a & b (lower left)
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

path_in   = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_2/revision_natcomm_final/';

path_source = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% - LOAD DATA - %%
 
hfa_lh_0   = struct();  % stores HFA for the 0% LH Condition
hfa_lh_25  = struct();  % stores HFA for the 25% LH Condition
hfa_lh_75  = struct();  % stores HFA for the 75% LH Condition

nepochs    = 3;
count      = zeros(nepochs,1); % set counter for subjects who have frontal cortex data

nboot      = 1000; % number of bootstrap distributions

for Isub = 1:numel(subjects)
    
    % get indicated subject ID
    subID = subjects{Isub};
    
    % load HFA data in ROIs
    cd(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data'));

    fprintf('Load Data %d/%d\n', Isub, numel(subjects));
    
    % load data
    load([subID '_activedata_ROI_pca.mat']);
            
    % get relevant data
    tmpData = {roidata.hll, roidata.br};
    
    % loop over epochs 
    for Iepoch = 1:numel(tmpData)
           
        % increase counter
        count(Iepoch) = count(Iepoch) + 1;
        
        % #################
        % Frontal Cortex
        % #################

        % get HFA data for below and above median split
        if ~isempty(tmpData{Iepoch}.frontal.active)
                                              
            % get likelihoods of stopping
            idx_lh0  = find(tmpData{Iepoch}.frontal.active.trialinfo(:,13) == 0);
            idx_lh25 = find(tmpData{Iepoch}.frontal.active.trialinfo(:,13) == 25);
            idx_lh75 = find(tmpData{Iepoch}.frontal.active.trialinfo(:,13) == 75);
            
            % ----------------------------------------------------------
            % Bootstrapping to Account for Trial Differences
            % ----------------------------------------------------------
                        
            mintrl = min([length(idx_lh0), length(idx_lh25), length(idx_lh75)]);
            
            tmp_lh0  = NaN(nboot, mintrl, numel(tmpData{Iepoch}.frontal.active.label), length(tmpData{Iepoch}.frontal.active.time));
            tmp_lh25 = NaN(nboot, mintrl, numel(tmpData{Iepoch}.frontal.active.label), length(tmpData{Iepoch}.frontal.active.time));
            tmp_lh75 = NaN(nboot, mintrl, numel(tmpData{Iepoch}.frontal.active.label), length(tmpData{Iepoch}.frontal.active.time));
            
            for Iboot = 1:nboot

                % temporary data for 0% LH
                rndidx = randsample(length(idx_lh0), mintrl);
                tmp_lh0(Iboot,:,:,:)  = tmpData{Iepoch}.frontal.active.trial(idx_lh0(rndidx),:,:); clear rndidx
               
                % temporary data for 25% LH
                rndidx = randsample(length(idx_lh25), mintrl);
                tmp_lh25(Iboot,:,:,:) = tmpData{Iepoch}.frontal.active.trial(idx_lh25(rndidx),:,:);

                % temporary data for 75% LH
                rndidx = randsample(length(idx_lh75), mintrl);
                tmp_lh75(Iboot,:,:,:) = tmpData{Iepoch}.frontal.active.trial(idx_lh75(rndidx),:,:);

            end
            
            if Iepoch == 1
                hfa_lh_0.hll.frontal(count(Iepoch),:)    = squeeze(mean(mean(mean(tmp_lh0,1),2),3))';
                hfa_lh_25.hll.frontal(count(Iepoch),:)   = squeeze(mean(mean(mean(tmp_lh25,1),2),3))';
                hfa_lh_75.hll.frontal(count(Iepoch),:)   = squeeze(mean(mean(mean(tmp_lh75,1),2),3))';
            elseif Iepoch == 2
                hfa_lh_0.br.frontal(count(Iepoch),:)     = squeeze(mean(mean(mean(tmp_lh0,1),2),3))';
                hfa_lh_25.br.frontal(count(Iepoch),:)    = squeeze(mean(mean(mean(tmp_lh25,1),2),3))';
                hfa_lh_75.br.frontal(count(Iepoch),:)    = squeeze(mean(mean(mean(tmp_lh75,1),2),3))';
            end
            
        end
        
        % #################
        % Motor Cortex
        % #################

        % get HFA data for below and above median split
        if ~isempty(tmpData{Iepoch}.motor.active)
                                  
            % remove trials where acc = 0
            idx_false = find(tmpData{Iepoch}.motor.active.trialinfo(:,9) == 0);
            tmpData{Iepoch}.motor.active.trial(idx_false,:,:)     = [];
            tmpData{Iepoch}.motor.active.trialinfo(idx_false,:,:) = [];
            
            % get likelihoods of stopping
            idx_lh0  = find(tmpData{Iepoch}.motor.active.trialinfo(:,13) == 0);
            idx_lh25 = find(tmpData{Iepoch}.motor.active.trialinfo(:,13) == 25);
            idx_lh75 = find(tmpData{Iepoch}.motor.active.trialinfo(:,13) == 75);
            
            % ----------------------------------------------------------
            % Bootstrapping to Account for Trial Differences
            % ----------------------------------------------------------
            
            nboot = 1000;
            
            mintrl = min([length(idx_lh0), length(idx_lh25), length(idx_lh75)]);
            
            tmp_lh0  = NaN(nboot, mintrl, numel(tmpData{Iepoch}.motor.active.label), length(tmpData{Iepoch}.motor.active.time));
            tmp_lh25 = NaN(nboot, mintrl, numel(tmpData{Iepoch}.motor.active.label), length(tmpData{Iepoch}.motor.active.time));
            tmp_lh75 = NaN(nboot, mintrl, numel(tmpData{Iepoch}.motor.active.label), length(tmpData{Iepoch}.motor.active.time));
            
            for Iboot = 1:nboot

                % temporary data for 0% LH
                rndidx = randsample(length(idx_lh0), mintrl);
                tmp_lh0(Iboot,:,:,:)  = tmpData{Iepoch}.motor.active.trial(idx_lh0(rndidx),:,:); clear rndidx
               
                % temporary data for 25% LH
                rndidx = randsample(length(idx_lh25), mintrl);
                tmp_lh25(Iboot,:,:,:) = tmpData{Iepoch}.motor.active.trial(idx_lh25(rndidx),:,:);

                % temporary data for 75% LH
                rndidx = randsample(length(idx_lh75), mintrl);
                tmp_lh75(Iboot,:,:,:) = tmpData{Iepoch}.motor.active.trial(idx_lh75(rndidx),:,:);

            end
            
            if Iepoch == 1
                hfa_lh_0.hll.motor(count(Iepoch),:)    = squeeze(mean(mean(mean(tmp_lh0,1),2),3))';
                hfa_lh_25.hll.motor(count(Iepoch),:)   = squeeze(mean(mean(mean(tmp_lh25,1),2),3))';
                hfa_lh_75.hll.motor(count(Iepoch),:)   = squeeze(mean(mean(mean(tmp_lh75,1),2),3))';
            elseif Iepoch == 2
                hfa_lh_0.br.motor(count(Iepoch),:)     = squeeze(mean(mean(mean(tmp_lh0,1),2),3))';
                hfa_lh_25.br.motor(count(Iepoch),:)    = squeeze(mean(mean(mean(tmp_lh25,1),2),3))';
                hfa_lh_75.br.motor(count(Iepoch),:)    = squeeze(mean(mean(mean(tmp_lh75,1),2),3))';
            end
            
        end        
        
    end % Iepoch
    
end % Isub

%% - SAVE SOURCE DATA - %%

% go into results folder
mkdir(path_source, 'Fig_2a_b_lower_left');
cd(fullfile(path_source, 'Fig_2a_b_lower_left'));

time = tmpData{1}.frontal.active.time;

% -- FRONTAL CORTEX 0% LH -- %%
tmp = NaN(size(hfa_lh_0.hll.frontal,1)+1, size(hfa_lh_0.hll.frontal,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = hfa_lh_0.hll.frontal;

table_frontal.lh0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_2a_hfa_frontal_lh0.xlsx';
writetable(table_frontal.lh0,filename, 'Sheet', 1);

% -- FRONTAL CORTEX 25% LH -- %%
tmp = NaN(size(hfa_lh_25.hll.frontal,1)+1, size(hfa_lh_25.hll.frontal,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = hfa_lh_25.hll.frontal;

table_frontal.lh25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_2a_hfa_frontal_lh25.xlsx';
writetable(table_frontal.lh25,filename, 'Sheet', 1);

% -- FRONTAL CORTEX 75% LH -- %%
tmp = NaN(size(hfa_lh_75.hll.frontal,1)+1, size(hfa_lh_75.hll.frontal,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = hfa_lh_75.hll.frontal;

table_frontal.lh75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_2a_hfa_frontal_lh75.xlsx';
writetable(table_frontal.lh75,filename, 'Sheet', 1);

% -- MOTOR CORTEX 0% LH -- %%
tmp = NaN(size(hfa_lh_0.hll.motor,1)+1, size(hfa_lh_0.hll.motor,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = hfa_lh_0.hll.motor;

table_motor.lh0   = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_2a_hfa_motor_lh0.xlsx';
writetable(table_motor.lh0,filename, 'Sheet', 1);

% -- MOTOR CORTEX 25% LH -- %%
tmp = NaN(size(hfa_lh_25.hll.motor,1)+1, size(hfa_lh_25.hll.motor,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = hfa_lh_25.hll.motor;

table_motor.lh25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_2a_hfa_motor_lh25.xlsx';
writetable(table_motor.lh25,filename, 'Sheet', 1);

% -- MOTOR CORTEX 75% LH -- %%
tmp = NaN(size(hfa_lh_75.hll.motor,1)+1, size(hfa_lh_75.hll.motor,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = hfa_lh_75.hll.motor;

table_motor.lh75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_2a_hfa_motor_lh75.xlsx';
writetable(table_motor.lh75,filename, 'Sheet', 1);

%% - REMOVE ROWS WITH NO DATA - %%

hfa_lh_0.hll.frontal(find(all(hfa_lh_0.hll.frontal == 0, 2)),:)       = [];
hfa_lh_25.hll.frontal(find(all(hfa_lh_25.hll.frontal == 0, 2)),:)     = [];
hfa_lh_75.hll.frontal(find(all(hfa_lh_75.hll.frontal == 0, 2)),:)     = [];

hfa_lh_0.hll.motor(find(all(hfa_lh_0.hll.motor == 0, 2)),:)           = [];
hfa_lh_25.hll.motor(find(all(hfa_lh_25.hll.motor == 0, 2)),:)         = [];
hfa_lh_75.hll.motor(find(all(hfa_lh_75.hll.motor == 0, 2)),:)         = [];

hfa_lh_0.br.frontal(find(all(hfa_lh_0.br.frontal == 0, 2)),:)         = [];
hfa_lh_25.br.frontal(find(all(hfa_lh_25.br.frontal == 0, 2)),:)       = [];
hfa_lh_75.br.frontal(find(all(hfa_lh_75.br.frontal == 0, 2)),:)       = [];

hfa_lh_0.br.motor(find(all(hfa_lh_0.br.motor == 0, 2)),:)             = [];
hfa_lh_25.br.motor(find(all(hfa_lh_25.br.motor == 0, 2)),:)           = [];
hfa_lh_75.br.motor(find(all(hfa_lh_75.br.motor == 0, 2)),:)           = [];

%% - GET TIME INFO - %%

% get time information
time  = {tmpData{1}.frontal.active.time, tmpData{2}.frontal.active.time};

%% - CLUSTER-BASED PERMUTATION STATS - %%

data_lh0  = struct();
data_lh25 = struct();
data_lh75 = struct();

% Frontal data
data_lh0.frontal   = {hfa_lh_0.hll.frontal, hfa_lh_0.br.frontal};
data_lh25.frontal  = {hfa_lh_25.hll.frontal, hfa_lh_25.br.frontal};
data_lh75.frontal  = {hfa_lh_75.hll.frontal, hfa_lh_75.br.frontal};

% Motor data
data_lh0.motor     = {hfa_lh_0.hll.motor, hfa_lh_0.br.motor};
data_lh25.motor    = {hfa_lh_25.hll.motor, hfa_lh_25.br.motor};
data_lh75.motor    = {hfa_lh_75.hll.motor, hfa_lh_75.br.motor};

% Get data of interest
q       = 'Statistic: Frontal = 1; Motor = 2: ';
ans_roi = input(q);

% Get segment of interest 
q       = 'Which segment: HLL = 1; BR = 2: ';
ans_seg = input(q);

if ans_roi == 1
    
    tmp_1 = data_lh0.frontal{ans_seg};
    tmp_2 = data_lh25.frontal{ans_seg};
    tmp_3 = data_lh75.frontal{ans_seg};
    label = 'Frontal';
        
elseif ans_roi == 2
    
    tmp_1 = data_lh0.motor{ans_seg};
    tmp_2 = data_lh25.motor{ans_seg};
    tmp_3 = data_lh75.motor{ans_seg};
    label = 'Motor';

end

statsdata_1  = cell(size(tmp_1,1),1);
statsdata_2  = cell(size(tmp_2,1),1);
statsdata_3  = cell(size(tmp_3,1),1);

for Isub = 1:numel(statsdata_1)
    
    statsdata_1{Isub}.avg      = tmp_1(Isub,:);
    statsdata_1{Isub}.time     = time{ans_seg};
    statsdata_1{Isub}.dimord   = 'chan_time';
    statsdata_1{Isub}.label    = {label};
    statsdata_1{Isub}.cfg      = [];
    
    statsdata_2{Isub}.avg      = tmp_2(Isub,:);
    statsdata_2{Isub}.time     = time{ans_seg};
    statsdata_2{Isub}.dimord   = 'chan_time';
    statsdata_2{Isub}.label    = {label};
    statsdata_2{Isub}.cfg      = [];
    
    statsdata_3{Isub}.avg      = tmp_3(Isub,:);
    statsdata_3{Isub}.time     = time{ans_seg};
    statsdata_3{Isub}.dimord   = 'chan_time';
    statsdata_3{Isub}.label    = {label};
    statsdata_3{Isub}.cfg      = [];
    
    
end
    
% --------------------------------------------------
% Compute F-Statistics for Main Effect of Prediction
% --------------------------------------------------

cfg                     = [];
cfg.channel             = 'all';
cfg.tail                = 1;
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesFunivariate';
cfg.alpha               = 0.05;
cfg.correctm            = 'cluster';
cfg.numrandomization    = 10000;  
cfg.spmversion          = 'spm12';
cfg.computestat         = 'yes'; 
cfg.computecritval      = 'yes'; 
cfg.computeprob         = 'yes'; 

% make design matrix
subj = numel(statsdata_1);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
for i = 1:subj
design(1,subj+subj+i) = i;
end
design(2,1:subj)             = 1;
design(2,subj+1:2*subj)      = 2;
design(2,subj+subj+1:3*subj) = 3;

% clarify independent and dependent variable
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

% do the stats
stat = ft_timelockstatistics(cfg, statsdata_1{:}, statsdata_2{:}, statsdata_3{:});

% --------------------------------------------------
% Post-Hoc Pairwise Comparisons
% --------------------------------------------------

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = 'all';
cfg.tail                = 0;
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.alpha               = 0.05;
cfg.correctm            = 'cluster';
cfg.correcttail         = 'prob';
cfg.numrandomization    = 10000;  
cfg.spmversion          = 'spm12';
cfg.computestat         = 'yes'; 
cfg.computecritval      = 'yes'; 
cfg.computeprob         = 'yes'; 

subj = numel(statsdata_1);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

stat_0vs25   = ft_timelockstatistics(cfg, statsdata_2{:}, statsdata_1{:});
stat_0vs75   = ft_timelockstatistics(cfg, statsdata_3{:}, statsdata_1{:});
stat_25vs75  = ft_timelockstatistics(cfg, statsdata_3{:}, statsdata_2{:});

%% - PLOT DATA - %%

% plot data
h = shadedErrorBar(time{ans_seg}, nanmean(tmp_1), nanstd(tmp_1) / ...
    sqrt(size(tmp_1,1)));
h.mainLine.Color = '#008450'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#008450'; h.patch.EdgeColor = '#008450'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#008450'; h.edge(2).Color = '#008450';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

hold on

h = shadedErrorBar(time{ans_seg}, nanmean(tmp_2), nanstd(tmp_2) / ...
    sqrt(size(tmp_2,1)));
h.mainLine.Color = '#EFB700'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#EFB700'; h.patch.EdgeColor = '#EFB700'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#EFB700'; h.edge(2).Color = '#EFB700';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

hold on

h = shadedErrorBar(time{ans_seg}, nanmean(tmp_3), nanstd(tmp_3) / ...
    sqrt(size(tmp_3,1)));
h.mainLine.Color = '#B81D13'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#B81D13'; h.patch.EdgeColor = '#B81D13'; h.patch.FaceAlpha = 0.5;
h.edge(1).Color = '#B81D13'; h.edge(2).Color = '#B81D13';
h.edge(1).LineWidth = 0.1; h.edge(2).LineWidth = 0.1;

ylabel('HFA [z]');
if ans_seg == 1
    xlabel('Time to HLL');
elseif ans_seg == 2
    xlabel('Time to BR');
end

% plot significant timepoints
a = ylim;

island = bwconncomp(stat.mask);
for jj = 1:island.NumObjects
    h = line(time{ans_seg}(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.1, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    b = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 4);
    b.Color = 'k';
    hold on
    b = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 4);
    b.Color = 'k';
end
clear island

hold on

island = bwconncomp(stat_0vs75.mask);
for jj = 1:island.NumObjects
    h = line(time{ans_seg}(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.3, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    b = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 4);
    b.Color = '#008450';
    hold on
    b = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 4);
    b.Color = '#B81D13';
end
clear island

island = bwconncomp(stat_0vs25.mask);
for jj = 1:island.NumObjects
    h = line(time{ans_seg}(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.7, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    b1 = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 4);
    b1.Color = '#008450';
    hold on
    b1 = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 4);
    b1.Color = '#EFB700';
end
clear island

island = bwconncomp(stat_25vs75.mask);
for jj = 1:island.NumObjects
    h = line(time{ans_seg}(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.5, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    b2 = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 4);
    b2.Color = '#EFB700';
    hold on
    b2 = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 4);
    b2.Color = '#B81D13';
end
clear island

a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [time{ans_seg}(1) time{ans_seg}(end)], ...
    'Ylim', [a(1)-0.1 a(2)+0.1]); box off; hold on;
set(gcf, 'Position', [491 320 201 175]);


