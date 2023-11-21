%% Figure S5 a & b

%%

clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_data = '/Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/HFA/ExpVariance_Approach/HFA_state_space_dynamics';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_S5/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% Load Data

% select condition of interest
selcond = input('Select Data: (1) HLL, (2) BR: ');

if selcond == 1
    load(fullfile(path_data, 'state_space_data_hll_all_elecs.mat'));
else
    load(fullfile(path_data, 'state_space_data_br_all_elecs.mat'));
end
   
% load single subject to get time
load /Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/OSL24/Data/HFA/PCA/PCA_ROI_Data/OSL24_activedata_ROI_pca.mat

%% Define Colors

color = {'#008450', '#EFB700', '#B81D13'};

%% Plot grand average in low dimensional space

% full time
if selcond == 1
    time = roidata.hll.frontal.active.time;
else
    time = roidata.br.frontal.active.time;
end

% get time zero
t0 = nearest(time,0);

for Isub = 1:numel(subjects)
    
    % -----------------------
    % Frontal Cortex
    % -----------------------

    if isfield(sp_data{Isub}, 'frontal')
        
        % euclidean distance between 0 and 25%
        data_frontal.state.lh_0{Isub}.avg     = squeeze(mean(sp_data{Isub}.frontal.lowDdata.lh_0(1,:,:),2))';
        data_frontal.state.lh_0{Isub}.time    = time;
        data_frontal.state.lh_0{Isub}.dimord  = 'chan_time';
        data_frontal.state.lh_0{Isub}.label   = {'PC1'};
        data_frontal.state.lh_0{Isub}.cfg     = [];
           
        % euclidean distance between 0 and 75%
        data_frontal.state.lh_25{Isub}.avg     = squeeze(mean(sp_data{Isub}.frontal.lowDdata.lh_25(1,:,:),2))';
        data_frontal.state.lh_25{Isub}.time    = time;
        data_frontal.state.lh_25{Isub}.dimord  = 'chan_time';
        data_frontal.state.lh_25{Isub}.label   = {'PC1'};
        data_frontal.state.lh_25{Isub}.cfg     = [];
        
        % euclidean distance between 25 and 75%
        data_frontal.state.lh_75{Isub}.avg     = squeeze(mean(sp_data{Isub}.frontal.lowDdata.lh_75(1,:,:),2))';
        data_frontal.state.lh_75{Isub}.time    = time;
        data_frontal.state.lh_75{Isub}.dimord  = 'chan_time';
        data_frontal.state.lh_75{Isub}.label   = {'PC1'};
        data_frontal.state.lh_75{Isub}.cfg     = [];
                
    end
    
    % -----------------------
    % Motor Cortex
    % -----------------------

    if isfield(sp_data{Isub}, 'motor')
        
        % euclidean distance between 0 and 25%
        data_motor.state.lh_0{Isub}.avg     = squeeze(mean(sp_data{Isub}.motor.lowDdata.lh_0(1,:,:),2))';
        data_motor.state.lh_0{Isub}.time    = time;
        data_motor.state.lh_0{Isub}.dimord  = 'chan_time';
        data_motor.state.lh_0{Isub}.label   = {'PC1'};
        data_motor.state.lh_0{Isub}.cfg     = [];
           
        % euclidean distance between 0 and 75%
        data_motor.state.lh_25{Isub}.avg     = squeeze(mean(sp_data{Isub}.motor.lowDdata.lh_25(1,:,:),2))';
        data_motor.state.lh_25{Isub}.time    = time;
        data_motor.state.lh_25{Isub}.dimord  = 'chan_time';
        data_motor.state.lh_25{Isub}.label   = {'PC1'};
        data_motor.state.lh_25{Isub}.cfg     = [];
        
        % euclidean distance between 25 and 75%
        data_motor.state.lh_75{Isub}.avg     = squeeze(mean(sp_data{Isub}.motor.lowDdata.lh_75(1,:,:),2))';
        data_motor.state.lh_75{Isub}.time    = time;
        data_motor.state.lh_75{Isub}.dimord  = 'chan_time';
        data_motor.state.lh_75{Isub}.label   = {'PC1'};
        data_motor.state.lh_75{Isub}.cfg     = [];
                
    end
        
end


% remove empty cells
data_frontal.state.lh_0     = data_frontal.state.lh_0(~cellfun(@isempty, data_frontal.state.lh_0));
data_frontal.state.lh_25    = data_frontal.state.lh_25(~cellfun(@isempty, data_frontal.state.lh_25));
data_frontal.state.lh_75    = data_frontal.state.lh_75(~cellfun(@isempty, data_frontal.state.lh_75));

data_motor.state.lh_0       = data_motor.state.lh_0(~cellfun(@isempty, data_motor.state.lh_0));
data_motor.state.lh_25      = data_motor.state.lh_25(~cellfun(@isempty, data_motor.state.lh_25));
data_motor.state.lh_75      = data_motor.state.lh_75(~cellfun(@isempty, data_motor.state.lh_75));

%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_S5');
cd(fullfile(path_source, 'Fig_S5'));

sourcedata = struct;

time = data_frontal.state.lh_0{1}.time;

for Isub = 1:numel(data_frontal.state.lh_0)
    
    if isempty(data_frontal.state.lh_0{Isub})
        
        sourcedata.lh_0(Isub,:)  = NaN(1, length(time));
        sourcedata.lh_25(Isub,:) = NaN(1, length(time));
        sourcedata.lh_75(Isub,:) = NaN(1, length(time));

    else
        
        sourcedata.lh_0(Isub,:)  = data_frontal.state.lh_0{Isub}.avg(1,:);
        sourcedata.lh_25(Isub,:) = data_frontal.state.lh_25{Isub}.avg(1,:);
        sourcedata.lh_75(Isub,:) = data_frontal.state.lh_75{Isub}.avg(1,:);

    end        
        
end%Isub

% - STORE 0% LH - %
tmp = NaN(size(sourcedata.lh_0,1)+1, size(sourcedata.lh_0,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_frontal.lh_0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S5_frontal_PCA_lh0.xlsx';
writetable(table_frontal.lh_0,filename, 'Sheet', 1);

% - STORE 25% LH - %
tmp = NaN(size(sourcedata.lh_25,1)+1, size(sourcedata.lh_25,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_25;

table_frontal.lh_25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S5_frontal_PCA_lh25.xlsx';
writetable(table_frontal.lh_25,filename, 'Sheet', 1);

% - STORE 75% LH - %
tmp = NaN(size(sourcedata.lh_75,1)+1, size(sourcedata.lh_75,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_75;

table_frontal.lh_75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S5_frontal_PCA_lh75.xlsx';
writetable(table_frontal.lh_75,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_S5');
cd(fullfile(path_source, 'Fig_S5'));

sourcedata = struct;

time = data_motor.state.lh_0{1}.time;

for Isub = 1:numel(data_motor.state.lh_0)
    
    if isempty(data_motor.state.lh_0{Isub})
        
        sourcedata.lh_0(Isub,:)  = NaN(1, length(time));
        sourcedata.lh_25(Isub,:) = NaN(1, length(time));
        sourcedata.lh_75(Isub,:) = NaN(1, length(time));

    else
        
        sourcedata.lh_0(Isub,:)  = data_motor.state.lh_0{Isub}.avg(1,:);
        sourcedata.lh_25(Isub,:) = data_motor.state.lh_25{Isub}.avg(1,:);
        sourcedata.lh_75(Isub,:) = data_motor.state.lh_75{Isub}.avg(1,:);

    end        
        
end%Isub

% - STORE 0% LH - %
tmp = NaN(size(sourcedata.lh_0,1)+1, size(sourcedata.lh_0,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:14)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_motor.lh_0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S5_motor_PCA_lh0.xlsx';
writetable(table_motor.lh_0,filename, 'Sheet', 1);

% - STORE 25% LH - %
tmp = NaN(size(sourcedata.lh_25,1)+1, size(sourcedata.lh_25,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:14)';
tmp(2:end,2:end)  = sourcedata.lh_25;

table_motor.lh_25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S5_motor_PCA_lh25.xlsx';
writetable(table_motor.lh_25,filename, 'Sheet', 1);

% - STORE 75% LH - %
tmp = NaN(size(sourcedata.lh_75,1)+1, size(sourcedata.lh_75,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:14)';
tmp(2:end,2:end)  = sourcedata.lh_75;

table_motor.lh_75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S5_motor_PCA_lh75.xlsx';
writetable(table_motor.lh_75,filename, 'Sheet', 1);

%% Cluster-Based Permutation Stats 

% 1. First one-way ANOVA is performed (ft_statfun_depsamplesFunivariate)
% 2. If 1 is significant, post-hoc follow up contrasts are performed

% Get data of interest
q       = 'Statistic: Frontal = 1; Motor = 2: ';
ans_roi = input(q);

statsdata = struct();

if ans_roi == 1
    
    statsdata.lh_0   = data_frontal.state.lh_0;
    statsdata.lh_25  = data_frontal.state.lh_25;
    statsdata.lh_75  = data_frontal.state.lh_75;

elseif ans_roi == 2
    
    statsdata.lh_0   = data_motor.state.lh_0;
    statsdata.lh_25  = data_motor.state.lh_25;
    statsdata.lh_75  = data_motor.state.lh_75;

end
    
% Run stats

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = 'all';
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
subj = numel(statsdata.lh_0);
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

fstat = ft_timelockstatistics(cfg, statsdata.lh_0{:}, statsdata.lh_25{:}, statsdata.lh_75{:});    

% check positive clusters

if isfield(fstat, 'posclusters')

    for Icluster = 1:numel(fstat.posclusters)

        fprintf('%%---------------------------------------------------\n');
        fprintf('%% POSITIVE CLUSTERS---------------------------------\n');
        fprintf('%%---------------------------------------------------\n');

        fprintf('\n');

        fprintf('p-value cluster %d = %.3f\n', Icluster, fstat.posclusters(Icluster).prob);
        fprintf('clustermass cluster %d = %.3f\n', Icluster, fstat.posclusters(Icluster).clusterstat);

        begeffect = fstat.time...
            (find(fstat.posclusterslabelmat == Icluster, 1, 'first'));
        endeffect = fstat.time...
            (find(fstat.posclusterslabelmat == Icluster, 1, 'last'));

        fprintf('cluster %d from %.3f to %.3f ms\n', Icluster, begeffect, endeffect);

        fprintf('\n');

    end

end

% %% Plot F-Stats
% 
% plot(fstat.time, fstat.stat, 'k', 'LineWidth', 2);
% xlim([fstat.time(1) fstat.time(end)])
% ylabel('Stats [F]', 'FontSize', 15);
% xlabel('Time to BR');
% set(gca, 'fontsize', 13, 'linewidth', 3);
% box off
% set(gcf, 'Position', [639 347 328 249]);

%% Perform Follow-Up Comparisons

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

subj = numel(statsdata.lh_0);
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

stat = cell(3,1);
stat{1}  = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_0{:});
stat{2}  = ft_timelockstatistics(cfg, statsdata.lh_25{:}, statsdata.lh_0{:});
stat{3}  = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_25{:});

%% Plot Data

% get current grand averages

cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter      = 'avg';

tmpGA_lh0     = ft_timelockgrandaverage(cfg, statsdata.lh_0{:});
tmpGA_lh25    = ft_timelockgrandaverage(cfg, statsdata.lh_25{:});
tmpGA_lh75    = ft_timelockgrandaverage(cfg, statsdata.lh_75{:});

% plot data

a = shadedErrorBar(tmpGA_lh0.time, squeeze(mean(tmpGA_lh0.individual(:,1,:))), ...
    std(squeeze(tmpGA_lh0.individual(:,1,:))) / sqrt(size(tmpGA_lh0.individual,1)));
a.mainLine.Color = color{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{1}; a.patch.EdgeColor = color{1}; a.patch.FaceAlpha = 0.5;
hold on
a = shadedErrorBar(tmpGA_lh25.time, squeeze(mean(tmpGA_lh25.individual(:,1,:))), ...
    std(squeeze(tmpGA_lh25.individual(:,1,:))) / sqrt(size(tmpGA_lh25.individual,1)));
a.mainLine.Color = color{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{2}; a.patch.EdgeColor = color{2}; a.patch.FaceAlpha = 0.5;
hold on
a = shadedErrorBar(tmpGA_lh75.time, squeeze(mean(tmpGA_lh75.individual(:,1,:))), ...
    std(squeeze(tmpGA_lh75.individual(:,1,:))) / sqrt(size(tmpGA_lh75.individual,1)));
a.mainLine.Color = color{3}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{3}; a.patch.EdgeColor = color{3}; a.patch.FaceAlpha = 0.5;

ylabel('Population Activity');

if selcond == 1  
    xlabel('Time to HLL');
else
    xlabel('Time to BR');
end

% plot significant timepoints
a = ylim;

hold on

% mask f-stats
island = bwconncomp(fstat.mask);
for jj = 1:island.NumObjects
    h = line(time(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.5, 1, numel(island.PixelIdxList{jj})), 'linewidth', 5);
    plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', 'k', 'LineWidth', 5);
    hold on
    plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', 'k', 'LineWidth', 5);
end 
clear island

if isfield(fstat, 'posclusters') && fstat.posclusters(1).prob < 0.05

    % mask 0 vs. 75%
    island = bwconncomp(stat{1}.mask);
    for jj = 1:island.NumObjects
        h = line(time(island.PixelIdxList{jj}), ...
            repmat(a(2)+1.5, 1, numel(island.PixelIdxList{jj})), 'linewidth', 5);
        plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#008450', 'LineWidth', 5);
        hold on
        plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#B81D13', 'LineWidth', 5);
    end
    clear island

    % mask 0 vs. 25%
    island = bwconncomp(stat{2}.mask);
    for jj = 1:island.NumObjects
        h = line(time(island.PixelIdxList{jj}), ...
            repmat(a(2)+2.5, 1, numel(island.PixelIdxList{jj})), 'linewidth', 5);
        plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#008450', 'LineWidth', 5);
        hold on
        plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#EFB700', 'LineWidth', 5);
    end
    clear island

    % mask 25 vs. 75%
    island = bwconncomp(stat{3}.mask);
    for jj = 1:island.NumObjects
        h = line(time(island.PixelIdxList{jj}), ...
            repmat(a(2)+3.5, 1, numel(island.PixelIdxList{jj})), 'linewidth', 5);
        plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#EFB700', 'LineWidth', 5);
        hold on
        plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#B81D13', 'LineWidth', 5);
    end
    clear island
    
end

a = ylim;

set(gca, 'linewidth', 1, 'FontSize',13, 'xlim', [time(1) time(end)]); box off; hold on;
set(gcf, 'Position', [471 330 486 325]);

% save data
cd(path_out);

if ans_roi == 1
    print(gcf,'PC1_PFC_Locked_HLL.pdf','-dpdf','-r400');  
else
    print(gcf,'PC1_M1_Locked_HLL.pdf','-dpdf','-r400');  
end

