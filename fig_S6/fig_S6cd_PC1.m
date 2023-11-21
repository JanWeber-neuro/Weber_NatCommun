%% Figure S6c & d (left panel)
% this script rungs additional analysis required for the plots of figure
% S6c & d (left panel)

%%

clear
ft_defaults

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_data = '/Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/HFA/ExpVariance_Approach/HFA_state_space_dynamics';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_5/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% - LOAD DATA - %%

load(fullfile(path_data, 'state_space_data_br_all_elecs.mat'));

% load single subject to get time
load /Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/OSL24/Data/HFA/PCA/PCA_ROI_Data/OSL24_activedata_ROI_pca.mat

%% - EXTRACT DATA - %%

% full time
time = roidata.br.frontal.active.time;

fsample       = 512;

% get time zero
t0 = nearest(time,0);

for Isub = 1:numel(subjects)
    
    % -----------------------
    % Frontal Cortex
    % -----------------------

    if isfield(sp_data{Isub}, 'frontal')
        
        % PCs 1-3 for 0% 
        data_frontal.state.lh_0{Isub}.avg     = squeeze(mean(sp_data{Isub}.frontal.lowDdata.lh_0,2));
        data_frontal.state.lh_0{Isub}.time    = time;
        data_frontal.state.lh_0{Isub}.dimord  = 'chan_time';
        data_frontal.state.lh_0{Isub}.label   = {'PC1', 'PC2', 'PC3'};
        data_frontal.state.lh_0{Isub}.cfg     = [];
           
        % PCs 1-3 for 25% 
        data_frontal.state.lh_25{Isub}.avg     = squeeze(mean(sp_data{Isub}.frontal.lowDdata.lh_25,2));
        data_frontal.state.lh_25{Isub}.time    = time;
        data_frontal.state.lh_25{Isub}.dimord  = 'chan_time';
        data_frontal.state.lh_25{Isub}.label   = {'PC1', 'PC2', 'PC3'};
        data_frontal.state.lh_25{Isub}.cfg     = [];
        
        % PCs 1-3 for 75% 
        data_frontal.state.lh_75{Isub}.avg     = squeeze(mean(sp_data{Isub}.frontal.lowDdata.lh_75,2));
        data_frontal.state.lh_75{Isub}.time    = time;
        data_frontal.state.lh_75{Isub}.dimord  = 'chan_time';
        data_frontal.state.lh_75{Isub}.label   = {'PC1', 'PC2', 'PC3'};
        data_frontal.state.lh_75{Isub}.cfg     = [];
                
    end
    
    % -----------------------
    % Motor Cortex
    % -----------------------

    if isfield(sp_data{Isub}, 'motor')
        
        % PCs 1-3 for 0% 
        data_motor.state.lh_0{Isub}.avg     = squeeze(mean(sp_data{Isub}.motor.lowDdata.lh_0,2));
        data_motor.state.lh_0{Isub}.time    = time;
        data_motor.state.lh_0{Isub}.dimord  = 'chan_time';
        data_motor.state.lh_0{Isub}.label   = {'PC1', 'PC2', 'PC3'};
        data_motor.state.lh_0{Isub}.cfg     = [];
           
        % PCs 1-3 for 25% 
        data_motor.state.lh_25{Isub}.avg     = squeeze(mean(sp_data{Isub}.motor.lowDdata.lh_25,2));
        data_motor.state.lh_25{Isub}.time    = time;
        data_motor.state.lh_25{Isub}.dimord  = 'chan_time';
        data_motor.state.lh_25{Isub}.label   = {'PC1', 'PC2', 'PC3'};
        data_motor.state.lh_25{Isub}.cfg     = [];
        
        % PCs 1-3 for 75% 
        data_motor.state.lh_75{Isub}.avg     = squeeze(mean(sp_data{Isub}.motor.lowDdata.lh_75,2));
        data_motor.state.lh_75{Isub}.time    = time;
        data_motor.state.lh_75{Isub}.dimord  = 'chan_time';
        data_motor.state.lh_75{Isub}.label   = {'PC1', 'PC2', 'PC3'};
        data_motor.state.lh_75{Isub}.cfg     = [];
                
    end
    
    % -----------------------
    % Frontal % Motor Cortex
    % -----------------------

    if isfield(sp_data{Isub}, 'frontal') && isfield(sp_data{Isub}, 'motor')
        
        % difference PC1 75% - 0%
        data_frontal.state.diff{Isub}.avg     = data_frontal.state.lh_75{Isub}.avg - ...
            data_frontal.state.lh_0{Isub}.avg;
        data_frontal.state.diff{Isub}.time    = time;
        data_frontal.state.diff{Isub}.dimord  = 'chan_time';
        data_frontal.state.diff{Isub}.label   = {'PC1', 'PC2', 'PC3'};
        data_frontal.state.diff{Isub}.cfg     = [];
        
        % difference PC1 75% - 0%
        data_motor.state.diff{Isub}.avg     = data_motor.state.lh_75{Isub}.avg - ...
            data_motor.state.lh_0{Isub}.avg;
        data_motor.state.diff{Isub}.time    = time;
        data_motor.state.diff{Isub}.dimord  = 'chan_time';
        data_motor.state.diff{Isub}.label   = {'PC1', 'PC2', 'PC3'};
        data_motor.state.diff{Isub}.cfg     = [];
                        
    end    
        
end

%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_5c_left');
cd(fullfile(path_source, 'Fig_5c_left'));

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
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_frontal.lh_0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5c_frontal_PCA_lh0.xlsx';
writetable(table_frontal.lh_0,filename, 'Sheet', 1);

% - STORE 25% LH - %
tmp = NaN(size(sourcedata.lh_25,1)+1, size(sourcedata.lh_25,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.lh_25;

table_frontal.lh_25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5c_frontal_PCA_lh25.xlsx';
writetable(table_frontal.lh_25,filename, 'Sheet', 1);

% - STORE 75% LH - %
tmp = NaN(size(sourcedata.lh_75,1)+1, size(sourcedata.lh_75,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:18)';
tmp(2:end,2:end)  = sourcedata.lh_75;

table_frontal.lh_75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5c_frontal_PCA_lh75.xlsx';
writetable(table_frontal.lh_75,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_5d_left');
cd(fullfile(path_source, 'Fig_5d_left'));

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
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_motor.lh_0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5d_motor_PCA_lh0.xlsx';
writetable(table_motor.lh_0,filename, 'Sheet', 1);

% - STORE 25% LH - %
tmp = NaN(size(sourcedata.lh_25,1)+1, size(sourcedata.lh_25,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_25;

table_motor.lh_25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5d_motor_PCA_lh25.xlsx';
writetable(table_motor.lh_25,filename, 'Sheet', 1);

% - STORE 75% LH - %
tmp = NaN(size(sourcedata.lh_75,1)+1, size(sourcedata.lh_75,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_75;

table_motor.lh_75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_5d_motor_PCA_lh75.xlsx';
writetable(table_motor.lh_75,filename, 'Sheet', 1);

%%

% remove empty cells
data_frontal.state.lh_0     = data_frontal.state.lh_0(~cellfun(@isempty, data_frontal.state.lh_0));
data_frontal.state.lh_25    = data_frontal.state.lh_25(~cellfun(@isempty, data_frontal.state.lh_25));
data_frontal.state.lh_75    = data_frontal.state.lh_75(~cellfun(@isempty, data_frontal.state.lh_75));
data_frontal.state.diff     = data_frontal.state.diff(~cellfun(@isempty, data_frontal.state.diff));

data_motor.state.lh_0       = data_motor.state.lh_0(~cellfun(@isempty, data_motor.state.lh_0));
data_motor.state.lh_25      = data_motor.state.lh_25(~cellfun(@isempty, data_motor.state.lh_25));
data_motor.state.lh_75      = data_motor.state.lh_75(~cellfun(@isempty, data_motor.state.lh_75));
data_motor.state.diff       = data_motor.state.diff(~cellfun(@isempty, data_motor.state.diff));

% get grand average
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter      = 'avg';

GA_frontal = struct();
GA_motor   = struct();
GA_mtl     = struct();

% frontal grand average
GA_frontal.lowD.lh_0   = ft_timelockgrandaverage(cfg, data_frontal.state.lh_0{:});
GA_frontal.lowD.lh_25  = ft_timelockgrandaverage(cfg, data_frontal.state.lh_25{:});
GA_frontal.lowD.lh_75  = ft_timelockgrandaverage(cfg, data_frontal.state.lh_75{:});

% motor grand average
GA_motor.lowD.lh_0     = ft_timelockgrandaverage(cfg, data_motor.state.lh_0{:});
GA_motor.lowD.lh_25    = ft_timelockgrandaverage(cfg, data_motor.state.lh_25{:});
GA_motor.lowD.lh_75    = ft_timelockgrandaverage(cfg, data_motor.state.lh_75{:});

%% - CLUSTER STATS ON PC1 - %%

% Get data of interest
q       = 'Statistic: Frontal = 1; Motor = 2: ';
selroi = input(q);

statsdata = struct();

if selroi == 1
    
    statsdata.lh_0   = data_frontal.state.lh_0;
    statsdata.lh_25  = data_frontal.state.lh_25;
    statsdata.lh_75  = data_frontal.state.lh_75;

elseif selroi == 2
    
    statsdata.lh_0   = data_motor.state.lh_0;
    statsdata.lh_25  = data_motor.state.lh_25;
    statsdata.lh_75  = data_motor.state.lh_75;

end
    

% ################################
% F-Statistic
% ################################

cfg                     = [];
cfg.channel             = 'PC1';
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

% extract start point of cluster
if isfield(fstat, 'posclusters')
    
    for Icluster = 1:numel(fstat.posclusters)
        
        fprintf('\n');
        
        fprintf('p-value cluster %d = %.3f\n', Icluster, fstat.posclusters(Icluster).prob);
        
        begeffect = fstat.time(find(fstat.posclusterslabelmat == Icluster, 1, 'first'));
        endeffect = fstat.time(find(fstat.posclusterslabelmat == Icluster, 1, 'last'));
        
        fprintf('Effect cluster %d from %.3f to %.3f ms\n', Icluster, begeffect, endeffect);
        
        fprintf('\n');

    end
    
end

% ################################
% Post-Hoc Pairwise Comparisons
% ################################

cfg                     = [];
cfg.channel             = 'PC1';
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

stat_0vs25   = ft_timelockstatistics(cfg, statsdata.lh_25{:}, statsdata.lh_0{:});
stat_0vs75   = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_0{:});
stat_25vs75  = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_25{:});

% compute cohen's d as effect size
cfg            = [];
cfg.channel    = 'PC1';
cfg.method     = 'analytic';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

d_0vs25        = ft_timelockstatistics(cfg, statsdata.lh_25{:}, statsdata.lh_0{:});
d_0vs75        = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_0{:});
d_25vs75       = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_25{:});

% get grand average
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter      = 'avg';

% frontal grand average
tmp_lh0   = ft_timelockgrandaverage(cfg, statsdata.lh_0{:});
tmp_lh25  = ft_timelockgrandaverage(cfg, statsdata.lh_25{:});
tmp_lh75  = ft_timelockgrandaverage(cfg, statsdata.lh_75{:});

% colors for plotting
color = {'#008450', '#EFB700', '#B81D13'};

figure;

% plot first PC frontal for different conditions
a = shadedErrorBar(tmp_lh0.time, squeeze(mean(tmp_lh0.individual(:,1,:))), ...
    std(squeeze(tmp_lh0.individual(:,1,:))) / sqrt(size(tmp_lh0.individual,1)));
a.mainLine.Color = color{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{1}; a.patch.EdgeColor = color{1}; a.patch.FaceAlpha = 0.5;
hold on
a = shadedErrorBar(tmp_lh25.time, squeeze(mean(tmp_lh25.individual(:,1,:))), ...
    std(squeeze(tmp_lh25.individual(:,1,:))) / sqrt(size(tmp_lh25.individual,1)));
a.mainLine.Color = color{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{2}; a.patch.EdgeColor = color{2}; a.patch.FaceAlpha = 0.5;
hold on
a = shadedErrorBar(tmp_lh75.time, squeeze(mean(tmp_lh75.individual(:,1,:))), ...
    std(squeeze(tmp_lh75.individual(:,1,:))) / sqrt(size(tmp_lh75.individual,1)));
a.mainLine.Color = color{3}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = color{3}; a.patch.EdgeColor = color{3}; a.patch.FaceAlpha = 0.5;

% plot significant timepoints
a = ylim;

island = bwconncomp(fstat.mask);
for jj = 1:island.NumObjects
    h = line(tmp_lh0.time(island.PixelIdxList{jj}), ...
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
    h = line(tmp_lh25.time(island.PixelIdxList{jj}), ...
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
    h = line(tmp_lh75.time(island.PixelIdxList{jj}), ...
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
    h = line(tmp_lh75.time(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.5, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    b2 = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 4);
    b2.Color = '#EFB700';
    hold on
    b2 = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 4);
    b2.Color = '#B81D13';
end
clear island

if selroi == 1
    ylim([-3.5 7]);
end

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [tmp_lh0.time(1) tmp_lh0.time(end)]); box off; hold on;

set(gcf, 'Position', [530 333 348 277]);

allposthoc = {stat_0vs75, stat_0vs25, stat_25vs75};
allcohens  = {d_0vs75, d_0vs25, d_25vs75};
posthocidx = {'0 vs 75', '0 vs 25', '25 vs 75'};

% check positive clusters
for Ipair = 1:numel(allposthoc)

    % extract start point of cluster
    if isfield(allposthoc{Ipair}, 'posclusters')

        for Icluster = 1:numel(allposthoc{Ipair}.posclusters)

            fprintf('%%---------------------------------------------------\n');
            fprintf('%% POSITIVE CLUSTERS---------------------------------\n');
            fprintf('%%---------------------------------------------------\n');
            
            fprintf('\n');

            fprintf('p-value cluster %d = %.3f\n', Icluster, allposthoc{Ipair}.posclusters(Icluster).prob);
            fprintf('clustermass cluster %d = %.3f\n', Icluster, allposthoc{Ipair}.posclusters(Icluster).clusterstat);
            
            begeffect = allposthoc{Ipair}.time...
                (find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'first'));
            endeffect = allposthoc{Ipair}.time...
                (find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'last'));

            fprintf('Effect (%s): cluster %d from %.3f to %.3f ms\n', posthocidx{Ipair}, Icluster, begeffect, endeffect);

            fprintf('Cohens d cluster %d = %.3f\n', Icluster, mean(allcohens{Ipair}.cohensd(...
                find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'first'):...
                find(allposthoc{Ipair}.posclusterslabelmat == Icluster, 1, 'last'))));
            
            fprintf('\n');

        end

    end
    
end
    
% check negative clusters
for Ipair = 1:numel(allposthoc)

    % extract start point of cluster
    if isfield(allposthoc{Ipair}, 'negclusters')

        for Icluster = 1:numel(allposthoc{Ipair}.negclusters)
            
            fprintf('%%---------------------------------------------------\n');
            fprintf('%% NEGATIVE CLUSTERS---------------------------------\n');
            fprintf('%%---------------------------------------------------\n');
            
            fprintf('\n');

            fprintf('p-value cluster %d = %.3f\n', Icluster, allposthoc{Ipair}.negclusters(Icluster).prob);
            fprintf('clustermass cluster %d = %.3f\n', Icluster, allposthoc{Ipair}.negclusters(Icluster).clusterstat);
            
            begeffect = allposthoc{Ipair}.time...
                (find(allposthoc{Ipair}.negclusterslabelmat == Icluster, 1, 'first'));
            endeffect = allposthoc{Ipair}.time...
                (find(allposthoc{Ipair}.negclusterslabelmat == Icluster, 1, 'last'));

            fprintf('Effect (%s): cluster %d from %.3f to %.3f ms\n', posthocidx{Ipair}, Icluster, begeffect, endeffect);

            fprintf('Cohens d cluster %d = %.3f\n', Icluster, mean(allcohens{Ipair}.cohensd(...
                find(allposthoc{Ipair}.negclusterslabelmat == Icluster, 1, 'first'):...
                find(allposthoc{Ipair}.negclusterslabelmat == Icluster, 1, 'last'))));
            
            fprintf('\n');

        end

    end
    
end
    
cd(path_out);

if selroi == 1
    print(gcf,'fig_5_frontal_pc1.pdf','-dpdf','-r400');    
else
    print(gcf,'fig_5_motor_pc1.pdf','-dpdf','-r400');    
end


%% - COMPUTE INTERACTION EFFECT FRONTAL VS. MOTOR - %%

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

subj = numel(data_frontal.state.diff);
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

stat_int  = ft_timelockstatistics(cfg, data_frontal.state.diff{:}, data_motor.state.diff{:});

% compute cohen's d as effect size
cfg            = [];
cfg.method     = 'analytic';
cfg.parameter  = 'avg';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

cohensd        = ft_timelockstatistics(cfg, data_frontal.state.diff{:}, data_motor.state.diff{:});

% extract start point of cluster
if isfield(stat_int, 'posclusters')

    for Icluster = 1:numel(stat_int.posclusters)

        fprintf('\n');

        fprintf('p-value cluster %d = %.3f\n', Icluster, stat_int.posclusters(Icluster).prob);
        fprintf('clustermass cluster %d = %.3f\n', Icluster, stat_int.posclusters(Icluster).clusterstat);

        begeffect = stat_int.time...
            (find(stat_int.posclusterslabelmat == Icluster, 1, 'first'));
        endeffect = stat_int.time...
            (find(stat_int.posclusterslabelmat == Icluster, 1, 'last'));

        fprintf('Interaction Effect: cluster %d from %.3f to %.3f ms\n', Icluster, begeffect, endeffect);

        fprintf('Cohens d cluster %d = %.3f\n', Icluster, mean(cohensd.cohensd(...
            find(stat_int.posclusterslabelmat == Icluster, 1, 'first'):...
            find(stat_int.posclusterslabelmat == Icluster, 1, 'last'))));

        fprintf('\n');

    end

end

% plot data
plot(stat_int.time, stat_int.stat, 'k', 'LineWidth', 2);

hold on

% shade grey where significant
sigpnts = stat_int.time(find(stat_int.mask));
x_lim   = xlim;
y_lim   = [0 stat_int.critval(2)]; % based on stats in frontal cortex

if ~isempty(sigpnts)
    pgon = polyshape([sigpnts(1) sigpnts(1) (sigpnts(end) - sigpnts(1)) (sigpnts(end) - sigpnts(1))], ...
        [y_lim(end) y_lim(1) y_lim(1) y_lim(end)]);
    h = plot(pgon);
    h.FaceAlpha = 0.3;
    h.FaceColor = [0.5 0.5 0.5];
    h.EdgeColor = [0.5 0.5 0.5];
    h.EdgeAlpha = 0.3;
end

ylabel('Stats [t]', 'FontSize', 13);
xlabel('Time to BR [s]', 'FontSize', 13);
set(gca, 'fontsize', 13, 'linewidth', 3);
box off
set(gcf, 'Position', [552 396 205 97]);

%% - PLOT FIRST PCS AGAINST EACH OTHER - %%

figure;

% get timing at button release (t0)
trange = 1:nearest(GA_frontal.lowD.lh_0.time, 0.3);

% extract low dimensional features for 0% condition
lowfrontal_lh0 = squeeze(GA_frontal.lowD.lh_0.individual(:,1,trange));
lowmotor_lh0 = squeeze(GA_motor.lowD.lh_0.individual(:,1,trange));

% extract low dimensional features for 25% condition
lowfrontal_lh25 = squeeze(GA_frontal.lowD.lh_25.individual(:,1,trange));
lowmotor_lh25 = squeeze(GA_motor.lowD.lh_25.individual(:,1,trange));

% extract low dimensional features for 75% condition
lowfrontal_lh75 = squeeze(GA_frontal.lowD.lh_75.individual(:,1,trange));
lowmotor_lh75 = squeeze(GA_motor.lowD.lh_75.individual(:,1,trange));

% ----------
% Plot 0%
% ----------
plot(smoothdata(mean(lowfrontal_lh0), 'gaussian', round(fsample*0.05)), ...
      smoothdata(mean(lowmotor_lh0), 'gaussian', round(fsample*0.05)), ...
      'Color', color{1}, 'LineWidth', 1);
   
hold on
x1 = smoothdata(mean(lowfrontal_lh0), 'gaussian', round(fsample*0.05));
x2 = smoothdata(mean(lowmotor_lh0), 'gaussian', round(fsample*0.05));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{1}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{1}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);   
hold on

% ----------
% Plot 25%
% ----------
plot(smoothdata(mean(lowfrontal_lh25), 'gaussian', round(fsample*0.05)), ...
      smoothdata(mean(lowmotor_lh25), 'gaussian', round(fsample*0.05)), ...
      'Color', color{2}, 'LineWidth', 1);
   
hold on
x1 = smoothdata(mean(lowfrontal_lh25), 'gaussian', round(fsample*0.05));
x2 = smoothdata(mean(lowmotor_lh25), 'gaussian', round(fsample*0.05));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{2}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{1}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);   

hold on


% ----------
% Plot 75%
% ----------
plot(smoothdata(mean(lowfrontal_lh75), 'gaussian', round(fsample*0.05)), ...
      smoothdata(mean(lowmotor_lh75), 'gaussian', round(fsample*0.05)), ...
      'Color', color{3}, 'LineWidth', 1);
   
hold on
x1 = smoothdata(mean(lowfrontal_lh75), 'gaussian', round(fsample*0.05));
x2 = smoothdata(mean(lowmotor_lh75), 'gaussian', round(fsample*0.05));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{3}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{3}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);   

set(gca, 'linewidth', 3, 'fontsize', 15);
xlabel('Frontal');
ylabel('Motor');
xlim([-4 6]);
box off

set(gcf, 'Position', [530 333 348 277]);

view(0,90)

print(gcf,'fig_5_frontal_motor_interaction.pdf','-dpdf','-r400');    

%% - PLOT FRONTAL CORTEX PC1 TO PC2 - %%

% define smoothing parameter for plotting
smoothwin = 0.05;

trange = 1:nearest(time, 0.3);

figure;

% ----------
% Plot 0%
% ----------
plot(smoothdata(mean(squeeze(GA_frontal.lowD.lh_0.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin)), ...
      smoothdata(mean(squeeze(GA_frontal.lowD.lh_0.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin)), ...
      'Color', color{1}, 'LineWidth', 1);

hold on

x1 = smoothdata(mean(squeeze(GA_frontal.lowD.lh_0.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin));
x2 = smoothdata(mean(squeeze(GA_frontal.lowD.lh_0.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{1}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{1}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);

% ----------
% Plot 25%
% ----------
plot(smoothdata(mean(squeeze(GA_frontal.lowD.lh_25.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin)), ...
      smoothdata(mean(squeeze(GA_frontal.lowD.lh_25.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin)), ...
      'Color', color{2}, 'LineWidth', 1);

hold on

x1 = smoothdata(mean(squeeze(GA_frontal.lowD.lh_25.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin));
x2 = smoothdata(mean(squeeze(GA_frontal.lowD.lh_25.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{2}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{2}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);

hold on

% ----------
% Plot 75%
% ----------
plot(smoothdata(mean(squeeze(GA_frontal.lowD.lh_75.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin)), ...
      smoothdata(mean(squeeze(GA_frontal.lowD.lh_75.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin)), ...
      'Color', color{3}, 'LineWidth', 1);

hold on

x1 = smoothdata(mean(squeeze(GA_frontal.lowD.lh_75.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin));
x2 = smoothdata(mean(squeeze(GA_frontal.lowD.lh_75.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{3}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{3}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);

set(gca, 'linewidth', 3, 'fontsize', 13);
xlabel('Dim. 1');
ylabel('Dim. 2');
box off

set(gcf, 'Position', [552 342 205 151]);

xlim([-5 5]);
ylim([-1 1]);

print(gcf,'fig_5_frontal_PC1_PC2.pdf','-dpdf','-r400');    

%% - PLOT MOTOR CORTEX PC1 TO PC2 - %%

% define smoothing parameter for plotting
smoothwin = 0.05;

trange = 1:nearest(time, 0.3);

figure;

% ----------
% Plot 0%
% ----------
plot(smoothdata(mean(squeeze(GA_motor.lowD.lh_0.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin)), ...
      smoothdata(mean(squeeze(GA_motor.lowD.lh_0.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin)), ...
      'Color', color{1}, 'LineWidth', 1);

hold on

x1 = smoothdata(mean(squeeze(GA_motor.lowD.lh_0.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin));
x2 = smoothdata(mean(squeeze(GA_motor.lowD.lh_0.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{1}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{1}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);

% ----------
% Plot 25%
% ----------
plot(smoothdata(mean(squeeze(GA_motor.lowD.lh_25.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin)), ...
      smoothdata(mean(squeeze(GA_motor.lowD.lh_25.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin)), ...
      'Color', color{2}, 'LineWidth', 1);

hold on

x1 = smoothdata(mean(squeeze(GA_motor.lowD.lh_25.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin));
x2 = smoothdata(mean(squeeze(GA_motor.lowD.lh_25.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{2}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{2}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);

hold on

% ----------
% Plot 75%
% ----------
plot(smoothdata(mean(squeeze(GA_motor.lowD.lh_75.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin)), ...
      smoothdata(mean(squeeze(GA_motor.lowD.lh_75.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin)), ...
      'Color', color{3}, 'LineWidth', 1);

hold on

x1 = smoothdata(mean(squeeze(GA_motor.lowD.lh_75.individual(:,1,trange))), 'gaussian', round(fsample*smoothwin));
x2 = smoothdata(mean(squeeze(GA_motor.lowD.lh_75.individual(:,2,trange))), 'gaussian', round(fsample*smoothwin));
% mark start
scatter(x1(1), x2(1), 50, 'filled', ...
    'MarkerFaceColor', color{3}, 'MarkerFaceAlpha', 1);
% mark end
scatter(x1(t0), x2(t0), 150, ...
    'MarkerEdgeColor', color{3}, 'MarkerFaceAlpha', 0.7, 'Marker', '+', 'LineWidth', 1);

set(gca, 'linewidth', 3, 'fontsize', 13);
xlabel('Dim. 1');
ylabel('Dim. 2');
box off

set(gcf, 'Position', [552 342 205 151]);

% save data
print(gcf,'fig_5_motor_PC1_PC2.pdf','-dpdf','-r400');    


