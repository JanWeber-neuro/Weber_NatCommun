%% Figure 3c & d
% this script rungs additional analysis required for the plots of figure
% 3c & 3d.

%%

clear 
ft_defaults

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_3/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/ % add functions

%% - COLORMAP - %%

cmap = flipud(cbrewer('div', 'Spectral', 2000, 'PCHIP'));

%% - LOAD ROI INFORMATION - %%

% electrode information containing electrode index for ROIs
load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/HFA/ExpVariance_Approach/TaskInformative_Electrodes/ROIs_active_elecs.mat

%% - LOAD DATA - %%

all_tfr = struct();

for subji = 1:numel(subjects)
    
    % get subject ID
    subID = subjects{subji};
    
    fprintf('Load TFR Data %s\n', subID);   
    
    % get ROI information and extract electrode indices
    
    % frontal ROI
    if isfield(ROIs(subji).hll, 'frontal')
        if isfield(ROIs(subji).hll.frontal, 'active')
            elecs.frontal = ROIs(subji).hll.frontal.active.elecidx; % electrode index
        else
        elecs.frontal = [];
        end
    else
        elecs.frontal = [];
    end
    
    % motor ROI
    if isfield(ROIs(subji).hll, 'motor')
        if isfield(ROIs(subji).hll.motor, 'active')
            elecs.motor = ROIs(subji).hll.motor.active.elecidx;
        else
        elecs.motor = [];
        end
    else
        elecs.motor = [];
    end
    
    % go into folder and load data
    cd([path_in filesep subID filesep 'Data' filesep 'TFR']);
    
    % TFR for HLL (z-scored using baseline bootstrapping)
    load([subID '_HLL_TFR.mat']);
    
    % rename for simplicity
    tfr_hll = freq_task_hll; clear freq_task_hll
    
    % get trial information
           
    idx_0  = find(tfr_hll.trialinfo(:,13) == 0);
    idx_25 = find(tfr_hll.trialinfo(:,13) == 25);
    idx_75 = find(tfr_hll.trialinfo(:,13) == 75);    
    
    % check whether only one channel is task-informative in a region (get
    % matrix of 1 x Frequency x Time (average across elecs and trial)
    if ~isempty(elecs.frontal)
        if size(elecs.frontal,2) ~= 1
            tfr_hll.powspctrm_frontal_lh0   = nanmean(squeeze(nanmean(tfr_hll.powspctrm(idx_0,elecs.frontal,:,:),1)),1);
            tfr_hll.powspctrm_frontal_lh25  = nanmean(squeeze(nanmean(tfr_hll.powspctrm(idx_25,elecs.frontal,:,:),1)),1);
            tfr_hll.powspctrm_frontal_lh75  = nanmean(squeeze(nanmean(tfr_hll.powspctrm(idx_75,elecs.frontal,:,:),1)),1);
        else
            tfr_hll.powspctrm_frontal_lh0   = nanmean(squeeze(tfr_hll.powspctrm(idx_0,elecs.frontal,:,:)),1);
            tfr_hll.powspctrm_frontal_lh25  = nanmean(squeeze(tfr_hll.powspctrm(idx_25,elecs.frontal,:,:)),1);
            tfr_hll.powspctrm_frontal_lh75  = nanmean(squeeze(tfr_hll.powspctrm(idx_75,elecs.frontal,:,:)),1);
        end
    end
        
    if ~isempty(elecs.motor)
        if size(elecs.motor,2) ~= 1
            tfr_hll.powspctrm_motor_lh0     = nanmean(squeeze(nanmean(tfr_hll.powspctrm(idx_0,elecs.motor,:,:),1)),1);
            tfr_hll.powspctrm_motor_lh25    = nanmean(squeeze(nanmean(tfr_hll.powspctrm(idx_25,elecs.motor,:,:),1)),1);
            tfr_hll.powspctrm_motor_lh75    = nanmean(squeeze(nanmean(tfr_hll.powspctrm(idx_75,elecs.motor,:,:),1)),1);
        else
            tfr_hll.powspctrm_motor_lh0     = nanmean(squeeze(tfr_hll.powspctrm(idx_0,elecs.motor,:,:)),1);
            tfr_hll.powspctrm_motor_lh25    = nanmean(squeeze(tfr_hll.powspctrm(idx_25,elecs.motor,:,:)),1);
            tfr_hll.powspctrm_motor_lh75    = nanmean(squeeze(tfr_hll.powspctrm(idx_75,elecs.motor,:,:)),1);
        end
    end
    
    % average data across channels
    cfg = [];
    cfg.avgoverrpt  = 'yes';
    cfg.avgoverchan = 'yes';
    cfg.nanmean     = 'yes';
        
    all_tfr.hll{subji} = ft_selectdata(cfg, tfr_hll);
    all_tfr.hll{subji}.label = {'ECoG'};

end

% frequency to use (from 2Hz to end)
freq2use  = nearest(all_tfr.hll{1}.freq, 2);
% time2 to use ± 0.5 sec. from HLL
time2use  = nearest(all_tfr.hll{1}.time, -0.5):nearest(all_tfr.hll{1}.time, 0.5);

%% - EXTRACT DATA PER ROI - %%

data_frontal.lh_0  = cell(numel(subjects),1);
data_frontal.lh_25 = cell(numel(subjects),1);
data_frontal.lh_75 = cell(numel(subjects),1);

data_motor.lh_0    = cell(numel(subjects),1);
data_motor.lh_25   = cell(numel(subjects),1);
data_motor.lh_75   = cell(numel(subjects),1);

for subji = 1:numel(subjects)
    
    if isfield(all_tfr.hll{subji}, 'powspctrm_frontal_lh0')
        
        data_frontal.lh_0{subji}.label      = {'Frontal'};
        data_frontal.lh_0{subji}.freq       = all_tfr.hll{subji}.freq(freq2use:end);
        data_frontal.lh_0{subji}.time       = all_tfr.hll{subji}.time(time2use);
        data_frontal.lh_0{subji}.powspctrm  = all_tfr.hll{subji}.powspctrm_frontal_lh0(1,freq2use:end,time2use);
        data_frontal.lh_0{subji}.dimord     = all_tfr.hll{subji}.dimord;
        data_frontal.lh_0{subji}.cfg        = [];
    
        data_frontal.lh_25{subji}.label      = {'Frontal'};
        data_frontal.lh_25{subji}.freq       = all_tfr.hll{subji}.freq(freq2use:end);
        data_frontal.lh_25{subji}.time       = all_tfr.hll{subji}.time(time2use);
        data_frontal.lh_25{subji}.powspctrm  = all_tfr.hll{subji}.powspctrm_frontal_lh25(1,freq2use:end,time2use);
        data_frontal.lh_25{subji}.dimord     = all_tfr.hll{subji}.dimord;
        data_frontal.lh_25{subji}.cfg        = [];
    
        data_frontal.lh_75{subji}.label      = {'Frontal'};
        data_frontal.lh_75{subji}.freq       = all_tfr.hll{subji}.freq(freq2use:end);
        data_frontal.lh_75{subji}.time       = all_tfr.hll{subji}.time(time2use);
        data_frontal.lh_75{subji}.powspctrm  = all_tfr.hll{subji}.powspctrm_frontal_lh75(1,freq2use:end,time2use);
        data_frontal.lh_75{subji}.dimord     = all_tfr.hll{subji}.dimord;
        data_frontal.lh_75{subji}.cfg        = [];
        
    end
    
    if isfield(all_tfr.hll{subji}, 'powspctrm_motor_lh0')
        
        data_motor.lh_0{subji}.label      = {'motor'};
        data_motor.lh_0{subji}.freq       = all_tfr.hll{subji}.freq(freq2use:end);
        data_motor.lh_0{subji}.time       = all_tfr.hll{subji}.time(time2use);
        data_motor.lh_0{subji}.powspctrm  = all_tfr.hll{subji}.powspctrm_motor_lh0(1,freq2use:end,time2use);
        data_motor.lh_0{subji}.dimord     = all_tfr.hll{subji}.dimord;
        data_motor.lh_0{subji}.cfg        = [];
    
        data_motor.lh_25{subji}.label      = {'motor'};
        data_motor.lh_25{subji}.freq       = all_tfr.hll{subji}.freq(freq2use:end);
        data_motor.lh_25{subji}.time       = all_tfr.hll{subji}.time(time2use);
        data_motor.lh_25{subji}.powspctrm  = all_tfr.hll{subji}.powspctrm_motor_lh25(1,freq2use:end,time2use);
        data_motor.lh_25{subji}.dimord     = all_tfr.hll{subji}.dimord;
        data_motor.lh_25{subji}.cfg        = [];
    
        data_motor.lh_75{subji}.label      = {'motor'};
        data_motor.lh_75{subji}.freq       = all_tfr.hll{subji}.freq(freq2use:end);
        data_motor.lh_75{subji}.time       = all_tfr.hll{subji}.time(time2use);
        data_motor.lh_75{subji}.powspctrm  = all_tfr.hll{subji}.powspctrm_motor_lh75(1,freq2use:end,time2use);
        data_motor.lh_75{subji}.dimord     = all_tfr.hll{subji}.dimord;
        data_motor.lh_75{subji}.cfg        = [];
        
    end
    
end

% remove empty cells
data_frontal.lh_0  = data_frontal.lh_0(~cellfun(@isempty, data_frontal.lh_0));
data_frontal.lh_25 = data_frontal.lh_25(~cellfun(@isempty, data_frontal.lh_25));
data_frontal.lh_75 = data_frontal.lh_75(~cellfun(@isempty, data_frontal.lh_75));

data_motor.lh_0  = data_motor.lh_0(~cellfun(@isempty, data_motor.lh_0));
data_motor.lh_25 = data_motor.lh_25(~cellfun(@isempty, data_motor.lh_25));
data_motor.lh_75 = data_motor.lh_75(~cellfun(@isempty, data_motor.lh_75));

%% - GRAND AVERAGE - %%

cfg = [];
cfg.parameter         = 'powspctrm';

GA_hll.frontal.lh_0   = ft_freqgrandaverage(cfg, data_frontal.lh_0{:});
GA_hll.frontal.lh_25  = ft_freqgrandaverage(cfg, data_frontal.lh_25{:});
GA_hll.frontal.lh_75  = ft_freqgrandaverage(cfg, data_frontal.lh_75{:});

GA_hll.motor.lh_0     = ft_freqgrandaverage(cfg, data_motor.lh_0{:});
GA_hll.motor.lh_25    = ft_freqgrandaverage(cfg, data_motor.lh_25{:});
GA_hll.motor.lh_75    = ft_freqgrandaverage(cfg, data_motor.lh_75{:});

%% - CLUSTER F-STATS - %%

% Get data of interest
q       = 'Statistic: Frontal = 1; Motor = 2: ';
selroi  = input(q);

if selroi == 1
    
    statsdata  = data_frontal;
    
elseif selroi == 2
    
    statsdata  = data_motor;
        
end
    
% Run stats

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = [-0.5 0.5];
cfg.tail                = 1;
cfg.parameter           = 'powspctrm';
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

stat = ft_freqstatistics(cfg, statsdata.lh_0{:}, statsdata.lh_25{:}, statsdata.lh_75{:});    
    
% Plot F-Stats

figure;

load /Users/janweber/Documents/MATLAB/general_functions/colormaps/lajolla/lajolla.mat

cmap_fstats = lajolla;

% cmap_fstats = cbrewer('seq', 'YlOrRd', 1000,  'PCHIP');

% plot data
contourf(stat.time, stat.freq, squeeze(stat.stat), 150, 'LineColor', 'none'); colormap(cmap_fstats);
hold on
contour(stat.time, stat.freq, squeeze(stat.mask),1,'linecolor','k', 'LineWidth', 2);
maxstat = max(stat.stat, [], 'all');
set(gca, 'Clim', [0 7], 'linewidth', 0.5, 'fontsize', 13);
set(gca,'ytick',[4 8 16 32 64 128],'yscale','log');
xlabel('Time to HLL [s]');
ylabel('Frequency [Hz]');
set(gca,'xticklabel',{[]}, 'yticklabel', {[]});
box off
set(gcf, 'Position', [533 328 247 185]);

cd(path_out);

if selroi == 1
    print(gcf,'fig_3_tfr_fstat_frontal.png','-dpng','-r500');
else
    print(gcf,'fig_3_tfr_fstat_motor.png','-dpng','-r500');
end

% plot the colormap and save
figure;
set(gca, 'Clim', [0 7], 'linewidth', 2, 'fontsize', 13);
y = colorbar;
colormap(cmap_fstats);

print(gcf,'fig_3_tfr_fstat_colorbar.pdf','-dpdf','-r500');

%% Perform Follow-Up Comparisons

cfg                     = [];
cfg.channel             = 'all';
cfg.latency             = [-0.5 0.5];
cfg.tail                = 0;
cfg.parameter           = 'powspctrm';
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

stat_0vs25  = ft_freqstatistics(cfg, statsdata.lh_25{:}, statsdata.lh_0{:});

stat_0vs75  = ft_freqstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_0{:});

stat_25vs75 = ft_freqstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_25{:});


%% - PLOT CHANGES PER FREQUENCY BAND - %%

% color for conditons
color = {'#008450', '#EFB700', '#B81D13'};

% get grand average and keep individuals
cfg = [];
cfg.keepindividual    = 'yes';
cfg.parameter         = 'powspctrm';

GA_hll.lh_0   = ft_freqgrandaverage(cfg, statsdata.lh_0{:});
GA_hll.lh_25  = ft_freqgrandaverage(cfg, statsdata.lh_25{:});
GA_hll.lh_75  = ft_freqgrandaverage(cfg, statsdata.lh_75{:});


% extract significant frequencies from f-statistic
[idx,~]  = find(squeeze(stat.mask));
sigfreqs = stat.freq(unique(idx));

figure;

a = plot(GA_hll.lh_0.freq, smoothdata(squeeze(mean(mean(GA_hll.lh_0.powspctrm(:,1,:,:),1),4))', 'movmean', 5), ...
    'Color', color{1}, 'LineWidth', 4);
hold on
a = plot(GA_hll.lh_25.freq, smoothdata(squeeze(mean(mean(GA_hll.lh_25.powspctrm(:,1,:,:),1),4))', 'movmean', 5), ...
    'Color', color{2}, 'LineWidth', 4);
hold on
a = plot(GA_hll.lh_75.freq, smoothdata(squeeze(mean(mean(GA_hll.lh_75.powspctrm(:,1,:,:),1),4))', 'movmean', 5), ...
    'Color', color{3}, 'LineWidth', 4);

ylabel('Power [z]');
xlabel('Frequency [Hz]');
set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [GA_hll.lh_0.freq(1) GA_hll.lh_0.freq(end)]); 
box off; 

if selroi == 1

    % find transition
    x = find(diff(diff(sigfreqs)) > 10);

    % mask desynch
    scatter(sigfreqs(1:x+1), ones(1, x+1)*-3, 70, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2);

    hold on

    % mask synch
    scatter(sigfreqs(x+2:end), ones(1, length(sigfreqs(x+2:end)))*2, 70, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2);

    set(gca, 'xscale', 'log', 'xtick', [4 8 16 32 64 128], 'ylim', [-3.2 2.2]);

else
    
    % mask synch
    scatter(sigfreqs, ones(1, length(sigfreqs))*4, 70, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
        'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.2);

    set(gca, 'xscale', 'log', 'xtick', [4 8 16 32 64 128], 'ylim', [-3 4.2]);
    
end


set(gcf, 'Position', [479 351 199 184]);

if selroi == 1
    print(gcf,'fig_3_powchange_frontal.pdf','-dpdf','-r500');
else
    print(gcf,'fig_3_powchange_motor.pdf','-dpdf','-r500');
end

%% - SAVE SOURCE DATA - %%

% go into results folder
mkdir(path_source, 'Fig_3c_d_right');
cd(fullfile(path_source, 'Fig_3c_d_right'));

freq = GA_hll.lh_0.freq;

tmp = [freq; smoothdata(squeeze(mean(mean(GA_hll.lh_0.powspctrm(:,1,:,:),1),4))', 'movmean', 5); ...
    smoothdata(squeeze(mean(mean(GA_hll.lh_25.powspctrm(:,1,:,:),1),4))', 'movmean', 5); ...
    smoothdata(squeeze(mean(mean(GA_hll.lh_75.powspctrm(:,1,:,:),1),4))', 'movmean', 5)];

tfr_table = array2table(tmp); clear tmp

% write data into excelsheet
if selroi == 1
    filename = 'tfr_pfc.xlsx';
    writetable(tfr_table,filename, 'Sheet', 1);
elseif selroi == 2
    filename = 'tfr_motor.xlsx';
    writetable(tfr_table,filename, 'Sheet', 1);
end



