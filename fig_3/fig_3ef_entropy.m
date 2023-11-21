%% Figure 3e & f (right panel)
% this script rungs additional analysis required for the plots of figure
% 3e & 3f.

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

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/sample_entropy/

%% - LOAD ELECTRODE INFO - %%

load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/HFA/ExpVariance_Approach/TaskInformative_Electrodes/ROIs_active_elecs.mat

%% - LOAD RAW DATA - %%

fsample    = 512; % sampling freq

all_Data   = struct(); % structure where data is stored

for Isub = 1:numel(subjects)
    
    subID = subjects{Isub};
        
    fprintf('%s\n', subID);
    
    % load ROI Data
    cd([path_in filesep subID filesep 'Data']);
    
    load([subID '_clean_segmented.mat']);
    
    % get indices for conditions
    idx_lh0  = find(data_clean.hll.trialinfo(:,13) == 0);
    idx_lh25 = find(data_clean.hll.trialinfo(:,13) == 25);
    idx_lh75 = find(data_clean.hll.trialinfo(:,13) == 75);
        
    alltrialidx = {idx_lh0, idx_lh25, idx_lh75};

    % get only task informative electrodes
    
    elec    = struct();
    roidata = struct();
    
    if isfield(ROIs(Isub).hll, 'frontal') && isfield(ROIs(Isub).hll.frontal, 'active')
        elec.frontal            = ROIs(Isub).hll.frontal.active.elecidx;
        
        roidata.frontal.trial   = data_clean.hll.trial(:,elec.frontal,:);
        roidata.frontal.time    = data_clean.hll.time;
        roidata.frontal.dimord  = data_clean.hll.dimord;
        roidata.frontal.label   = data_clean.hll.label(elec.frontal);
        roidata.frontal.cfg     = [];
        
    end

    if isfield(ROIs(Isub).hll, 'motor') && isfield(ROIs(Isub).hll.motor, 'active')
        elec.motor              = ROIs(Isub).hll.motor.active.elecidx;
        roidata.motor.trial     = data_clean.hll.trial(:,elec.motor,:);
        roidata.motor.time      = data_clean.hll.time;
        roidata.motor.dimord    = data_clean.hll.dimord;
        roidata.motor.label     = data_clean.hll.label(elec.motor);
        roidata.motor.cfg       = [];
        
    end

    % store data
    if isfield(roidata, 'frontal') || isfield(roidata, 'motor')
        
        all_Data(Isub).trialidx     = alltrialidx;
        all_Data(Isub).roidata      = roidata;
        
    end

end

%% - COMPUTE ENTROPY - %%

winsize        = round(fsample*0.1); % window size (100ms window)

stepsize       = round(fsample*0.02); % step size (20ms stepsize)

timewin        = all_Data(5).roidata.frontal.time;

smoothwin      = 0.005; % smoothing window

ncondition     = 3; % number of conditions

neural_entropy = struct(); % structure where entropy over time is stored
    
for Isub = 1:numel(subjects)
    
    fprintf('Compute entropy %s\n', subjects{Isub});
    
    if ~isempty(all_Data(Isub).roidata)

        % get data, trialmatrix and index for active electrodes
        roidata         = all_Data(Isub).roidata;
        alltrialidx     = all_Data(Isub).trialidx;

        for Icond = 1:ncondition

            % ---------------------------
            % Frontal Cortex
            % ---------------------------

            if isfield(roidata, 'frontal')

                % current data
                tmpData_task = roidata.frontal.trial(alltrialidx{Icond},:,:);

                % overlapping samples (Nx2) to use 
                samples2use = [1:stepsize:length(timewin)-winsize; ...
                    1+winsize:stepsize:length(timewin)]';            

                % initialize entropy matrix
                entropy_task  = NaN(size(tmpData_task,1), size(tmpData_task,2), size(samples2use,1));

                for Ichan = 1:size(entropy_task,2) %loop over channel

                    for Itrial = 1:size(entropy_task,1) % over trials

                        for Isample = 1:size(samples2use,1) % over samples 

                            tmpsamples = samples2use(Isample,:); % current samples to use

                            % compute entropy 
                            entropy_task(Itrial,Ichan,Isample) = sampen(squeeze(tmpData_task(Itrial,Ichan,tmpsamples(1):tmpsamples(2))), 2, 0.2, 'chebychev');

                        end

                        % compute entropy 
                        entropy_task(Itrial,Ichan,:) = smoothdata(squeeze(entropy_task(Itrial,Ichan,:))', 'movmean', [round(fsample*smoothwin) round(fsample*smoothwin)]);

                    end

                end

                % store things in structure
                if Icond == 1
                    
                    neural_entropy.frontal.lh_0{Isub}.entropy        = mean(squeeze(nanmean(entropy_task,2))); 
                    neural_entropy.frontal.lh_0{Isub}.time           = roidata.frontal.time(1:stepsize:end-winsize);              
                    neural_entropy.frontal.lh_0{Isub}.label          = {'Frontal'};
                    neural_entropy.frontal.lh_0{Isub}.dimord         = 'chan_time';
                    neural_entropy.frontal.lh_0{Isub}.cfg            = [];

                elseif Icond == 2
                    
                    neural_entropy.frontal.lh_25{Isub}.entropy        = mean(squeeze(nanmean(entropy_task,2))); 
                    neural_entropy.frontal.lh_25{Isub}.time           = roidata.frontal.time(1:stepsize:end-winsize);              
                    neural_entropy.frontal.lh_25{Isub}.label          = {'Frontal'};
                    neural_entropy.frontal.lh_25{Isub}.dimord         = 'chan_time';
                    neural_entropy.frontal.lh_25{Isub}.cfg            = [];            

                elseif Icond == 3
                    
                    neural_entropy.frontal.lh_75{Isub}.entropy        = mean(squeeze(nanmean(entropy_task,2))); 
                    neural_entropy.frontal.lh_75{Isub}.time           = roidata.frontal.time(1:stepsize:end-winsize);              
                    neural_entropy.frontal.lh_75{Isub}.label          = {'Frontal'};
                    neural_entropy.frontal.lh_75{Isub}.dimord         = 'chan_time';
                    neural_entropy.frontal.lh_75{Isub}.cfg            = [];            
                
                end

                clear tmpData_task

            end

            % ---------------------------
            % Motor Cortex
            % ---------------------------

            if isfield(roidata, 'motor')

                % current data
                tmpData_task = roidata.motor.trial(alltrialidx{Icond},:,:);

                % overlapping samples (Nx2) to use 
                samples2use = [1:stepsize:length(timewin)-winsize; ...
                    1+winsize:stepsize:length(timewin)]';            

                % initialize entropy matrix
                entropy_task  = NaN(size(tmpData_task,1), size(tmpData_task,2), size(samples2use,1));

                for Ichan = 1:size(entropy_task,2) %loop over channel

                    for Itrial = 1:size(entropy_task,1) % over trials

                        for Isample = 1:size(samples2use,1) % over samples 

                            tmpsamples = samples2use(Isample,:); % current samples to use

                            % compute entropy 
                            entropy_task(Itrial,Ichan,Isample) = sampen(squeeze(tmpData_task(Itrial,Ichan,tmpsamples(1):tmpsamples(2))), 2, 0.2, 'chebychev');

                        end

                        % compute entropy 
                        entropy_task(Itrial,Ichan,:) = smoothdata(squeeze(entropy_task(Itrial,Ichan,:))', 'movmean', [round(fsample*smoothwin) round(fsample*smoothwin)]);

                    end

                end

                % store things in structure
                if Icond == 1
                    
                    neural_entropy.motor.lh_0{Isub}.entropy        = mean(squeeze(nanmean(entropy_task,2))); 
                    neural_entropy.motor.lh_0{Isub}.time           = roidata.motor.time(1:stepsize:end-winsize);              
                    neural_entropy.motor.lh_0{Isub}.label          = {'motor'};
                    neural_entropy.motor.lh_0{Isub}.dimord         = 'chan_time';
                    neural_entropy.motor.lh_0{Isub}.cfg            = [];

                elseif Icond == 2
                    
                    neural_entropy.motor.lh_25{Isub}.entropy        = mean(squeeze(nanmean(entropy_task,2))); 
                    neural_entropy.motor.lh_25{Isub}.time           = roidata.motor.time(1:stepsize:end-winsize);              
                    neural_entropy.motor.lh_25{Isub}.label          = {'motor'};
                    neural_entropy.motor.lh_25{Isub}.dimord         = 'chan_time';
                    neural_entropy.motor.lh_25{Isub}.cfg            = [];            

                elseif Icond == 3
                    
                    neural_entropy.motor.lh_75{Isub}.entropy        = mean(squeeze(nanmean(entropy_task,2))); 
                    neural_entropy.motor.lh_75{Isub}.time           = roidata.motor.time(1:stepsize:end-winsize);              
                    neural_entropy.motor.lh_75{Isub}.label          = {'motor'};
                    neural_entropy.motor.lh_75{Isub}.dimord         = 'chan_time';
                    neural_entropy.motor.lh_75{Isub}.cfg            = [];            
                
                end

                clear tmpData_task

            end
            
        end %Icond
        
        % -------------------------------------------------
        % Frontal & Motor Cortex for Interaction Testing
        % -------------------------------------------------

        if isfield(roidata, 'frontal') && isfield(roidata, 'motor')

            % frontal cortex
            neural_entropy.frontal.diff{Isub}.entropy        = neural_entropy.frontal.lh_75{Isub}.entropy - ...
                neural_entropy.frontal.lh_0{Isub}.entropy;
            neural_entropy.frontal.diff{Isub}.time           = roidata.frontal.time(1:stepsize:end-winsize);              
            neural_entropy.frontal.diff{Isub}.label          = {'Dummy'};
            neural_entropy.frontal.diff{Isub}.dimord         = 'chan_time';
            neural_entropy.frontal.diff{Isub}.cfg            = [];

            % motor cortex
            neural_entropy.motor.diff{Isub}.entropy        = neural_entropy.motor.lh_75{Isub}.entropy - ...
                neural_entropy.motor.lh_0{Isub}.entropy;
            neural_entropy.motor.diff{Isub}.time           = roidata.motor.time(1:stepsize:end-winsize);              
            neural_entropy.motor.diff{Isub}.label          = {'Dummy'};
            neural_entropy.motor.diff{Isub}.dimord         = 'chan_time';
            neural_entropy.motor.diff{Isub}.cfg            = [];

            clear tmpData_task

        end            

        clear roidata alltrialidx active_elecs
        
    end

end %Isub
    
%% - REMOVE EMPTY CELLS - %%

% frontal
neural_entropy.frontal.lh_0  = neural_entropy.frontal.lh_0(~cellfun(@isempty, neural_entropy.frontal.lh_0));
neural_entropy.frontal.lh_25 = neural_entropy.frontal.lh_25(~cellfun(@isempty, neural_entropy.frontal.lh_25));
neural_entropy.frontal.lh_75 = neural_entropy.frontal.lh_75(~cellfun(@isempty, neural_entropy.frontal.lh_75));
neural_entropy.frontal.diff  = neural_entropy.frontal.diff(~cellfun(@isempty, neural_entropy.frontal.diff));
 
% motor
neural_entropy.motor.lh_0    = neural_entropy.motor.lh_0(~cellfun(@isempty, neural_entropy.motor.lh_0));
neural_entropy.motor.lh_25   = neural_entropy.motor.lh_25(~cellfun(@isempty, neural_entropy.motor.lh_25));
neural_entropy.motor.lh_75   = neural_entropy.motor.lh_75(~cellfun(@isempty, neural_entropy.motor.lh_75));
neural_entropy.motor.diff    = neural_entropy.motor.diff(~cellfun(@isempty, neural_entropy.motor.diff));

%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% - FRONTAL CORTEX - %

sourcedata_frontal = struct;

for Isub = 1:numel(neural_entropy.frontal.lh_0)
    
    if isempty(neural_entropy.frontal.lh_0{Isub})
        
        sourcedata_frontal.frontal.lh_0(Isub,:)  = NaN(1, size(samples2use,1));
        sourcedata_frontal.frontal.lh_25(Isub,:) = NaN(1, size(samples2use,1));
        sourcedata_frontal.frontal.lh_75(Isub,:) = NaN(1, size(samples2use,1));

    else
        
        sourcedata_frontal.frontal.lh_0(Isub,:) = neural_entropy.frontal.lh_0{Isub}.entropy;
        sourcedata_frontal.frontal.lh_25(Isub,:) = neural_entropy.frontal.lh_25{Isub}.entropy;
        sourcedata_frontal.frontal.lh_75(Isub,:) = neural_entropy.frontal.lh_75{Isub}.entropy;

    end        
        
end%Isub

% go into results folder
mkdir(path_source, 'Fig_3e_f_right');
cd(fullfile(path_source, 'Fig_3e_f_right'));

% get time window of interest
time = timewin(samples2use(:,1));

% -- FRONTAL CORTEX 0% LH -- %%
tmp = NaN(size(sourcedata_frontal.frontal.lh_0,1)+1, size(sourcedata_frontal.frontal.lh_0,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata_frontal.frontal.lh_0,1))';
tmp(2:end,2:end)  = sourcedata_frontal.frontal.lh_0;

table_frontal.lh0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_3f_sampEn_frontal_lh0.xlsx';
writetable(table_frontal.lh0,filename, 'Sheet', 1);

% -- FRONTAL CORTEX 25% LH -- %%
tmp = NaN(size(sourcedata_frontal.frontal.lh_25,1)+1, size(sourcedata_frontal.frontal.lh_25,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata_frontal.frontal.lh_25,1))';
tmp(2:end,2:end)  = sourcedata_frontal.frontal.lh_25;

table_frontal.lh25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_3f_sampEn_frontal_lh25.xlsx';
writetable(table_frontal.lh25,filename, 'Sheet', 1);

% -- FRONTAL CORTEX 75% LH -- %%
tmp = NaN(size(sourcedata_frontal.frontal.lh_75,1)+1, size(sourcedata_frontal.frontal.lh_75,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata_frontal.frontal.lh_75,1))';
tmp(2:end,2:end)  = sourcedata_frontal.frontal.lh_75;

table_frontal.lh75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_3f_sampEn_frontal_lh75.xlsx';
writetable(table_frontal.lh75,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% - MOTOR CORTEX - %

sourcedata_motor = struct;

for Isub = 1:numel(neural_entropy.motor.lh_0)
    
    if isempty(neural_entropy.motor.lh_0{Isub})
        
        sourcedata_motor.motor.lh_0(Isub,:)  = NaN(1, size(samples2use,1));
        sourcedata_motor.motor.lh_25(Isub,:) = NaN(1, size(samples2use,1));
        sourcedata_motor.motor.lh_75(Isub,:) = NaN(1, size(samples2use,1));

    else
        
        sourcedata_motor.motor.lh_0(Isub,:) = neural_entropy.motor.lh_0{Isub}.entropy;
        sourcedata_motor.motor.lh_25(Isub,:) = neural_entropy.motor.lh_25{Isub}.entropy;
        sourcedata_motor.motor.lh_75(Isub,:) = neural_entropy.motor.lh_75{Isub}.entropy;

    end        
        
end%Isub

% get time window of interest
time = timewin(samples2use(:,1));

% -- MOTOR CORTEX 0% LH -- %%
tmp = NaN(size(sourcedata_motor.motor.lh_0,1)+1, size(sourcedata_motor.motor.lh_0,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata_motor.motor.lh_0,1))';
tmp(2:end,2:end)  = sourcedata_motor.motor.lh_0;

table_motor.lh0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_3f_sampEn_motor_lh0.xlsx';
writetable(table_motor.lh0,filename, 'Sheet', 1);

% -- MOTOR CORTEX 25% LH -- %%
tmp = NaN(size(sourcedata_motor.motor.lh_25,1)+1, size(sourcedata_motor.motor.lh_25,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata_motor.motor.lh_25,1))';
tmp(2:end,2:end)  = sourcedata_motor.motor.lh_25;

table_motor.lh25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_3f_sampEn_motor_lh25.xlsx';
writetable(table_motor.lh25,filename, 'Sheet', 1);

% -- MOTOR CORTEX 75% LH -- %%
tmp = NaN(size(sourcedata_motor.motor.lh_75,1)+1, size(sourcedata_motor.motor.lh_75,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata_motor.motor.lh_75,1))';
tmp(2:end,2:end)  = sourcedata_motor.motor.lh_75;

table_motor.lh75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_3f_sampEn_motor_lh75.xlsx';
writetable(table_motor.lh75,filename, 'Sheet', 1);

%% - MERGE DATA - %%

data_frontal = {neural_entropy.frontal.lh_0, neural_entropy.frontal.lh_25, neural_entropy.frontal.lh_75};
data_motor   = {neural_entropy.motor.lh_0, neural_entropy.motor.lh_25, neural_entropy.motor.lh_75};

%% - CLUSTER F-STATS - %%

% Get data of interest
q       = 'Statistic: Frontal = 1; Motor = 2: ';
selroi = input(q);

statsdata = struct();

if selroi == 1
    
    statsdata.lh_0   = data_frontal{1};
    statsdata.lh_25  = data_frontal{2};
    statsdata.lh_75  = data_frontal{3};

elseif selroi == 2
    
    statsdata.lh_0   = data_motor{1};
    statsdata.lh_25  = data_motor{2};
    statsdata.lh_75  = data_motor{3};
        
end
    
% Run stats

cfg                     = [];
cfg.channel             = 'all';
cfg.tail                = 1;
cfg.parameter           = 'entropy';
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

stat = ft_timelockstatistics(cfg, statsdata.lh_0{:}, statsdata.lh_25{:}, statsdata.lh_75{:});    

% extract start point of cluster
if isfield(stat, 'posclusters')
    
    for Icluster = 1:numel(stat.posclusters)
        
        fprintf('\n');
        
        fprintf('p-value cluster %d = %.3f\n', Icluster, stat.posclusters(Icluster).prob);
        
        begeffect = stat.time(find(stat.posclusterslabelmat == Icluster, 1, 'first'));
        endeffect = stat.time(find(stat.posclusterslabelmat == Icluster, 1, 'last'));
        
        fprintf('Effect cluster %d from %.3f to %.3f ms\n', Icluster, begeffect, endeffect);
        
        fprintf('\n');

    end
    
end
    
%% - PLOT F-STATS - %%

% plot data
smoothwin_plot = 0.002;
plot(stat.time, smoothdata(stat.stat, 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'k', 'LineWidth', 2);

hold on

% shade grey where significant
sigpnts = stat.time(find(stat.mask));
x_lim   = xlim;
y_lim   = [0 6.5]; % based on stats in frontal cortex

if ~isempty(sigpnts)
    pgon = polyshape([sigpnts(1) sigpnts(1) (sigpnts(end) - sigpnts(1)) (sigpnts(end) - sigpnts(1))], ...
        [y_lim(end) y_lim(1) y_lim(1) y_lim(end)]);
    h = plot(pgon);
    h.FaceAlpha = 0.3;
    h.FaceColor = [0.5 0.5 0.5];
    h.EdgeColor = [0.5 0.5 0.5];
    h.EdgeAlpha = 0.3;
end

xlim([x_lim(1) x_lim(end)])
ylim([y_lim(1) y_lim(end)])
ylabel('Stats [F]', 'FontSize', 13);
xlabel('Time to HLL [s]', 'FontSize', 13);
set(gca, 'fontsize', 13, 'linewidth', 3);
box off
set(gcf, 'Position', [552 396 205 97]);

cd(path_out);

% save figure
if selroi == 1
    print(gcf,'fig_3_frontal_fstat_sampEn.pdf','-dpdf','-r400');
else
    print(gcf,'fig_3_motor_fstat_sampEn.pdf','-dpdf','-r400');
end
    
%% - FOLLOW-UP CLUSTER-BASED STATS - %%

cfg                     = [];
cfg.channel             = 'all';
cfg.tail                = 0;
cfg.parameter           = 'entropy';
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
cfg.parameter  = 'entropy';
cfg.method     = 'analytic';
cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar       = 1;
cfg.design     = design;
cfg.uvar       = 1;
cfg.ivar       = 2;

d_0vs25        = ft_timelockstatistics(cfg, statsdata.lh_25{:}, statsdata.lh_0{:});
d_0vs75        = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_0{:});
d_25vs75       = ft_timelockstatistics(cfg, statsdata.lh_75{:}, statsdata.lh_25{:});

% get current grand averages

cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter      = 'entropy';

tmpGA_lh0          = ft_timelockgrandaverage(cfg, statsdata.lh_0{:});
tmpGA_lh25         = ft_timelockgrandaverage(cfg, statsdata.lh_25{:});
tmpGA_lh75         = ft_timelockgrandaverage(cfg, statsdata.lh_75{:});

% plot data
smoothwin_plot = 0.002;
plot(tmpGA_lh0.time, smoothdata(squeeze(nanmean(tmpGA_lh0.individual)), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#008450', 'LineWidth', 4);
hold on
plot(tmpGA_lh25.time, smoothdata(squeeze(nanmean(tmpGA_lh25.individual)), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#EFB700', 'LineWidth', 4);
hold on
plot(tmpGA_lh75.time, smoothdata(squeeze(nanmean(tmpGA_lh75.individual)), 'movmean', [round(fsample*smoothwin_plot) round(fsample*smoothwin_plot)]), 'Color', '#B81D13', 'LineWidth', 4);

% xlabel('Time to HLL [s]');
ylabel('Sample entropy');

% plot significant timepoints
a = ylim;

hold on

island = bwconncomp(stat.mask);
for jj = 1:island.NumObjects
    h = line(tmpGA_lh0.time(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.007, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    b = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 4);
    b.Color = 'k';
    hold on
    b = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 4);
    b.Color = 'k';
end
clear island


island = bwconncomp(stat_0vs75.mask);
for jj = 1:island.NumObjects
    h = line(tmpGA_lh0.time(island.PixelIdxList{jj}), ...
        repmat(a(2) + 0.005, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#008450', 'LineWidth', 4);
    hold on
    plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#B81D13', 'LineWidth', 4);
end
clear island

island = bwconncomp(stat_0vs25.mask);
for jj = 1:island.NumObjects
    h = line(tmpGA_lh25.time(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.015, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#008450', 'LineWidth', 4);
    hold on
    plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#EFB700', 'LineWidth', 4);
end
clear island

island = bwconncomp(stat_25vs75.mask);
for jj = 1:island.NumObjects
    h = line(tmpGA_lh75.time(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.025, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'Color', '#EFB700', 'LineWidth', 4);
    hold on
    plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'Color', '#B81D13', 'LineWidth', 4);
end
clear island

a = ylim;

set(gca, 'linewidth', 3, 'FontSize',13, 'xlim', [tmpGA_lh0.time(1) tmpGA_lh0.time(end)], ...
    'Ylim', [a(1)-0.001 a(2)+0.001]); box off; hold on;

set(gcf, 'Position', [503 336 240 162]);

allposthoc = {stat_0vs75, stat_0vs25, stat_25vs75};
allcohens  = {d_0vs75, d_0vs25, d_25vs75};
posthocidx = {'0 vs 75', '0 vs 25', '25 vs 75'};

for Ipair = 1:numel(allposthoc)

    % extract start point of cluster
    if isfield(allposthoc{Ipair}, 'posclusters')

        for Icluster = 1:numel(allposthoc{Ipair}.posclusters)

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
    
%% - SAVE FIGURE - %%

cd(path_out);
if selroi == 1
    print(gcf,'fig_3_frontal_sampEn.pdf','-dpdf','-r400');
else
    print(gcf,'fig_3_motor_sampEn.pdf','-dpdf','-r400');
end





