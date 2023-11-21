%% Figure 4c
% this script rungs additional analysis required for the plots of figure
% 4c

%%

clear
ft_defaults

%% - SUBJECTS - %%
 
% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_4/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

%% - LOAD DATA - %%

load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/HFA/ExpVariance_Approach/Peak_Triggered_Average/ByCondition/GA_FreqDomain.mat

% condition names
conditions = fieldnames(allData_freq);

% number of levels for factor prediction
numlevel   = fieldnames(allData_freq.probe);

disp('Loading Data...')

% load any data to extract frequency range 
load /Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/OSL13/Data/HFA/Peak_Triggered_Avg/OSL13_peakavg_freqdomain_pca_active_bycondition.mat
freq  = peakdatafreq.bp{1}.frontal.inactive.freq; 
label = peakdatafreq.probe{1}.frontal.inactive.label(1);

%% - CREATE FT DUMMY - %%

data_frontal = struct();

data_motor   = struct();

%% - EXTRACT DATA FRONTAL CORTEX - %%

count = zeros(2,1);

disp('Extract Data Frontal Cortex');

for Isub = 1:numel(subjects)
        
    %% - HLL - %%
    
    if ~all(isnan(allData_freq.hll.lh_0.frontal.active{Isub}), 'all')
        
        count(1,1) = count(1,1) + 1;
        
        % make sure the dimension is correct (channel x freq)
        dim = size(allData_freq.hll.lh_0.frontal.active{Isub},1);
        
        if dim(1) == length(freq) % if first dimension is freq, permute to channel x freq
            
            allData_freq.hll.lh_0.frontal.active{Isub}  = permute(allData_freq.hll.lh_0.frontal.active{Isub}, [2 1]);
            allData_freq.hll.lh_25.frontal.active{Isub} = permute(allData_freq.hll.lh_25.frontal.active{Isub}, [2 1]);
            allData_freq.hll.lh_75.frontal.active{Isub} = permute(allData_freq.hll.lh_75.frontal.active{Isub}, [2 1]);
            
        end
        
        % average across channel within a ROI
        if size(allData_freq.hll.lh_0.frontal.active{Isub},1) == 1 % check if data has more than one channel (if not, do not take the average)
            data_frontal.hll_lh0{count(1)}.powspctrm = allData_freq.hll.lh_0.frontal.active{Isub};
        else
            data_frontal.hll_lh0{count(1)}.powspctrm  = nanmean(allData_freq.hll.lh_0.frontal.active{Isub});
        end
        
        if size(allData_freq.hll.lh_25.frontal.active{Isub},1) == 1
            data_frontal.hll_lh25{count(1)}.powspctrm = allData_freq.hll.lh_25.frontal.active{Isub};
        else
            data_frontal.hll_lh25{count(1)}.powspctrm  = nanmean(allData_freq.hll.lh_25.frontal.active{Isub});
        end
    
        if size(allData_freq.hll.lh_75.frontal.active{Isub},1) == 1
            data_frontal.hll_lh75{count(1)}.powspctrm = allData_freq.hll.lh_75.frontal.active{Isub};
        else
            data_frontal.hll_lh75{count(1)}.powspctrm  = nanmean(allData_freq.hll.lh_75.frontal.active{Isub});
        end
        
        % fake FT structure
        data_frontal.hll_lh0{count(1)}.dimord  = 'chan_freq';
        data_frontal.hll_lh0{count(1)}.freq    = freq;
        data_frontal.hll_lh0{count(1)}.label   = label;
        data_frontal.hll_lh0{count(1)}.cfg     = [];

        data_frontal.hll_lh25{count(1)}.dimord = 'chan_freq';
        data_frontal.hll_lh25{count(1)}.freq   = freq;
        data_frontal.hll_lh25{count(1)}.label  = label;
        data_frontal.hll_lh25{count(1)}.cfg    = [];

        data_frontal.hll_lh75{count(1)}.dimord = 'chan_freq';
        data_frontal.hll_lh75{count(1)}.freq   = freq;
        data_frontal.hll_lh75{count(1)}.label  = label;
        data_frontal.hll_lh75{count(1)}.cfg    = [];        
        
        % get the difference between 0% & 75%
        data_frontal.hll_diff{count(1)}.powspctrm = data_frontal.hll_lh0{count(1)}.powspctrm - ...
            data_frontal.hll_lh75{count(1)}.powspctrm;
        data_frontal.hll_diff{count(1)}.dimord = 'chan_freq';
        data_frontal.hll_diff{count(1)}.freq   = freq;
        data_frontal.hll_diff{count(1)}.label  = label;
        data_frontal.hll_diff{count(1)}.cfg    = [];        

    end    
    
    %% - BR - %%
    
    if ~all(isnan(allData_freq.br.lh_0.frontal.active{Isub}), 'all')
        
        count(2,1) = count(2,1) + 1;
        
        % make sure the dimension is correct (channel x freq)
        dim = size(allData_freq.br.lh_0.frontal.active{Isub},1);
        
        if dim(1) == length(freq) % if first dimension is freq, permute to channel x freq
            
            allData_freq.br.lh_0.frontal.active{Isub}  = permute(allData_freq.br.lh_0.frontal.active{Isub}, [2 1]);
            allData_freq.br.lh_25.frontal.active{Isub} = permute(allData_freq.br.lh_25.frontal.active{Isub}, [2 1]);
            allData_freq.br.lh_75.frontal.active{Isub} = permute(allData_freq.br.lh_75.frontal.active{Isub}, [2 1]);
            
        end
        
        % average across channel within a ROI
        if size(allData_freq.br.lh_0.frontal.active{Isub},1) == 1 % check if data has more than one channel (if not, do not take the average)
            data_frontal.br_lh0{count(2)}.powspctrm = allData_freq.br.lh_0.frontal.active{Isub};
        else
            data_frontal.br_lh0{count(2)}.powspctrm  = nanmean(allData_freq.br.lh_0.frontal.active{Isub});
        end
        
        if size(allData_freq.br.lh_25.frontal.active{Isub},1) == 1
            data_frontal.br_lh25{count(2)}.powspctrm = allData_freq.br.lh_25.frontal.active{Isub};
        else
            data_frontal.br_lh25{count(2)}.powspctrm  = nanmean(allData_freq.br.lh_25.frontal.active{Isub});
        end
    
        if size(allData_freq.br.lh_75.frontal.active{Isub},1) == 1
            data_frontal.br_lh75{count(2)}.powspctrm = allData_freq.br.lh_75.frontal.active{Isub};
        else
            data_frontal.br_lh75{count(2)}.powspctrm  = nanmean(allData_freq.br.lh_75.frontal.active{Isub});
        end
        
        % fake FT structure
        data_frontal.br_lh0{count(2)}.dimord  = 'chan_freq';
        data_frontal.br_lh0{count(2)}.freq    = freq;
        data_frontal.br_lh0{count(2)}.label   = label;
        data_frontal.br_lh0{count(2)}.cfg     = [];

        data_frontal.br_lh25{count(2)}.dimord = 'chan_freq';
        data_frontal.br_lh25{count(2)}.freq   = freq;
        data_frontal.br_lh25{count(2)}.label  = label;
        data_frontal.br_lh25{count(2)}.cfg    = [];

        data_frontal.br_lh75{count(2)}.dimord = 'chan_freq';
        data_frontal.br_lh75{count(2)}.freq   = freq;
        data_frontal.br_lh75{count(2)}.label  = label;
        data_frontal.br_lh75{count(2)}.cfg    = []; 
        
        % get the difference between 0% & 75%
        data_frontal.br_diff{count(2)}.powspctrm = data_frontal.br_lh0{count(2)}.powspctrm - ...
            data_frontal.br_lh75{count(2)}.powspctrm;
        data_frontal.br_diff{count(2)}.dimord = 'chan_freq';
        data_frontal.br_diff{count(2)}.freq   = freq;
        data_frontal.br_diff{count(2)}.label  = label;
        data_frontal.br_diff{count(2)}.cfg    = [];        
        
    end    
    
end

%% - EXTRACT DATA MOTOR CORTEX - %%

count = zeros(2,1);

disp('Extract Data Motor Cortex');

for Isub = 1:numel(subjects)
        
    %% - HLL - %%

    if ~all(isnan(allData_freq.hll.lh_0.motor.active{Isub}), 'all')
        
        count(1,1) = count(1,1) + 1;
        
        % make sure the dimension is correct (channel x freq)
        dim = size(allData_freq.hll.lh_0.motor.active{Isub},1);
        
        if dim(1) == length(freq) % if first dimension is freq, permute to channel x freq
            
            allData_freq.hll.lh_0.motor.active{Isub}  = permute(allData_freq.hll.lh_0.motor.active{Isub}, [2 1]);
            allData_freq.hll.lh_25.motor.active{Isub} = permute(allData_freq.hll.lh_25.motor.active{Isub}, [2 1]);
            allData_freq.hll.lh_75.motor.active{Isub} = permute(allData_freq.hll.lh_75.motor.active{Isub}, [2 1]);
            
        end
        
        % average across channel within a ROI
        if size(allData_freq.hll.lh_0.motor.active{Isub},1) == 1 % check if data has more than one channel (if not, do not take the average)
            data_motor.hll_lh0{count(1)}.powspctrm = allData_freq.hll.lh_0.motor.active{Isub};
        else
            data_motor.hll_lh0{count(1)}.powspctrm  = nanmean(allData_freq.hll.lh_0.motor.active{Isub});
        end
        
        if size(allData_freq.hll.lh_25.motor.active{Isub},1) == 1
            data_motor.hll_lh25{count(1)}.powspctrm = allData_freq.hll.lh_25.motor.active{Isub};
        else
            data_motor.hll_lh25{count(1)}.powspctrm  = nanmean(allData_freq.hll.lh_25.motor.active{Isub});
        end
    
        if size(allData_freq.hll.lh_75.motor.active{Isub},1) == 1
            data_motor.hll_lh75{count(1)}.powspctrm = allData_freq.hll.lh_75.motor.active{Isub};
        else
            data_motor.hll_lh75{count(1)}.powspctrm  = nanmean(allData_freq.hll.lh_75.motor.active{Isub});
        end
        
        % fake FT structure
        data_motor.hll_lh0{count(1)}.dimord  = 'chan_freq';
        data_motor.hll_lh0{count(1)}.freq    = freq;
        data_motor.hll_lh0{count(1)}.label   = label;        
        data_motor.hll_lh0{count(1)}.cfg     = [];

        data_motor.hll_lh25{count(1)}.dimord = 'chan_freq';
        data_motor.hll_lh25{count(1)}.freq   = freq;
        data_motor.hll_lh25{count(1)}.label  = label;
        data_motor.hll_lh25{count(1)}.cfg    = [];

        data_motor.hll_lh75{count(1)}.dimord = 'chan_freq';
        data_motor.hll_lh75{count(1)}.freq   = freq;
        data_motor.hll_lh75{count(1)}.label  = label;
        data_motor.hll_lh75{count(1)}.cfg    = [];
        
        % get the difference between 0% & 75%
        data_motor.hll_diff{count(1)}.powspctrm = data_motor.hll_lh0{count(1)}.powspctrm - ...
            data_motor.hll_lh75{count(1)}.powspctrm;
        data_motor.hll_diff{count(1)}.dimord = 'chan_freq';
        data_motor.hll_diff{count(1)}.freq   = freq;
        data_motor.hll_diff{count(1)}.label  = label;
        data_motor.hll_diff{count(1)}.cfg    = [];        
        
    end    
    
    %% - BR - %%
    
    if ~all(isnan(allData_freq.br.lh_0.motor.active{Isub}), 'all')
        
        count(2,1) = count(2,1) + 1;
        
        % make sure the dimension is correct (channel x freq)
        dim = size(allData_freq.br.lh_0.motor.active{Isub},1);
        
        if dim(1) == length(freq) % if first dimension is freq, permute to channel x freq
            
            allData_freq.br.lh_0.motor.active{Isub}  = permute(allData_freq.br.lh_0.motor.active{Isub}, [2 1]);
            allData_freq.br.lh_25.motor.active{Isub} = permute(allData_freq.br.lh_25.motor.active{Isub}, [2 1]);
            allData_freq.br.lh_75.motor.active{Isub} = permute(allData_freq.br.lh_75.motor.active{Isub}, [2 1]);
            
        end
        
        % average across channel within a ROI
        if size(allData_freq.br.lh_0.motor.active{Isub},1) == 1 % check if data has more than one channel (if not, do not take the average)
            data_motor.br_lh0{count(2)}.powspctrm = allData_freq.br.lh_0.motor.active{Isub};
        else
            data_motor.br_lh0{count(2)}.powspctrm  = nanmean(allData_freq.br.lh_0.motor.active{Isub});
        end
        
        if size(allData_freq.br.lh_25.motor.active{Isub},1) == 1
            data_motor.br_lh25{count(2)}.powspctrm = allData_freq.br.lh_25.motor.active{Isub};
        else
            data_motor.br_lh25{count(2)}.powspctrm  = nanmean(allData_freq.br.lh_25.motor.active{Isub});
        end
    
        if size(allData_freq.br.lh_75.motor.active{Isub},1) == 1
            data_motor.br_lh75{count(2)}.powspctrm = allData_freq.br.lh_75.motor.active{Isub};
        else
            data_motor.br_lh75{count(2)}.powspctrm  = nanmean(allData_freq.br.lh_75.motor.active{Isub});
        end
        
        % fake FT structure
        data_motor.br_lh0{count(2)}.dimord  = 'chan_freq';
        data_motor.br_lh0{count(2)}.freq    = freq;
        data_motor.br_lh0{count(2)}.label    = label;
        data_motor.br_lh0{count(2)}.cfg     = [];

        data_motor.br_lh25{count(2)}.dimord = 'chan_freq';
        data_motor.br_lh25{count(2)}.freq   = freq;
        data_motor.br_lh25{count(2)}.label  = label;
        data_motor.br_lh25{count(2)}.cfg    = [];

        data_motor.br_lh75{count(2)}.dimord = 'chan_freq';
        data_motor.br_lh75{count(2)}.freq   = freq;
        data_motor.br_lh75{count(2)}.label  = label;
        data_motor.br_lh75{count(2)}.cfg    = [];
            
        % get the difference between 0% & 75%
        data_motor.br_diff{count(2)}.powspctrm = data_motor.br_lh0{count(2)}.powspctrm - ...
            data_motor.br_lh75{count(2)}.powspctrm;
        data_motor.br_diff{count(2)}.dimord = 'chan_freq';
        data_motor.br_diff{count(2)}.freq   = freq;
        data_motor.br_diff{count(2)}.label  = label;
        data_motor.br_diff{count(2)}.cfg    = [];        
        
    end   

end

%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% - FRONTAL CORTEX - %

sourcedata_frontal = struct;

for Isub = 1:numel(data_frontal.hll_lh0)
    
    if isempty(data_frontal.hll_lh0{Isub})
        
        sourcedata_frontal.lh_0(Isub,:)  = NaN(1, length(freq));
        sourcedata_frontal.lh_25(Isub,:) = NaN(1, length(freq));
        sourcedata_frontal.lh_75(Isub,:) = NaN(1, length(freq));

    else
        
        sourcedata_frontal.lh_0(Isub,:)  = data_frontal.hll_lh0{Isub}.powspctrm;
        sourcedata_frontal.lh_25(Isub,:) = data_frontal.hll_lh25{Isub}.powspctrm;
        sourcedata_frontal.lh_75(Isub,:) = data_frontal.hll_lh75{Isub}.powspctrm;

    end        
        
end%Isub

% go into results folder
mkdir(path_source, 'Fig_4c');
cd(fullfile(path_source, 'Fig_4c'));

% -- FRONTAL CORTEX 0% LH -- %%
tmp = NaN(size(sourcedata_frontal.lh_0,1)+1, size(sourcedata_frontal.lh_0,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:16)';
tmp(2:end,2:end)  = sourcedata_frontal.lh_0;

table_frontal.lh0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_4c_pta_frontal_lh0.xlsx';
writetable(table_frontal.lh0,filename, 'Sheet', 1);

% -- FRONTAL CORTEX 25% LH -- %%
tmp = NaN(size(sourcedata_frontal.lh_25,1)+1, size(sourcedata_frontal.lh_25,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:16)';
tmp(2:end,2:end)  = sourcedata_frontal.lh_25;

table_frontal.lh25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_4c_pta_frontal_lh25.xlsx';
writetable(table_frontal.lh25,filename, 'Sheet', 1);

% -- FRONTAL CORTEX 75% LH -- %%
tmp = NaN(size(sourcedata_frontal.lh_75,1)+1, size(sourcedata_frontal.lh_75,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:16)';
tmp(2:end,2:end)  = sourcedata_frontal.lh_75;

table_frontal.lh75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_4c_pta_frontal_lh75.xlsx';
writetable(table_frontal.lh75,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% - MOTOR CORTEX - %

sourcedata_motor = struct;

for Isub = 1:numel(data_motor.hll_lh0)
    
    if isempty(data_motor.hll_lh0{Isub})
        
        sourcedata_motor.lh_0(Isub,:)  = NaN(1, length(freq));
        sourcedata_motor.lh_25(Isub,:) = NaN(1, length(freq));
        sourcedata_motor.lh_75(Isub,:) = NaN(1, length(freq));

    else
        
        sourcedata_motor.lh_0(Isub,:)  = data_motor.hll_lh0{Isub}.powspctrm;
        sourcedata_motor.lh_25(Isub,:) = data_motor.hll_lh25{Isub}.powspctrm;
        sourcedata_motor.lh_75(Isub,:) = data_motor.hll_lh75{Isub}.powspctrm;

    end        
        
end%Isub

% -- MOTOR CORTEX 0% LH -- %%
tmp = NaN(size(sourcedata_motor.lh_0,1)+1, size(sourcedata_motor.lh_0,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:11)';
tmp(2:end,2:end)  = sourcedata_motor.lh_0;

table_motor.lh0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_supp_fig_4a_pta_motor_lh0.xlsx';
writetable(table_motor.lh0,filename, 'Sheet', 1);

% -- MOTOR CORTEX 25% LH -- %%
tmp = NaN(size(sourcedata_motor.lh_25,1)+1, size(sourcedata_motor.lh_25,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:11)';
tmp(2:end,2:end)  = sourcedata_motor.lh_25;

table_motor.lh25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_supp_fig_4a_pta_motor_lh25.xlsx';
writetable(table_motor.lh25,filename, 'Sheet', 1);

% -- MOTOR CORTEX 75% LH -- %%
tmp = NaN(size(sourcedata_motor.lh_75,1)+1, size(sourcedata_motor.lh_75,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:11)';
tmp(2:end,2:end)  = sourcedata_motor.lh_75;

table_motor.lh75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_supp_fig_4a_pta_motor_lh75.xlsx';
writetable(table_motor.lh75,filename, 'Sheet', 1);

%% - RUN STATS - %%

% Get data of interest
q       = 'Statistic: Frontal = 1; Motor = 2: ';
selroi = input(q);

if selroi == 1
    
    statsdata = data_frontal;
    
elseif selroi == 2
    
    statsdata = data_motor;
    
end
    
q           = 'Statistic: HLL = 1; BR = 2: ';
ans_predict = input(q);

if ans_predict == 1
        
    data_lh0  = statsdata.hll_lh0;
    data_lh25 = statsdata.hll_lh25;
    data_lh75 = statsdata.hll_lh75;
    
    text_2 = 'HLL';
    
elseif ans_predict == 2
    
    data_lh0  = statsdata.br_lh0;
    data_lh25 = statsdata.br_lh25;
    data_lh75 = statsdata.br_lh75;
    
    text_2 = 'BR';
    
end


% Run stats

cfg                     = [];
cfg.channel             = 'all';
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
subj = numel(data_lh0);
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

stat = ft_freqstatistics(cfg, data_lh0{:}, data_lh25{:}, data_lh75{:});    

%% - COMPUTE PAIRWISE COMPARISONS - %%

cfg                     = [];
cfg.channel             = 'all';
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

subj = numel(data_lh0);
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

stat_0vs25  = ft_freqstatistics(cfg, data_lh25{:}, data_lh0{:});

stat_0vs75  = ft_freqstatistics(cfg, data_lh75{:}, data_lh0{:});

stat_25vs75 = ft_freqstatistics(cfg, data_lh75{:}, data_lh25{:});

%% Grand Average

figure;

color = {'#008450', '#EFB700', '#B81D13'};

cfg = [];
cfg.keepindividual = 'yes';

data_0  = ft_freqgrandaverage(cfg, data_lh0{:});
data_25 = ft_freqgrandaverage(cfg, data_lh25{:});
data_75 = ft_freqgrandaverage(cfg, data_lh75{:});

% -------------
% Plot Data
% -------------

plot(data_0.freq, mean(squeeze(data_0.powspctrm)), 'Color', color{1}, 'LineWidth', 4);
hold on
plot(data_25.freq, mean(squeeze(data_25.powspctrm)), 'Color', color{2}, 'LineWidth', 4);
hold on
plot(data_25.freq, mean(squeeze(data_75.powspctrm)), 'Color', color{3}, 'LineWidth', 4);

hold on

% ----------
% Mask Stats
% ----------

% plot significant timepoints
a = ylim;

island = bwconncomp(stat.mask);
for jj = 1:island.NumObjects
    h = line(stat.freq(island.PixelIdxList{jj}), ...
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
    h = line(stat.freq(island.PixelIdxList{jj}), ...
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
    h = line(stat.freq(island.PixelIdxList{jj}), ...
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
    h = line(stat.freq(island.PixelIdxList{jj}), ...
        repmat(a(2)+0.5, 1, numel(island.PixelIdxList{jj})), 'linewidth', 4);
    b2 = plot(h.XData(1:round(end/2)), h.YData(1:round(end/2)), 'LineWidth', 4);
    b2.Color = '#EFB700';
    hold on
    b2 = plot(h.XData(round(end/2):end), h.YData(round(end/2):end), 'LineWidth', 4);
    b2.Color = '#B81D13';
end
clear island

a = ylim;

% make plot look nice
xlabel('Frequency [Hz]');
ylabel('Power')
set(gca, 'linewidth', 3, 'fontsize', 13, 'xlim', [1 15]);
box off
set(gcf, 'Position', [489 361 270 201]);

%% - SAVE FIGURES - %%

% save figures
if selroi == 1
    cd(path_out);
    print(gcf,'fig_4_pta_psd_frontal.pdf','-dpdf','-r400');
else
    cd(path_out_supps);
    print(gcf,'fig_4_pta_psd_motor.pdf','-dpdf','-r400');
end




