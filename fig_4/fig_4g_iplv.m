%% Figure 4g
% this script rungs additional analysis required for the plots of figure
% 4g.

%%

clear
ft_defaults

%% 

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_4/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

path_in  = '/Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/PLV/';

%% - LOAD DATA - %%

% Select the data: Differs in how the PLV surrogate distribution was derived

data     = 'imag_plv_trialshuffling.mat';
datapeak = 'peakfreq_plv_trialshuffling.mat';

load(fullfile(path_in, data));

% load peak information
load(fullfile(path_in, datapeak));

%% - LOAD ROI INFORMATION - %%

load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/Electrode_Information_ROI/elec_info_ROI.mat

%% - CREATE FREQUENCY VECTOR - %%

minfreq = 3;
maxfreq = 32;

centerfreq = logspace(log10(minfreq), log10(maxfreq), (maxfreq - minfreq)*2^(1/8));
foi = [centerfreq - (centerfreq/4); centerfreq + (centerfreq/4)]';

%% - BOOTSTRAP PARAMS - %%

nBoot = 1000; % number of bootstrap iterations

%% - EXTRACT DATA - %%

% preallocate memory
iplv_f_motor     = struct();

% extract data per subject

for Isub = 1:numel(imagPLV)
        
    fprintf('%d/%d\n', Isub, numel(imagPLV));
    
    if isempty(imagPLV{Isub})
        continue;
    end
    
    % get trial indices
    idx_lh0    = imagPLV{Isub}.idx_lh0;
    idx_lh25   = imagPLV{Isub}.idx_lh25;
    idx_lh75   = imagPLV{Isub}.idx_lh75;

    % ------------------------------------
    % Frontal & Motor
    % ------------------------------------
    
    if isfield(imagPLV{Isub}, 'f_motor')
                
        % obtain bootstrap distribution from surrogate data 
        % (channel x trial x freq)
        bootnormdata = NaN(size(imagPLV{Isub}.f_motor.true,2), size(imagPLV{Isub}.f_motor.true,3), ...
            length(centerfreq));
        
        for Ichan = 1:size(imagPLV{Isub}.f_motor.true,2) % loop over channels
            
            if isnan(peak_plv_f_motor{Isub}(Ichan))
                continue;
            end
            
            tmpData = squeeze(imagPLV{Isub}.f_motor.true(:,Ichan,:,:)); % current true data
            tmpData_surro = squeeze(imagPLV{Isub}.f_motor.surro(:,1,Ichan,:)); % current surrogate data

            % get bootstrap distribution
            bootDist = NaN(nBoot, length(centerfreq));
            for Iboot = 1:nBoot
                x = randperm(size(tmpData_surro,1));
                % take the mean over as many trials as in real data
                bootDist(Iboot,:) = mean(tmpData_surro(x(1:size(tmpData,1)),:));
            end
            
            % get mean and SD of bootstrap distribution
            mu = mean(bootDist);
            sigma = std(bootDist);
            
            % z-score data with respect to bootstrapped distribution
            y = (tmpData - mu) ./ sigma;
            
            bootnormdata(Ichan,:,:) = y;
                
        end

        clear count
        
        % get mean for conditions
        iplv_f_motor.all{Isub}.data       = squeeze(nanmean(nanmean(bootnormdata,1),2))';
        iplv_f_motor.all{Isub}.freq       = centerfreq;
        iplv_f_motor.all{Isub}.dimord     = 'chan_freq';
        iplv_f_motor.all{Isub}.label      = {'mtl_motor'};
        iplv_f_motor.all{Isub}.cfg        = [];

        iplv_f_motor.lh_0{Isub}.data      = squeeze(nanmean(nanmean(bootnormdata(:,idx_lh0,:),1),2))';
        iplv_f_motor.lh_0{Isub}.freq      = centerfreq;
        iplv_f_motor.lh_0{Isub}.dimord    = 'chan_freq';
        iplv_f_motor.lh_0{Isub}.label     = {'f_motor'};
        iplv_f_motor.lh_0{Isub}.cfg       = [];
        
        iplv_f_motor.lh_25{Isub}.data     = squeeze(nanmean(nanmean(bootnormdata(:,idx_lh25,:),1),2))';
        iplv_f_motor.lh_25{Isub}.freq     = centerfreq;
        iplv_f_motor.lh_25{Isub}.dimord   = 'chan_freq';
        iplv_f_motor.lh_25{Isub}.label    = {'f_motor'};
        iplv_f_motor.lh_25{Isub}.cfg      = [];
        
        iplv_f_motor.lh_75{Isub}.data     = squeeze(nanmean(nanmean(bootnormdata(:,idx_lh75,:),1),2))';
        iplv_f_motor.lh_75{Isub}.freq     = centerfreq;
        iplv_f_motor.lh_75{Isub}.dimord   = 'chan_freq';
        iplv_f_motor.lh_75{Isub}.label    = {'f_motor'};
        iplv_f_motor.lh_75{Isub}.cfg      = [];
 
    end     
    
end

%% - SAVE SOURCE DATA - %%

sourcedata = struct;

freq = iplv_f_motor.lh_0{2}.freq;

for Isub = 1:numel(iplv_f_motor.lh_0)
    
    if isempty(iplv_f_motor.lh_0{Isub})
        
        sourcedata.lh_0(Isub,:)  = NaN(1, length(freq));
        sourcedata.lh_25(Isub,:) = NaN(1, length(freq));
        sourcedata.lh_75(Isub,:) = NaN(1, length(freq));

    else
        
        sourcedata.lh_0(Isub,:)  = iplv_f_motor.lh_0{Isub}.data;
        sourcedata.lh_25(Isub,:) = iplv_f_motor.lh_25{Isub}.data;
        sourcedata.lh_75(Isub,:) = iplv_f_motor.lh_75{Isub}.data;

    end        
        
end%Isub

% go into results folder
mkdir(path_source, 'Fig_4g');
cd(fullfile(path_source, 'Fig_4g'));

% -- 0% LH -- %%
tmp = NaN(size(sourcedata.lh_0,1)+1, size(sourcedata.lh_0,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_0;

table_frontal.lh0 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_4g_iplv_lh0.xlsx';
writetable(table_frontal.lh0,filename, 'Sheet', 1);

% -- 25% LH -- %%
tmp = NaN(size(sourcedata.lh_25,1)+1, size(sourcedata.lh_25,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_25;

table_frontal.lh25 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_4g_iplv_lh25.xlsx';
writetable(table_frontal.lh25,filename, 'Sheet', 1);

% -- 75% LH -- %%
tmp = NaN(size(sourcedata.lh_75,1)+1, size(sourcedata.lh_75,2)+1);
tmp(1,2:end)      = freq;
tmp(2:end,1)      = (1:17)';
tmp(2:end,2:end)  = sourcedata.lh_75;

table_frontal.lh75 = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_4g_iplv_lh75.xlsx';
writetable(table_frontal.lh75,filename, 'Sheet', 1);

%% - REMOVE EMPTY CELLS - %%

% frontal & motor
iplv_f_motor.all         = iplv_f_motor.all(~cellfun(@isempty, iplv_f_motor.all));
iplv_f_motor.lh_0        = iplv_f_motor.lh_0(~cellfun(@isempty, iplv_f_motor.lh_0));
iplv_f_motor.lh_25       = iplv_f_motor.lh_25(~cellfun(@isempty, iplv_f_motor.lh_25));
iplv_f_motor.lh_75       = iplv_f_motor.lh_75(~cellfun(@isempty, iplv_f_motor.lh_75));

%% - GET GRAND AVERAGE - %%

cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter      = 'data';

%---- Grand Average for Frontal & Motor ----&

GA_f_motor   = struct();

% get condition-wise average
GA_f_motor.all    = ft_freqgrandaverage(cfg, iplv_f_motor.all{:});
GA_f_motor.lh0    = ft_freqgrandaverage(cfg, iplv_f_motor.lh_0{:});
GA_f_motor.lh25   = ft_freqgrandaverage(cfg, iplv_f_motor.lh_25{:});
GA_f_motor.lh75   = ft_freqgrandaverage(cfg, iplv_f_motor.lh_75{:});

% plot grand average 

figure;
h = shadedErrorBar(GA_f_motor.all.freq, squeeze(nanmean(GA_f_motor.all.data)), ...
    nanstd(squeeze(GA_f_motor.all.data))/sqrt(size(GA_f_motor.all.data,1)));
h.mainLine.Color = 'r'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = 'r'; h.patch.EdgeColor = 'r'; h.patch.FaceAlpha = 0.5;
ylabel('iPLV [z]'); xlabel('Frequency [Hz]'); xlim([centerfreq(1) centerfreq(end)]);
set(gca, 'linewidth', 3, 'fontsize', 13, 'xtick', [4 8 16 32], 'xscale', 'log');
box off

% save figures
cd(path_out);
print(gcf,'fig_4_iplv_f_motor_all.pdf','-dpdf','-r400');    

figure;
h = shadedErrorBar(GA_f_motor.lh0.freq, squeeze(nanmean(GA_f_motor.lh0.data)), ...
    nanstd(squeeze(GA_f_motor.lh0.data))/sqrt(size(GA_f_motor.lh0.data,1)));
h.mainLine.Color = '#008450'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#008450'; h.patch.EdgeColor = '#008450'; h.patch.FaceAlpha = 0.5;
hold on
h = shadedErrorBar(GA_f_motor.lh25.freq, squeeze(nanmean(GA_f_motor.lh25.data)), ...
    nanstd(squeeze(GA_f_motor.lh25.data))/sqrt(size(GA_f_motor.lh25.data,1)));
h.mainLine.Color = '#EFB700'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#EFB700'; h.patch.EdgeColor = '#EFB700'; h.patch.FaceAlpha = 0.5;
hold on
h = shadedErrorBar(GA_f_motor.lh75.freq, squeeze(nanmean(GA_f_motor.lh75.data)), ...
    nanstd(squeeze(GA_f_motor.lh75.data))/sqrt(size(GA_f_motor.lh75.data,1)));
h.mainLine.Color = '#B81D13'; h.mainLine.LineWidth = 1;
h.patch.FaceColor = '#B81D13'; h.patch.EdgeColor = '#B81D13'; h.patch.FaceAlpha = 0.5;
ylabel('normalized iPLV [z]'); xlabel('Frequency [Hz]'); xlim([centerfreq(1) centerfreq(end)]);
set(gca, 'linewidth', 3, 'fontsize', 13, 'xtick', [4 8 16 32], 'xscale', 'log');

% save figures
cd(path_out);
print(gcf,'fig_4_iplv_f_motor_percond.pdf','-dpdf','-r400');    


