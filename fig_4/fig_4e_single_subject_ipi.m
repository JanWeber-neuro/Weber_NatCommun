%% Figure 4e
% this script rungs additional analysis required for the plots of figure
% 4e.

%%

clear
ft_defaults

rng('default');

%% - PATH SETTINGS - %%

path_in = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/OSL24/Data/HFA/PCA/PCA_ROI_Data/';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/Figure_4/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/ % add functions

%% - LOAD DATA - %%
           
load(fullfile(path_in, 'OSL24_activedata_ROI_pca.mat'));

tmpData  = roidata.hll;

fsample  = 512;

%% - GET CONDITION INDEX - %%

if ~isempty(tmpData.frontal.inactive)

    idx_l0  = find(tmpData.frontal.inactive.trialinfo(:,13) == 0);
    idx_l25 = find(tmpData.frontal.inactive.trialinfo(:,13) == 25);
    idx_l75 = find(tmpData.frontal.inactive.trialinfo(:,13) == 75);

elseif ~isempty(tmpData.motor.inactive)

    idx_l0  = find(tmpData.motor.inactive.trialinfo(:,13) == 0);
    idx_l25 = find(tmpData.motor.inactive.trialinfo(:,13) == 25);
    idx_l75 = find(tmpData.motor.inactive.trialinfo(:,13) == 75);

elseif ~isempty(tmpData.mtl.inactive)

    idx_l0  = find(tmpData.mtl.inactive.trialinfo(:,13) == 0);
    idx_l25 = find(tmpData.mtl.inactive.trialinfo(:,13) == 25);
    idx_l75 = find(tmpData.mtl.inactive.trialinfo(:,13) == 75);

end

alltrialidx = {idx_l0, idx_l25, idx_l75}; % all trials concatenated

% get minimum number of trials within all conditions
min_numtrl  = min(cellfun('size',alltrialidx,1));

peakdistance  = struct(); 
peakamplitude = struct();

% loop over trials with different likelihood of stop
for jj = 1:numel(alltrialidx)

    %% Extract Frontal Inter-Peak Interval

    if ~isempty(tmpData.frontal.active)

        % get time for epoch
        time  = tmpData.frontal.active.time;

        % get data for current likelihood of stop
        dat   = tmpData.frontal.active.trial...
            (alltrialidx{jj},:,:);

        % make sure dimensions are correct
        if size(dat,1) == length(time)
            dat = permute(dat, [2 1]);
        end

        % randomly subsample the data according to the minimum
        % number of trials in one condition to not bias kurtosis
        subsample_trls = randperm(size(dat,1));
        dat = dat(subsample_trls(1:min_numtrl),:,:);

        % loop over channel
        for Ichan = 1:size(dat,2)

            % storage variable for peak interval
            peakinterval = [];
            % storage variable for peak amplitudes
            p_amplitude  = [];

            % loop over trials (restrict trial number to minimum
            % trial number across trials to not bias kurtosis)
            for triali = 1:min_numtrl

                % current data for the trial
                trialdata = squeeze(dat(triali,Ichan,:));

                [val, idx] = findpeaks(trialdata);

                if ~isempty(idx)

                    tmp = fsample ./ diff(idx); % get frequency between peaks
                    tmp(isoutlier(tmp)) = [];   % remove outliers
                    peakinterval = vertcat(peakinterval, tmp); 
                    val(isoutlier(val)) = [];
                    p_amplitude  = vertcat(p_amplitude, val);

                else

                    peakinterval = vertcat(peakinterval, NaN); 
                    p_amplitude  = vertcat(p_amplitude, NaN);

                end

            end % triali

            peakinterval(isoutlier(peakinterval)) = [];
            peakdistance.frontal{jj}{Ichan}  = peakinterval;

            peakamplitude.frontal{jj}{Ichan} = p_amplitude;

        end % chani

    else % if no data for ROI

        peakdistance.frontal{jj}  = NaN;
        peakamplitude.frontal{jj} = NaN;

    end

end % jj

%% - PLOT EXAMPLE - %%

chan2use = 1;

[f,xi] = ksdensity(peakdistance.frontal{1}{chan2use});
a = plot(xi, f,'LineWidth', 3); a.Color(4) = 0.8; a.Color = '#008450';
[~, idx] = max(f); hold on
a = xline(xi(idx), 'LineWidth', 4, 'LineStyle', '--');
a.Color = '#008450';
hold on
[f,xi] = ksdensity(peakdistance.frontal{2}{chan2use});
a = plot(xi, f,'LineWidth', 3); a.Color(4) = 0.8; a.Color = '#EFB700';
[~, idx] = max(f); hold on
a = xline(xi(idx), 'LineWidth', 4, 'LineStyle', '--');
a.Color = '#EFB700';
hold on
[f,xi] = ksdensity(peakdistance.frontal{3}{chan2use});
a = plot(xi, f,'LineWidth', 3); a.Color(4) = 0.8; a.Color = '#B81D13';
[~, idx] = max(f); hold on
a = xline(xi(idx), 'LineWidth', 4, 'LineStyle', '--');
a.Color = '#B81D13';

set(gca, 'fontsize', 13, 'linewidth', 3);
xlabel('IPI [Hz]');
ylabel('Probability');
box off
set(gcf, 'Position', [510 356 291 230]);

%% - SAVE SOURCE DATA - %%

% go into results folder
mkdir(path_source, 'Fig_4e');
cd(fullfile(path_source, 'Fig_4e'));

tmp = array2table(peakdistance.frontal{1}{chan2use}); 
filename = 'ipi_lh0.xlsx';
writetable(tmp,filename, 'Sheet', 1); clear tmp

tmp = array2table(peakdistance.frontal{2}{chan2use}); 
filename = 'ipi_lh25.xlsx';
writetable(tmp,filename, 'Sheet', 1); clear tmp

tmp = array2table(peakdistance.frontal{3}{chan2use}); 
filename = 'ipi_lh75.xlsx';
writetable(tmp,filename, 'Sheet', 1); clear tmp

%% - SAVE FIGURE - %%

% save figures
cd(path_out);
print(gcf,'fig_4_ipi_distribution.pdf','-dpdf','-r400');    
