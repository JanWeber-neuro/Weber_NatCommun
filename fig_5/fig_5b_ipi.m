%% Figure 5b
% this script rungs additional analysis required for the plots of figure
% 5b.

%% - HOUSE KEEPING - %%

clear
close all;

ft_defaults;

%% -SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path.in  = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';
path.fun = '/Users/janweber/Documents/MATLAB/general_functions/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath(path.fun); % add functions

%% - GENERAL SETTINGS - %%

nPad      = 2;     % zero paddings in seconds for PSD analysis

fsample   = 512;   % sampling frequency

%% - LOAD ELECTRODE LOCATION - %%

load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/Electrode_Information_ROI/elec_info_ROI.mat

%% - COMPUTE PC DIMENSION WITH STRONGEST THETA POWER - %%

% initialize peak components

avgipi.frontal = cell(numel(subjects),1);
avgipi.motor   = cell(numel(subjects),1);

for Isub = 1:numel(subjects)
    
    % subject ID
    subID = subjects{Isub};
    
    fprintf('%% ----------------------------------------------------------- %%\n');
    
    fprintf('%% Computation running for %s\n', subID);

    fprintf('%% ----------------------------------------------------------- %%\n');

    
    if isfile(fullfile(path.in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']))
        
        load(fullfile(path.in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']));

        % #########################
        % FRONTAL CORTEX
        % #########################
        
        if isfield(all_data, 'frontal')
            
            % ------------------------- 
            % Bring data into FT format
            % -------------------------

            ftdata = struct();
            ftdata.trial  = all_data.frontal.trial;
            ftdata.label  = sprintfc('PC%d',1:size(ftdata.trial,2));
            ftdata.dimord = all_data.frontal.dimord; % dimensions
            ftdata.time   = all_data.frontal.time;
            ftdata.cfg    = [];
            
            % ------------------------- 
            % Compute IRASA
            % ------------------------- 
            
            PadData = struct();              

            for chani = 1:numel(ftdata.label)

                % zero padding
                tmptrial  = squeeze(ftdata.trial(:,chani,:));

                % make FT structure
                PadData.trial(:,chani,:)   = [zeros(size(tmptrial,1),fsample*nPad) tmptrial zeros(size(tmptrial,1),fsample*nPad)];
                PadData.time               = 0:1/fsample:1/fsample*(size(PadData.trial,3)-1);
                PadData.dimord             = 'rpt_chan_time';
                PadData.label(chani)       = ftdata.label(chani);

            end

            % actual analysis

            cfg               = [];
            cfg.foilim        = [2 13];
            cfg.channel       = PadData.label;
            cfg.taper         = 'hanning';
            cfg.pad           = 'nextpow2';
            cfg.keeptrials    = 'yes';
            cfg.method        = 'irasa';
            powspec.frac      = ft_freqanalysis(cfg, PadData); % fractal part
            cfg.method        = 'mtmfft';
            powspec.orig      = ft_freqanalysis(cfg, PadData); % original spectrum

            clear PadData

            % correct for 1/f by subtracting fractal from original spectrum
            cfg               = [];
            cfg.parameter     = 'powspctrm';
            cfg.operation     = 'x2-x1';
            powspec.osci      = ft_math(cfg, powspec.frac, powspec.orig);
            powspec.osci. ...
                powspctrm(powspec.osci.powspctrm < 0) = 0;

            % ------------------------- 
            % Identify Peak 
            % ------------------------- 
            
            % loop over PCs
            
            peakamp = NaN(size(powspec.osci.powspctrm,2),1);
            
            for Icomp = 1:size(powspec.osci.powspctrm,2)
            
                % find peaks
                [pks,locs] = findpeaks(squeeze(mean(powspec.osci.powspctrm(:,Icomp,:))));
                
                peakamp(Icomp) = max(pks);
                
            end
            
            % find PC with max peak
            [~, peakPC] = max(peakamp); clear peakamp
            
            
            % admin
            if ~isempty(peakPC)
                                                
                % ------------------------------
                % Compute IPI
                % ------------------------------
                
                % current data (trial x time)
                tmpData = squeeze(all_data.frontal.trial(:,peakPC,:));
                
                ipi = NaN(size(tmpData,1),1);
                
                for Itrial = 1:size(tmpData)
                    
                    trialdata = tmpData(Itrial,:);
                    
                    [val, idx] = findpeaks(trialdata);

                    if ~isempty(idx)

                        tmp = fsample ./ diff(idx); % get frequency between peaks
                        tmp(isoutlier(tmp)) = [];   % remove outliers
                        ipi(Itrial) = mean(tmp);
                        
                    end
                    
                end %Itrial
                
                avgipi.frontal{Isub} = ipi;
                                
            end
                        
        end
        
        % #########################
        % MOTOR CORTEX
        % #########################
        
        if isfield(all_data, 'motor')
            
            % ------------------------- 
            % Bring data into FT format
            % -------------------------

            ftdata = struct();
            ftdata.trial  = all_data.motor.trial;
            ftdata.label  = sprintfc('PC%d',1:size(ftdata.trial,2));
            ftdata.dimord = all_data.motor.dimord; % dimensions
            ftdata.time   = all_data.motor.time;
            ftdata.cfg    = [];
            
            % ------------------------- 
            % Compute IRASA
            % ------------------------- 
            
            PadData = struct();              

            for chani = 1:numel(ftdata.label)

                % zero padding
                tmptrial  = squeeze(ftdata.trial(:,chani,:));

                % make FT structure
                PadData.trial(:,chani,:)   = [zeros(size(tmptrial,1),fsample*nPad) tmptrial zeros(size(tmptrial,1),fsample*nPad)];
                PadData.time               = 0:1/fsample:1/fsample*(size(PadData.trial,3)-1);
                PadData.dimord             = 'rpt_chan_time';
                PadData.label(chani)       = ftdata.label(chani);

            end

            % actual analysis

            cfg               = [];
            cfg.foilim        = [2 13];
            cfg.channel       = PadData.label;
            cfg.taper         = 'hanning';
            cfg.pad           = 'nextpow2';
            cfg.keeptrials    = 'yes';
            cfg.method        = 'irasa';
            powspec.frac      = ft_freqanalysis(cfg, PadData); % fractal part
            cfg.method        = 'mtmfft';
            powspec.orig      = ft_freqanalysis(cfg, PadData); % original spectrum

            clear PadData

            % correct for 1/f by subtracting fractal from original spectrum
            cfg               = [];
            cfg.parameter     = 'powspctrm';
            cfg.operation     = 'x2-x1';
            powspec.osci      = ft_math(cfg, powspec.frac, powspec.orig);
            powspec.osci. ...
                powspctrm(powspec.osci.powspctrm < 0) = 0;

            % ------------------------- 
            % Identify Peak 
            % ------------------------- 
            
            % loop over PCs
            
            peakamp = NaN(size(powspec.osci.powspctrm,2),1);
            
            for Icomp = 1:size(powspec.osci.powspctrm,2)
            
                % find peaks
                [pks,locs] = findpeaks(squeeze(mean(powspec.osci.powspctrm(:,Icomp,:))));
                
                peakamp(Icomp) = max(pks);
                
            end
            
            % find PC with max peak
            [~, peakPC] = max(peakamp); clear peakamp
            
            
            % admin
            if ~isempty(peakPC)
                                                
                % ------------------------------
                % Compute IPI
                % ------------------------------
                
                % current data (trial x time)
                tmpData = squeeze(all_data.motor.trial(:,peakPC,:));
                
                ipi = NaN(size(tmpData,1),1);
                
                for Itrial = 1:size(tmpData)
                    
                    trialdata = tmpData(Itrial,:);
                    
                    [val, idx] = findpeaks(trialdata);

                    if ~isempty(idx)

                        tmp = fsample ./ diff(idx); % get frequency between peaks
                        tmp(isoutlier(tmp)) = [];   % remove outliers
                        ipi(Itrial) = mean(tmp);
                        
                    end
                    
                end %Itrial
                
                avgipi.motor{Isub} = ipi;
                                
            end
                        
        end
        
    end
    
end %Isub

%% - REMOVE EMPTY CELLS - %%

avgipi.frontal = avgipi.frontal(~cellfun(@isempty, avgipi.frontal));
avgipi.motor   = avgipi.motor(~cellfun(@isempty, avgipi.motor));

%% - COMPUTE GRAND AVERAGE IPI - %%

GA_frontal = cat(1, avgipi.frontal{:});
GA_frontal(isoutlier(GA_frontal, 'percentiles', [0 95])) = []; % remove outliers

GA_motor = cat(1, avgipi.motor{:});
GA_motor(isoutlier(GA_motor, 'percentiles', [0 95])) = []; % remove outliers

figure;
a = histfit(GA_frontal);
a(1).FaceColor = 'r';
a(1).EdgeColor = a(1).FaceColor;
a(1).FaceAlpha = 1;
a(1).EdgeAlpha = 1;
a(1).LineWidth = 1;
a(2).Color     = 'k';

xlabel('IPI [Hz]');
ylabel('Count');

set(gca, 'linew', 0.5, 'fontsize', 13); box off;
set(gcf, 'Position', [524 394 230 193]);
            
cd(path.out);
print(gcf,'fig_6_histogram_thetadim_frontal_ipi.pdf','-dpdf','-r400');

figure;
a = histfit(GA_motor);
a(1).FaceColor = 'b';
a(1).EdgeColor = a(1).FaceColor;
a(1).FaceAlpha = 1;
a(1).EdgeAlpha = 1;
a(1).LineWidth = 1;
a(2).Color     = 'k';

xlabel('IPI [Hz]');
ylabel('Count');

set(gca, 'linew', 0.5, 'fontsize', 13); box off;
set(gcf, 'Position', [524 394 230 193]);
            
cd(path.out);
print(gcf,'fig_6_histogram_thetadim_motor_ipi.pdf','-dpdf','-r400');
                 
%% - SAVE SOURCE DATA - %%

% go into results folder
mkdir(path_source, 'Fig_6b');
cd(fullfile(path_source, 'Fig_6b'));

% - STORE DATA - %
table_ipi_frontal = array2table(GA_frontal); 
table_ipi_motor   = array2table(GA_motor); 

% write data into excelsheet
filename = 'fig_6b_thetadim_ipi_pfc.xlsx';
writetable(table_ipi_frontal,filename, 'Sheet', 1);

filename = 'fig_6b_thetadim_ipi_motor.xlsx';
writetable(table_ipi_motor,filename, 'Sheet', 1);

       
            
            