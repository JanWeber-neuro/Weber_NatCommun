%% Figure 5d
% this script rungs additional analysis required for the plots of figure
% 5d.

%%

clear
close all;
ft_defaults;

%% -SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_5/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% - GENERAL SETTINGS - %%

nPad      = 2;     % zero paddings in seconds for PSD analysis

fsample   = 512;   % sampling frequency

smoothwin = 0.025; % window for smoothing decoding timeseries

%% - LOAD DATA - %%

zval    = cell(numel(subjects),1);
rts     = cell(numel(subjects),1);
trlinfo = cell(numel(subjects),1);

for Isub = 1:numel(subjects)
    
    % subject ID
    subID = subjects{Isub};
    
    fprintf('%% ----------------------------------------------------------- %%\n');
    
    fprintf('%% Computation running for %s\n', subID);

    fprintf('%% ----------------------------------------------------------- %%\n');
    
    if isfile(fullfile(path_in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']))
        
        load(fullfile(path_in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']));

        % #########################
        % FRONTAL CORTEX
        % #########################
        
        if isfield(all_data, 'frontal') && isfield(all_data, 'motor')
            
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
                
                all_data.frontal.peakdim = peakPC;
                
            else
                
                all_data.frontal.peakdim = [];
                
            end
                    
            % #########################
            % MOTOR CORTEX
            % #########################

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

                all_data.motor.peakdim = peakPC;

            else

                all_data.motor.peakdim = [];

            end
                       
            % ------------------------- 
            % Compute power correlation
            % ------------------------- 

            if ~isempty(all_data.frontal.peakdim) && ~isempty(all_data.motor.peakdim)

               
                % get reaction times
                rts{Isub} = all_data.frontal.trialinfo(:,11);

                % get trial type info
                trlinfo{Isub} = all_data.frontal.trialinfo(:,13);

                % get data for both regions
                
                tidx = nearest(all_data.frontal.time, -0.5):nearest(all_data.frontal.time, 0.2);
                
                tmp_frontal = squeeze(all_data.frontal.trial(:, all_data.frontal.peakdim, tidx));

                tmp_motor   = squeeze(all_data.motor.trial(:, all_data.frontal.peakdim, tidx));

                % compute power correlation on a trial-by-trial basis

                rho = NaN(size(tmp_frontal, 1),1);

                for Itrial = 1:size(tmp_frontal, 1)

                    rho(Itrial) = corr(tmp_frontal(Itrial,:)', tmp_motor(Itrial,:)'); 

                end

                % create H0 by randomly cutting the data

                niter   = 1000;

                permrho = NaN(niter, size(tmp_frontal,1));

                for Iperm = 1:niter

                    for Itrial = 1:size(tmp_frontal,1)

                        % random time cut
                        rndcut = randsample(size(tmp_frontal,2),1);
                        
                        while rndcut == size(tmp_frontal,2)
                            rndcut = randsample(size(tmp_frontal,2),1);
                        end
                            
                        % shuffle data from frontal cortex
                        shuf_frontal = [tmp_frontal(Itrial,rndcut:end) tmp_frontal(Itrial,1:rndcut-1)];

                        r = corr(shuf_frontal', tmp_motor(Itrial,:)');

                        permrho(Iperm,Itrial) = r;
                        
                        clear rndcut

                    end

                end

                % z-transform rho values based on H0

                zval{Isub} = (rho - mean(permrho)') ./ std(permrho)';

            end
                
        end
        
    end
    
end %Isub

%% - SAVE SOURCE DATA Z-TRANSFORMED CONNECTIVITY - %%

% go into results folder
mkdir(path_source, 'Fig_6d');
cd(fullfile(path_source, 'Fig_6d'));

% - STORE DATA - %
table_zval = array2table(cellfun(@mean, zval)); clear tmp

% write data into excelsheet
filename = 'fig_6d_connectivity.xlsx';
writetable(table_zval,filename, 'Sheet', 1);

%% - REMOVE EMPTY CELLS - %%

zval    = zval(~cellfun(@isempty, zval));
rts     = rts(~cellfun(@isempty, rts));
trlinfo = trlinfo(~cellfun(@isempty, trlinfo));

%% - CONCATENATE ALL TRIALS - %%

all_zval = cat(1, zval{:}); % trial x 1
all_rts  = cat(1, rts{:});

%% - COMPUTE THE AVERAGE Z-VALUE - %%

close all;

c2use = {'#008450', '#EFB700', '#B81D13'};

% get average per subject 
avg_zval_all = cellfun(@mean, zval);

% # iterations for bootstrap distribution
nboot = 1000;

avg_zval_0   = NaN(numel(zval),1); 
avg_zval_25  = NaN(numel(zval),1); 
avg_zval_75  = NaN(numel(zval),1); 

for Isub = 1:numel(zval)
    
    % # number of trials within 75% condition
    min_trl = sum(trlinfo{Isub} == 75);

    bootdist_0  = NaN(nboot,1);
    bootdist_25 = NaN(nboot,1);

    for Iboot = 1:nboot

        % get shuffle index for 0%
        shufidx = randsample(find(trlinfo{Isub} == 0), min_trl);

        bootdist_0(Iboot) = mean(zval{Isub}(shufidx)); clear shufidx

        % get shuffle index for 25%
        shufidx = randsample(find(trlinfo{Isub} == 25), min_trl);

        bootdist_25(Iboot) = mean(zval{Isub}(shufidx)); clear shufidx

    end

    % get averages
    avg_zval_0(Isub)  = mean(bootdist_0);
    avg_zval_25(Isub) = mean(bootdist_25);
    avg_zval_75(Isub) = mean(zval{Isub}(trlinfo{Isub} == 75));
    
end %Isub
 
%% plot average z-value and test against 0

c2use = {'k'};

figure;

a = jw_scatterplot({avg_zval_all}, 0, 0, {'k'});
a.MarkerFaceColor = c2use{1};
a.MarkerEdgeColor = c2use{1};
a.MarkerEdgeAlpha = 0.5;
a.SizeData        = 80;

xlim([0.6 1.5]); ylim([-2 2]); yline(0, 'LineWidth', 2);
xlabel('Pow Corr.'); xticks('');
ylabel('rho [z]');
set(gca, 'linew', 2, 'fontsize', 13); box off;
set(gcf, 'Position', [547 457 227 197]);

p = signrank(avg_zval_all);
fprintf('Wilcoxon Test vs. 0; p = %.5f\n', p);

cd(path_out);
print(gcf,'fig_6_frontal_motor_thetadim_powcorr_zval.pdf','-dpdf','-r400');    

%% plot average z-value per condition

figure;

c2use = {'#008450', '#EFB700', '#B81D13'};

data2use = {avg_zval_0, avg_zval_25, avg_zval_75};

for ii = 1:3
    
    a = bar(ii, mean(data2use{ii}));
    a.FaceColor = c2use{ii};
    a.EdgeColor = c2use{ii};
    a.FaceAlpha = 0.7;
    a.EdgeAlpha = 0.7;
    a.LineWidth = 3;
    
    hold on
    
    vertical_errorbars(ii, 0.1, mean(data2use{ii}), std(data2use{ii}) / sqrt(length(data2use{ii})), c2use{ii}, 2)    

end
    
xticks(1:3);
xticklabels({'0%', '25%', '75%'}); xtickangle(45);
ylabel('rho');

set(gca, 'fontsize', 13, 'linewidth', 2);

box off

set(gcf, 'Position', [511 383 258 184]);

cd(path_out);
print(gcf,'fig_6_frontal_motor_thetadim_powcorr_zval_per_condition.pdf','-dpdf','-r400');    

%% plot histogram of all z values

addpath /Users/janweber/Documents/MATLAB/general_functions/colormaps/

c2use = 'k';

figure;
a = histfit(all_zval);
a(1).FaceColor = c2use;
a(1).EdgeColor = a(1).FaceColor;
a(1).FaceAlpha = 1;
a(1).EdgeAlpha = 1;
a(1).LineWidth = 3;

xlabel('rho [zval]');
ylabel('Count');
xline(0, 'LineStyle', '--', 'LineWidth', 1.5);

set(gca, 'linew', 2, 'fontsize', 13); box off;
set(gcf, 'Position', [547 457 227 197]);

cd(path_out);
print(gcf,'fig_6_frontal_motor_thetadim_powcorr_zval_hist.pdf','-dpdf','-r400');    

%% - SAVE SOURCE DATA Z-TRANSFORMED CONNECTIVITY - %%

% go into results folder
mkdir(path_source, 'Fig_6d');
cd(fullfile(path_source, 'Fig_6d'));

% - STORE DATA - %
table_zval = array2table(all_zval); clear tmp

% write data into excelsheet
filename = 'fig_6d_connectivity_histogram.xlsx';
writetable(table_zval,filename, 'Sheet', 1);
