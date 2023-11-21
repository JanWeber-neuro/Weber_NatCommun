%% Figure S12

%%

clear
close all;
ft_defaults;

%% -SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in  = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_S9/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

% add functions
addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% - GENERAL SETTINGS - %%

nPad      = 2;     % zero paddings in seconds for PSD analysis

fsample   = 512;   % sampling frequency

smoothwin = 0.025; % window for smoothing decoding timeseries

%% - LOAD DATA - %%

% preallocate data for data with theta dimension
ThetaDim  = struct();

% preallocate memory for 1/f corrected PSDs
PSD       = struct();

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
        
        if isfield(all_data, 'frontal')                                                
            
            % pass data
            ThetaDim.frontal{Isub} = all_data.frontal;
                        
        end
        
        % #########################
        % MOTOR CORTEX
        % #########################
        
        if isfield(all_data, 'motor')
            
            % pass data
            ThetaDim.motor{Isub} = all_data.motor;

            clear all_data
            
        end
            
    end
    
end %Isub

%% - EXTRACT PEAK LATENCIES OF DECODING ACCURACY - %%

% preallocate memory
peakdata = struct();

% set time window
time = all_data.frontal.time;
idx  = nearest(time,-0.5):nearest(time,0.3);

% - [EXTRACT PEAK LATENCY FOR CONTEXT DECODING PFC] - %

peakdata.frontal.context = NaN(18,1);
peakdata.frontal.action  = NaN(18,1);
peakdata.motor.action    = NaN(18,1);

for Isub = 1:numel(ThetaDim.frontal)

    if ~isempty(ThetaDim.frontal{Isub}) && ~isempty(ThetaDim.frontal{Isub}.stimdim) && ...
            ~isempty(ThetaDim.frontal{Isub}.actdim)
        
        count = count + 1;

        % get current data
        tmp = smoothdata(ThetaDim.frontal{Isub}.DecAcc_stim(...
            ThetaDim.frontal{Isub}.stimdim,idx), 'gaussian', [round(smoothwin*fsample) round(smoothwin*fsample)]);

        % find peak decoding latency
        [pks, locs] = findpeaks(tmp, time(idx));
        [~,x] = max(pks);
        
        % save latency
        peakdata.frontal.context(Isub) = locs(x); clear tmp
   
    end

end

% - [EXTRACT PEAK LATENCY FOR ACTION DECODING PFC] - %

count = 0;

for Isub = 1:numel(ThetaDim.frontal)

    if ~isempty(ThetaDim.frontal{Isub}) && ~isempty(ThetaDim.frontal{Isub}.actdim)

        count = count + 1;

        % get current data
        tmp = smoothdata(ThetaDim.frontal{Isub}.DecAcc_act(...
            ThetaDim.frontal{Isub}.actdim,idx), 'gaussian', [round(smoothwin*fsample) round(smoothwin*fsample)]);

        % find peak decoding latency
        [pks, locs] = findpeaks(tmp, time(idx));
        [~,x] = max(pks);
        
        % save latency
        peakdata.frontal.action(Isub) = locs(x); clear tmp
   
    end

end

% - [EXTRACT PEAK LATENCY FOR ACTION DECODING MOTOR CORTEX] - %

count = 0;

for Isub = 1:numel(ThetaDim.motor)

    if ~isempty(ThetaDim.motor{Isub}) && ~isempty(ThetaDim.motor{Isub}.actdim)

        count = count + 1;

        % get current data
        tmp = smoothdata(ThetaDim.motor{Isub}.DecAcc_act(...
            ThetaDim.motor{Isub}.actdim,idx), 'gaussian', [round(smoothwin*fsample) round(smoothwin*fsample)]);

        % find peak decoding latency
        [pks, locs] = findpeaks(tmp, time(idx));
        [~,x] = max(pks);
        
        % save latency
        peakdata.motor.action(Isub) = locs(x); clear tmp
   
    end

end

% descriptives
fprintf('peak context integration PFC = %.3f ± %.3f sec.\n', nanmean(peakdata.frontal.context), nanstd(peakdata.frontal.context));
fprintf('peak action encoding PFC     = %.3f ± %.3f sec.\n', nanmean(peakdata.frontal.action), nanstd(peakdata.frontal.action));
fprintf('peak action encoding M1      = %.3f ± %.3f sec.\n', nanmean(peakdata.motor.action), nanstd(peakdata.motor.action));

%% - SAVE SOURCE DATA - %%

% go into results folder
mkdir(path_source, 'Fig_S10');
cd(fullfile(path_source, 'Fig_S10'));

table_data = array2table([peakdata.frontal.context peakdata.frontal.action peakdata.motor.action], ...
    'VariableNames', {'PFC_Context', 'PFC_Action', 'Motor_Action'});

% write data into excelsheet
filename = 'fig_S10.xlsx';
writetable(table_data,filename, 'Sheet', 1);

%% - REMOVE NAN FOR PLOTTING AND STATISTICS - %%

peakdata.frontal.context(isnan(peakdata.frontal.context)) = [];
peakdata.frontal.action(isnan(peakdata.frontal.action)) = [];
peakdata.motor.action(isnan(peakdata.motor.action)) = [];

%% - PLOT PEAKS USING BAR PLOTS - %%

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/colormaps/

figure;
cmap = viridis(100);
col2use = [cmap(1,:); cmap(20,:); cmap(80,:)];

data2plot = {peakdata.frontal.context, peakdata.frontal.action, ...
    peakdata.motor.action};

for jj = 1:numel(data2plot)

    a = bar(jj, mean(data2plot{jj}));
    a.FaceColor = col2use(jj,:);
    a.EdgeColor = col2use(jj,:);
    a.FaceAlpha = 0.5;
    a.EdgeAlpha = 0.5;

    hold on

    vertical_errorbars(jj, 0.1, mean(data2plot{jj}), std(data2plot{jj}) / sqrt(length(data2plot{jj})), col2use(jj,:), 2)    

    jw_scatterplot({data2plot{jj}}, 0, jj-1+0.4, {col2use(jj,:)});

end

set(gca, 'fontsize', 13);

xticks('');

box off

% compute statistics 

statsdata = [peakdata.frontal.context peakdata.frontal.action peakdata.motor.action];
groupidx  = [ones(1,length(peakdata.frontal.context)) ...
             ones(1,length(peakdata.frontal.action))*2 ...
             ones(1,length(peakdata.motor.action))*3];

p = kruskalwallis(statsdata, groupidx, 'off');

fprintf('Kruskal Wallis: Main effect p = %.3f\n', p);

% go to output folder where figure should be stored
cd(path_out)

print(gcf,'peak_decoding_latencies_context_action_barplot.pdf','-dpdf','-r400');    

%% - PLOT PEAKS USING RAINCLOUD PLOTS - %%

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/colormaps/
addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

cmap = viridis(100);
col2use = [cmap(1,:); cmap(20,:); cmap(80,:)];

data2plot = {peakdata.frontal.context', peakdata.frontal.action', ...
    peakdata.motor.action'}';

% color 2 use for plotting
cb = {cmap(1,:), cmap(20,:), cmap(80,:)};

figure; 
h = rm_raincloud(data2plot, [1 0 0], 0, 'ks', 0.15);

% admin
for Icond = 1:numel(data2plot)
    
    % change density plot color
    h.p{Icond,1}.FaceColor        = cb{Icond};
    h.p{Icond,1}.EdgeColor        = cb{Icond};
    h.p{Icond,1}.FaceAlpha        = 0.8;

    % change color for indiv. dots
    h.s{Icond,1}.MarkerFaceColor  = cb{Icond};
    h.s{Icond,1}.MarkerFaceAlpha  = 0.8;
    h.s{Icond,1}.SizeData         = 40;

    % change mean dots
    h.m(Icond).MarkerFaceColor    = cb{Icond};
    h.m(Icond).SizeData           = 80;
    
end

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];
    
% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks('');
xline(0, 'LineStyle', '--');
set(gca, 'fontsize', 13);
box off

view([0 90]);

% compute statistics 

statsdata = [peakdata.frontal.context peakdata.frontal.action peakdata.motor.action];
groupidx  = [ones(1,length(peakdata.frontal.context)) ...
             ones(1,length(peakdata.frontal.action))*2 ...
             ones(1,length(peakdata.motor.action))*3];

p = kruskalwallis(statsdata, groupidx, 'off');

fprintf('Kruskal Wallis: Main effect p = %.3f\n', p);

% go to output folder where figure should be stored
cd(path_out)

print(gcf,'peak_decoding_latencies_context_action_raincloud.pdf','-dpdf','-r400');    
