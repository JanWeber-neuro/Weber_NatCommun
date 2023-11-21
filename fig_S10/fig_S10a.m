%% Figure S10a

%%

clear
ft_defaults

%% SUBJECTS

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% PATH SETTINGS

path_data = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_S8/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

%% LOAD DATA

% preallocate memory
truedata   = struct();
surrodata  = struct();

fsample   = 512;

smoothwin = 0.05; % window for smoothing

for Isub = 1:numel(subjects)
    
    subID = subjects{Isub};
    
    if isfile(fullfile(path_data, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']))
        
        load(fullfile(path_data, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']));
        
        % --------------------------------
        % Frontal Cortex
        % --------------------------------
        
        if isfield(all_data, 'frontal') && ~isempty(all_data.frontal.actdim) && ~isempty(all_data.frontal.stimdim)
            
            truedata.frontal{Isub}.avg          = smoothdata(all_data.frontal.DecAcc_stim(all_data.frontal.actdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
            truedata.frontal{Isub}.time         = all_data.frontal.time;
            truedata.frontal{Isub}.dimord       = 'chan_time';
            truedata.frontal{Isub}.label        = {'Frontal'};
            truedata.frontal{Isub}.cfg          = [];
            
            surrodata.frontal{Isub}.avg    = ones(1, length(all_data.frontal.time))*0.33;
            surrodata.frontal{Isub}.time   = all_data.frontal.time;
            surrodata.frontal{Isub}.dimord = 'chan_time';
            surrodata.frontal{Isub}.label  = {'Frontal'};
            surrodata.frontal{Isub}.cfg    = [];
        
        end
        
        % --------------------------------
        % Motor Cortex
        % --------------------------------
        
        if isfield(all_data, 'motor') && ~isempty(all_data.motor.actdim) && ~isempty(all_data.motor.stimdim)
            
            truedata.motor{Isub}.avg          = smoothdata(all_data.motor.DecAcc_stim(all_data.motor.actdim,:), 'movmean', [round(fsample*smoothwin)/2, round(fsample*smoothwin)/2]);
            truedata.motor{Isub}.time         = all_data.motor.time;
            truedata.motor{Isub}.dimord       = 'chan_time';
            truedata.motor{Isub}.label        = {'motor'};
            truedata.motor{Isub}.cfg          = [];
            
            surrodata.motor{Isub}.avg    = ones(1, length(all_data.motor.time))*0.33;
            surrodata.motor{Isub}.time   = all_data.motor.time;
            surrodata.motor{Isub}.dimord = 'chan_time';
            surrodata.motor{Isub}.label  = {'motor'};
            surrodata.motor{Isub}.cfg    = [];
        
        end
                            
    end
    
end %Isub
        
%% REMOVE EMPTY CELLS

truedata.frontal   = truedata.frontal(~cellfun(@isempty, truedata.frontal));
surrodata.frontal  = surrodata.frontal(~cellfun(@isempty, surrodata.frontal));

truedata.motor     = truedata.motor(~cellfun(@isempty, truedata.motor));
surrodata.motor    = surrodata.motor(~cellfun(@isempty, surrodata.motor));

%% - SAVE SOURCE DATA FRONTAL CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_S8a');
cd(fullfile(path_source, 'Fig_S8a'));

sourcedata = [];

time = truedata.frontal{1}.time;

for Isub = 1:numel(truedata.frontal)
    
    if isempty(truedata.frontal{Isub})
        
        sourcedata(Isub,:)  = NaN(1, length(time));

    else
        
        sourcedata(Isub,:)  = truedata.frontal{Isub}.avg(1,:);

    end        
        
end%Isub

% - STORE 0% LH - %
tmp = NaN(size(sourcedata,1)+1, size(sourcedata,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata,1))';
tmp(2:end,2:end)  = sourcedata;

table_frontal = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S8a_frontal.xlsx';
writetable(table_frontal,filename, 'Sheet', 1);

%% - SAVE SOURCE DATA MOTOR CORTEX - %%

% go into results folder
mkdir(path_source, 'Fig_S8a');
cd(fullfile(path_source, 'Fig_S8a'));

sourcedata = [];

time = truedata.motor{1}.time;

for Isub = 1:numel(truedata.motor)
    
    if isempty(truedata.motor{Isub})
        
        sourcedata(Isub,:)  = NaN(1, length(time));

    else
        
        sourcedata(Isub,:)  = truedata.motor{Isub}.avg(1,:);

    end        
        
end%Isub

% - STORE 0% LH - %
tmp = NaN(size(sourcedata,1)+1, size(sourcedata,2)+1);
tmp(1,2:end)      = time;
tmp(2:end,1)      = (1:size(sourcedata,1))';
tmp(2:end,2:end)  = sourcedata;

table_motor = array2table(tmp); clear tmp

% write data into excelsheet
filename = 'fig_S8a_motor.xlsx';
writetable(table_motor,filename, 'Sheet', 1);

%% STATS

% region of interest names
roinames = fieldnames(truedata);

for Iroi = 1:numel(roinames)

    % current ROI
    tmpROI = roinames{Iroi};

    % pass data
    statsdata_act   = truedata.(tmpROI);
    statssurro_act  = surrodata.(tmpROI);
        
    % compute stats
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
    
    subj = numel(statsdata_act);
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
    
    stat_act = ft_timelockstatistics(cfg, statssurro_act{:}, statsdata_act{:});
    
    % compute cohen's d as effect size
    cfg            = [];
    cfg.parameter  = 'avg';
    cfg.method     = 'analytic';
    cfg.statistic  = 'cohensd'; % see FT_STATFUN_COHENSD
    cfg.ivar       = 1;
    cfg.design     = design;
    cfg.uvar       = 1;
    cfg.ivar       = 2;
    
    d_act         = ft_timelockstatistics(cfg, statssurro_act{:}, statsdata_act{:});
    
    % check positive clusters
    
    if isfield(stat_act, 'posclusters')
    
        for Icluster = 1:numel(stat_act.posclusters)
    
            fprintf('%%---------------------------------------------------\n');
            fprintf('%% POSITIVE CLUSTERS---------------------------------\n');
            fprintf('%%---------------------------------------------------\n');
    
            fprintf('\n');
    
            fprintf('p-value cluster %d = %.3f\n', Icluster, stat_act.posclusters(Icluster).prob);
            fprintf('clustermass cluster %d = %.3f\n', Icluster, stat_act.posclusters(Icluster).clusterstat);
    
            begeffect = stat_act.time...
                (find(stat_act.posclusterslabelmat == Icluster, 1, 'first'));
            endeffect = stat_act.time...
                (find(stat_act.posclusterslabelmat == Icluster, 1, 'last'));
    
            fprintf('Action Decoding Effect: cluster %d from %.3f to %.3f ms\n', Icluster, begeffect, endeffect);
    
            fprintf('Cohens d cluster %d = %.3f\n', Icluster, nanmean(d_act.cohensd(...
                find(stat_act.posclusterslabelmat == Icluster, 1, 'first'):...
                find(stat_act.posclusterslabelmat == Icluster, 1, 'last'))));
    
            fprintf('\n');
    
        end
    
    end
    
    %% GRAND AVERAGE
    
    % get grand average
    cfg = [];
    cfg.keepindividual = 'yes';
    
    GA_Act   = ft_timelockgrandaverage(cfg, statsdata_act{:});
    GA_Surro   = ft_timelockgrandaverage(cfg, statssurro_act{:});
    
    % plot data
    time = GA_Act.time;
    
    figure;
    a = shadedErrorBar(time, squeeze(mean(GA_Act.individual)), ...
        std(squeeze(GA_Act.individual)) / sqrt(size(GA_Act.individual,1)));
    a.mainLine.Color = 'k'; a.mainLine.LineWidth = 3;
    a.patch.FaceColor = 'k'; a.patch.EdgeColor =' k'; a.patch.FaceAlpha = 0.5;
    
    hold on
    
    a = shadedErrorBar(time, squeeze(mean(GA_Surro.individual)), ...
        std(squeeze(GA_Surro.individual)) / sqrt(size(GA_Surro.individual,1)));
    a.mainLine.Color = 'r'; a.mainLine.LineWidth = 3;
    a.patch.FaceColor = 'r'; a.patch.EdgeColor =' r'; a.patch.FaceAlpha = 0.5;
    
    
    xlabel('Time to BR [s]');
    ylabel('Classifier Accuracy');
    yline(0.33, 'LineWidth', 3, 'LineStyle', '--');
    set(gca, 'fontsize', 13, 'xlim', [time(1) time(end)], 'ylim', [0.27 0.47]); 
    box off
    
    if isfield(stat_act, 'posclusters')
        
        ncluster = length(stat_act.posclusters);
    
        for Icluster = 1:ncluster    
    
            if stat_act.posclusters(Icluster).prob < 0.025
    
                idx = find(stat_act.posclusterslabelmat == Icluster);
    
                % shade grey where significant
                sigpnts = stat_act.time(idx);
    
                hold on
    
                scatter(sigpnts, ones(1, length(sigpnts))*0.31, 30, 'filled',...
                    'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.2);
                clear sigpnts
    
            end
    
        end
    
    end
    
    set(gcf, 'Position', [596 417 229 198]);

    fname = ['decoding_context_from_action_' tmpROI '.pdf'];

    cd(path_out);
    
    print(gcf,fname,'-dpdf','-r400');    

end



