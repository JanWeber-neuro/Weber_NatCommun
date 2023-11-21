%% Figure 3e & f (left panel)
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

path_in  = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_3/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% - LOAD TASK-ENCODING ELECS - %%

load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/HFA/ExpVariance_Approach/TaskInformative_Electrodes/ROIs_active_elecs.mat

%% - LOAD DATA - %%

niterations = 1000; % # iterations per bootstrap

% initialize structure for aperiodic features
ap_slope  = struct();
ap_offset = struct();

for Isub = 1:numel(subjects)
    
    % get subject ID
    subID = subjects{Isub};

    fprintf('%%-----------------------------------------------------------\n');
    fprintf('Computation for %s\n', subID);
    fprintf('%%-----------------------------------------------------------\n');

    % load PSD
    load(fullfile(path_in, subID, 'Data', 'PSD', [subID '_PSD.mat']));
      
    % get the trial information
    load(fullfile(path_in, subID, 'Data', [subID '_Trial_Info_Clean.mat']))
    
    % get index for stop trials or trials that could not be defined
    idx = find(isnan(trl_mod_table.TrlType) | trl_mod_table.TrlType == 2);
    
    trl_mod_table(idx,:)= [];
    
    % get trial indices
    alltrlidx = {find(trl_mod_table.LikehoodStop == 0), ...
                 find(trl_mod_table.LikehoodStop == 25), ...
                 find(trl_mod_table.LikehoodStop == 75)};
     
    % get minimum # trials across conditions
    mintrl = min([numel(alltrlidx{1}) numel(alltrlidx{2}) ...
        numel(alltrlidx{3})]);
         
    % --------------------------------------------------------------------
    % Compute Aperiodic Features
    % --------------------------------------------------------------------

    % --------------------------
    % Frontal Cortex
    % --------------------------
    
    if isfield(ROIs(Isub).hll, 'frontal') && isfield(ROIs(Isub).hll.frontal, 'active')

        % get frontal electrodes that are task-informative
        elec_frontal = ROIs(Isub).hll.frontal.active.elecidx;
        
        % select data of interest
        cfg = [];
        cfg.channel = powspec.hll.frac.label(elec_frontal);
        
        data        = ft_selectdata(cfg, powspec.hll.frac); clear elec_frontal
        
        % --------------------------------------
        % compute aperiodic fit in log-log space
        % --------------------------------------
            
        % preallocate memory
        tmpSlope  = NaN(numel(data.label), mintrl);
        tmpOffset = NaN(numel(data.label), mintrl);
       
        for Ichan = 1:numel(data.label)
           
            for Itrial = 1:size(data.powspctrm,1)

                % do the fitting
                tmpData = squeeze(data.powspctrm(Itrial,Ichan,:))';

                % get frequency vector
                tmpFreq = data.freq;

                % only fit between 30-45 Hz
                fitrange = [30 45];
                idxfreq = nearest(tmpFreq,fitrange(1)):nearest(tmpFreq,fitrange(2));

                % compute linear fit
                b = polyfit(log10(tmpFreq(idxfreq)), log10(tmpData(idxfreq)), 1);

                % store regression parameters

                tmpSlope(Ichan,Itrial)  = b(1);
                tmpOffset(Ichan,Itrial) = b(2);

            end %Itrial

        end %Ichan
          
        % ----------------
        % Store ID Frontal
        % ----------------
        ap_slope.frontal.ID(Isub)    = str2double(subID(4:end));
        ap_offset.frontal.ID(Isub)   = str2double(subID(4:end));     
        
        % ----------------
        % Store Data 0%
        % ----------------
        
        % get a bootstrap distributon for the slope
        bootslope = NaN(niterations,size(tmpSlope,1),mintrl);
        for Ichan = 1:size(tmpSlope,1)
            for Iboot = 1:niterations
                rndidx = datasample(alltrlidx{1}, mintrl);
                bootslope(Iboot,Ichan,:) = tmpSlope(Ichan,rndidx);
            end
        end
        
        % get a bootstrap distributon for the offset
        bootoffset = NaN(niterations,size(tmpOffset,1),mintrl);
        for Ichan = 1:size(tmpOffset,1)
            for Iboot = 1:niterations
                rndidx = datasample(alltrlidx{1}, mintrl);
                bootoffset(Iboot,Ichan,:) = tmpOffset(Ichan,rndidx);
            end
        end
        
        ap_slope.frontal.lh0{Isub}   = mean(bootslope, 'all', 'omitnan'); % average across channels
        ap_offset.frontal.lh0{Isub}  = mean(bootoffset, 'all', 'omitnan');

        clear bootslope bootoffset
        
        % ----------------
        % Store Data 25%
        % ----------------
        
        % get a bootstrap distributon for the slope
        bootslope = NaN(niterations,size(tmpSlope,1),mintrl);
        for Ichan = 1:size(tmpSlope,1)
            for Iboot = 1:niterations
                rndidx = datasample(alltrlidx{2}, mintrl);
                bootslope(Iboot,Ichan,:) = tmpSlope(Ichan,rndidx);
            end
        end
        
        % get a bootstrap distributon for the offset
        bootoffset = NaN(niterations,size(tmpOffset,1),mintrl);
        for Ichan = 1:size(tmpOffset,1)
            for Iboot = 1:niterations
                rndidx = datasample(alltrlidx{2}, mintrl);
                bootoffset(Iboot,Ichan,:) = tmpOffset(Ichan,rndidx);
            end
        end
        
        ap_slope.frontal.lh25{Isub}   = mean(bootslope, 'all', 'omitnan'); % average across channels
        ap_offset.frontal.lh25{Isub}  = mean(bootoffset, 'all', 'omitnan');

        clear bootslope bootoffset
        
        % ----------------
        % Store Data 75%
        % ----------------
        ap_slope.frontal.lh75{Isub}  = mean(tmpSlope(:,alltrlidx{3})); % average across channels
        ap_offset.frontal.lh75{Isub} = mean(tmpOffset(:,alltrlidx{3}));
        
        clear data
        
    end
        
    % --------------------------
    % Motor Cortex
    % --------------------------    
    
    if isfield(ROIs(Isub).hll, 'motor') && isfield(ROIs(Isub).hll.motor, 'active')

        % get motor electrodes that are task-informative
        elec_motor = ROIs(Isub).hll.motor.active.elecidx;
        
        % select data of interest
        cfg = [];
        cfg.channel = powspec.hll.frac.label(elec_motor);
        
        data        = ft_selectdata(cfg, powspec.hll.frac); clear elec_motor
        
        % --------------------------------------
        % compute aperiodic fit in log-log space
        % --------------------------------------
            
        % preallocate memory
        tmpSlope  = NaN(numel(data.label), mintrl);
        tmpOffset = NaN(numel(data.label), mintrl);
       
        for Ichan = 1:numel(data.label)
           
            for Itrial = 1:size(data.powspctrm,1)

                % do the fitting
                tmpData = squeeze(data.powspctrm(Itrial,Ichan,:))';

                % get frequency vector
                tmpFreq = data.freq;

                % only fit between 30-45 Hz
                fitrange = [30 45];
                idxfreq = nearest(tmpFreq,fitrange(1)):nearest(tmpFreq,fitrange(2));

                % compute linear fit
                b = polyfit(log10(tmpFreq(idxfreq)), log10(tmpData(idxfreq)), 1);

                % store regression parameters

                tmpSlope(Ichan,Itrial)  = b(1);
                tmpOffset(Ichan,Itrial) = b(2);

            end %Itrial

        end %Ichan
          
        % ----------------
        % Store ID motor
        % ----------------
        ap_slope.motor.ID(Isub)    = str2double(subID(4:end));
        ap_offset.motor.ID(Isub)   = str2double(subID(4:end));     
        
        % ----------------
        % Store Data 0%
        % ----------------
        
        % get a bootstrap distributon for the slope
        bootslope = NaN(niterations,size(tmpSlope,1),mintrl);
        for Ichan = 1:size(tmpSlope,1)
            for Iboot = 1:niterations
                rndidx = datasample(alltrlidx{1}, mintrl);
                bootslope(Iboot,Ichan,:) = tmpSlope(Ichan,rndidx);
            end
        end
        
        % get a bootstrap distributon for the offset
        bootoffset = NaN(niterations,size(tmpOffset,1),mintrl);
        for Ichan = 1:size(tmpOffset,1)
            for Iboot = 1:niterations
                rndidx = datasample(alltrlidx{1}, mintrl);
                bootoffset(Iboot,Ichan,:) = tmpOffset(Ichan,rndidx);
            end
        end
        
        ap_slope.motor.lh0{Isub}   = mean(bootslope, 'all', 'omitnan'); % average across channels
        ap_offset.motor.lh0{Isub}  = mean(bootoffset, 'all', 'omitnan');

        clear bootslope bootoffset
        
        % ----------------
        % Store Data 25%
        % ----------------
        
        % get a bootstrap distributon for the slope
        bootslope = NaN(niterations,size(tmpSlope,1),mintrl);
        for Ichan = 1:size(tmpSlope,1)
            for Iboot = 1:niterations
                rndidx = datasample(alltrlidx{2}, mintrl);
                bootslope(Iboot,Ichan,:) = tmpSlope(Ichan,rndidx);
            end
        end
        
        % get a bootstrap distributon for the offset
        bootoffset = NaN(niterations,size(tmpOffset,1),mintrl);
        for Ichan = 1:size(tmpOffset,1)
            for Iboot = 1:niterations
                rndidx = datasample(alltrlidx{2}, mintrl);
                bootoffset(Iboot,Ichan,:) = tmpOffset(Ichan,rndidx);
            end
        end
        
        ap_slope.motor.lh25{Isub}   = mean(bootslope, 'all', 'omitnan'); % average across channels
        ap_offset.motor.lh25{Isub}  = mean(bootoffset, 'all', 'omitnan');

        clear bootslope bootoffset
        
        % ----------------
        % Store Data 75%
        % ----------------
        ap_slope.motor.lh75{Isub}  = mean(tmpSlope(:,alltrlidx{3})); % average across channels
        ap_offset.motor.lh75{Isub} = mean(tmpOffset(:,alltrlidx{3}));
        
        clear data
        
    end
    
end %Isub
        
%% - REMOVE EMPTY CELLS - %%

% -----------------
% Frontal Cortex
% -----------------

ap_slope.frontal.ID    = ap_slope.frontal.ID(find(ap_slope.frontal.ID ~= 0));
ap_offset.frontal.ID   = ap_offset.frontal.ID(find(ap_offset.frontal.ID ~= 0));

ap_slope.frontal.lh0   = cellfun(@mean, ap_slope.frontal.lh0(~cellfun(@isempty, ap_slope.frontal.lh0)));
ap_slope.frontal.lh25  = cellfun(@mean, ap_slope.frontal.lh25(~cellfun(@isempty, ap_slope.frontal.lh25)));
ap_slope.frontal.lh75  = cellfun(@mean, ap_slope.frontal.lh75(~cellfun(@isempty, ap_slope.frontal.lh75)));

ap_offset.frontal.lh0  = cellfun(@mean, ap_offset.frontal.lh0(~cellfun(@isempty, ap_offset.frontal.lh0)));
ap_offset.frontal.lh25 = cellfun(@mean, ap_offset.frontal.lh25(~cellfun(@isempty, ap_offset.frontal.lh25)));
ap_offset.frontal.lh75 = cellfun(@mean, ap_offset.frontal.lh75(~cellfun(@isempty, ap_offset.frontal.lh75)));

% -----------------
% Motor Cortex
% -----------------

ap_slope.motor.ID      = ap_slope.motor.ID(find(ap_slope.motor.ID ~= 0));
ap_offset.motor.ID     = ap_offset.motor.ID(find(ap_offset.motor.ID ~= 0));

ap_slope.motor.lh0     = cellfun(@mean, ap_slope.motor.lh0(~cellfun(@isempty, ap_slope.motor.lh0)));
ap_slope.motor.lh25    = cellfun(@mean, ap_slope.motor.lh25(~cellfun(@isempty, ap_slope.motor.lh25)));
ap_slope.motor.lh75    = cellfun(@mean, ap_slope.motor.lh75(~cellfun(@isempty, ap_slope.motor.lh75)));

ap_offset.motor.lh0    = cellfun(@mean, ap_offset.motor.lh0(~cellfun(@isempty, ap_offset.motor.lh0)));
ap_offset.motor.lh25   = cellfun(@mean, ap_offset.motor.lh25(~cellfun(@isempty, ap_offset.motor.lh25)));
ap_offset.motor.lh75   = cellfun(@mean, ap_offset.motor.lh75(~cellfun(@isempty, ap_offset.motor.lh75)));

%% Plot Data

selroi = input('Which regions should be plotted? Frontal = 1, Motor = 2: ');

if selroi == 1
    
    data2plot = {ap_slope.frontal.lh0', ap_slope.frontal.lh25', ap_slope.frontal.lh75'}';
    
elseif selroi == 2
    
    data2plot = {ap_slope.motor.lh0', ap_slope.motor.lh25', ap_slope.motor.lh75'}';
        
end
    

cb =  {'#008450', '#EFB700', '#B81D13'};

figure(2); 
h = rm_raincloud(data2plot, [1 0 0], 0, 'ks');

% admin
for Icond = 1:3
    
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
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('LL Stop'); % x and y axis are flipped in the function
xlabel('1/f Slope');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off

set(gcf, 'Position', [544 272 186 189]);

%% - SAVE FIGURE - %%

cd(path_out);
if selroi == 1
    print(gcf,'fig_3_frontal_aperiodic_slope.pdf','-dpdf','-r400');
else
    print(gcf,'fig_3_motor_aperiodic_slope.pdf','-dpdf','-r400');
end

%% Write into Excel

% get only subject number without "OSL"
xx  = regexp(subjects,'\d+(\.)?(\d+)?','match');
subjectnumber =str2double([xx{:}])';      

statsdata = [];

% include frontal data
statsdata(:,4) = vertcat(ap_slope.frontal.lh0', ...
                         ap_slope.frontal.lh25', ...
                         ap_slope.frontal.lh75', ...
                         ...
                         ap_slope.motor.lh0', ...
                         ap_slope.motor.lh25', ...
                         ap_slope.motor.lh75');


statsdata(1:end,1) = vertcat(repmat(ap_slope.frontal.ID, 1, 3)', ...
    repmat(ap_slope.motor.ID, 1, 3)');

% 1 = Frontal, 2 = Motor
statsdata(1:end,2) = vertcat(ones(numel(ap_slope.frontal.lh0)*3, 1), ...
                                   repmat(3, numel(ap_slope.motor.lh0)*3, 1));

% document factor conditions (0, 25, 75% likelihood of Stop) 
statsdata(1:end,3) = vertcat(...
    [zeros(length(ap_slope.frontal.ID),1); ...
    repmat(25, length(ap_slope.frontal.ID), 1); ...
    repmat(75, length(ap_slope.frontal.ID), 1)], ...
    ...
    [zeros(length(ap_slope.motor.ID),1); ...
    repmat(25, length(ap_slope.motor.ID), 1); ...
    repmat(75, length(ap_slope.motor.ID), 1)]);

%% Write Data into Table

% go into results folder
mkdir(path_source, 'Fig_3e_f_left');
cd(fullfile(path_source, 'Fig_3e_f_left'));

statstable  = array2table(statsdata, 'VariableNames', {'ID', 'ROI', 'Prediction', 'Slope'});

writetable(statstable,filename, 'Sheet', 'Slope', 'WriteRowNames', true);


