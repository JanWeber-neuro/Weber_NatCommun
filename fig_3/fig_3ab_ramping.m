%% Figure 3a & b (right panel)
% this script rungs additional analysis required for the plots of figure
% 3a & 3b.

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

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

%% - PARAMS FOR BOOTSTRAPPING - %%

nboot = 1000;

%% - FIT THE DATA - %%
    
rampest = struct();

for Isub = 1:numel(subjects)
        
    warning off
    
    % get indicated subject ID
    subID = subjects{Isub};
    
    % load HFA data
    load(fullfile(path_in, subID, 'Data', 'HFA', 'PCA', 'PCA_ROI_Data', [subID '_activedata_ROI_pca.mat']))
    
    fprintf('ID %d/%d\n', Isub, numel(subjects));
    
    % ------------------------------------
    % Extract Trial Information
    % ------------------------------------

    % extract trial information to get index of LHs
    if ~isempty(roidata.hll.frontal.inactive)
        trlinfo = roidata.hll.frontal.inactive.trialinfo;
    elseif ~isempty(roidata.hll.motor.inactive)
        trlinfo = roidata.hll.motor.inactive.trialinfo;
    end
    
    % get indices
    idx_0  = find(trlinfo(:,13) == 0);  % 0% LH
    idx_25 = find(trlinfo(:,13) == 25); % 25% LH
    idx_75 = find(trlinfo(:,13) == 75); % 75% LH
    
    % combine
    alltrlidx = {idx_0, idx_25, idx_75};
    
    % minimum # trials
    min_numtrl = min([length(idx_0) length(idx_25) length(idx_75)]);
    
    % merge data
    MergedData = {roidata.hll};
        
    % Loop over epochs
    for Iepoch = 1:numel(MergedData)
        
        % Loop over conditions
        
        for Icond = 1:numel(alltrlidx)
        
        % #####################
        % Frontal Cortex
        % #####################

            if ~isempty(MergedData{Iepoch}.frontal.active)

                % get time
                time = MergedData{Iepoch}.frontal.active.time;
                timetoevent = nearest(time, -0.5):nearest(time, 0); 

                % current data
                tmpdata = MergedData{Iepoch}.frontal.active.trial(:,:,timetoevent);
                
                % preallocate memory for slope
                slope = NaN(nboot, size(tmpdata,2), min_numtrl);
                
                for Iboot = 1:nboot

                    % get random trial index
                    rndidx = randsample(numel(alltrlidx{Icond}), min_numtrl);

                    for Ichan = 1:size(tmpdata,2) % loop over channels

                        for Itrial = 1:min_numtrl % loop over trials

                            % current trial for slope estimation and smooth data using linear regression
                            tmptrial = squeeze(tmpdata(alltrlidx{Icond}(rndidx(Itrial)), Ichan, :))';

                            % get beta estimates (intercept & slope)
                            b_est = robustfit(time(timetoevent), tmptrial);

                            slope(Iboot,Ichan,Itrial) = b_est(2);

                        end

                    end

                    % bring slope estimates into structure
                    rampest.hll.frontal{Icond}{Isub} = mean(slope, 'all', 'omitnan');
                    
                end

            else

                % if no data obtained
                rampest.hll.frontal{Icond}{Isub} = NaN;

            end
            
            % #####################
            % Motor Cortex
            % #####################

            if ~isempty(MergedData{Iepoch}.motor.active)

                % get time
                time = MergedData{Iepoch}.motor.active.time;
                timetoevent = nearest(time, -0.5):nearest(time, 0); 

                % current data
                tmpdata = MergedData{Iepoch}.motor.active.trial(:,:,timetoevent);
                
                % preallocate memory for slope
                slope = NaN(nboot, size(tmpdata,2), min_numtrl);
                
                for Iboot = 1:nboot

                    % get random trial index
                    rndidx = randsample(numel(alltrlidx{Icond}), min_numtrl);

                    for Ichan = 1:size(tmpdata,2) % loop over channels

                        for Itrial = 1:min_numtrl % loop over trials

                            % current trial for slope estimation and smooth data using linear regression
                            tmptrial = squeeze(tmpdata(alltrlidx{Icond}(rndidx(Itrial)), Ichan, :))';

                            % get beta estimates (intercept & slope)
                            b_est = robustfit(time(timetoevent), tmptrial);

                            slope(Iboot,Ichan,Itrial) = b_est(2);

                        end

                    end

                    % bring slope estimates into structure
                    rampest.hll.motor{Icond}{Isub} = mean(slope, 'all', 'omitnan');

                end

            else
                
                % if no data obtained
                rampest.hll.motor{Icond}{Isub} = NaN;

            end    

        end % Ipredict
        
    end % Iepoch
    
end %Isub
           
%% - EXTRACT AVERAGE PER FRONTAL CORTEX - %%

ga_frontal_lh0  = struct();
ga_frontal_lh25 = struct();
ga_frontal_lh75 = struct();

%% - FRONTAL HLL - %%  

% Extract slope
ga_frontal_lh0.hll.slope   = cellfun(@mean, rampest.hll.frontal{1});
ga_frontal_lh25.hll.slope  = cellfun(@mean, rampest.hll.frontal{2});
ga_frontal_lh75.hll.slope  = cellfun(@mean, rampest.hll.frontal{3});

%% - EXTRACT AVERAGE PER MOTOR CORTEX - %%

ga_motor_lh0    = struct();
ga_motor_lh25   = struct();
ga_motor_lh75   = struct();

%% - MOTOR HLL - %%  

ga_motor_lh0.hll.slope     = cellfun(@mean, rampest.hll.motor{1});
ga_motor_lh25.hll.slope    = cellfun(@mean, rampest.hll.motor{2});
ga_motor_lh75.hll.slope    = cellfun(@mean, rampest.hll.motor{3});

%% Plot Data

cd(path_out);

MergeData_hll = {ga_frontal_lh0.hll.slope, ga_frontal_lh25.hll.slope, ga_frontal_lh75.hll.slope; ...
    ga_motor_lh0.hll.slope, ga_motor_lh25.hll.slope, ga_motor_lh75.hll.slope};

cb =  {'#008450', '#EFB700', '#B81D13'};

for Iroi = 1:size(MergeData_hll,1)
    
    figure(Iroi);
    
    h = rm_raincloud({MergeData_hll{Iroi,1}, MergeData_hll{Iroi,2}, MergeData_hll{Iroi,3}}', [1 0 0], 0, 'ks', 1.5);
    
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
    xlabel('Ramping Slope')
    set(gca, 'linewidth', 3, 'fontsize', 13);
    box off
    set(gcf, 'Position', [527 269 180 191]);

    if Iroi == 1
        print(gcf,'fig_3_frontal_slope_grouplevel_HLL.pdf','-dpdf','-r400');    
    else
        print(gcf,'fig_3_motor_slope_grouplevel_HLL.pdf','-dpdf','-r400');    
    end

end

%% Extract Data for Stats in Python (Mean)

% go into results folder
mkdir(path_source, 'Fig_3a_b_right');
cd(fullfile(path_source, 'Fig_3a_b_right'));

% get only subject number without "OSL"
xx  = regexp(subjects,'\d+(\.)?(\d+)?','match');
subjectnumber =str2double([xx{:}])';      

statsdata = struct();

%% Slope for Until HLL

% bring data together
statsdata.hll(:,4) = horzcat(ga_frontal_lh0.hll.slope, ...
                             ga_frontal_lh25.hll.slope, ...
                             ga_frontal_lh75.hll.slope, ...
                             ga_motor_lh0.hll.slope, ...
                             ga_motor_lh25.hll.slope, ...
                             ga_motor_lh75.hll.slope)';

statsdata.hll(1:end,1) = repmat(subjectnumber, 6, 1);

% 1 = Frontal, 2 = Motor
statsdata.hll(1:end,2) = vertcat(ones(size(statsdata.hll, 1) / 2, 1), ...
                                   repmat(3, size(statsdata.hll, 1) / 2, 1));

% document factor conditions (0, 25, 75% likelihood of Stop)                               
statsdata.hll(1:end,3) = repmat(vertcat(zeros(size(statsdata.hll, 1) / 6, 1), ...
                                   repmat(25, size(statsdata.hll, 1) / 6, 1), ... 
                                   repmat(75, size(statsdata.hll, 1) / 6, 1)), 2 , 1);

%% Write Data into Table

statstable_hll  = array2table(statsdata.hll, 'VariableNames', {'ID', 'ROI', 'Prediction', 'Slope'});

% write data into excelsheet
filename = 'fig_3a_b_right.xlsx';

writetable(statstable_hll,filename, 'Sheet', 'HLL', 'WriteRowNames', true);



        