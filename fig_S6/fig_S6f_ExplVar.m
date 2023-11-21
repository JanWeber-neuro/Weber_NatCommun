%% Figure S6f
% this script rungs additional analysis required for the plots of figure
% S6f

%%

clear 
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_5/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/

%% Load Electrode Information

load /Volumes/IMPECOG/iEEG_Data_JW/Results/GroupLevel/Raw_Data/Electrode_Information_ROI/elec_info_ROI.mat

%%

% preallocate memory for variance explained by first 3 PCs relative to 100%
ExplVar_frontal = NaN(numel(subjects),1);
ExplVar_motor   = NaN(numel(subjects),1);

% preallocate memory for storing all PCs
PCs_frontal     = cell(numel(subjects),1);
PCs_motor       = cell(numel(subjects),1);

for Isub = 1:numel(subjects)
    
    % get subject ID
    subID = subjects{Isub};
    
    % check if folder exists and if there is data
    if isfolder(fullfile(path_in, subID, 'Data', 'State_Space'))

        load(fullfile(path_in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']));

        % -----------------------------------------
        % Frontal Cortex
        % -----------------------------------------

        if isfield(all_data, 'frontal')

            x = cumsum(all_data.frontal.pcainfo.explained(1:3)); % get cumulative sum for first 3 components
            ExplVar_frontal(Isub) = x(end) / 100; clear x 
            PCs_frontal{Isub}     = all_data.frontal.pcainfo.explained;
        
        end
        
        % -----------------------------------------
        % Motor Cortex
        % -----------------------------------------

        if isfield(all_data, 'motor')

            x = cumsum(all_data.motor.pcainfo.explained(1:3)); % get cumulative sum for first 3 components
            ExplVar_motor(Isub) = x(end) / 100; clear x 
            PCs_motor{Isub}     = all_data.motor.pcainfo.explained;
        
        end

    end
    
end %Isub

%% Plot Relative Percent of First 3 Components

figure(1);
jw_scatterplot({ExplVar_frontal, ExplVar_motor}, 1, 0, {'r', 'b'});

set(gca, 'linewidth', 3, 'fontsize', 13);
ylabel({'rel. contribition', 'first 3 PCs'});
xticks(1:2);
xticklabels({'Frontal', 'Motor'});
box off

set(gcf, 'Position', [515 450 243 173]);

%% Plot the Cumulative Sum of first 5 Components

nPCs = 5;

% Frontal
GA_frontal = NaN(numel(PCs_frontal),nPCs);

for Isub = 1:size(GA_frontal,1)
    
    if ~isempty(PCs_frontal{Isub}) && numel(PCs_frontal{Isub}) > nPCs-1
        
        GA_frontal(Isub,:) = cumsum(PCs_frontal{Isub}(1:nPCs)');
        
    end
    
end
   
% Motor
GA_motor   = NaN(numel(PCs_motor),nPCs);

for Isub = 1:size(GA_motor,1)
    
    if ~isempty(PCs_motor{Isub}) &&  numel(PCs_motor{Isub}) > nPCs-1
        
        GA_motor(Isub,:) = cumsum(PCs_motor{Isub}(1:nPCs)');
        
    end
    
end
       
figure(2);
a = shadedErrorBar(1:nPCs, nanmean(GA_frontal), nanstd(GA_frontal) / sqrt(size(GA_frontal,1) - sum(isnan(GA_frontal(:,1)))));
a.mainLine.Color = 'r'; a.mainLine.LineWidth = 2; 
a.patch.FaceColor = 'r'; a.patch.EdgeColor = 'r'; a.patch.FaceAlpha = 0.5;
a.mainLine.Marker = 'o'; a.mainLine.MarkerFaceColor = 'r'; a.mainLine.MarkerEdgeColor = 'r';

hold on

a = shadedErrorBar(1:nPCs, nanmean(GA_motor), nanstd(GA_motor) / sqrt(size(GA_motor,1) - sum(isnan(GA_motor(:,1)))));
a.mainLine.Color = 'b'; a.mainLine.LineWidth = 2;
a.patch.FaceColor = 'b'; a.patch.EdgeColor = 'b'; a.patch.FaceAlpha = 0.5;
a.mainLine.Marker = 'o'; a.mainLine.MarkerFaceColor = 'b'; a.mainLine.MarkerEdgeColor = 'b';

set(gca, 'linewidth', 3, 'fontsize', 13);
ylabel('Cumulative Sum');
xlim([1 nPCs]);
xlabel('PCs');
xticks(1:nPCs);
box off

set(gcf, 'Position', [515 450 243 173]);

%% Fit the Slope of the Explained Variance Trend 

% --------------------------

% a steeper slope indicates high variance explained by first components
% whereas a flatter slope indicates that later components still add fair
% amount of variance

% higher offset indicates that first PC explains more variance

% --------------------------

% Frontal
slope_frontal  = NaN(numel(PCs_frontal),1);
offset_frontal = NaN(numel(PCs_frontal),1);

for Isub = 1:numel(PCs_frontal)
    
    if ~isempty(PCs_frontal{Isub})
        % number of PCs
        nPCs = 1:length(PCs_frontal{Isub});

        % fit the data
        b = polyfit(nPCs, PCs_frontal{Isub}, 1);

        slope_frontal(Isub)  = b(1); % store slope
        offset_frontal(Isub) = b(2); % store offset
    end
    
end

% Motor
slope_motor  = NaN(numel(PCs_motor),1);
offset_motor = NaN(numel(PCs_motor),1);

for Isub = 1:numel(PCs_motor)
    
    if ~isempty(PCs_motor{Isub})
        nPCs = 1:length(PCs_motor{Isub});

        b = polyfit(nPCs, PCs_motor{Isub}, 1);

        slope_motor(Isub)  = b(1);
        offset_motor(Isub) = b(2);
    end

end

% get grand average fit

datafit_frontal = NaN(numel(PCs_frontal), length(nPCs));
datafit_motor   = NaN(numel(PCs_motor), length(nPCs));

for Isub = 1:size(datafit_frontal,1)
    datafit_frontal(Isub,:) = offset_frontal(Isub) + slope_frontal(Isub)*nPCs;
end
    
for Isub = 1:size(datafit_motor,1)
    datafit_motor(Isub,:) = offset_motor(Isub) + slope_motor(Isub)*nPCs;
end

% plot data
figure(3);

plot(nPCs, mean(datafit_frontal, 'omitnan'), 'r', 'LineWidth', 4)
hold on
plot(nPCs, mean(datafit_motor, 'omitnan'), 'b', 'LineWidth', 4)

set(gca, 'linewidth', 3, 'fontsize', 13);
xlim([1 length(datafit_frontal)]);
ylabel('Model Fit');
xlabel('PCs');
box off

set(gcf, 'Position', [515 450 243 173]);

%% Save source data

% go into results folder
mkdir(path_source, 'Fig_5f');
cd(fullfile(path_source, 'Fig_5f'));

% get only subject number without "OSL"
xx  = regexp(subjects,'\d+(\.)?(\d+)?','match');
subjectnumber =str2double([xx{:}])';      

statsdata = [];

% include relative explained variance of cumulative sum of first 3 components
statsdata(:,3) = vertcat(ExplVar_frontal, ...
                         ExplVar_motor);
          
% include fitted slope 
statsdata(:,4) = vertcat(slope_frontal, ...
                         slope_motor);
   
                     
% include fitted offset 
statsdata(:,5) = vertcat(offset_frontal, ...
                         offset_motor);

statsdata(1:end,1) = repmat(subjectnumber, 2, 1);

% 1 = Frontal, 2 = Motor, 3 = MTL
statsdata(1:end,2) = vertcat(ones(size(statsdata, 1) / 2, 1), ...
    repmat(3, size(statsdata, 1) / 2, 1));
                               
% save data
statstable = array2table(statsdata, 'VariableNames', {'ID', 'ROI', 'rel_explvar', 'slope', 'offset'});

% write data into excelsheet
filename = 'PC_ExplVar.xlsx';

writetable(statstable,filename);

%% Save source data

% go into results folder
mkdir(path_source, 'Fig_5f');
cd(fullfile(path_source, 'Fig_5f'));

tmp = array2table(GA_frontal);
filename = 'PC_ExplVar_pfc_PC1to5.xlsx';
writetable(tmp,filename);

clear tmp

tmp = array2table(GA_motor);
filename = 'PC_ExplVar_motor_PC1to5.xlsx';
writetable(tmp,filename);


