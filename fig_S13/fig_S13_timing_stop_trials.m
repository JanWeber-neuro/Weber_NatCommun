%% Figure S13

%%

clear
close all;

%% Set path

path_main = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/';
cd(path_main);

path_out = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm/Figures/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/raincloud_plots/

%% Subjects

subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', ...
    'OSL27', 'OSL28', 'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', ...
    'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Loop over subjects and extract reaction time and accuracy

% initilize outcome variables reaction time (RT) and accuracy (ACC)
SSD = [];

for Isub = 1:numel(subjects)
    
    subID = subjects{Isub};

    load(fullfile(path_main, subID, 'Data', [subID '_Trial_Info.mat']));
    
    SSD = cat(1, SSD, trl_mod_table.SSD);

%     SSD(find(SSD < 0.1)) = NaN;

end

SSD(logical(isnan(SSD))) = [];
SSD(SSD > 1) = [];
SSD(SSD < 0.4) = [];

SSD = SSD - 0.56;

a = histogram(SSD, 'Normalization', 'probability');
a.FaceColor = 'k';
a.FaceAlpha = 1;

hold on

[f,xi] = ksdensity(SSD);

plot(xi, f*0.012, 'r', 'LineWidth', 2)

ylabel('probability');
xlabel('Time [s]');
box off

set(gca, 'fontsize', 13);

% go to output folder and save figure
cd(path_out);

print(gcf,'timing_stop_trials.pdf','-dpdf','-r400');   

%% SAVE SOURCE DATA

% go into results folder
mkdir(path_source, 'Fig_S13');
cd(fullfile(path_source, 'Fig_S13'));

table_data = array2table(SSD);

% write data into excelsheet
filename = 'fig_S13_SSD.xlsx';
writetable(table_data,filename, 'Sheet', 1);
