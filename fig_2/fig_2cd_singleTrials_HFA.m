%% Figure 2c & d

%%

clear
ft_defaults

%%
rng('default');

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/ % add functions

%% - LOAD EXAMPLE DATA - %%
       
path_data = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects/OSL24/Data/HFA/PCA/PCA_ROI_Data/';
load(fullfile(path_data, 'OSL24_activedata_ROI_pca.mat'));

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_2/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

%% - MERGE DATA - %%

tmpData    = {roidata.hll};

%% - GET CONDITION INDEX - %%

idx_l0  = find(roidata.hll.frontal.inactive.trialinfo(:,13) == 0);
idx_l25 = find(roidata.hll.frontal.inactive.trialinfo(:,13) == 25);
idx_l75 = find(roidata.hll.frontal.inactive.trialinfo(:,13) == 75);

alltrialidx = {idx_l0, idx_l25, idx_l75}; % all trials concatenated

%% - PLOT FRONTAL SINGLE TRIALS - %%

cd(path_out);

% get time for epoch
time  = roidata.hll.frontal.active.time;                

% get data with LH 0
data_lh0  = roidata.hll.frontal.active.trial(alltrialidx{1},:,:);

% get data with LH 25
data_lh25 = roidata.hll.frontal.active.trial(alltrialidx{2},:,:);

% get data with LH 75
data_lh75 = roidata.hll.frontal.active.trial(alltrialidx{3},:,:);

figure;

a = plot(time, squeeze(data_lh0(15,5,:)), 'LineWidth', 3);
a.Color = '#008450'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

print(gcf,'fig_2_single_trl_frontal_0.pdf','-dpdf','-r400');

figure;

a = plot(time, squeeze(data_lh25(12,5,:)), 'LineWidth', 3);
a.Color = '#EFB700'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

print(gcf,'fig_2_single_trl_frontal_25.pdf','-dpdf','-r400');

figure;

a = plot(time, squeeze(data_lh75(7,5,:)), 'LineWidth', 3);
a.Color = '#B81D13'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

print(gcf,'fig_2_single_trl_frontal_75.pdf','-dpdf','-r400');

% save source data
tmp = [squeeze(data_lh0(15,5,:)) ...
    squeeze(data_lh25(12,5,:)) ...
    squeeze(data_lh75(7,5,:))];

% go into results folder
mkdir(path_source, 'fig_2c');
cd(fullfile(path_source, 'fig_2c'));

tmptable = array2table(tmp);

% write data into excelsheet
filename = 'fig_2c_pfc_single_trials.xlsx';

writetable(tmptable,filename);

%% - PLOT MOTOR SINGLE TRIALS - %%

% get time for epoch
time  = roidata.hll.motor.active.time;                

% get data with LH 0
data_lh0  = roidata.hll.motor.active.trial(alltrialidx{1},:,:);

% get data with LH 25
data_lh25 = roidata.hll.motor.active.trial(alltrialidx{2},:,:);

% get data with LH 75
data_lh75 = roidata.hll.motor.active.trial(alltrialidx{3},:,:);

figure;

a = plot(time, squeeze(data_lh0(2,6,:)), 'LineWidth', 3);
a.Color = '#008450'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

print(gcf,'fig_2_single_trl_motor_0.pdf','-dpdf','-r400');

figure;

a = plot(time, squeeze(data_lh25(22,5,:)), 'LineWidth', 3);
a.Color = '#EFB700'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

print(gcf,'fig_2_single_trl_motor_25.pdf','-dpdf','-r400');

figure;

a = plot(time, squeeze(data_lh75(13,5,:)), 'LineWidth', 3);
a.Color = '#B81D13'; ylim([-8 25]);
% make plot look nice
xticks([-0.4 -0.2 0 0.2]);
xticklabels({'-0.4', '-0.2', '0', '0.2'});
set(gca, 'LineWidth', 3, 'FontSize', 13); box off
xlim([time(1) time(end)]);
yline(0, 'LineStyle', '--'); yline(10, 'LineStyle', '--'); yline(20, 'LineStyle', '--');

% remove ticks
h = gca; 
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];

set(gcf, 'Position', [574 436 215 136]);

print(gcf,'fig_2_single_trl_motor_75.pdf','-dpdf','-r400');

% save source data
tmp = [squeeze(data_lh0(2,6,:)) ...
    squeeze(data_lh25(22,5,:)) ...
    squeeze(data_lh75(13,5,:))];

% go into results folder
mkdir(path_source, 'fig_2c');
cd(fullfile(path_source, 'fig_2c'));

tmptable = array2table(tmp);

% write data into excelsheet
filename = 'fig_2c_m1_single_trials.xlsx';

writetable(tmptable,filename);
