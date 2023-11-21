%% - PLOT FIGURES - %%
% This script plots all figures found in the manuscript "Ramping dynamics and theta oscillations 
% reflect dissociable signatures during rule-guided human behavior". 
% All necessary data can be found within the subfolder "source_data".

% Data exceptions
%
%   Fig. 1a: This figure is an illustration and contains no data.
%
%   Fig. 1b: This figure is a schematic illustration and contains no data.
%
%   Fig. 1e: This figure contains data about electrode implanation
%            trajectories that would compromise the anonymity of the
%            patient.
%
%   Fig. 2a & b (right panel): This figure contains data about electrode implanation
%            trajectories that would compromise the anonymity of the
%            patient.
% 
%   Fig. 2e & f (left panel): This figure contains 3D data.
%
%   Fig. 3c & d (left panel): This figure contains 3D data.
%
%   Fig. 4a/d: These figures only illustrate the analytical approach.
%
%   Fig. 4i: This figure contains data about electrode implanation
%            trajectories that would compromise the anonymity of the
%            patient.
%
%   Fig. S2d: This figure depicts a schematic illustration and contains no
%             data
%
%   Fig. S4: This figure only shows single electrode data to illustrate the
%            "blindness" of the ANOVA.
%
%   Fig. S6h: This figure contains 3D data.
%
% Following plotting, figures were exported to Adobe Illustrator for
% resizing/positioning to tidy the complete figures.
% Consequently, some dimensions/colours/fonts/labels may vary. 
%
% Dependencies: 
%   raincloud plot
%   Fieldtrip
%   cbrewer & viridis colormap
%   shadedErrorBar

clear
close all;
clc;

%% - PATH SETTINGS - %%

% define data directory
dir_repo = pwd;
addpath(genpath(fullfile(dir_repo, 'subfunctions')));

% setup path dependencies
path = setup_path;

%% - FIGURE 1 - %%

addpath(genpath(fullfile(dir_repo, 'fig_1')));

plot_fig_1c(path);

plot_fig_1d(path);

%% - FIGURE 2 - %%

addpath(genpath(fullfile(dir_repo, 'fig_2')));

plot_fig_2ab_GrandAvg_HFA(path);

plot_fig_2ab_PeakMetrics_HFA(path); 

plot_fig_2cd_singleTrials_HFA(path);

plot_fig_2e_neurobehav_corr(path);

plot_fig_2f_neurobehav_corr(path);

%% - FIGURE 3 - %%

addpath(genpath(fullfile(dir_repo, 'fig_3')));

plot_fig_3ab_ramping(path);

plot_fig_3cd_tfr(path);

plot_fig_3ef_slope(path);

plot_fig_3ef_entropy(path);

%% - FIGURE 4 - %%

addpath(genpath(fullfile(dir_repo, 'fig_4')));

plot_fig_4c_psd(path);

plot_fig_4e_ipi(path);

plot_fig_4f_ipi(path);

plot_fig_4g_iplv(path);

plot_fig_4h_psi(path);

%% - FIGURE 5 - %%

addpath(genpath(fullfile(dir_repo, 'fig_5')));

plot_fig_5a_pta(path);

plot_fig_5b_ipi(path);

plot_fig_5c_decoding(path);

plot_fig_5d_connectivity(path);

%% - SUPPLEMENTARY FIGURE 1 - %%

addpath(genpath(fullfile(dir_repo, 'fig_S1')));

plot_fig_S1_RTDist(path);

%% - SUPPLEMENTARY FIGURE 2  - %%

addpath(genpath(fullfile(dir_repo, 'fig_S2')));

plot_fig_S2ab_sdt(path);

plot_fig_S2c_lba(path);

%% - SUPPLEMENTARY FIGURE 3  - %%

addpath(genpath(fullfile(dir_repo, 'fig_S3')));

plot_fig_S3(path);

%% - SUPPLEMENTARY FIGURE 5  - %%

addpath(genpath(fullfile(dir_repo, 'fig_S5')));

plot_fig_S5a(path);

plot_fig_S5b(path);

%% - SUPPLEMENTARY FIGURE 6  - %%

addpath(genpath(fullfile(dir_repo, 'fig_S6')));

plot_fig_S6ab_MDD(path);

plot_fig_S6ab_StateTrans(path);

plot_fig_S6cd_PC1(path);

plot_fig_S6e_inset(path);

plot_fig_S6f_ExplVar(path);

plot_fig_S6g_decoding(path);

%% - SUPPLEMENTARY FIGURE 7  - %%

addpath(genpath(fullfile(dir_repo, 'fig_S7')));

plot_fig_S7_PC1_HLL(path);

%% - SUPPLEMENTARY FIGURE 8 - %%

addpath(genpath(fullfile(dir_repo, 'fig_S8')));

plot_fig_S8_PC1_to_PC5(path);

%% - SUPPLEMENTARY FIGURE 9 - %%

addpath(genpath(fullfile(dir_repo, 'fig_S9')));

plot_fig_S9_action_decoding_HLL(path);

%% - SUPPLEMENTARY FIGURE 10  - %%

addpath(genpath(fullfile(dir_repo, 'fig_S10')));
 
plot_fig_S10a(path);

plot_fig_S10b(path);

%% - SUPPLEMENTARY FIGURE 11  - %%
 
addpath(genpath(fullfile(dir_repo, 'fig_S11')));

plot_fig_S11_decoding(path);

%% - SUPPLEMENTARY FIGURE 12  - %%

addpath(genpath(fullfile(dir_repo, 'fig_S12')));

plot_fig_S12_decoding_latencies(path);

%% - SUPPLEMENTARY FIGURE 13  - %%

addpath(genpath(fullfile(dir_repo, 'fig_S13')));

plot_fig_S13_timing_stop_trials(path);

%% - SUPPLEMENTARY FIGURE 14 - %%

addpath(genpath(fullfile(dir_repo, 'fig_S14')));

plot_fig_S14_pca_shuffling(path);
