%% Figure S6e (inset)

%%

clear
ft_defaults

close all;

%% - SUBJECTS - %%

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% - PATH SETTINGS - %%

path_in   = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out  = '/Users/Janweber 1/Desktop/PhD_HelfrichLab/Project_1_Predictions_ECoG/Results/Figures/submission_natcomm_v3/Figure_5/revision_natcomm_final/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

%% - COMPUTATION START - %%

% PC of interest

pc2use  = 1; % PC1

% preallocate memory for zvalue, RTs and trial information

zval    = cell(numel(subjects),1);
rts     = cell(numel(subjects),1);
trlinfo = cell(numel(subjects),1);

for Isub = 1:numel(subjects)
    
    % get subject ID
    subID = subjects{Isub};
    
    % print subject ID
    fprintf('Running # %d/%d\n', Isub, numel(subjects));
    
    % check if folder exists and if there is data
    if isfolder(fullfile(path_in, subID, 'Data', 'State_Space'))
               
        % load data
        load(fullfile(path_in, subID, 'Data', 'State_Space', [subID '_stimcoding_all_elecs_BR_w_surro.mat']));
        
    end

    % make sure both frontal and motor cortex contain data
    if isfield(all_data, 'frontal') && isfield(all_data, 'motor')
        
        % get time index for computation
        tidx = nearest(all_data.frontal.time, -0.5):nearest(all_data.frontal.time, 0.2);

        % get reaction times
        rts{Isub} = all_data.frontal.trialinfo(:,11);
        
        % get trial type info
        trlinfo{Isub} = all_data.frontal.trialinfo(:,13);
        
        % get data for both regions
        
        tmp_frontal = squeeze(all_data.frontal.trial(:, pc2use, tidx));
        
        tmp_motor   = squeeze(all_data.motor.trial(:, pc2use, tidx));
        
        % compute power correlation on a trial-by-trial basis
        
        rho = NaN(size(tmp_frontal, 1),1);
        
        for Itrial = 1:size(tmp_frontal, 1)
            
            rho(Itrial) = corr(tmp_frontal(Itrial,:)', tmp_motor(Itrial,:)'); 
            
        end
        
        % create H0 by randomly cutting the data
        
        niter   = 1000;
        
        % preallocate memory (permutation x trial)
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
                              
                % compute correlation
                r = corr(shuf_frontal', tmp_motor(Itrial,:)');
                
                permrho(Iperm,Itrial) = r;
                
            end %Itrial
            
        end %Iperm
           
        % z-transform rho values based on H0
        
        zval{Isub} = (rho - mean(permrho)') ./ std(permrho)';
        
    end
    
end % Isub
        
%% - REMOVE EMPTY CELLS - %%

zval    = zval(~cellfun(@isempty, zval));
rts     = rts(~cellfun(@isempty, rts));
trlinfo = trlinfo(~cellfun(@isempty, trlinfo));

%% - CONCATENATE ALL TRIALS - %%

all_zval = cat(1, zval{:}); % trial x 1
all_rts  = cat(1, rts{:});

%% - COMPUTE THE AVERAGE Z-VALUE - %%

close all;

% addpath /Users/janweber/Documents/MATLAB/general_functions/colormaps/
% 
% cm_data = viridis(100);
% 
% c2use   = cm_data(70,:);

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
set(gcf, 'Position', [547 443 157 211]);

p = signrank(avg_zval_all);
fprintf('Wilcoxon Test vs. 0; p = %.5f\n', p);

%% plot average z-value per condition

figure;

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

    hold on
    
    scatter(ones(size(data2use{ii}))*ii, data2use{ii}, 'filled', 'MarkerFaceColor', ...
        c2use{ii}, 'jitter','on', 'jitterAmount',0.04);

end
    
xticks(1:3);
xticklabels({'0%', '25%', '75%'}); xtickangle(45);
ylabel('rho');

set(gca, 'fontsize', 13, 'linewidth', 2);

box off

set(gcf, 'Position', [511 383 258 184]);

% - SAVE FIGURES - %

cd(path_out);
print(gcf,'fig_5_inset_powcorr_pc1.pdf','-dpdf','-r400');    

%% plot histogram of all z values

cm_data=viridis(100);

c2use   = cm_data(70,:);

figure;
a = histfit(all_zval);
a(1).FaceColor = c2use;
a(1).EdgeColor = a(1).FaceColor;
a(1).FaceAlpha = 0.7;
a(1).EdgeAlpha = 0.7;
a(1).LineWidth = 3;

xlabel('rho [zval]');
ylabel('Count');
xline(0, 'LineStyle', '--', 'LineWidth', 1.5);

set(gca, 'linew', 2, 'fontsize', 13); box off;
set(gcf, 'Position', [514 442 172 138]);

% - SAVE FIGURES - %

cd(path_out);
print(gcf,'fig_5_hist_rho_powcorr.pdf','-dpdf','-r400');    


%% Plot Correlation z-value and RT

idxout = isoutlier(all_rts, 'median') | ...
    isoutlier(all_zval, 'median');

all_rts(idxout)  = [];
all_zval(idxout) = [];

figure;
scatter(all_rts,all_zval, 50, 'filled', 'MarkerFaceColor', c2use, 'MarkerEdgeColor', c2use, ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);

xlabel('RT [s]');
ylabel('rho [z]');
set(gca, 'linew', 2, 'fontsize', 13); box off;

set(gcf, 'Position', [511 429 194 117]);


a           = lsline;
a.LineWidth = 2;
a.Color     = 'k';

[r,p] = corr(all_rts,all_zval, 'rows', 'complete');

fprintf('Correlation RT ~ Power; r = %.3f, p = %.3f\n', r, p);

% compute linear regression
mdl = fitlm(all_zval, all_rts);

anova(mdl, 'summary')
fprintf('\nAdjusted Rsquared = %.4f\n', mdl.Rsquared.Adjusted);

% - SAVE FIGURES - %

cd(path_out);
print(gcf,'fig_5_corr_rt_powcorr_rho.pdf','-dpdf','-r400');    

%% Save Data

% -------------------------------------------------
% Save Data into Excel Sheet for Stats in Python
% -------------------------------------------------

% go into results folder
mkdir(path_source, 'Fig_5e_inset');
cd(fullfile(path_source, 'Fig_5e_inset'));

statsdata = [];

% include relative explained variance of cumulative sum of first 3 components
statsdata(:,3) = vertcat(avg_zval_0, avg_zval_25, avg_zval_75);
                     
% give subject an ID                    
statsdata(1:end,1) = vertcat(repmat(1:length(avg_zval_0), 1, 3)');

statsdata(1:end,2) = vertcat(ones(1, length(avg_zval_0))', ...
    ones(1, length(avg_zval_25))'*2, ones(1, length(avg_zval_75))'*3);
                                                    
% save data
statstable = array2table(statsdata, 'VariableNames', {'ID', 'Prediction', 'rho'});

% write data into excelsheet
filename = 'Fig_5e_inset.xlsx';

writetable(statstable,filename);


