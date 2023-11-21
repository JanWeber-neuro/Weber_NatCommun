%% Figure S14

clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_out    = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm/Figures/';

path_source  = '/Volumes/IMPECOG/iEEG_Data_JW/Key_Code_Revision_NatComm_Final/Data/SourceData/';

path_elec   = '/Volumes/IMPECOG/elecs_and_info';

addpath('/Volumes/IMPECOG/iEEG_Data_JW/functions/');

%% Initialize Variables

% active data (active electrodes)
activedata.hll.trial   = [];
activedata.hll.elecs   = []; 
activedata.hll.subID   = [];

%%

for subji = 1:numel(subjects)
    
    %% Step 1: Load Data and Electrode Files
    
    % get indicated subject ID

    subID = subjects{subji};

    % go into HFA folder
    cd([path_in filesep subID filesep 'Data']);

    if ~isfolder('HFA') % if folder does not exist

        continue;

    else

        cd([path_in filesep subID filesep 'Data' filesep 'HFA']);

        if exist([subID '_TaskActive_ExpVariance.mat'], 'file') == 2

            fprintf('Load Data %s\n', subID);

            % load explained variance data

            load([subID '_TaskActive_ExpVariance.mat']);
            
            data = anovadata; clear anovadata;
            
        else 

            warning('Data for the current subject is not available');

        end

    end
    
    % load MNI electrode file

    cd([path_elec filesep subID filesep 'Electrodes']);

    elecfile = dir('*elec_mni_frv.mat');
    
    load(elecfile.name);
       
    %% Extract MNI coordinates for active channels
 
    % hll
    chanidx = data.hll.w2.active_elecs.chanidx;
    data.hll.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.hll.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.hll.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.hll.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
        
    %% Extract trials from task active electrodes

    if ~isempty(data.hll.w2.idxchan_sigthrestime)
        activedata.hll.trial   = cat(1, activedata.hll.trial, data.hll.w2.fstat(data.hll.w2.idxchan_sigthrestime,:));
        activedata.hll.elecs   = horzcat(activedata.hll.elecs, data.hll.w2.active_elecs.mni_elec.label);
        activedata.hll.subID   = horzcat(activedata.hll.subID, repmat({subID}, 1, numel(data.hll.w2.idxchan_sigthrestime)));
    end

end

%% Get Time Information Into Structure

activedata.hll.time   = data.hll.w2.time;

%% Perform PCA

% go to output folder to store figures
cd(path_out);

[pcainfo.hll.coeff,~,~,~,pcainfo.hll.explained,~] = pca(activedata.hll.trial');

% determine significant PCs using permutation testing

nperm         = 1000;
permexplained = zeros(nperm,1);

checkprocess = linspace(10,100,10);

for permi = 1:nperm
    
    if ismember(permi/nperm*100,checkprocess)
        fprintf('hll process %.0f%%\n', permi/nperm*100);
    end
    % random data by randomizing data
    randdat = reshape(activedata.hll.trial(randperm(numel(activedata.hll.trial))),size(activedata.hll.trial));
    [coeff,~,~,~,randexplained,~]  = pca(randdat');

    permexplained(permi,1) = max(randexplained);
    
end

figure;
h = histogram(permexplained);
h.FaceColor = 'k';
h.FaceAlpha = 1;

ylabel('Count');
xlabel('% explained variance');

box off
set(gca, 'fontsize', 13);

pcainfo.hll.pcbychance = mean(permexplained);

xline(pcainfo.hll.pcbychance, 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);

print(gcf,'Histogram_Chance_Level_PCA.pdf','-dpdf','-r400');    

% plot permutation distribution and explained variance by last sig. PC

figure;
h = histogram(permexplained);
h.FaceColor = 'k';
h.FaceAlpha = 1;

ylabel('Count');
xlabel('% explained variance');

box off
set(gca, 'fontsize', 13);

% mean variance explained by chance
pcainfo.hll.pcbychance = mean(permexplained);

% get significant PCs
sigpcs = find(pcainfo.hll.explained > pcainfo.hll.pcbychance);

xline(pcainfo.hll.pcbychance, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');

xline(pcainfo.hll.explained(sigpcs(end)), 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);

xline(pcainfo.hll.explained(sigpcs(end)+1), 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);

xlim([0.8 1.5]);

print(gcf,'Histogram_Chance_Level_PCA_2.pdf','-dpdf','-r400');    

% save source data
mkdir(path_source, 'Fig_S14');
cd(fullfile(path_source, 'Fig_S14'));

table_data = array2table(pcainfo.hll.explained);

% write data into excelsheet
filename = 'fig_S14_PCA_ev.xlsx';
writetable(table_data,filename, 'Sheet', 1);

table_data = array2table(permexplained);

% write data into excelsheet
filename = 'fig_S14_PCA_permexplained.xlsx';
writetable(table_data,filename, 'Sheet', 1);

table_data = array2table(pcainfo.hll.pcbychance);

% write data into excelsheet
filename = 'fig_S14_PCA_chancelevel.xlsx';
writetable(table_data,filename, 'Sheet', 1);

table_data = array2table(sigpcs);

% write data into excelsheet
filename = 'fig_S14_PCA_sigpcs.xlsx';
writetable(table_data,filename, 'Sheet', 1);

% select channels loading high on significant components

% find significant PCs

figure;
plot(cumsum(pcainfo.hll.explained(sigpcs)), 'k.', 'MarkerSize', 20);
hold on
plot(cumsum(pcainfo.hll.explained(sigpcs)), 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

xlabel('Principal Component');
ylabel('Cumulative %EV');

box off
set(gca, 'fontsize', 13);

print(gcf,'Cumulative_Variance_Explained_PCA.pdf','-dpdf','-r400');    

x = cumsum(pcainfo.hll.explained(sigpcs));

fprintf('number of sig. components = %d, explained var = %.3f\n', sigpcs(end), x(end));

% save as source data
table_data = array2table(cumsum(pcainfo.hll.explained(sigpcs)));

% write data into excelsheet
filename = 'fig_S14_PCA_cumsum_ev.xlsx';
writetable(table_data,filename, 'Sheet', 1);

% find channels exceeding 75th percentile for the sig components

pcchannels = cell(length(sigpcs),1);

for compi = 1:length(sigpcs)
    
    % find percentile for respective component
    
    tmppercentile = prctile(pcainfo.hll.coeff(:,compi), 75);
    
    pcchannels{compi} = find(pcainfo.hll.coeff(:,compi) > tmppercentile);
    
end
    
% get unique channels

pcainfo.hll.loadingchans = unique(reshape([pcchannels{:}], [], 1));

%% Plot Time Series of PCs

addpath /Users/'Janweber 1'/Documents/MATLAB/general_functions/colormaps/

cmap = viridis(3);

% HLL
figure;

for jj = 1:size(cmap,1)
    plot(activedata.hll.time, pcainfo.hll.coeff(:,jj)' * activedata.hll.trial, 'Color', cmap(jj,:), 'LineWidth', 3);
    hold on
end

xlim([activedata.hll.time(1) activedata.hll.time(end)]);
set(gca, 'fontsize', 13);
xlabel('Time [s]');
ylabel('Activity [au]');
legend('PC1', 'PC2', 'PC3', 'Location', 'best');

box off                

print(gcf,'PCA_traces_F_series.pdf','-dpdf','-r400');    
