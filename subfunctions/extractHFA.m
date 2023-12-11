function extractHFA 

%% Extract High-Frequency Activity

% 1. Extract high gamma in 10Hz bins and average later
% 2. Baseline correction using a bootstrap z-score procedure 

%% 

clear
ft_defaults

%% Paths

addpath('/Volumes/IMPECOG/iEEG_Data_JW/functions/');

path_in  = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

%% Subjects

subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL18', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Features

hg_vec   = 75:10:145;      % center frequencies
win      = 5;              % freq window around center
niter    = 1000;           % number of bootstrap iterations
bl_start = 0.3;            % baseline from -200 - 0ms prior to probe stimulus onset
bl_end   = 0.5;

%% Loop over subjects

for subji = 1:numel(subjects)
    
    subID = subjects{subji};
    
    fprintf('Computing HFA for %s\n', subID)
    
    % go into subject folder
    
    cd([path_in filesep subID filesep 'Data']);
    
    if exist([subID '_clean.mat'], 'file') == 2
        
        load([subID '_clean.mat']);

        load([subID '_Trial_Info_clean.mat']);

        load([subID '_Fs.mat']);

        data_clean.fsample = fsample;
        
    else 
        
        continue;
        
    end
    
    % assign memory
    hgtrace = zeros([length(hg_vec) size(data_clean.trial,1) length(data_clean.label) length(data_clean.time)]);

    for h = 1:length(hg_vec)

        % filter the data in the range first
        cfg              = [];
        cfg.bpfilter     = 'yes';
        cfg.bpfreq       = [hg_vec(h)-win hg_vec(h)+win];
        cfg.hilbert      = 'abs';
        filt             = ft_preprocessing(cfg, data_clean);

        % zbaseline correct
        filt             = rh_zbaseline_erp(filt, bl_start, bl_end, niter);

        hgtrace(h,:,:,:) = filt.ztrial;
        clear filt

    end
    
    % and edit filt to save high freq activity
    hfa       = data_clean;
    hfa.trial = hgtrace; 
    hfa.cfg   = [];
    
    clear hgtrace   
    
    % go into output folder
    
    cd([path_in filesep subID filesep 'Data' filesep 'HFA']);
        
    savename = [subID '_HFA.mat'];
    
    save(savename, 'hfa', '-v7.3');
    
    % save the averaged HFA file as welll
    hfa.trial = squeeze(nanmean(hfa.trial,1)); % average over HFA bins

    savename = [subID '_HFAavgoverbins.mat'];
    
    save(savename, 'hfa', '-v7.3');
    
end


end