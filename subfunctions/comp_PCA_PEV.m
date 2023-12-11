%% Principal Component Analysis on F-statistic obtained from sliding window ANOVA

%%

% Step 1: Load Data obtained from sliding window ANOVA 

% Step 2: Extract MNI coordinates for active channels for later group-level analyses

% Step 3: Extract trials for task-active electrodes and concatenate along subjects

% Step 4: Perform PCA on subject-concatenated data (channel x time)

% Step 5: PCA on randomly permuted data to get threshold for significant
% components (in terms of explained variance)

% Step 6: Extract highly loading channels: A loading channel is defined as
% a channel that exceeds the 75th percentile of all positive r values as the cutoff criterion 
% (important: cutoff is not restricted to only the first PC, but to PCs
% that explain significantly more variance as compared to an empirical null
% distribution)

% Step 7: Store data

%%

clear
ft_defaults

%% Subjects

% ID 18 missing because of LP at 50Hz
subjects = {'OSL13', 'OSL14', 'OSL17', 'OSL21', 'OSL24', 'OSL26', 'OSL27', 'OSL28', ...
    'OSL29', 'OSL30', 'OSL31', 'OSL32', 'OSL34', 'OSL35', 'OSL36', 'OSL38', 'OSL39', 'OSL40'};

%% Paths

path_in     = '/Volumes/IMPECOG/iEEG_Data_JW/SingleSubjects';

path_elec   = '/Volumes/IMPECOG/elecs_and_info';

addpath('/Volumes/IMPECOG/iEEG_Data_JW/functions/');

%% Initialize Variables

% active data (active electrodes)
activedata.probe.trial = []; activedata.bp.trial    = [];
activedata.hll.trial   = []; activedata.br.trial    = []; 
activedata.fb.trial    = [];

activedata.probe.elecs = []; activedata.bp.elecs = [];
activedata.hll.elecs   = []; activedata.br.elecs = [];
activedata.fb.elecs    = [];

activedata.probe.subID = []; activedata.bp.subID = [];
activedata.hll.subID   = []; activedata.br.subID = [];
activedata.fb.subID    = [];

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
    
    % probe
    chanidx = data.probe.w2.active_elecs.chanidx;
    data.probe.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.probe.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.probe.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.probe.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
    
    % bp
    chanidx = data.bp.w2.active_elecs.chanidx;
    data.bp.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.bp.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.bp.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.bp.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
 
    % hll
    chanidx = data.hll.w2.active_elecs.chanidx;
    data.hll.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.hll.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.hll.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.hll.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
    
    % br
    chanidx = data.br.w2.active_elecs.chanidx;
    data.br.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.br.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.br.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.br.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
    
    % fb
    chanidx = data.fb.w2.active_elecs.chanidx;
    data.fb.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.fb.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.fb.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.fb.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
    
    %% Extract trials from task active electrodes

    if ~isempty(data.probe.w2.idxchan_sigthrestime)
        activedata.probe.trial = cat(1, activedata.probe.trial, data.probe.w2.fstat(data.probe.w2.idxchan_sigthrestime,:));
        activedata.probe.elecs = horzcat(activedata.probe.elecs, data.probe.w2.active_elecs.mni_elec.label);
        activedata.probe.subID = horzcat(activedata.probe.subID, repmat({subID}, 1, numel(data.probe.w2.idxchan_sigthrestime)));
    end
    if ~isempty(data.bp.w2.idxchan_sigthrestime)
        activedata.bp.trial    = cat(1, activedata.bp.trial, data.bp.w2.fstat(data.bp.w2.idxchan_sigthrestime,:));
        activedata.bp.elecs    = horzcat(activedata.bp.elecs, data.bp.w2.active_elecs.mni_elec.label);
        activedata.bp.subID    = horzcat(activedata.bp.subID, repmat({subID}, 1, numel(data.bp.w2.idxchan_sigthrestime)));
    end
    if ~isempty(data.hll.w2.idxchan_sigthrestime)
        activedata.hll.trial   = cat(1, activedata.hll.trial, data.hll.w2.fstat(data.hll.w2.idxchan_sigthrestime,:));
        activedata.hll.elecs   = horzcat(activedata.hll.elecs, data.hll.w2.active_elecs.mni_elec.label);
        activedata.hll.subID   = horzcat(activedata.hll.subID, repmat({subID}, 1, numel(data.hll.w2.idxchan_sigthrestime)));
    end
    if ~isempty(data.br.w2.idxchan_sigthrestime)
        activedata.br.trial    = cat(1, activedata.br.trial, data.br.w2.fstat(data.br.w2.idxchan_sigthrestime,:));
        activedata.br.elecs    = horzcat(activedata.br.elecs, data.br.w2.active_elecs.mni_elec.label);
        activedata.br.subID    = horzcat(activedata.br.subID, repmat({subID}, 1, numel(data.br.w2.idxchan_sigthrestime)));
    end
    if ~isempty(data.fb.w2.idxchan_sigthrestime)
        activedata.fb.trial    = cat(1, activedata.fb.trial, data.fb.w2.fstat(data.fb.w2.idxchan_sigthrestime,:));
        activedata.fb.elecs    = horzcat(activedata.fb.elecs, data.fb.w2.active_elecs.mni_elec.label);
        activedata.fb.subID    = horzcat(activedata.fb.subID, repmat({subID}, 1, numel(data.fb.w2.idxchan_sigthrestime)));
    end

end

%% Get Time Information Into Structure

activedata.probe.time = data.probe.w2.time;
activedata.bp.time    = data.bp.w2.time;
activedata.hll.time   = data.hll.w2.time;
activedata.br.time    = data.br.w2.time;
activedata.fb.time    = data.fb.w2.time;

%% Perform PCA

%% Probe

[pcainfo.probe.coeff,~,~,~,pcainfo.probe.explained,~] = pca(activedata.probe.trial');

% determine significant PCs using permutation testing

nperm         = 1000;
permexplained = zeros(nperm,1);

checkprocess = linspace(10,100,10);

for permi = 1:nperm
    
    if ismember(permi/nperm*100,checkprocess)
        fprintf('probe process %.0f%%\n', permi/nperm*100);
    end
    % random data by randomizing data
    randdat = reshape(activedata.probe.trial(randperm(numel(activedata.probe.trial))),size(activedata.probe.trial));
    [coeff,~,~,~,randexplained,~]  = pca(randdat');

    permexplained(permi,1) = max(randexplained);
    
end

pcainfo.probe.pcbychance = mean(permexplained);

% select channels loading high on significant components

% find significant PCs

sigpcs = find(pcainfo.probe.explained > pcainfo.probe.pcbychance);

% find channels exceeding 75th percentile for the sig components

pcchannels = cell(length(sigpcs),1);

for compi = 1:length(sigpcs)
    
    % find percentile for respective component
    
    tmppercentile = prctile(pcainfo.probe.coeff(:,compi), 75);
    
    pcchannels{compi} = find(pcainfo.probe.coeff(:,compi) > tmppercentile);
    
end
    
% get unique channels

pcainfo.probe.loadingchans = unique(reshape([pcchannels{:}], [], 1));

%% BP

[pcainfo.bp.coeff,~,~,~,pcainfo.bp.explained,~] = pca(activedata.bp.trial');

% determine significant PCs using permutation testing

nperm         = 1000;
permexplained = zeros(nperm,1);

checkprocess = linspace(10,100,10);

for permi = 1:nperm
    
    if ismember(permi/nperm*100,checkprocess)
        fprintf('bp process %.0f%%\n', permi/nperm*100);
    end
    % random data by randomizing data
    randdat = reshape(activedata.bp.trial(randperm(numel(activedata.bp.trial))),size(activedata.bp.trial));
    [coeff,~,~,~,randexplained,~]  = pca(randdat');

    permexplained(permi,1) = max(randexplained);
    
end

pcainfo.bp.pcbychance = mean(permexplained);

% select channels loading high on significant components

% find significant PCs

sigpcs = find(pcainfo.bp.explained > pcainfo.bp.pcbychance);

% find channels exceeding 75th percentile for the sig components

pcchannels = cell(length(sigpcs),1);

for compi = 1:length(sigpcs)
    
    % find percentile for respective component
    
    tmppercentile = prctile(pcainfo.bp.coeff(:,compi), 75);
    
    pcchannels{compi} = find(pcainfo.bp.coeff(:,compi) > tmppercentile);
    
end
    
% get unique channels

pcainfo.bp.loadingchans = unique(reshape([pcchannels{:}], [], 1));

%% HLL

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

pcainfo.hll.pcbychance = mean(permexplained);

% select channels loading high on significant components

% find significant PCs

sigpcs = find(pcainfo.hll.explained > pcainfo.hll.pcbychance);

% find channels exceeding 75th percentile for the sig components

pcchannels = cell(length(sigpcs),1);

for compi = 1:length(sigpcs)
    
    % find percentile for respective component
    
    tmppercentile = prctile(pcainfo.hll.coeff(:,compi), 75);
    
    pcchannels{compi} = find(pcainfo.hll.coeff(:,compi) > tmppercentile);
    
end
    
% get unique channels

pcainfo.hll.loadingchans = unique(reshape([pcchannels{:}], [], 1));


%% BR

[pcainfo.br.coeff,~,~,~,pcainfo.br.explained,~] = pca(activedata.br.trial');

% determine significant PCs using permutation testing

nperm         = 1000;
permexplained = zeros(nperm,1);

checkprocess = linspace(10,100,10);

for permi = 1:nperm
    
    if ismember(permi/nperm*100,checkprocess)
        fprintf('br process %.0f%%\n', permi/nperm*100);
    end
    % random data by randomizing data
    randdat = reshape(activedata.br.trial(randperm(numel(activedata.br.trial))),size(activedata.br.trial));
    [coeff,~,~,~,randexplained,~]  = pca(randdat');

    permexplained(permi,1) = max(randexplained);
    
end

pcainfo.br.pcbychance = mean(permexplained);

% select channels loading high on significant components

% find significant PCs

sigpcs = find(pcainfo.br.explained > pcainfo.br.pcbychance);

% find channels exceeding 75th percentile for the sig components

pcchannels = cell(length(sigpcs),1);

for compi = 1:length(sigpcs)
    
    % find percentile for respective component
    
    tmppercentile = prctile(pcainfo.br.coeff(:,compi), 75);
    
    pcchannels{compi} = find(pcainfo.br.coeff(:,compi) > tmppercentile);
    
end
    
% get unique channels

pcainfo.br.loadingchans = unique(reshape([pcchannels{:}], [], 1));

%% FB

[pcainfo.fb.coeff,~,~,~,pcainfo.fb.explained,~] = pca(activedata.fb.trial');

% determine significant PCs using permutation testing

nperm         = 1000;
permexplained = zeros(nperm,1);

checkprocess = linspace(10,100,10);

for permi = 1:nperm
    
    if ismember(permi/nperm*100,checkprocess)
        fprintf('fb process %.0f%%\n', permi/nperm*100);
    end
    % random data by randomizing data
    randdat = reshape(activedata.fb.trial(randperm(numel(activedata.fb.trial))),size(activedata.fb.trial));
    [coeff,~,~,~,randexplained,~]  = pca(randdat');

    permexplained(permi,1) = max(randexplained);
    
end

pcainfo.fb.pcbychance = mean(permexplained);

% select channels loading high on significant components

% find significant PCs

sigpcs = find(pcainfo.fb.explained > pcainfo.fb.pcbychance);

% find channels exceeding 75th percentile for the sig components

pcchannels = cell(length(sigpcs),1);

for compi = 1:length(sigpcs)
    
    % find percentile for respective component
    
    tmppercentile = prctile(pcainfo.fb.coeff(:,compi), 75);
    
    pcchannels{compi} = find(pcainfo.fb.coeff(:,compi) > tmppercentile);
    
end
    
% get unique channels

pcainfo.fb.loadingchans = unique(reshape([pcchannels{:}], [], 1));

%% Store Highly Loading Channels for each Participant

activedata.probe.pca    = pcainfo.probe;
activedata.bp.pca       = pcainfo.bp;
activedata.hll.pca      = pcainfo.hll;
activedata.br.pca       = pcainfo.br;
activedata.fb.pca       = pcainfo.fb;

mergedata = {activedata.probe, activedata.bp, activedata.hll, activedata.br, activedata.fb};

channelstore = cell(numel(subjects),numel(mergedata));

for condi = 1:numel(mergedata)
    
    counter      = zeros(numel(subjects),1); % counter to store multiple channels in channelstore

    for chani = 1:numel(mergedata{condi}.pca.loadingchans)
        
        % get index for loading channels
        tmpchanidx = mergedata{condi}.pca.loadingchans(chani);
        
        % get subject to which the channels corresponds to
        tmpsubject_idx = find(strcmp(subjects, mergedata{condi}.subID{tmpchanidx}));
                
        counter(tmpsubject_idx) = counter(tmpsubject_idx) + 1;
        
        % store channel
        
        channelstore{tmpsubject_idx,condi}.elecs(counter(tmpsubject_idx)) = mergedata{condi}.elecs(tmpchanidx);
        
    end
    
end

clear mergedata
clear data

%% Store Data in Individual Data from the Sliding-Window ANOVA 
        
for subji = 1:numel(subjects)
    
    
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
    
    % probe
    chanidx = data.probe.w2.active_elecs.chanidx;
    data.probe.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.probe.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.probe.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.probe.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
    
    % bp
    chanidx = data.bp.w2.active_elecs.chanidx;
    data.bp.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.bp.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.bp.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.bp.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
 
    % hll
    chanidx = data.hll.w2.active_elecs.chanidx;
    data.hll.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.hll.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.hll.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.hll.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
    
    % br
    chanidx = data.br.w2.active_elecs.chanidx;
    data.br.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.br.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.br.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.br.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
    
    % fb
    chanidx = data.fb.w2.active_elecs.chanidx;
    data.fb.w2.active_elecs.mni_elec.label     = elec_mni_frv.label(chanidx)';
    data.fb.w2.active_elecs.mni_elec.chanpos   = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        data.fb.w2.active_elecs.mni_elec.chantype  = elec_mni_frv.chantype(chanidx);
    end
    data.fb.w2.active_elecs.mni_elec.elecpos   = elec_mni_frv.elecpos(chanidx,:);
        
    
    %% Get significant electrodes as estimated via PCA
    
    MergedData = {data.probe, data.bp, data.hll, data.br, data.fb};
    IdxData  = {'Probe', 'BP', 'HLL', 'BR', 'FB'}; 
    
    for condi = 1:numel(MergedData)
                
        % compare labels
        xx      = cellfun(@(c)strcmp(c,MergedData{condi}.w2.active_elecs.mni_elec.label), ...
            channelstore{subji,condi}.elecs,'UniformOutput',false);
        % concatenate cell arrays
        xx_mod = vertcat(xx{:});
        
        % extract channel labels that correspond to channels in PCA
        [~,chanidx] = find(xx_mod);
        
        % store labels in data structure for PCA
                                                
        MergedData{condi}.pca = MergedData{condi}.w2;
        
        if ~isempty(MergedData{condi}.pca.active_elecs.chanidx)
            MergedData{condi}.pca.active_elecs.chanidx = MergedData{condi}.pca.active_elecs.chanidx(chanidx);
        end
        if ~isempty(MergedData{condi}.pca.idxchan_sigthrestime)
            MergedData{condi}.pca.idxchan_sigthrestime = MergedData{condi}.pca.idxchan_sigthrestime(chanidx);
        end
        % index electrode in acpc format
        if isfield(MergedData{condi}.pca.active_elecs.elec, 'label')
            MergedData{condi}.pca.active_elecs.elec.label    = MergedData{condi}.pca.active_elecs.elec.label(chanidx);
        end
        if isfield(MergedData{condi}.pca.active_elecs.elec, 'chanpos')
            MergedData{condi}.pca.active_elecs.elec.chanpos  = MergedData{condi}.pca.active_elecs.elec.chanpos(chanidx,:);
        end
        if isfield(MergedData{condi}.pca.active_elecs.elec, 'chantype')
            MergedData{condi}.pca.active_elecs.elec.chantype = MergedData{condi}.pca.active_elecs.elec.chantype(chanidx);
        end
        if isfield(MergedData{condi}.pca.active_elecs.elec, 'elecpos')
            MergedData{condi}.pca.active_elecs.elec.elecpos  = MergedData{condi}.pca.active_elecs.elec.elecpos(chanidx,:);
        end
        
        % index electrode in MNI format
        if isfield(MergedData{condi}.pca.active_elecs.mni_elec, 'label')
            MergedData{condi}.pca.active_elecs.mni_elec.label    = MergedData{condi}.pca.active_elecs.mni_elec.label(chanidx);
        end
        if isfield(MergedData{condi}.pca.active_elecs.mni_elec, 'chanpos')
            MergedData{condi}.pca.active_elecs.mni_elec.chanpos  = MergedData{condi}.pca.active_elecs.mni_elec.chanpos(chanidx,:);
        end
        if isfield(MergedData{condi}.pca.active_elecs.mni_elec, 'chantype')
            MergedData{condi}.pca.active_elecs.mni_elec.chantype = MergedData{condi}.pca.active_elecs.mni_elec.chantype(chanidx);
        end
        if isfield(MergedData{condi}.pca.active_elecs.mni_elec, 'elecpos')
            MergedData{condi}.pca.active_elecs.mni_elec.elecpos  = MergedData{condi}.pca.active_elecs.mni_elec.elecpos(chanidx,:);
        end
        
        % remove fields
        MergedData{condi}.pca = rmfield(MergedData{condi}.pca, 't_sig_threstime');
        MergedData{condi}     = rmfield(MergedData{condi}, 'w2');
        
        if strcmp(IdxData{condi}, 'Probe')
            
            pcadata.probe = MergedData{condi};
                        
        elseif strcmp(IdxData{condi}, 'BP')
            
            pcadata.bp    = MergedData{condi};
            
        elseif strcmp(IdxData{condi}, 'HLL')
            
            pcadata.hll   = MergedData{condi};
            
        elseif strcmp(IdxData{condi}, 'BR')
            
            pcadata.br    = MergedData{condi};
            
        elseif strcmp(IdxData{condi}, 'FB')
            
            pcadata.fb    = MergedData{condi};
            
        end
        
    end
    
    % save data 
    savingname = [path_in filesep subID filesep 'Data' filesep 'HFA' ...
        filesep strcat(subID, '_TaskActive_PCAexpvar')];
    
    save(savingname, 'pcadata', '-v7.3');
    
    clear pcadata

end
   
        
