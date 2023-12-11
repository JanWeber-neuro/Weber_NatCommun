%% Extract task informative electrodes

%%

% Step 1: Load Data and Electrode Files

% Step 2: Segment the data to all conditions (probe, button press, hit lower limit, button release)

% Step 3: Perform sliding window ANOVA

% Step 4: Find active electrodes that exceeded critical threshold

% Step 5: Extract Index and Electrode Positions for Active Elecs (for plotting purposes)


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

%% Features

fsample     = 512;   % sampling freq
smoothwind  = 0.05;  % smoothingwindow for sliding explained variance

%% 

% initiliaze data structure
hfa_all        = struct();

% initialize structure containing info which electrodes are considered as active
AllActiveElecs = struct();

% initialize structure containing all trial information after segmentation
AllSampleInfo  = struct();

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

        if exist([subID '_HFAavgoverbins.mat'], 'file') == 2

            fprintf('Load Data %s\n', subID);

            % load high frequency activity

            load([subID '_HFAavgoverbins.mat']);

            % go into data folder

            cd([path_in filesep subID filesep 'Data']);

            load([subID '_Trial_Info_clean.mat']);

            hfa.fsample = fsample;

        else 

            warning('Data for the current subject is not available');

        end

    end

    %% Load electrode and surface files

    disp('Loading electrode files...');

    % load electrode file

    cd([path_elec filesep subID filesep 'Electrodes']);

    elecfile = dir('*elec_acpc*.mat');

    if numel(elecfile) > 1 % if more than one electrode file exists

        elecfile = elecfile(1); % take the first file

    end

    elec = load(elecfile.name);

    if isfield(elec, 'elec_acpc_f')

        elec = elec.elec_acpc_f;

    else

        elec = elec.elec_acpc;

    end  

    if strcmp(subID, 'OSL32')

        for ii = 1:numel(elec.label)

            tmplabel  = elec.label{ii}; % complete label

            labelname = tmplabel(isstrprop(tmplabel,'alpha')); % label name

            labelnum  = regexp(tmplabel,'\d*','Match'); % label number

            elec.label{ii} = [labelname(2) labelnum{1}];

        end

    elseif strcmp(subID, 'OSL35') || strcmp(subID, 'OSL39')


        for ii = 1:numel(elec.label)

            tmplabel  = elec.label{ii}; % complete label

            labelname = tmplabel(isstrprop(tmplabel,'alpha')); % label name

            labelname = labelname(4:5);

            labelnum  = regexp(tmplabel,'\d*','Match'); % label number

            if str2double(labelnum{1}(1)) == 0

                mod_label = [labelname labelnum{1}(2)]; % combine label and number again

            else

                mod_label = [labelname labelnum{1}];

            end

            elec.label{ii}   = mod_label;

        end

    end


    % load surface files

    cd([path_elec filesep subID filesep 'Surfaces']);

    surffile = dir('*cortex_lh.mat'); % left hemisphere

    load(surffile.name);

    surffile = dir('*cortex_rh.mat'); % right hemisphere

    load(surffile.name);
    
    %% Modify channel labels to match with electrode file

    mod_label = cell(numel(hfa.label),1);

    for chani = 1:numel(hfa.label)

        tmplabel  = hfa.label{chani}; % complete label

        labelname = tmplabel(isstrprop(tmplabel,'alpha')); % label name

        labelnum  = regexp(tmplabel,'\d*','Match'); % label number

        % exceptions where labels are weird
        if ~strcmp(subID, 'OSL26') &&  ~strcmp(subID, 'OSL31') &&  ~strcmp(subID, 'OSL32') ...
                &&  ~strcmp(subID, 'OSL34') &&  ~strcmp(subID, 'OSL38') &&  ~strcmp(subID, 'OSL40')

            if str2double(labelnum{1}(1)) == 0

                mod_label{chani} = [labelname labelnum{1}(2)]; % combine label and number again

            else

                mod_label{chani} = [labelname labelnum{1}];

            end

        else

            if str2double(labelnum{1}(1)) == 0

                mod_label{chani} = [labelname(2) labelnum{1}(2)];

            else

                mod_label{chani} = [labelname(2) labelnum{1}];

            end

        end

    end
    
    %% Step 2: Re-Epoch Data

    disp('Re-Epoching Data...');

    %% Segment Data to Onset = Probe Onset

    trl = NaN(size(hfa.sampleinfo,1),3);

    % segmentation: from -0.5 pre probe onset to onset of the bar
    for triali = 1:size(hfa.trialinfo,1)

        curr_onset = hfa.sampleinfo(triali,1); % current onset = baseline
%         trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + min(trl_mod_table.ProbeOnset) - fsample/2 ...
%             curr_onset + abs(trl_mod_table.DistTrigger(triali)) + min(trl_mod_table.BarOnset) -fsample/2];
        % 273 samples = 0.53 sec. = common onset of probe along subjects
        % 691 samples = 1.3496 sec. (probe + 0.8164) = common onset of bar onset along subjects
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + 273 - fsample/2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + 691 -fsample/2];

    end
    
    AllSampleInfo(subji).trl_probe = trl;
    
    clear trl

    %% Segment Data to Onset = Button Press

    trl = NaN(size(hfa.sampleinfo,1),3);

    % segmentation: from -0.5 pre button press to hit of lower limit
    for triali = 1:size(hfa.trialinfo,1)

        curr_onset = hfa.sampleinfo(triali,1); % current onset = baseline
%         trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonPress(triali) - fsample/2 ...
%             curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonPress(triali) + ... 
%             min(trl_mod_table.HitLowLim - trl_mod_table.ButtonPress) -fsample/2];
        
        % take 286 instead of "min(trl_mod_table.HitLowLim -
        % trl_mod_table.ButtonPress)" as this is not completely consistent
        % across subjects
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonPress(triali) - fsample/2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonPress(triali) + ... 
            286 -fsample/2];

    end

    AllSampleInfo(subji).trl_bp = trl;
    
    clear trl
    
    %% Segment Data to Onset = Hit Lower Limit

    trl = NaN(size(hfa.sampleinfo,1),3);

    % segmentation: from -0.5 pre HLL to 0.3 sec post-HLL
    for triali = 1:size(hfa.trialinfo,1)

        curr_onset = hfa.sampleinfo(triali,1); % current onset = baseline
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.HitLowLim(triali) - fsample/2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.HitLowLim(triali) + round(fsample * 0.3) -fsample/2];

    end
    
    AllSampleInfo(subji).trl_hll = trl;
    
    clear trl

    %% Segment Data to Onset = Button Release (BL)

    trl = NaN(size(hfa.sampleinfo,1),3);

    % segmentation: from -0.5 pre BL to 0.3 sec post-BL
    for triali = 1:size(hfa.trialinfo,1)

        curr_onset = hfa.sampleinfo(triali,1); % current onset = baseline  
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonRelease(triali) - fsample/2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.ButtonRelease(triali) + round(fsample * 0.3) -fsample/2];

    end 
    
    AllSampleInfo(subji).trl_br = trl;
    
    clear trl
    
    %% Segment Data to Feedback (FB)
    
    trl = NaN(size(hfa.sampleinfo,1),3);
    
    % segmentation: from -0.5 pre FB to 1 sec post-FB
    for triali = 1:size(hfa.trialinfo,1)
        
        curr_onset = hfa.sampleinfo(triali,1); % current onset = baseline
        trl(triali,:) = [curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.feedbackOnset(triali) - fsample/2 ...
            curr_onset + abs(trl_mod_table.DistTrigger(triali)) + trl_mod_table.feedbackOnset(triali) + fsample -fsample/2];
                
    end
    
    AllSampleInfo(subji).trl_fb = trl;
    
    clear trl
    
    %% Remove Stop Trials
    
    % get index for stop trials to remove
    IdxStopTrl = find(trl_mod_table.TrlType == 2 | isnan(trl_mod_table.TrlType)); % trialtype 2 = stop trial

    % remove stop trials from table
    tmptable               = trl_mod_table;
    tmptable(IdxStopTrl,:) = [];
    
    AllSampleInfo(subji).trl_probe(IdxStopTrl,:,:)  = [];
    AllSampleInfo(subji).trl_bp(IdxStopTrl,:,:)     = [];
    AllSampleInfo(subji).trl_hll(IdxStopTrl,:,:)    = [];
    AllSampleInfo(subji).trl_br(IdxStopTrl,:,:)     = [];
    AllSampleInfo(subji).trl_fb(IdxStopTrl,:,:)     = [];
    
    % redefine trial to probe onset
    cfg                     = [];
    cfg.trl                 = AllSampleInfo(subji).trl_probe;
    hfa_mod.probe           = ft_redefinetrial(cfg, hfa);
    hfa_mod.probe.mod_label = mod_label;  
    
    % redefine trial to button press
    cfg                     = [];
    cfg.trl                 = AllSampleInfo(subji).trl_bp;
    hfa_mod.bp              = ft_redefinetrial(cfg, hfa);
    hfa_mod.bp.mod_label    = mod_label;  
   
    % redefine trial to hit lower limit
    cfg                     = [];
    cfg.trl                 = AllSampleInfo(subji).trl_hll;
    hfa_mod.hll             = ft_redefinetrial(cfg, hfa);
    hfa_mod.hll.mod_label   = mod_label;  
    
    % redefine trial to button release
    cfg                     = [];
    cfg.trl                 = AllSampleInfo(subji).trl_br;
    hfa_mod.br              = ft_redefinetrial(cfg, hfa);
    hfa_mod.br.mod_label    = mod_label;  
    
    % redefine trial to feedback
    cfg                     = [];
    cfg.trl                 = AllSampleInfo(subji).trl_fb;
    hfa_mod.fb              = ft_redefinetrial(cfg, hfa);
    hfa_mod.fb.mod_label    = mod_label;  
    
    clear hfa
    
    %% Step 3: Sliding Window ANOVA
    
    % loop over conditions
    IdxData  = {'Probe', 'BP', 'HLL', 'BR', 'FB'}; 
    AllData  = {hfa_mod.probe, hfa_mod.bp, hfa_mod.hll, hfa_mod.br, hfa_mod.fb};
    
    % critical f-value for thresholding
    CriticalF_Prediction = finv(0.95,2,size(tmptable,1)); % critical f-value for factor prediction (df1 = 2 (0%, 25%, 75% = 3 - 1), df2 = # trials)
    min_pc               = 0.1;   % 10% of consecutive samples (Helfrich - Neuron)
    
    factors  = {'Prediction'}; 
        
    %% Loop over conditions
    
    for condi = 1:numel(IdxData)
        
        % initialize storage for significant timepoints per channel
        sigtimepoints = cell(numel(hfa_mod.probe.label),1);

        % time range 
        trange   = [AllData{condi}.time(1) + smoothwind AllData{condi}.time(end) - smoothwind];
        twin     = round(smoothwind * fsample); % define smoothing window in samples
        
        % initialize omega square structure
        
        w2.time         = AllData{condi}.time(nearest(AllData{condi}.time,trange(1)):nearest(AllData{condi}.time,trange(2)));
        w2.label        = AllData{condi}.label;
        w2.trial        = NaN(length(w2.label), length(w2.time)); % initialize tiemcourse for exp. variance
        w2.pval         = NaN(length(w2.label), length(w2.time)); % initialize pval timecourse
        w2.fstat        = NaN(length(w2.label), length(w2.time)); % initialize f statistic timecourse
        w2.dimord       = 'rpt_chan_time';
        w2.trialinfo    = factors;
        w2.condition    = IdxData(condi);

        % loop across channels
        for chani = 1:length(w2.label)
            
            fprintf('Subject %d/%d, Channel %d/%d, Condition %d/%d\n', subji, numel(subjects), chani, length(AllData{condi}.label), condi, numel(IdxData));
            
            t = 1; % counter for timepoint
        
            for tidx = nearest(AllData{condi}.time,trange(1)):nearest(AllData{condi}.time,trange(2))

                % do the test with factor prediction
                [p, table] = anovan(squeeze(nanmean(AllData{condi}.trial(:,chani,tidx-twin:tidx+twin),3)), ...
                    {tmptable.LikehoodStop}, 'model', 'linear', 'sstype',2, ...
                    'varnames', {'Prediction'}, 'display', 'off');

                % calc w2 [ω2 = (SSeffect – (dfeffect)(MSerror)) / MSerror + SStotal]
                w2.trial(chani,t) = (table{2,2} - (table{2,3} * table{3,5})) / (table{3,5} + table{4,2});

                % get significance
                w2.pval(chani,t)  = p; 
                w2.fstat(chani,t) = table{2,6};
                          
                % trial keeping
                clear p table
                t = t+1;

            end %tidx
            
            % thresholding
            fstat = w2.fstat(chani,:);
            
            fstat(fstat < CriticalF_Prediction) = 0;

            % search for cluster for factor prediction

            island_predict = bwconncomp(fstat);
            
            sigtime_predict = [];
            for jj = 1:island_predict.NumObjects % loop through non-zero traces
                
                % check if any continous non-zero element exists and check whether the length exceeds the predefined minimum 
                if condi == 1 || condi == 5 && any(w2.time(island_predict.PixelIdxList{jj})) > 0
                    
                    % make sure that if locked to probe onset or button press any
                    % activity after event of interest (probe, FB) is above thres
                    if numel(island_predict.PixelIdxList{jj}) > ceil(sum(w2.time > 0) * min_pc)
                        
                        sigtime_predict = cat(1, sigtime_predict, island_predict.PixelIdxList{jj});
                        
                    end
                    
                elseif condi == 2 || condi == 3 || condi == 4
                    
                    if numel(island_predict.PixelIdxList{jj}) > ceil(length(w2.time) * min_pc)
                        
                        sigtime_predict = cat(1, sigtime_predict, island_predict.PixelIdxList{jj});
                        
                    end
                    
                end
                
            end
    
            sigtimepoints{chani,1} = sigtime_predict;
            
            clear fstat island_predict
                       
        end % chani
        
        w2.t_sig_threstime = sigtimepoints;
        AllData{condi}.w2  = w2;
                
        clear w2 sigtimepoints
               
    end % condi
    
    %% Pass data to hfa_mod

    hfa_mod.probe = AllData{strcmp(IdxData, 'Probe')};
    hfa_mod.bp    = AllData{strcmp(IdxData, 'BP')};
    hfa_mod.hll   = AllData{strcmp(IdxData, 'HLL')};
    hfa_mod.br    = AllData{strcmp(IdxData, 'BR')};
    hfa_mod.fb    = AllData{strcmp(IdxData, 'FB')};
    clear AllData
    
    %% Step 4: Find all active electrodes
    
    % probe
    hfa_mod.probe.w2.idxchan_sigthrestime = find(~cellfun(@isempty, hfa_mod.probe.w2.t_sig_threstime)); % thresholded timepoints using a predefined min time
    
    % button press
    hfa_mod.bp.w2.idxchan_sigthrestime    = find(~cellfun(@isempty, hfa_mod.bp.w2.t_sig_threstime));    
    
    % hit lower limt
    hfa_mod.hll.w2.idxchan_sigthrestime   = find(~cellfun(@isempty, hfa_mod.hll.w2.t_sig_threstime));
    
    % button release
    hfa_mod.br.w2.idxchan_sigthrestime    = find(~cellfun(@isempty, hfa_mod.br.w2.t_sig_threstime)); 
    
    % feedback
    hfa_mod.fb.w2.idxchan_sigthrestime    = find(~cellfun(@isempty, hfa_mod.fb.w2.t_sig_threstime));  
    
    %% Bring everything together 
   
    hfa_all(subji).data = hfa_mod;
    
    clear hfa_mod
    
    %% Step 5: Extract Electrode Positions for Active Elecs 
        
    % loop over conditions
    IdxData  = {'Probe', 'BP', 'HLL', 'BR', 'FB'}; 
    AllData  = {hfa_all(subji).data.probe, hfa_all(subji).data.bp, hfa_all(subji).data.hll, hfa_all(subji).data.br, hfa_all(subji).data.fb};
   
    for condi = 1:numel(AllData)
        
        fprintf('Extracting Elec Location %d. Condition\n', condi);
        
        % current condition
        tmpdata = AllData{condi};
        
        if ~isempty(tmpdata.w2.idxchan_sigthrestime)
            
            % modified electrode
            mod_elec = elec;

            for chani = 1:numel(tmpdata.w2.idxchan_sigthrestime)

                % find corresponding label in elec structure (needed for
                % plotting on cortical surface)
                tmplabel = tmpdata.mod_label{tmpdata.w2.idxchan_sigthrestime(chani)};
                
                AllActiveElecs.chanidx(chani) = find(strcmpi(elec.label, tmplabel));

                if isfield(elec, 'label')
                    AllActiveElecs.elec.label(chani)       = elec.label(logical(strcmpi(elec.label, tmplabel)));
                end
                if isfield(elec, 'chanpos')
                    AllActiveElecs.elec.chanpos(chani,:)   = elec.chanpos(logical(strcmpi(elec.label, tmplabel)),:);
                end
                if isfield(elec, 'tra')
                    AllActiveElecs.elec.tra(chani)         = elec.tra(logical(strcmpi(elec.label, tmplabel)), ...
                        logical(strcmpi(elec.label, tmplabel)));
                end
                if isfield(elec, 'chantype')
                    AllActiveElecs.elec.chantype(chani)    = elec.chantype(logical(strcmpi(elec.label, tmplabel)));
                end
                if isfield(elec, 'elecpos')
                    AllActiveElecs.elec.elecpos(chani,:)   = elec.elecpos(logical(strcmpi(elec.label, tmplabel)),:);
                end
                if isfield(elec, 'chanside')
                    AllActiveElecs.elec.chanside(chani)    = elec.chanside(logical(strcmpi(elec.label, tmplabel)));
                end

            end

        else

             AllActiveElecs.chanidx = [];
             AllActiveElecs.elec    = [];
             
        end

        if strcmp(IdxData{condi}, 'Probe')
            
            hfa_all(subji).data.probe.w2.active_elecs = AllActiveElecs;
                        
        elseif strcmp(IdxData{condi}, 'BP')
            
            hfa_all(subji).data.bp.w2.active_elecs    = AllActiveElecs;
            
        elseif strcmp(IdxData{condi}, 'HLL')
            
            hfa_all(subji).data.hll.w2.active_elecs   = AllActiveElecs;
            
        elseif strcmp(IdxData{condi}, 'BR')
            
            hfa_all(subji).data.br.w2.active_elecs    = AllActiveElecs;
            
        elseif strcmp(IdxData{condi}, 'FB')
            
            hfa_all(subji).data.fb.w2.active_elecs    = AllActiveElecs;
            
        end
        
        clear AllActiveElecs
                
    end % condi
    
    %% Load MNI electrode file

    cd([path_elec filesep subID filesep 'Electrodes']);

    elecfile = dir('*elec_mni_frv.mat');
    
    load(elecfile.name);
       
    % Extract MNI coordinates for active channels
    
    % probe
    chanidx = hfa_all(subji).data.probe.w2.active_elecs.chanidx;
    hfa_all(subji).data.probe.w2.active_elecs.mni_elec.label = elec_mni_frv.label(chanidx)';
    hfa_all(subji).data.probe.w2.active_elecs.mni_elec.chanpos = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        hfa_all(subji).data.probe.w2.active_elecs.mni_elec.chantype = elec_mni_frv.chantype(chanidx);
    end
    hfa_all(subji).data.probe.w2.active_elecs.mni_elec.elecpos = elec_mni_frv.elecpos(chanidx,:);
    
    % bp
    chanidx = hfa_all(subji).data.bp.w2.active_elecs.chanidx;
    hfa_all(subji).data.bp.w2.active_elecs.mni_elec.label = elec_mni_frv.label(chanidx)';
    hfa_all(subji).data.bp.w2.active_elecs.mni_elec.chanpos = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        hfa_all(subji).data.bp.w2.active_elecs.mni_elec.chantype = elec_mni_frv.chantype(chanidx);
    end
    hfa_all(subji).data.bp.w2.active_elecs.mni_elec.elecpos = elec_mni_frv.elecpos(chanidx,:);
 
    % hll
    chanidx = hfa_all(subji).data.hll.w2.active_elecs.chanidx;
    hfa_all(subji).data.hll.w2.active_elecs.mni_elec.label = elec_mni_frv.label(chanidx)';
    hfa_all(subji).data.hll.w2.active_elecs.mni_elec.chanpos = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        hfa_all(subji).data.hll.w2.active_elecs.mni_elec.chantype = elec_mni_frv.chantype(chanidx);
    end
    hfa_all(subji).data.hll.w2.active_elecs.mni_elec.elecpos = elec_mni_frv.elecpos(chanidx,:);
    
    % br
    chanidx = hfa_all(subji).data.br.w2.active_elecs.chanidx;
    hfa_all(subji).data.br.w2.active_elecs.mni_elec.label = elec_mni_frv.label(chanidx)';
    hfa_all(subji).data.br.w2.active_elecs.mni_elec.chanpos = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        hfa_all(subji).data.br.w2.active_elecs.mni_elec.chantype = elec_mni_frv.chantype(chanidx);
    end
    hfa_all(subji).data.br.w2.active_elecs.mni_elec.elecpos = elec_mni_frv.elecpos(chanidx,:);
    
    % fb
    chanidx = hfa_all(subji).data.fb.w2.active_elecs.chanidx;
    hfa_all(subji).data.fb.w2.active_elecs.mni_elec.label = elec_mni_frv.label(chanidx)';
    hfa_all(subji).data.fb.w2.active_elecs.mni_elec.chanpos = elec_mni_frv.chanpos(chanidx,:);
    if isfield(elec_mni_frv, 'chantype')
        hfa_all(subji).data.fb.w2.active_elecs.mni_elec.chantype = elec_mni_frv.chantype(chanidx);
    end
    hfa_all(subji).data.fb.w2.active_elecs.mni_elec.elecpos = elec_mni_frv.elecpos(chanidx,:);
        
    %% Add surface file

    hfa_all(subji).cortex_rh = cortex_rh; clear cortex_rh
    hfa_all(subji).cortex_lh = cortex_lh; clear cortex_lh
    
    anovadata = hfa_all(subji).data;
    
    savingname = [path_in filesep subID filesep 'Data' filesep 'HFA' ...
        filesep strcat(subID, '_TaskActive_ExpVariance')];
    
    save(savingname, 'anovadata', '-v7.3');
    
    clear anovadata
    
end % subji


