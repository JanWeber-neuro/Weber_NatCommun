function [all_plv, all_iplv] = channelwise_plv(data_1, data_2, foi, nseed, Niterations, ...
    fsample, padding, time2use, method)

% #########################################################################
% ---- INPUT ---- %

% 1. data_1         = seed data
% 2. data_2         = data that is used to compute PLV to seed data
% 3. foi            = frequency of interest
% 4. nseed          = number of seed electrodes
% 5. Niterations    = number of iterations for permutation testing
% 6. fsample        = sampling frequency
% 7. padding        = "1" or "0" to indicate whether extra data was added
%                   for filtering purposes. The added "wings" will be removed after
%                   filtering.
% 8. time2use       = relevant timepoints to be used after filtering
% 9. method         = 1) "cutdata": cuts data from the seed electrode at a
%                        random point and puts it in the reverse order (benefit: does not destroy
%                        the temporal structure of the data)
%                     2) "trialshuffling": shuffles trials from seed and
%                     target randomly. Benefit: does not destroy the
%                     temporal structure of the data at all, but might be a
%                     bit conservative if the trials are similar.
% #########################################################################

% #########################################################################
% ---- OUTPUT ---- %

% 1. all_plv  = structure that contains different fields, such as the raw
% PLVs, the surrogate PLVs and the z-normalized PLVs based on the surrogate
% distribution

% 2. all_iplv = same as "all_plv", just for the imaginary PLV

% #########################################################################


%%

% assign memory to output variables
all_plv  = struct();
all_iplv = struct();

%%

for Ifreq = 1:size(foi,1) % loop over frequencies

    fprintf('progress: %.d%%\n', round(Ifreq / size(foi,1)*100))
    
    % current frequency range for filtering
    tmpFreq = foi(Ifreq,:);

    % loop over seed channel
    for Ichan_one = 1:nseed

        for Ichan_two = 1:size(data_2,2)

            tmpTrial_1 = squeeze(data_1(:,Ichan_one,:));   % seed data (output station)
            tmpTrial_2 = squeeze(data_2(:,Ichan_two,:));   % sender data 

            % check dimension (make sure it is trials x time)
            if size(tmpTrial_1,1) == size(data_1,3)
                tmpTrial_1 = tmpTrial_1';
            end
            if size(tmpTrial_2,1) == size(data_2,3)
                tmpTrial_2 = tmpTrial_2';
            end

            % filter data
            tmpFilt_1 = ft_preproc_bandpassfilter(tmpTrial_1, ...
            fsample, tmpFreq, [], [], [], [], 'hanning'); 

            tmpFilt_2 = ft_preproc_bandpassfilter(tmpTrial_2, ...
            fsample, tmpFreq, [], [], [], [], 'hanning');

            % apply hilbert transform (computation along columns)
            tmpHilbert_1 = hilbert(tmpFilt_1')';
            tmpHilbert_2 = hilbert(tmpFilt_2')';
            
            %------- optional -------%
            % cut-off the wings if some padding was used to prevent edge artifacts
            if padding
                tmpHilbert_1 = tmpHilbert_1(:, time2use(1):time2use(end));
                tmpHilbert_2 = tmpHilbert_2(:, time2use(1):time2use(end));
            end

            % loop over trials 
            for Itrial = 1:size(tmpHilbert_1,1)

                % compute phase locking value
                [plv,iplv] = rh_plv(tmpHilbert_1(Itrial,:), tmpHilbert_2(Itrial,:));

                % store phase locking values
                all_plv.tmp_plv(Ichan_one, Ichan_two, Itrial, Ifreq)      = plv;
                all_iplv.tmp_imagplv(Ichan_one, Ichan_two, Itrial, Ifreq) = iplv;

            end

            % -----------------------------------
            % create null distribution
            % -----------------------------------

            if strcmp(method, 'trialshuffle')

                flag = 0;
                while ~flag % loop until there is no overlap in trials

                    % generate N random trial numbers 
                    rndTrials_1 = randsample(size(data_1,1), Niterations, true);
                    rndTrials_2 = randsample(size(data_2,1), Niterations, true);

                    % find index where samples overlap
                    idxoverlap  = find((rndTrials_1 - rndTrials_2) == 0);

                    if isempty(idxoverlap)
                        flag = 1;
                    end

                end

                % create null distribution
                for Iperm = 1:Niterations
                    % compute phase locking value
                    [plv,iplv] = rh_plv(tmpHilbert_1(rndTrials_1(Iperm),:), tmpHilbert_2(rndTrials_2(Iperm),:));

                    % store phase locking values
                    all_plv.tmp_plv_surro(Iperm, Ichan_one, Ichan_two, Ifreq)      = plv;
                    all_iplv.tmp_imagplv_surro(Iperm, Ichan_one, Ichan_two, Ifreq) = iplv;

                end
                
            elseif strcmp(method, 'cutdata')
                
                rndTrials = randsample(size(data_1,1), Niterations, true);
                
                % create null distribution
                for Iperm = 1:Niterations

                    % get a random time point where the data is cut
                    rndidx = randsample(length(tmpHilbert_1),1);
                    tmpHilbert_1 = [tmpHilbert_1(:, rndidx:end) tmpHilbert_1(:, 1:rndidx-1)]; 

                    % compute phase locking value
                    [plv,iplv] = rh_plv(tmpHilbert_1(rndTrials(Iperm),:), tmpHilbert_2(rndTrials(Iperm),:));

                    % store phase locking values
                    all_plv.tmp_plv_surro(Iperm, Ichan_one, Ichan_two, Ifreq)      = plv;
                    all_iplv.tmp_imagplv_surro(Iperm, Ichan_one, Ichan_two, Ifreq) = iplv;


                end 
                
            end

        end %Ichan_two

    end %Ichan_one

end %Ifreq

% ---------------------------------------------------------------
% z-score true data based on surrogate distribution
% ---------------------------------------------------------------

% preallocate memory for z-normalized PLV
all_plv.plv_z   = NaN(size(all_plv.tmp_plv));
all_iplv.iplv_z = NaN(size(all_iplv.tmp_imagplv));

for Ichan_one = 1:size(all_plv.tmp_plv_surro,2) % loop over channel pairs

    for Ichan_two = 1:size(all_plv.tmp_plv_surro,3)

        for Ifreq = 1:size(all_plv.tmp_plv_surro,4) % loop over frequency

            % mean and standard deviation for PLV surrogate distribution
            mu_plv     = mean(squeeze(all_plv.tmp_plv_surro(:,Ichan_one,Ichan_two,Ifreq)));
            sigma_plv  = std(squeeze(all_plv.tmp_plv_surro(:,Ichan_one,Ichan_two,Ifreq)));

            % mean and standard deviation for imaginary PLV surrogate distribution
            mu_iplv    = mean(squeeze(all_iplv.tmp_imagplv_surro(:,Ichan_one,Ichan_two,Ifreq)));
            sigma_iplv = std(squeeze(all_iplv.tmp_imagplv_surro(:,Ichan_one,Ichan_two,Ifreq)));

            for Itrial = 1:size(all_plv.tmp_plv,3) % z-score each trial

                % z-score PLV
                tmp   = squeeze(all_plv.tmp_plv(Ichan_one, Ichan_two, Itrial, Ifreq));
                all_plv.plv_z(Ichan_one,Ichan_two,Itrial,Ifreq) = (tmp - mu_plv) / sigma_plv; clear tmp

                % z-score ImagPLV
                tmp   = squeeze(all_iplv.tmp_imagplv(Ichan_one, Ichan_two, Itrial, Ifreq));
                all_iplv.iplv_z(Ichan_one,Ichan_two,Itrial,Ifreq) = (tmp - mu_iplv) / sigma_iplv; clear tmp

            end

        end %Ifreq

    end %Ichan_two

end %Ichan_one




