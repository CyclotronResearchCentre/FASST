function Do = crc_par_iica(Di)

% Pulse artefact rejection with a iICA method.
% ICA is applied multiple times on the same stretch of data:
% at each repetition the artefact-related component are automatically
% selected and removed from the data. The N-different-times corrected data
% are then averaged.
%_______________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by C. Phillips, 2011.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$


crcdef = crc_get_defaults('par');
prefiICA = crcdef.iicapref;

% Recover Peaks
Peaks = Di.CRC.EKGPeaks;

% Recover number of iterations
nit = Di.cache.nit;

% Initializaing the cleaned EEG datafile
fn_dat = fullfile(spm_str_manip(fnamedat(Di),'H'), ...
    [ prefiICA, spm_str_manip(fnamedat(Di),'t')]);

if exist(fullfile(path(Di), fn_dat),'file')
    % This check is because I want append data in the new file
    % I have to make sure that it does not exist yet.
    delete(fullfile(path(Di), fn_dat))
    disp(['!!! WARNING : EXISTING ',fn_dat,' FILE OVERWRITTEN !!!' ])
end

% create empty object
Do = clone(Di,fn_dat);


% Creating a mean BCG template
disp('................................')
disp('.... Creating BCG template .....')
disp('................................')

TEMPLATE = zeros(1,Di.nsamples);
for ipoint = 1:size(Peaks,2)
    % use 1/5 sec
    TEMPLATE(Peaks(ipoint):Peaks(ipoint)+Di.fsample/5) = 1;
end

% Split good EEG channel, to correct, and others
l_eeg = Di.cache.l2cor;
l_other = 1:Di.nchannels;
l_other(l_eeg) = [];

% Computing the bounds of the sweeps that will be processed 1 by 1
beatstep = ceil(3/2*length(l_eeg).^2);
% x-beat step to scan the data
beatsteplength = 2*beatstep;
% length of EEG data processed by steps (in number of bins)
% #points ~= 3*N.^2
if beatsteplength>Di.nsamples
    disp('******************************************')
    disp('**** WARNING: Not enough data points  ****')
    disp('****          to proceed with ICA !   ****')
    disp('**** Still calculating but can''t     ****')
    disp('**** guarantee the results...         ****')
    disp('******************************************')
    % number of steps needed and bounds of chunks
    nbeatstep = 1;
    CHKDTB = [1 Di.nsamples];
    bounds_mean = 1:Di.nsamples;
    bounds_keep = Di.nsamples;
else
    % number of steps needed and bounds of chunks
    nbeatstep = floor(size(TEMPLATE,2)/beatstep);
    CHKDTB = [((0:nbeatstep-1)*beatstep)'+1 ...
        ((0:nbeatstep-1)*beatstep)'+beatsteplength];
    CHKDTB(end) = Di.nsamples;
    bounds_mean = 1:beatsteplength/2;
    bounds_keep = (beatsteplength/2+1):beatsteplength ;
end

% Computing ICA
keepEEG = [];
meanEEG = [];
for jbeatstep = 1:nbeatstep
    data = Di(l_eeg, CHKDTB(jbeatstep,1):CHKDTB(jbeatstep,2));
    tmpdata = reshape( data, size(data,1), size(data,2));
    tmpdata = detrend(tmpdata','constant')'; % zero mean
    tmprank = rank(tmpdata(:,1:floor(size(tmpdata,2)/2)));
    AV_cleanedEEG = zeros(size(data,1),size(data,2), nit);

    % Creating a binary template
    % The main BCG artefact should occur within 500 ms after the
    % systole starts.
    tt = TEMPLATE(1,CHKDTB(jbeatstep,1):CHKDTB(jbeatstep,2));
    % '-2' because I get rid of the EOG and EKG
    epochlength = round(prctile(diff(Peaks),99)*500/1000);
    stag = 1;
    nstag = epochlength/stag; % nombre de pas pour couvrir RR
    tt = [tt;zeros(nstag-1,size(tt,2))];
    for istag = 2:nstag
        tt(istag,:) = [zeros(1,stag) tt(istag-1,1:end-stag)];
    end
    % used later on in the code, for the correlation bit
    d_tt = detrend(tt','constant');
    ssd_tt = sqrt(sum(d_tt.^2)');

    RRbound = max(diff(Peaks));
    RRmean = mean(diff(Peaks));

    for jnit = 1:nit % nit iterations to get a mean clean set
        clc
        disp('..............................................')
        disp(['.... File ' num2str(ifile) '/' ...
            num2str(size(files,1)) '; step ' ...
            num2str(jbeatstep) '/' num2str(nbeatstep) ...
            '; itération  '  num2str(jnit) ' ....'])
        disp('.... Computing the independent components ....')
        disp('..............................................')

        % if there is no rank reduction then, don't use it!
        % This avoids some large matrix multiplication
        if tmprank == size(tmpdata,1)
            [icaweights,icasphere] = runica( tmpdata, ...
                'lrate', 0.001 ,'verbose','off' );
        else
            [icaweights,icasphere] = runica( tmpdata, ...
                'lrate', 0.001, 'pca', tmprank ,'verbose','off' );
        end
        % use the detrended data.
        icaact    = (icaweights*icasphere)*tmpdata;
        invweights = pinv(icaweights*icasphere);

        % Choosing the relevant ic automatically
        disp('......................................................')
        disp('.... Automatically identifying the BCG components ....')
        disp('......................................................')

        % Initializing correlation variables
        CC = zeros(1,size(icaact,1));

        for icomp = 1:size(icaact,1)
            % Creating the binary sweep ttt,  to be compared with the
            % binary template tt

            dicaact = icaact(icomp,:);
            % as data are detrended, mean=0, just look at
            % the variances.
            QRSidentifier = find(dicaact > 2*std(dicaact));
            diff_QRS = diff(QRSidentifier);
            ica_QRS = diff_QRS(diff_QRS>60/180*Di.fsample);
            % number of point between RR if 180 bpm
            %HR = 60sec/min; 200 bpm; Di.Radcpoints/sec
            if ~isempty(ica_QRS)
                meanicaQRS = prctile(ica_QRS,50);
                % If this icomp is a bcg_ica, there should be
                % an event nearly at each beat.
                % We check that there is no gap of more than
                % 5 max(RR intervals) between the elements of ttt
                % and that the mean RR is ~=  mean RR+- 100 points
                if max(ica_QRS) < 5*RRbound && ...
                        (meanicaQRS < RRmean+100 && meanicaQRS > RRmean-100 )
                    % create ttt directly for the list of
                    % useful peak detections.
                    ttt = zeros(size(tt,2),1);
                    for ipoint = find(diff_QRS>60/180*Di.fsample)
                        ttt(QRSidentifier(ipoint):QRSidentifier(ipoint)+Di.fsample/5) = 1;
                    end

                    % make things much faster by avoiding
                    % corrcoef routine and calculating directly the
                    % coef of interest, for all the stags
                    % considered at once!
                    if any(ttt)
                        d_ttt = ttt-mean(ttt);
                        cc = (d_tt'*d_ttt)./(ssd_tt*sqrt(d_ttt'*d_ttt));
                        CC(icomp) = prctile(abs(cc),99);
                    end
                end
            end
        end
        icaoi_index = find(CC> mean(CC) + std(CC)) ;

        % skewness of ica_bcg should be high
        % calculate skewness only once
        sk_icaact = skewness(icaact(icaoi_index,:)');
        if any(sk_icaact< 0.33) %  ~= prctile(SK,5)
            disp('WARNING SOME ICAs were discarded')
            icaoi_index(sk_icaact< 0.33) = [];
            % 0.3 is the cutoff skewness computed empirically
        end

        [cleanedEEG] = compvar(data, icaact, ...
            invweights, setdiff(1:size(icaact,1),icaoi_index));
        cleanedEEG = detrend(cleanedEEG','constant');
        % zero mean, to avoid offset avery beatstep
        cleanedEEG = cleanedEEG';
        AV_cleanedEEG(:,:,jnit) = cleanedEEG;
    end % end jnit
    p_cleanedEEG = mean(AV_cleanedEEG,3);
    % add the other channels, in the same order !!!
    cleanedEEG = zeros(Di.nchannels,size(p_cleanedEEG,2));
    cleanedEEG(l_eeg,:) = p_cleanedEEG;
    cleanedEEG(l_other,:) = Di( l_other, ...
        CHKDTB(jbeatstep,1):CHKDTB(jbeatstep,2),1);

    % Mean, Keep, Write
    % tomean : the first half of the epoch;
    % keepEEG = the second half of an epoch
    % meanEEG = the mean of the present tomean and the previous keepEEG
    if diff(CHKDTB(jbeatstep,:)) < beatstep
        % if last cleanEEG smaller than beatstep
        tomeanEEG = cleanedEEG;
        meanEEG = [(tomeanEEG+keepEEG(:,1:size(tomeanEEG,2)))/2 ...
            keepEEG(:,size(tomeanEEG,2)+1:end)];
    elseif diff(CHKDTB(jbeatstep,:)) >= beatstep
        % if cleanEEG is larger than beatstep but possibly
        % smaller than 2*beatstep
        if jbeatstep == 1
            % if first pass, no keepEEG available, mean = tomean
            tomeanEEG = cleanedEEG(:,bounds_mean);
            meanEEG = tomeanEEG;
            if size(cleanedEEG,2) < beatsteplength
                % in case it is a very short EEG bit
                keepEEG = cleanedEEG(:,bounds_keep(1):end);
            else
                keepEEG = cleanedEEG(:,bounds_keep);
            end
        elseif jbeatstep > 1 && jbeatstep < nbeatstep
            % after 1st pass, mean tomean with keepEEG of
            % the previous pass
            meanEEG = (keepEEG+cleanedEEG(:,bounds_mean))/2;
            keepEEG = cleanedEEG(:,bounds_keep);
        elseif jbeatstep == nbeatstep
            % if last pass, mean tomean and previous keepEEG
            % and add the present (remaining) keepEEG
            bounds_keep = (beatstep+1):size(cleanedEEG,2);
            meanEEG = [(keepEEG+cleanedEEG(:,bounds_mean))/2 ...
                cleanedEEG(:,bounds_keep)];
        end
    end
    % detrending to avoid steps in the reconstructed EEG.
    % This possibly induces a change in spectrum a period =
    % beatstep (11532)/250 = 46 s , i.e. 0.021 Hz
    meanEEG = detrend(meanEEG','constant');
    meanEEG = meanEEG';
    % writing the first part of cleanedEEG == meanEEG
    %                 CLEANEEG.data.scale = ones(size(eeg.data, 1), 1);
    Do(:,CHKDTB(jbeatstep,1):CHKDTB(jbeatstep,2)) = meanEEG;
%     fwrite(fpd_clean, meanEEG, 'float32');
end % jbeatstep

return

