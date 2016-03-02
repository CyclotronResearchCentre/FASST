function crc_SP_detect(handles)

% Automatically identifies spindles in periods of NREM sleep, excluding 
% periods of arousals and artefacts.
% Runs only on pre-referenced data which have been fully scored using FASST.
%
% Detection based on amplitude criterion as in Molle et al. (2002)
% Threshold for spindle detection varies for each channel (percentile 90)
% The routine requires a visual confirmation of detected spindles
% Data are subsequently epoched and fed into a time frequency analysis
% (resulting files saved on the disk)
%
% Results are saved in D, under 2 different fields
% - in D.events:
% 555 : anterior spindles, referred to as Ant_SP
% 666 : posterior spindles, referred to as Post_SP
% 777 : undetermined spindles, referred to as Spindles
% - in D.CRC.spindles:
% bounds        starting and ending time of each spindle (in time points)
% duration      spindle duration (ms)
% amplitude     spindle amplitude peak to peak (microV)
% maxelectrode  electrode with the maximum amplitude
% index_slowsp  index of slow spindles
% index_fastsp  index of fast spindles
% sp_frequency  spindle  maximal frequency

%__________________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Written by P.Maquet (2005) and adapted for FASST by J.Schrouff (2010).
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$


%--------------------------------------------------------------------------
%-------------------- Loading file and checking the fields ----------------
%--------------------------------------------------------------------------

%Select and load file
if ~nargin
    file = spm_select(1, 'mat', 'Select EEG file','' ,pwd,'.*');
    Dor = crc_eeg_load(file);
else
    Dor = crc_eeg_load(handles.fname);
end
close all
crcdef = crc_get_defaults('sp');

%checking that the file was rereferenced: it is assumed that people know
%they have to re-reference to the mean of the mastoids.
disp('---- Checking reference----')
histD=history(Dor);
cont=0;
conti=1;
for i=1:size(histD,2)
    if strcmpi(histD(i).fun,'spm_eeg_montage')
        cont=1;
    end
end
if ~cont && ~handles.reref
    disp('WARNING: no re-referencing found!-- please re-reference using spm8 (spm_eeg_montage)')
    conti=spm_input('Continue?',1,'y/n',[1,0],1);
end
if ~conti
    return
end

%checking that the file was scored using FASST
try
    Dor.CRC;
    if ~isfield(Dor.CRC,'score')
        disp('File not scored using the FASST Toolbox')
        disp('Please score before using the spindle detection tool')
        return
    end
catch
   disp('File not scored using the FASST Toolbox')
   disp('Please score before using the spindle detection tool')
   return
end



%selecting the scorer
scorer=handles.scorer;
%extract periods of interest
if handles.analyse==3
    if isfield(Dor.CRC, 'score') && size(Dor.CRC.score,1)>4
        if ~isempty(handles.stagesp)
            [D, TPOI]=crc_extractSW(Dor,handles.analyse,handles.stagesp,scorer,0);
        else
            [D, TPOI]=crc_extractSW(Dor,handles.analyse,crcdef.stagesp,scorer,0);
        end
    end
elseif handles.analyse==2
    [D,TPOI]=crc_extractSW(Dor,handles.analyse,crcdef.stagesp,scorer,0);
else
    [D,cleanSW]=crc_extractSW(Dor,handles.analyse, crcdef.stagesp,scorer,0);
    %get the part of the vector comprised between the begin and the
    %start times edited by the user
    belbound=find(cleanSW>=handles.Begpts);
    upbound=find(cleanSW<=handles.Endpts);
    TPOI=intersect(cleanSW(belbound),cleanSW(upbound));
end
stoppt=min(max(TPOI),nsamples(D));
TPOI=TPOI(find(TPOI<=stoppt));


% -------------------------------------------------------------------------
% ------------------------ Getting the electrodes of interest -------------
%--------------------------------------------------------------------------
%Electrodes of reference are Fz, Pz and Cz in the 10-20 system
elecoi = crcdef.elecoi; nelecoi = size(elecoi,2);
load('CRC_electrodes.mat')
%Get the 2D theoretical position of the reference electrodes
posref=zeros(nelecoi,2);
for iref = 1:nelecoi
    isref=strcmpi(elecoi{iref},names);
    posref(iref,:)=pos(find(isref),:);
end

%Get the 2D theoretical position of the file electrodes
eegchan=find(strcmpi(crcdef.type,chantype(D)));
labchan=chanlabels(D);
labchan=labchan(eegchan);

poselec=zeros(size(labchan,2),2);
for ielec = 1:size(labchan,2)
    if ~crcdef.usetheor
        try
            poselec(ielec,:)=(coor2D(D,ielec))';
        catch
            iselec=strcmpi(labchan{ielec},names);
            poselec(ielec,:)=pos(find(iselec),:);
        end
    else
        iselec=strcmpi(labchan{ielec},names);
        poselec(ielec,:)=pos(find(iselec),:);
    end
end

%Computes the distance between the reference and the electrodes to find
%which electrodes best fit the reference electrodes.
distpos=ones(size(labchan,2),nelecoi);
indelecoi=zeros(nelecoi,1);
front=[];
par=[];
fronti=[];
pari=[];
for iref = 1:nelecoi
    refpos=repmat(posref(iref,:),size(labchan,2),1);
    distpos(:,iref)=((refpos(:,1)-poselec(:,1)).^2+(refpos(:,2)-poselec(:,2)).^2).^0.5;
    [dum1,posmin]=min(distpos(:,iref));
    indelecoi(iref,1)=eegchan(posmin);
    if posref (iref,2)>=0.5
        front=[front, posmin];
        fronti=[fronti,iref];
    elseif posref (iref,2)<0.5
        par=[par, posmin];
        pari=[pari, iref];
    end
end

%check that, amongst the considered reference electrodes, none is marked as
%bad
bad=badchannels(D);
fbad=intersect(front,bad);
pbad=intersect(par,bad);
nbad=[];
if ~isempty(fbad)
    nbad=chanlabels(D,fbad);
elseif ~isempty(pbad)
    nbad=[nbad chanlabels(D,fbad)];
end
if ~isempty(nbad)
    disp(['One of the reference eletrode is bad. Please, edit the',...
        'crc_def.sp.elecoi (in the crc_defaults) and replace the',...
        'electrode nearest to: ',nbad])
    error('One of the reference electrode is bad')
end

%--------------------------------------------------------------------------
%-----------------  RMS method (Molle, 2002)  -----------------------------
%--------------------------------------------------------------------------
% pack
disp('..............................................')
disp('.... Computing RMS and detecting spindles ....')
disp('..............................................')

%Filter data between 11 and 18 Hz (default) whilst avoiding out-of-memory 
%errors for long recordings at high sampling rates

if ~isempty(handles.highfc)&& ~isempty(handles.lowfc)
    fqcut=[handles.highfc, handles.lowfc];
else
    fqcut=[crcdef.highfc, crcdef.lowfc];
end

fhc = fqcut(1)/(D.fsample/3);
flc = fqcut(2)/(D.fsample/3);
order=4;
[b1,a1] = butter(order,[fhc, flc],'pass');

% out=daqmem;
% memsz=out.AvailPhys/5;
memsz =crc_get_defaults('mem.sp_maxmemload');
szd=D.nsamples*8;
numblocks=ceil(szd/memsz);
data_f=zeros(nelecoi,D.nsamples,numblocks);
scales=crc_scales(D,indelecoi);
h = waitbar(0,'Please wait...');
for ielec = 1:nelecoi
    if numblocks==1
        data_f(ielec,:,1)  = filtfilt(b1, a1,D(indelecoi(ielec),:)/scales(ielec));
    else
        disp('Data this size are not handled yet by spindle detection')
        disp('Please chunk your data into smaller files')
        return
    end
    string = ['Please wait... ' num2str(100*(ielec/nelecoi)) ' %'];
    waitbar((ielec)/nelecoi,h,string)
end    
close(h)
%--------------------------------------------------------------------------
%Compute the threshold of detection if not provided by user
if ~crcdef.threshold
    %get the periods of stage 2
    flagthresh=0;
    if handles.analyse==3 || handles.analyse==2
        if isfield(Dor.CRC, 'score') && size(Dor.CRC.score,1)>4
            [Dthresh, TPOIthresh]=crc_extractSW(Dor,handles.analyse,crcdef.stagethresh,scorer,0);
        else
            disp('Warning: threshold computed on all data instead of stage 2 only')
            Dthresh=D;
            TPOIthresh=TPOI;
            flagthresh=1;
        end
    else
        [Dthresh,clthresh]=crc_extractSW(D,handles.analyse, crcdef.stagethresh,scorer,0);
        %get the part of the vector comprised between the begin and the
        %start times edited by the user
        belbound=find(clthresh>=handles.Begpts);
        upbound=find(clthresh<=handles.Endpts);
        TPOIthresh=intersect(belbound,upbound);
    end
    stoppt=min(max(TPOIthresh),nsamples(Dthresh));
    TPOIthresh=TPOIthresh(find(TPOIthresh<=stoppt));
    %filter the stage 2 data
    if ~flagthresh
        memsz =crc_get_defaults('mem.sp_maxmemload');
        szd=Dthresh.nsamples*8;
        numblocks=ceil(szd/memsz);
        data_thresh=zeros(nelecoi,Dthresh.nsamples,numblocks);
        h = waitbar(0,'Please wait...');
        for ielec = 1:nelecoi
            if numblocks==1
                data_thresh(ielec,:,1)  = filtfilt(b1, a1,Dthresh(indelecoi(ielec),:)/scales(ielec));
            else
                disp('Data this size are not handled yet by spindle detection')
                disp('Please chunk your data into smaller files')
                return
            end
            string = ['Please wait... ' num2str(100*(ielec/nelecoi)) ' %'];
            waitbar((ielec)/nelecoi,h,string)
        end    
        close(h)
    else
        data_thresh=data_f;
    end
    % find spindle detection threshold - allow different threshold for
    % different channels
%     w = rectwin(round(Dthresh.fsample/10)); % a boxcar window of .1s
%     data_rms = zeros(nelecoi,Dthresh.nsamples);
    spinthres = zeros(nelecoi,1);
    for ielec = 1:nelecoi
%         conv_rms = (conv(sqrt(data_thresh(ielec,:).^2),w))./size(w,1);
%         data_rms(ielec,:) = conv_rms(round(length(w)/2):round(length(w)/2)+Dthresh.nsamples-1);
        spinthres(ielec) = prctile(data_thresh(ielec,TPOIthresh),crcdef.prcthresh);
        disp(['Threshold of detection for channel ',labchan{indelecoi(ielec)},...
            ' : ',num2str(spinthres(ielec))])
    end
    clear Dthresh TPOIthresh data_thresh
else
    if numel(crcdef.threshold)==nelecoi
        spinthres=crcdef.threshold;
    elseif numel(crcdef.threshold)==1
        spinthres=crcdef.threshold*ones(1,nelecoi);
    else
        disp('Error: the number of thresholds should correspond to number of channels of interest or to 1')
        return
    end
end
%find spindles using the previously computed thresholds: combines the
%spindles detected on the different channels in a long one if they overlap
w = rectwin(round(D.fsample/10)); % a boxcar window of .1s
data_rms = zeros(nelecoi,D.nsamples);
drms=data_rms;
for ielec = 1:nelecoi
    conv_rms = (conv(sqrt(data_f(ielec,:).^2),w))./size(w,1);
    data_rms(ielec,:) = conv_rms(round(length(w)/2):round(length(w)/2)+D.nsamples-1);
    drms(ielec,find(data_rms(ielec,:)<spinthres(ielec))) = 0;
    drms(ielec,find(data_rms(ielec,:)>spinthres(ielec))) = 1;
end
darms=ceil(sum(drms)/nelecoi);
difrms = [diff(darms,1,2) 0]; % 1s = points of a specific spindle; >1 identifies jump to next spindle
indpn = find(difrms > 0)+1;
SpEnds = find(difrms < 0);
if length(indpn) == length(SpEnds)
    indpn = [indpn;SpEnds]';
elseif length(indpn) > length(SpEnds)
    indpn = [indpn(1:end-1);SpEnds]';
else
    indpn = [indpn;SpEnds(2:end)]';
end
jumpspind = find(diff(indpn,1,2) > D.fsample*crcdef.lengthsp); % 400 ms given that the rms does not detect the first and last oscillations of the spindle
indpn=indpn(jumpspind,:);
jumpspind = find(diff(indpn(:,1)) > D.fsample*crcdef.succsp); % successive spindles should be separed by 1000 ms 
Spindlebounds1 = indpn(jumpspind,:); % onset - offset - duration
Spindlebounds1 = sort(Spindlebounds1);

%--------------------------------------------------------------------------

% ordering and finding amplitude, identify the spindle channel with maximal amplitude
% ptp amplitude and corresponding electrode (electrode number coresponding
% to D.channels, for spindles in the TPOI interval
Spindleboundsok = [];
ibok=1;
for ibound = 1:size(Spindlebounds1, 1)
    if ~isempty(intersect(Spindlebounds1(ibound,1), TPOI)) && ...   % do not consider spindles starting in REM or artifacted periods
           numel(intersect(Spindlebounds1(ibound,1):Spindlebounds1(ibound,2),TPOI))>=D.fsample*crcdef.lengthsp % consider spindles starting at least 400ms before an artifact
        Spindleboundsok(ibok,1:3) = [Spindlebounds1(ibound,1) ...
             Spindlebounds1(ibound,2) ...
             diff([Spindlebounds1(ibound,1)  Spindlebounds1(ibound,2)],1,2)*1000/D.fsample];
        valmax = max(data_f(:,Spindlebounds1(ibound,1):Spindlebounds1(ibound,2)),[],2);
        valmin = min(data_f(:,Spindlebounds1(ibound,1):Spindlebounds1(ibound,2)),[],2);
        MaxAmp = abs(valmax)+abs(valmin);
        [valelec,poselec] = max(MaxAmp); % max amplitude across electrodes of interest
        Spindleboundsok(ibok,4:5) = [valelec indelecoi(poselec)]; 
        ibok=ibok+1;
    end
end

%Check that the spindles could be epoched within considered time window
%(i.e. -300 to +1000 ms)
timeseg=[-300/1000*fsample(D), 1000/1000*fsample(D)];
lowtimeseg=Spindleboundsok(:,1)+timeseg(1);
uptimeseg=Spindleboundsok(:,2)+timeseg(2);
boundsok=intersect(find(lowtimeseg>=1), find(uptimeseg<=nsamples(D)));
Spindleboundsok=Spindleboundsok(boundsok,:);

numsp=size(Spindleboundsok,1);

if isempty(Spindleboundsok)
    disp('No spindles in this data set using these parameters')
    return
else
    disp(['Number of spindles detected: ',num2str(numsp)])
end
 
% save spindle triggers in D (will be overwritten if wavelet analysis)
OrigD = D;
if isfield(D.CRC,'spindles')
    D.CRC.spindles=[];
end
D.CRC.spindles.bounds     = Spindleboundsok(:,1:2);
D.CRC.spindles.duration   = Spindleboundsok(:,3);
D.CRC.spindles.amplitude  = Spindleboundsok(:,4);
D.CRC.spindles.maxelectrode   = strvcat(labchan(Spindleboundsok(:,5))); %removed indelecoi
D.CRC.spindles.good=ones(size(Spindleboundsok,1),1);

evtsp=struct('type',cell(1,numsp),...
    'value',zeros(1,numsp),...
    'time',zeros(1,numsp),...
    'duration',zeros(1,numsp),...
    'offset',zeros(1,numsp));
allevt=events(D);
evt_orig=[];
for i=1:size(allevt,2)
    if isempty(strfind(allevt(i).type,'SP')) && ...
            ~isempty(allevt(i).time)
       evt_orig=[evt_orig allevt(i)];
    end
end

for isp = 1:numsp
    evtsp(isp).type         = 'SP';
    evtsp(isp).value        = '777';
    evtsp(isp).time         = Spindleboundsok(isp,1)/D.fsample; %in seconds
    evtsp(isp).duration     = Spindleboundsok(isp,3)/1000; % in seconds
    evtsp(isp).offset       = 0;
end
allevt=[evt_orig, evtsp];
D=events(D,1,allevt);
if isfield(D.CRC,'goodevents') 
    if ntrials(D)==1
        D.CRC.goodevents=ones(1,size(allevt,2));
    else
        D.CRC.goodevents=ones(1,ntrials(D));
    end
end
save(D);


% -------------------------------------------------------------------------
%Performing wavelet analysis to determine the anterior or posterior
%character of the spindles

if isfield(handles,'wav') && handles.wav



disp('.............  Wavelet analysis ..............')
disp('..............................................')

% Epoching
S.D         = D;
S.bc        = 0;
S.pretrig   = -300;  % - pre-trigger time [in ms]
S.posttrig  = 1000; % - post-trigger time [in ms]
S.trialdef.conditionlabel = 'Spindles'; % - string label for the condition
S.trialdef.eventtype  = 'SP'; %    - string
S.trialdef.eventvalue  = '777';  %   - string, numeric or empty
S.reviewtrials = 0;  %   - review individual trials after selection
S.save         =0;    % - save trial definition

D = spm_eeg_epochs(S);

% Time frequency
S.D             = D;
S.channels      = 'EEG';
S.frequencies   = 11:16;
S.timewin       = [-250 1000];%time interval in ms
S.method        = 'morlet';
S.phase         = 0;
Dtf = spm_eeg_tf(S);

foi = find(Dtf.frequencies>=11 & Dtf.frequencies<=16);

%  keep Fz, Cz, Pz, average over time, rearrange (elecxfeq) x ispindles
curdata = squeeze(mean(Dtf(indelecoi,foi,Dtf.fsample/4:Dtf.fsample,:),3));
curdata3 = squeeze(sum(curdata,2));

%More power in the band of interest on the frontal or posterior channels?
if size(fronti,2)>1
    datf=mean(curdata3(fronti,:));
else
    datf=curdata3(fronti,:);
end
if size(pari,2)>1
    datp=mean(curdata3(pari,:));
else
    datp=curdata3(pari,:);
end
Ant_sp = find( datf> datp);
Post_sp = find(datf < datp);
undef_sp=find(datf == datp);

Fmap = repmat(Dtf.frequencies(foi),nelecoi,1);
MaxF = zeros(size(curdata,3),1);
for  icur = 1:size(curdata,3)
    dataoi = curdata(:,:,icur);
    [val,pos] = max(dataoi(:));
    MaxF(icur) = Fmap(pos);
end
    
% -------------------------------------------------------------------------
% save fast and slow spindles

D= OrigD;
D.CRC.spindles=[];
D.CRC.spindles.bounds     = Spindleboundsok(:,1:2);
D.CRC.spindles.duration   = Spindleboundsok(:,3);
D.CRC.spindles.amplitude  = Spindleboundsok(:,4);
D.CRC.spindles.maxelectrode   = strvcat(labchan(Spindleboundsok(:,5))); % removed indelecoi
D.CRC.spindles.index_antsp   = Ant_sp;
D.CRC.spindles.index_postsp   = Post_sp;
D.CRC.spindles.index_undefsp   = undef_sp;
D.CRC.spindles.sp_frequency   = MaxF;
D.CRC.spindles.good=ones(size(Spindleboundsok,1),1);
% save version number of routine
[v,r] = crc_fasst_utils('Ver',mfilename);
D.CRC.spindles.verSP = struct('v_nr',v,'rel',r);        

evtsp=struct('type',cell(1,numsp),...
    'value',zeros(1,numsp),...
    'time',zeros(1,numsp),...
    'duration',zeros(1,numsp),...
    'offset',zeros(1,numsp));

for isp = 1:numsp
    if find(Ant_sp == isp)
        evtsp(isp).type         = 'Ant-SP';
        evtsp(isp).value        = '555';
    elseif find(Post_sp == isp)
        evtsp(isp).type         = 'Post-SP';
        evtsp(isp).value        = '666';
    else
        evtsp(isp).type         = 'SP';
        evtsp(isp).value        = '777';
    end    
    evtsp(isp).time         = Spindleboundsok(isp,1)/D.fsample; %in seconds
    evtsp(isp).duration     = Spindleboundsok(isp,3)/1000; % in seconds
    evtsp(isp).offset       = 0;
end

allevt=[evt_orig, evtsp];
D=events(D,[],allevt);
if isfield(D.CRC,'goodevents') 
    if ntrials(D)==1
        D.CRC.goodevents=ones(1,size(allevt,2));
    else
        D.CRC.goodevents=ones(1,ntrials(D));
    end
end
save(D);
end

% -------------------------------------------------------------------------
% -------------------------- Display --------------------------------------
%--------------------------------------------------------------------------


% accept/reject
if handles.review && ~isempty(Spindleboundsok)
    disp('.......... Visually check spindles ...........')
    disp('..............................................')
    %get data into flag
    flags=struct('index',[], 'file',[], 'Dmeg', []);
    flags.Dmeg{1}=D;
    flags.file{1}=D.fname;

    %select by defaults channels corresponding to Fz and Pz
    flags.index=zeros(1,nelecoi);
    for i=1:nelecoi
        flags.index(i)= indchannel(D, labchan(indelecoi(i)) ); %1st electrode of interest
    end
    %display
    dis_selchan(flags)
elseif isempty(Spindleboundsok)
    disp('No spindles detected in this file:')
    disp('Check the artefact detection if this result is not plausible')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
