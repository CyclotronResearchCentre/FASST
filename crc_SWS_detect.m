function [D]= crc_SWS_detect(handles)

% Detects Slow Waves Sleep (and delta waves)on EEG data. The criteria are
% based on Massimini's but were softened to detect all waves. The data are 
% first loaded and stages 3 and 4 are extracted if the file is scored using 
% FASST. If the file is not scored, we supposed that it contains only stages
% 3 or 4. The data are then filtered (0.2-4Hz) if needed. The detection is
% achieved on 4 Regions Of Interest (frontal, left central, right central,
% centro-parietal, which can be either automatically created or manually) 
% and trajectory of each wave is then detected on the whole scalp. 
%
% The criteria used are magnitude criteria in microV ( minimum negative peak amplitude
% -40 for delta and -80 for SWS, minimum total magnitude: 75 for delta and
% 140 for SWS) and duration criteria in ms (duration of negative peak:
% between 250 and 1250 ms, time between up zero crossing and positive 
%peak : maximum 2000 ms). Another criterion was added on the slope between
%the negative and the positive peaks (criteria of  minimum percentile 90).
%
% The outputs are new structures in D.CRC.SW where SW is a structure with:
% - SW, an array of substructures with the fields , describing each SW
%   individually:
%       down        index of downward zero crossing
%       up          index of upward zero crossing
%       negmax      index of maximum negativity using Massimini criteria
%       negmax_tp   negmax in time points
%       posmax_tp   index of maximum positivity
%       posmax      posmax in time points
%       upstart     start of the upstate
%       upend       end of upstate
%       maxslope    maximum of slope index in the upswing
%       channels    for each wave, channel numbers respecting the criteria in E.data
%       electrodes  for each wave, channel names respecting the criteria in E.data
%       delays      for each wave, delay of the minimum in the SW, for each channel
%       uponset     for each wave, the onset of the up state in terms of scan
%       amplitude   peak to peak magnitude
%       code        code of the wave
%       neg_slope   maximum of slope between DZC and negmax
%       negmax_TR   negmax in TR, computed after the beginning of the scan
%       good        good (1) or bad (0) SW
% - origin_count, number of waves 'starting' at each electrode and total 
%   number of waves detected by each electrode.
% - DATA4ROI, structure with information about the ROIs used for the 
%   preliminary SW detection.
%
% If the user chose to review the detected waves, another structure will be
% created containing only 'accepted' waves: D.CRC.goodSW
%
% The code of the wave is 3'(number of the first electrode detecting it)' 
% forSWS and 4'(number of the first electrode detecting it)' for delta 
% waves. All fields are in seconds or in time points when said explicitly.
%
% The field D.events was also modified to contain the type 'delta'
% for a delta wave and 'SW' for a slow wave, the value equals to the code of 
% the wave, the time corresponding to the maximum power of the negative peak
% (in seconds) and the duration corresponding to the maximum delay.
%
% The possibility to review all waves is given and allows the user to
% reject some waves. The 'accepted' waves will be saved in D.CRC.goodSW,
% separately from the automatically detected waves.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This function was built with the help of Pierre Maquet and Christophe
% Phillips who have kindly shared routines and parts of code.
% Jessica Schrouff, 23/08/2007, CRC Ulg
% last modified 19/01/2009 for SPM8
% last modified 22/10/2009 for FASST toolbox
% last modified 11/05/2010 for chunked files and default electrodes
% positions
%__________________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by J. Schrouff & C. Phillips, 2009.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

global gindex
clc
close all
% loading files
%--------------------------------------------------------------------------
if ~nargin
    file = spm_select(1, 'mat', 'Select cleaned EEG file','',pwd,'.*');
    Dss = crc_eeg_load(file);
else
    Dss = crc_eeg_load(handles.fname);
end

crcdef = crc_get_defaults('swsd');
start = 0;
%selecting the scorer
scorer=handles.scorer;
%extract periods of interest
if handles.analyse==3
    if isfield(Dss.CRC, 'score') && size(Dss.CRC.score,1)>4
        if ~isempty(handles.stagesw)
            [Ds, cleanSW] = crc_extractSW(Dss,handles.analyse,handles.stagesw,scorer);
        else
            [Ds, cleanSW] = crc_extractSW(Dss,handles.analyse,crcdef.stagesw,scorer);
        end
    end
elseif handles.analyse==2
    [Ds,cleanSW]=crc_extractSW(Dss,handles.analyse,crcdef.stagesw,scorer);
    disp('Warning: assuming file contains only SWS stages')
else
    [Ds,cleanSW]=crc_extractSW(Dss,handles.analyse, crcdef.stagesw,scorer);
    %get the part of the vector comprised between the begin and the
    %start times edited by the user
    belbound=find(cleanSW>=handles.Begpts);
    upbound=find(cleanSW<=handles.Endpts);
    cleanSW=intersect(belbound,upbound);
    start=handles.Begpts/fsample(Ds); %relative start in s
end
csw=cleanSW/fsample(Dss);

%consider joint fMRI-EEG case
startscan=0;
%get time of scan start w.r.t EEG start 
a=Ds.events(:);
if handles.fmri && ~isempty(handles.TR) && ~isempty(handles.marker) && ...
        ~isempty(a)
    for i=1:size(a,2)
        if strcmpi(num2str(a(i).value),handles.marker)
            startscan=a(i).time;
            break
        end
    end
    if startscan==0
        disp('Marker not found in EEG events structure!!')
        cont=spm_input('Continue?',1,'y|n',[1 0],0);
        close gcf
        if cont==0
            return
        end
    end
end


%Checking montage
%------------------------------------------------------------------
%if no re-reference field in D.history, then warning and ask to continue or
%not

disp('---- Checking reference----')
if handles.reref==0
    disp('WARNING: no re-referencing found!-- please re-reference using spm8 (spm_eeg_montage)')
    cont=spm_input('Continue?',1,'y/n',[1,0],1);
else
    disp('re-referencing achieved....ok')
    cont=1;
end
if cont~=1
    return
end

%Checking units
%------------------------------------------------------------------
%if 'unknown' units for the EEG channels, supposing it is in µV and warning

disp('---- Checking units----')

bad=badchannels(Ds);
eeg=find(strcmpi(crcdef.type,chantype(Ds)));
D_channels_eeg=setdiff(eeg,bad);


%filtering data
%--------------------------------------------------------------------------
%Use of a lowpass filter (cutoff rate asked in gui) followed by a highpass
%filter. Filtering is achieved during detection to avoid writing new data 
%on the hard drive.
%For Slow Waves Sleep detection, we recommand lowpass filtering at 4Hz and 
%highpass filtering at 0.25.

if ~isempty(handles.highfc)&& ~isempty(handles.lowfc)
    fqcut=[handles.highfc, handles.lowfc];
else
    fqcut=[crcdef.highfc, crcdef.lowfc];
end
disp('-----filtering data-----')
args=[];
dd=0;
a=history(Ds);
for i=1:size(a,2)
    if strcmpi('spm_eeg_filter',a(i).fun)
        dd=dd+1;
        args=[args, a(i).args.filter];
    end
end
vv=0;
PHz=[];
if ~isempty(args)
    for i=1:size(args,2)
        [d1,smax]=max(size(args(i).PHz));
        if d1>1 && smax==1
            freqcut=(args(i).PHz)';
        else
            freqcut=args(i).PHz;
        end
        PHz=[PHz,freqcut];
    end
    fc=sort(PHz);
    if fc(1)<0.3 || fc(2)<10 || fc(2)>3
        vv=1;
    end
end

fhc = fqcut(1)/(fsample(Ds)/2);
flc = fqcut(2)/(fsample(Ds)/2);
order = crcdef.butterorder;
[b1,a1] = butter(order,fhc,'high');
[b2,a2] = butter(order,flc,'low');

%creating non overlapping regions of interest
%--------------------------------------------------------------------------

disp('---creating non overlapping ROIs---') 
DATA4ROI = struct('data',[],'channels',[],'rate',[],'roisel',[],'nameroi',[]);

auto=handles.roisel;
h = waitbar(0,'Please wait...');
if auto==1
    %automatic selection of ROI based on electrodes position in
    load CRC_electrodes.mat
    eeg_chan=find(strcmpi(crcdef.type,chantype(Ds)));
    if ~crcdef.usetheor
        try
            pos_eeg_chan=(coor2D(Ds,eeg_chan))';
        catch
            pos_eeg_chan=zeros(numel(eeg_chan),2);
            for ich=1:numel(eeg_chan)
                iselec=strcmpi(chanlabels(Ds,eeg_chan(ich)),names);
                pos_eeg_chan(ich,:)=pos(find(iselec),:);
            end
        end
    else
        pos_eeg_chan=zeros(numel(eeg_chan),2);
        for ich=1:numel(eeg_chan)
            iselec=strcmpi(chanlabels(Ds,eeg_chan(ich)),names);
            pos_eeg_chan(ich,:)=pos(find(iselec),:);
        end
    end
    ROI_centers=crcdef.ROI_centers;
    for i=1:size(ROI_centers,1)
        pos_cent=ones(size(pos_eeg_chan,1),2);
        pos_cent(:,1)=ROI_centers(i,1);
        pos_cent(:,2)=ROI_centers(i,2);
        dist_pos=((pos_cent(:,1)-pos_eeg_chan(:,1)).^2+(pos_cent(:,2)-pos_eeg_chan(:,2)).^2).^0.5;
        roi_chan=find(dist_pos<=0.1);
        tmp = zeros(1,size(Ds,2));
        div=1;
        
        for j=1:length(roi_chan)
            scales=crc_scales(Ds,roi_chan(j));
            importdat=Ds(roi_chan(j),:)/scales;            
            if vv~=1
                importdat = filtfilt(b1, a1, importdat);
                importdat = filtfilt(b2, a2, importdat);
            end
            tmp(1,:) = tmp(1,:)*(1-(1/div))+importdat*(1/div);
            div=div+1;            
        end
        DATA4ROI.data= [DATA4ROI.data ; tmp];
        string = ['Please wait... ' num2str(100*(i/size(ROI_centers,1))) ' %'];
        waitbar((i)/size(ROI_centers,1),h,string)
    end
elseif auto==0
    numroi=handles.numroi;
    name_roi=handles.name_roi;
    for iroi = 1:numroi %regions of interest
        tmp = zeros(1,size(Ds,2));
        div=1;
        sel_ROI=handles.sel_ROI;
        for iselroi=1:size(sel_ROI{iroi},2)
            scales=crc_scales(Ds,sel_ROI{iroi}(iselroi));
            importdat = Ds(sel_ROI{iroi}(iselroi),:)/scales;
            if vv~=1
                importdat = filtfilt(b1, a1, importdat);
                importdat = filtfilt(b2, a2, importdat);
            end
            tmp(1,:) = tmp(1,:)*(1-(1/div))+importdat*(1/div);
            div=div+1;
            
        end
        DATA4ROI.data= [DATA4ROI.data ; tmp];
        string = ['Please wait... ' num2str(100*(iroi/numroi)) ' %'];
        waitbar((iroi)/numroi,h,string)
    end
end
string = ['Please wait... ' num2str(100*1) ' %'];
waitbar(1,h,string)
close(h)

% Completing DATA4ROI structure
DATA4ROI.roisel=auto;
if auto==0
    DATA4ROI.nameroi=name_roi;
else
    DATA4ROI.nameroi=[];
end
DATA4ROI.channels = chanlabels(Ds,D_channels_eeg);
DATA4ROI.rate = fsample(Ds);
clear  tmp ichannel ielec iroi 
       
%detecting criteria over the ROIs
%--------------------------------------------------------------------------

disp(['---detecting criteria over the ROIs---'])

%structure for each wave

SW=struct('down', [], ...  %index of downward zero crossing
    'up', [], ...          % index of upward zero crossing
    'negmax', [], ...      % index of maximum negativity using Massimini criteria
    'negmax_tp', [], ...   % negmax in time points
    'posmax', [], ...      % index of maximum positivity
    'posmax_tp', [], ...   % posmax in time points
    'upstart', [], ...     % start of the upstate
    'upend', [], ...       % end of upstate
    'maxslope', [], ...    % maximum of slope index in the upswing
    'channels', [], ...    % for each wave, channel numbers respecting the criteria in E.data
    'electrodes', [], ...  % for each wave, channel names respecting the criteria in E.data
    'delays', [], ...      % for each wave, delay of the minimum in the SW, for each channel
    'uponset', [], ...     % for each wave, the onset of the up state in terms of scan
    'amplitude', [],...    % peak to peak magnitude
    'neg_slope', [],...    % maximum of negative slope
    'code',[],...          % code of the wave: 4 for delta and 3 for SWS followed by the number of the channel
    'negmax_TR',[],...     % negmax in TR, computed after the beginning of the scan
    'good',[]);            % assumes the wave is a 'good' event, can be modified during display

param1 = struct('SWlength',crcdef.SWlength, ...
                'SWmAmpl',crcdef.SWmAmpl);
countwaves=0;
            
for idataf = 1:size(DATA4ROI.data,1) % looping over the 4 ROIs
    disp(['ROI number ' num2str(idataf) ])
    [SW,countwaves,Ds] = find_SW(DATA4ROI.data(idataf,:),DATA4ROI,Ds, ...
                                    param1,countwaves,SW,crcdef);
    disp(['Total number of waves detected:' num2str(countwaves)])
end

if ~(size(SW,2)>1)
    disp('-- No Slow Waves detected in this data set --')
    D=Ds;
    return
end

% Removing data from DATA4ROI, as it can take too much space (several Mb).
DATA4ROI = rmfield(DATA4ROI,'data');

%detecting waves on each channel - establishing travel directory
%--------------------------------------------------------------------------
%detects waves on all electrodes and establishes the trajectory followed by
%the wave on the scalp.

disp(['----establishing traveling direction----'])

ONSETS= struct('channel',[],'onsets',[]);
SWS_count=0;
origin_count=pm_origin_count(Ds);
delta_count=0;
for jchannel = 1:size(origin_count,1)
    ONSETS(jchannel).channel = origin_count{jchannel,1};
end
evtsw=struct('type',[],'value',[],'duration',[],'time',[],'offset',[]);
orevtsw=evtsw;

%Define additional indexes to load in order to filter properly and
%therefore avoid windowing effects

[h1,t1]=impz(b1,a1);
[h2,t2]=impz(b2,a2);
addpts=max(size(t2,2),size(t1,1));
addpts=min(addpts,20000);
h = waitbar(0,'Please wait...');
badwav=[];
for iwav = 1:size(SW,2)                          % for each detected wave

    %     Mind that the index cannot be <0
    sampleindex =round(SW(iwav).negmax_tp-param1.SWlength(4)*fsample(Ds)/1000:SW(iwav).negmax_tp+param1.SWlength(4)*fsample(Ds)/1000);
    extended_idx=[sampleindex(1)-addpts:sampleindex(1)-1, sampleindex, sampleindex(end)+1:sampleindex(end)+addpts];
    markervect=[zeros(1,addpts) ones(1,length(sampleindex)) zeros(1,addpts)];
    
    tokeep = find(and(extended_idx > 0, extended_idx < nsamples(Ds)));

    extended_idx=extended_idx(tokeep);
    markervect=markervect(tokeep);
    scales=crc_scales(Ds,D_channels_eeg);
    scales=repmat(scales,1,numel(extended_idx));
    sample = Ds(D_channels_eeg,extended_idx)./scales;             % display 10 s
    if vv~=1 %filter data if necessary
        for zchan=1:size(sample,1)
            sample(zchan,:) = filtfilt(b1, a1, sample(zchan,:));
            sample(zchan,:) = filtfilt(b2, a2, sample(zchan,:));
        end
    end

    recupidx=find(markervect==1);
    sample=sample(:,recupidx);
    
    
    rds=sample(:,round(size(sample,2)/2)-20:round(size(sample,2)/2)+20);
    power=sum(rds.^2);  %time = max of power
    [ipower,jpower]=max(power);
    SW(iwav).negmax_tp=SW(iwav).negmax_tp+(jpower-round(size(rds,2)/2));
    SW(iwav).negmax=SW(iwav).negmax_tp*1000/fsample(Ds);
    [sampmin,sampos] = min(sample');
    [sortedsamp, samporder] = sort(sampos);          % sortedsamp = ordered delays; samporder = ranking order

    %     delete shallow and 'up' channels
    deriv = diff(sample(samporder,:)')';
    if SW(iwav).code == 100
        samp_delta =(sampmin(samporder) < param1.SWmAmpl(1));
        sampindex = samp_delta .* (mean(deriv(:,1:round(size(sample,2)/2))') < 0);
        delta_count=delta_count+1;% modified for delta waves
    elseif SW(iwav).code == 101
        samp_SWS = (sampmin(samporder) < param1.SWmAmpl(2));
        sampindex = samp_SWS .* (mean(deriv(:,1:round(size(sample,2)/2))') < 0);
        SWS_count=SWS_count+1;
    end
    if ~any(sampindex)
        badwav=[badwav,iwav];
        continue
    end
    channels_eeg = chanlabels(Ds,D_channels_eeg);
    sortdelays = sortedsamp(find(sampindex));
    sortchannel = samporder(find(sampindex));
    SW(iwav).channels = D_channels_eeg(sortchannel); %channels respecting magnitude and slope criteria, ordered by delays
    SW(iwav).electrodes = channels_eeg(sortchannel); % corresponding electrodes
    SW(iwav).delays = sortdelays; % increasing delays (of minimum) 
    
    %define new codes
    if SW(iwav).code == 101
        if length(num2str(SW(iwav).channels(1)))==1
            a=['30',num2str(SW(iwav).channels(1))];
        elseif length(num2str(SW(iwav).channels(1)))==2
            a=['3',num2str(SW(iwav).channels(1))];
        end
        evtsw(iwav).type='SW';
    elseif SW(iwav).code == 100
        if length(num2str(SW(iwav).channels(1)))==1
            a=['40',num2str(SW(iwav).channels(1))];
        elseif length(num2str(SW(iwav).channels(1)))>=2
            a=['4',num2str(SW(iwav).channels(1))];
        end
        evtsw(iwav).type='delta';
    end
    evtsw(iwav).value=a;
    SW(iwav).code=str2double(a);
    evtsw(iwav).time=SW(iwav).negmax_tp/fsample(Ds);
    evtsw(iwav).duration=max(sortdelays)/1000;
    orevtsw(iwav)=evtsw(iwav);
    if handles.analyse==3 || handles.analyse==2
        orevtsw(iwav).time=csw(SW(iwav).negmax_tp);
    elseif handles.analyse==1
        orevtsw(iwav).time=csw(SW(iwav).negmax_tp)+ start;
    end
          
    % count the SW channel start    
    for kelec = 1:size(origin_count,1)
        if strcmpi((origin_count{kelec,1}),char(SW(iwav).electrodes(1)))
            origin_count{kelec,2} = origin_count{kelec,2}+1;
        end
        for i=1:size(SW(iwav).electrodes,2)
            if strcmpi((origin_count{kelec,1}),char(SW(iwav).electrodes(i)))
                origin_count{kelec,3} = origin_count{kelec,3}+1;
            end
        end
    end
    
    %get time in TR if joint EEG-fMRI acquisition
    if handles.fmri && startscan~=0
        SW(iwav).negmax_TR=(SW(iwav).negmax-startscan*1000)/handles.TR;
    else
        SW(iwav).negmax_TR=NaN;
    end
    string = ['Please wait... ' num2str(100*(iwav/size(SW,2))) ' %'];
    waitbar((iwav)/size(SW,2),h,string)
    
end
string = ['Please wait... ' num2str(100*1) ' %'];
waitbar(1,h,string)

if ~isempty(badwav)
    disp('Warning: spikes on one channel were detected')
    for ibad=1:size(badwav,2)
        SW=[SW(1:badwav(ibad)-1),SW(badwav(ibad)+1:end)];
        badwav=badwav-1;
    end
end

%save file with data of interest
Ds.CRC.SW.SW = SW;
Ds.CRC.SW.origin_count = origin_count;
%do not duplicate waves if multiple detections on one file
allevt=events(Ds);
evt_orig=[];
for i=1:size(allevt,2)
    if isempty(strfind(allevt(i).type,'SW')) && ...
            isempty(strfind(allevt(i).type,'delta')) && ...
            ~isempty(allevt(i).time)
       evt_orig=[evt_orig allevt(i)];
    end
end
allevt=[evt_orig, evtsw];
Ds=events(Ds,1,allevt);
if isfield(Ds.CRC,'goodevents') 
    if ntrials(Ds)==1
        Ds.CRC.goodevents=ones(1,size(allevt,2));
    else
        Ds.CRC.goodevents=ones(1,ntrials(Ds));
    end
end
D=Ds;
save(D);


%compute the timing of waves w.r.t. original data file
for iwav = 1:size(SW,2)
        SW(iwav).down = csw(round(SW(iwav).down/1000*fsample(Ds)))*1000+start*1000;         %position of downward zero crossing
        SW(iwav).up = csw(round(SW(iwav).up/1000*fsample(Ds)))*1000+start*1000;             %position of upward zero crossing
        SW(iwav).negmax = csw(SW(iwav).negmax_tp)*1000+start*1000;                 %negative peak position in ms
        SW(iwav).negmax_tp= csw(SW(iwav).negmax_tp)*fsample(Ds)+start*fsample(Ds);     %negative peak position in time points
        SW(iwav).posmax = csw(SW(iwav).posmax_tp)*1000+start*1000;                 %positive peak position in ms
        SW(iwav).posmax_tp= csw(SW(iwav).posmax_tp)*fsample(Ds)+start*fsample(Ds);     %positive peak position in time points
        SW(iwav).upstart = csw(round(SW(iwav).upstart/1000*fsample(Ds)))*1000+start*1000;   %begin of upstate
        SW(iwav).upend = csw(round(SW(iwav).upend/1000*fsample(Ds)))*1000+start*1000;       %end of upstate
end
       
%save both files
Dss.CRC.SW.SW = SW;
Dss.CRC.SW.origin_count=origin_count;
% save version number of routine
[v,r] = crc_fasst_utils('Ver',mfilename);
Dss.CRC.SW.verSW = struct('v_nr',v,'rel',r);        

%do not duplicate waves if multiple detections on one file
allevt=events(Dss);
evt_orig=[];
for i=1:size(allevt,2)
    if isempty(strfind(allevt(i).type,'SW')) && ...
            isempty(strfind(allevt(i).type,'delta')) &&...
            ~isempty(allevt(i).time)
       evt_orig=[evt_orig allevt(i)];
    end
end
allevt=[evt_orig, orevtsw];
Dss=events(Dss,1,allevt);
if isfield(Dss.CRC,'goodevents') 
    if ntrials(Dss)==1
        Dss.CRC.goodevents=ones(1,size(allevt,2));
    else
        Dss.CRC.goodevents=ones(1,ntrials(Dss));
    end
end
D=Dss;
save(D);
close(h)

% visual check of all SWS
%--------------------------------------------------------------------------
%opens the main display to allow the user to review each wave, either by
%using the events (delta or SW) or the scrolling. Also allows the
%construction of delay or potential maps in the bottom right corner.

vc=handles.review;
% Display loop
if vc==1 && ~isempty(D.CRC.SW.SW)
    %add fileio directory of spm8 to have access to read_sens
    dspm=spm('dir');
    addpath([dspm,filesep,'external',filesep,'fileio'])
    elpos=handles.sensauto;
    if elpos==0
        if ~isempty(handles.sensfname)
            zebris_name=handles.sensfname;
        else
            zebris_name=spm_select(1, 'any', 'Select electrodes positioning file','' ,pwd,'.*');
        end
    else
        zebris_name='CRC_electrodes.mat';
    end
    cas=handles.maps;
    
    %get data into flag
    flags=struct('index',[], 'file',[], 'Dmeg', [],'delmap',[]);
    flags.Dmeg{1}=D;
    flags.file{1}=D.fname;
    flags.delmap=struct('elpos',elpos,'zebris_name',zebris_name,'cas',cas);
    %display
    dis_selchan(flags)
end

%save data
D=Ds;
save(D);

%--------------------------------------------------------------------------
%---- SUBFUNCTION TO DETECT WAVES CORRESPONDING TO MASSIMINI CRITERIA  ----
%--------------------------------------------------------------------------


function [SW,countwaves,E] = find_SW(F,DATA4ROI,E,param1,countwaves,SW,crcdef)


% find zero crossings

F(2,:)=sign(F(1,:)); %gives the sign of data
F(3,:)=[0 diff(F(1,:))]; % gives the differential of data
F(4,:) = sign(F(3,:)); %gives the sign of differential

DZC=find(diff(F(2,:)) ==-2);
UZC=find(diff(F(2,:)) ==2);

%criterium of maximum slope index percentile 90

MSI=zeros(1,size(F,2));
MSI_plot=MSI;
MSI_plot(find(F(3,:) > crc_percentile(F(3,:),crcdef.prcentile)))=100;
MSI=find(MSI_plot==100); 

for imsi=1:size(MSI,2)-1
    
    %find nearest MSI and DZC
    indiceDZC = find((DZC-MSI(imsi))<0);
    if ~isempty(indiceDZC)
        indiceDZC=indiceDZC(end);
        iDZC=DZC(indiceDZC);
    else
        iDZC=1;
    end
           
    %find posmin and valmin
    [valmin,indposmin]=min(F(1,iDZC:MSI(imsi)));
    posmin=iDZC+indposmin;
    negslope=min(F(3,iDZC:posmin));
    
    %find iUZC between iDZC and iDZC+SW_length(2)
    upperbound=size(F,2)-(iDZC+DATA4ROI.rate*param1.SWlength(2)/1000);
    if upperbound >0
        iUZC = find(diff(F(2,iDZC:iDZC+round(DATA4ROI.rate*param1.SWlength(2)/1000))) == 2)+iDZC; 
    else
        iUZC= find(diff(F(2,:)) ==2)+iDZC;
    end
    
    
    if ~isempty(iUZC)&& ~isempty(indiceDZC)
        iUZC=iUZC(1);
    
        
        %verification of Massimini criteria on length and magnitude of SWS
        %and delta waves
        
        %criterion on negative peak magnitude
        if ((iUZC-iDZC) <= (param1.SWlength(2)*DATA4ROI.rate/1000) && ...
                (param1.SWlength(1)*DATA4ROI.rate/1000) <=(iUZC-iDZC) )
            
            %negative peak magnitude
            if valmin <= param1.SWmAmpl(1)
                upperbound= size(F,2)-(iUZC+(DATA4ROI.rate*param1.SWlength(3)/1000));
                
                if upperbound >0
                    posmax = iUZC + find(diff(F(4,iUZC + find(F(1,iUZC:iUZC+round(DATA4ROI.rate*param1.SWlength(3)/1000)) > 0))) == -2); 
                else
                    posmax= iUZC +find(diff(F(4,iUZC+find(F(1,iUZC:size(F,2)-1 >0)))) ==-2);
                end
                
                if ~isempty(posmax)
                    start2end_up=[];
                    posmax=posmax(1);
                    valmax=F(1,posmax);
                    
                    
                 %criterion on peak to peak magnitude using mimimal
                 %criteria of -40 and 75 microV
                     if ((indiceDZC+1)< size(DZC,2)) && ((iUZC+(DATA4ROI.rate*param1.SWlength(3)/1000)) < size(F,2))

                         if DZC(indiceDZC+1)-iUZC < DATA4ROI.rate*param1.SWlength(3)/1000
                                start2end_up = iUZC + find(F(1,iUZC:MSI(imsi+1))>= 0.8*valmax); 
                         else
                                start2end_up = iUZC + find(F(1,iUZC:iUZC+round(DATA4ROI.rate*param1.SWlength(3)/1000))>= 0.8*valmax);
                         end
                     end

                     if ~isempty(start2end_up)

                         if (abs (valmax)+ abs (valmin)) >= param1.SWmAmpl(3)

                         %avoid the doubloons
                         
                          if  isempty(SW(end).negmax) || ...% nothing in SW.negmax (1st pass)
                                  ((all(abs(posmin - squeeze(cat(1,SW(:).negmax)./1000*fsample(E)))  > DATA4ROI.rate*param1.SWlength(3)/5/1000))&&...   % need 400ms between SW negativity
                                     (all(abs(posmax - squeeze(cat(1,SW(:).posmax)./1000*fsample(E)))  > DATA4ROI.rate*param1.SWlength(3)/5/1000)))     
                                
                                 %fill the SW structure being careful of
                                %delta waves (different code)
                                
                                countwaves =  countwaves+1;
                                                                
                                ms=max(F(3,posmin: posmax));
                               
                                if (valmin <= param1.SWmAmpl(2) && (abs (valmax)+ abs (valmin)) >= param1.SWmAmpl(4))
                                    code =101;
                                else 
                                    code =100;
                                end
                               
                                disp([ num2str(imsi)  ' detected possible SW - ' ...
                                    num2str(countwaves) ' SW kept in total;  at ' ...
                                    num2str(posmin) ' : ' ...
                                    num2str(valmin) ' microV; at '   ...
                                    num2str(posmax) ' : '  num2str(valmax) ' microV; '])                               
                                    
                                SW(countwaves).down = iDZC*1000/fsample(E);                  %position of downward zero crossing
                                SW(countwaves).up = iUZC*1000/fsample(E);                    %position of upward zero crossing
                                SW(countwaves).negmax = posmin*1000/fsample(E);              %negative peak position in ms
                                SW(countwaves).negmax_tp= posmin;                           %negative peak position in time points
                                SW(countwaves).posmax = posmax*1000/fsample(E);              %positive peak position in ms
                                SW(countwaves).posmax_tp= posmax;                           %positive peak position in time points
                                SW(countwaves).upstart = start2end_up(1)*1000/fsample(E);    %begin of upstate
                                SW(countwaves).upend = start2end_up(end)*1000/fsample(E);    %end of upstate
                                SW(countwaves).amplitude = (abs(valmax)+abs(valmin));       %total magnitude
                                SW(countwaves).uponset = SW(countwaves).negmax ;
                                SW(countwaves).maxslope = ms;
                                SW(countwaves).code= code;
                                SW(countwaves).neg_slope = negslope;
                                SW(countwaves).good=1;
                               
                            end
                        end
                    end
                end
            end
        end
    end
end

N_sw = size(SW,2);
if N_sw>1
    countwaves = N_sw;
else % check if at least one SW was found in 1st ROI(s)
    if ~isempty(SW.down)
        countwaves = N_sw;
    end
end
% countwaves= size(SW,2);

%--------------------------------------------------------------------------
%-----------  SUBFUNCTION TO INITIALIZE COUNT ON CHANNELS   ---------------
%--------------------------------------------------------------------------
function [origin_count] = pm_origin_count(data)
bad=badchannels(data);
eeg=find(strcmpi('EEG',chantype(data)));
data_channels_eeg=setdiff(eeg,bad);
all_names=chanlabels(data);
origin_count=all_names(data_channels_eeg)';
for i=1:size(data_channels_eeg,2)
    origin_count{i}=upper(deblank(origin_count{i}));
    origin_count(i,2)={0};
    origin_count(i,3)={0};
end
origin_count(:,2)={0};
origin_count(:,3)={0};



                