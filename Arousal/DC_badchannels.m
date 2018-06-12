function data = DC_badchannels(data, cfg)


% FORMAT DC_badchannels(args)
% Bad channels detection for whole sleep recordings. 
% Bad channels are detected by scoring window and are either:
%       * Noisy
%      or
%       * Flat
% Badchannels are saved in: 
% * C.CRC.DC.badchannels.chan_defaillant = flat channels saved in a matrix.
%                           the number of the row correspond to the
%                           EEG channel number, from 1 to the number of EEG channels 
% * C.CRC.DC.badchannels.chan_incoherent
% the thresholds are available in CRC_get_default('bc')    
%
% INPUT
%       .file   - data file (.mat files)
%__________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% function to detect badchannels by window of 20s or 30s according to the
% time used for scoring.
% Badchannels are recorded in: 
% * C.CRC.DC.badchannels.chan_defaillant
% * C.CRC.DC.badchannels.chan_incoherent
% the thresholds are available in CRC_get_default('bc')
% ----------------------

% Loading and parameters values 
% *****************************
winsize     =   cfg.winsize;
channels    =   cfg.channels;

fs          =   fsample(data);   
nfs         =   fs/20;   % new sample frequency
R           =   ceil(fs/nfs);
res_data    =   zeros(numel(channels),ceil(size(data,2)/R));
% decimate the original signal to a sampling frequency around 50Hz in order
% to have comparable results from different sampling frequency and just
% observe low transition comparison.
% RESAMPLING DATA HERE
try 
    res_data(1,:)    =   downsample(data(channels(1),:),R);
catch 
    res_data    =   zeros(numel(channels), floor(size(data,2)/R));
    res_data(1,:)    =   downsample(data(channels(1),:),R);
end
for ie = 2 : numel(channels)
        res_data(ie,:)    =   downsample(data(channels(ie),:),R);
end

nspl        =   ceil(nsamples(data)/R);
Time        =   nspl/nfs;
% ***************************************
% I added here a new criterium for NofW so that it resembles more the
% overall NofW; before there were to less windows
NofW_all    =   ceil(nsamples(data)/fs/winsize);
NofW        =   ceil(Time / winsize);
if (NofW ~= NofW_all)
    NofW = floor(Time / winsize);
end
% ***************************************
Etime       =   winsize*nfs;

% Thresholds: get them from default file: crc_defaults
% ***********
clear global crc_def;
tr  = crc_get_defaults('qc.bc');

% Obvious channels --------------------------------------------
tr_n    = tr.n;     % noisy channel threshold
tr_f1   = tr.f1;    % 1st threshold of flat channel

% Flat channels -----------------------------------------------
tr_f2   = tr.f2;    % 2nd threshold of flat channel
tr_tf   = tr.tf;    % duration of flat channel
tr_ampl = tr.ampl;  % amplitude depending on standard deviation

% initialization
chan_def = cell(NofW,1);%NaN(length(channels),NofW*winsize);
chan_incoh = cell(NofW,1);%NaN(length(channels),NofW*winsize);

% *************************************************************************
%                    Obviously bad channels detection
% *************************************************************************
h = waitbar(0,'badchannels detection');  
for w = 1 : NofW
    window = res_data(:,(w-1)*Etime +1 : min(w*Etime,nspl));
       
    % --- obvious noisy channel
    eeg_an    =   std(window,[],2); 
    if any(eeg_an >= tr_n)
       idc = find(eeg_an >= tr_n);
       chan_incoh{w} = channels(idc);
       
   % --- obvious flat channel
   elseif any(eeg_an <= tr_f1) 
       idc = find(eeg_an <= tr_f1);
       chan_def{w} = channels(idc); 
   end
    String  =  ['Obvious bad channels detection: ' num2str(w/NofW*100) ' %'];
    waitbar((w/NofW),h,String);
end

% *************************************************************************
%                Finer detection for noisy and flat channels
% *************************************************************************
for w = 1 : NofW
    cki = chan_incoh{w};
    ckn = chan_def{w};
    good_chan = setdiff(channels,union(cki,ckn)); % index of obvious bad channels = [cki and ckn] so we keep index of potential good channels
    [gum igood] = intersect(channels, good_chan);
    if ~isempty(igood) % wenn der Kanal eben nicht leer ist
       window     =   res_data(igood,(w-1)*Etime+1 : min(w*Etime,nspl));
       st_5epo    =   std(window,[],2); % standard deviation along each channel
       for c = 1:numel(igood)
            w_c = window(c,:);        
% ***** flat channel ***** (Devuyst)
            if std(w_c)<tr_f2
                def = abs(w_c);
                v = find(def <= tr_ampl*st_5epo(c)); 
                diff_def = diff(v);
                fin = unique([find(diff_def~=1) length(diff_def)]);
                deb = [1 fin(1:end-1)+1];
                duration_def = (fin - deb)/fs;
                if any(duration_def >= tr_tf)
                    chan_def{w} = union(chan_def{w},channels(igood(c)));  %artefact de défaillance des électrodes
                end  
            end
% ***** noisy channels *****
            other     =     igood~=igood(c); %all good channels except the one we are analyzing
            eeg_oth   =  	mean(window(other,:));
            R         =     corr(w_c(:),eeg_oth(:));%mean(abs((w_c-mean(w_c))./std(eeg_oth)));
            if (R <= 0.2)%tr_r)
               chan_incoh{w} = union(chan_incoh{w},channels(igood(c)));  %artefact de défaillance des électrodes
            end
       end
    end
    String  =  ['Finer bad channels detection: ' num2str(w/NofW*100) ' %'];
    waitbar((w/NofW),h,String);
end
close(h);
% if more than the half of bad channels, all epochs from the relative 30s-time window are removed...
ep = [];
for iw = 1 : NofW 
    cki = chan_incoh{iw};
    ckn = chan_def{iw};
    ckt = union(cki, ckn);
    if numel(ckt)>length(channels)/2
        ep = [ep (iw-1)*winsize+1 : iw*winsize];
    end
end
% save bad channels detected in the CRC.DC structure.
fprintf(1,'Bad channels detected per scoring window \n')
data.CRC.DC.shortartf.badchan = unique(ep);
data.CRC.DC.badchannels.chan_defaillant = chan_def;
data.CRC.DC.badchannels.chan_incoherent = chan_incoh;
save(data);


