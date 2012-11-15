function D = crc_eeg_rdata_edf(Fdata)
% converts EEG data from EDF- to SPM8-format
%
% FORMAT D = crc_eeg_rdata_edf(Fdata)
% 
% Fdata       - filename of edf-file
% D           - SPM8 compatible matlab structure for an M/EEG object.
%
%_______________________________________________________________________
% 
% crc_eeg_rdata_edf reads a continuous *.edf file, stores everything
% in struct D and saves struct D to mat-file. The data is stored separately
% in a 'dat-file', created from the edf file.
%
% Note:
% - data are assumed continuous.
% - data format is assumed int16.
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% But largely inspired from spm_rdata_bdf (from an old version of SPM) 
% which had been written by Stefan Kiebel, thanks mate !
% $Id$

%% Read in file information
if nargin<1
    Fdata = spm_select(1,'^.*\.[eE][dD][fF]$','select edf-file');
end

%% read data/header, load channel setup file
[fdata,header] = crc_eeg_readHeaderEDF(Fdata);

%Introduce a template with recognized electrodes.
Fchannels = 'CRC_electrodes.mat';
if ~exist(Fchannels,'file'),
    % Create electrode template if not there
    crc_electrodes;
end
CRC_el = load('CRC_electrodes.mat');

%% Create blank spm8 structure
% D = struct('type',[],'Nsamples',[],'Fsample',[],'timeOnset',[], ...
%     'trials',[], 'channels', struct, 'data', [], 'fname', [], 'path', [], ...
%     'sensors', [], 'fiducials', [], 'artifacts', [], 'transform', [], ...
%     'other', [], 'history', []);
D = meeg;
D = struct(D);

%% Set key information first: type, fname, path, Fsample, timeOnset

% Assume data are continuous...
D.type = 'continuous' ;

[pth,fname,ext] = fileparts(Fdata);
D.fname = [fname,'.mat'];
D.path = pth;
D.timeOnset = 0;
D.Fsample = header.samplerate(1);
fnamedat = spm_str_manip(fdata,'t');
Nchannels = header.channels;

%% Deal with the data
datatype = 'int16';
D.Nsamples = header.nrsamples(1)*header.records;

% Scaling issues: not sure how it works... 
% accounting now for real data unit, as suggested by Giovanni Piantoni.
rescale = ones(Nchannels,1);
for k = 1:Nchannels
    if strfind(header.physdime(k,:), 'mV')
        rescale(k) = 1000;
    else
        rescale(k) = 1;
    end
end
scl_slope = (header.physmax - header.physmin)./(header.digimax-header.digimin) .*rescale;

dat = file_array( ...
    fullfile(D.path, fnamedat), ...     % fname     - filename
    [Nchannels D.Nsamples],...          % dim       - dimensions
    spm_type(datatype), ...             % dtype     - datatype
    0, ...                              % offset    - offset into file (default = 0)
    scl_slope);                         % scl_slope - scalefactor
data = struct(...
    'fnamedat',fnamedat, ...
    'datatype',datatype, ... % Not sure if should indicate '-le' at the end
    'y',dat);
%     'scale',ones(Nchannels,1), ... % Initialised at 1 as scaling is already stored in file_array object
D.data = data;

%% Deal with channels
blankchan = struct('label',[],'bad', [], 'type', 'unknown', ...
    'X_plot2D', [],'Y_plot2D',[], 'units', 'unknown');

for ii=1:header.channels
    tempchan = blankchan;
    % Get chnnel label from header, i.e. EEG montage
    tempchan.label = deblank(header.channelname(ii,:));
    el_ind = strmatch(upper(tempchan.label),upper(CRC_el.names),'exact');
        % Index in crc_electrodes, or empty
    tempchan.bad = 0; % Not a bad channel by default

    % By default type is 'unknown', unless it is clearly some EOG/EMG/ECG
    % channel, or if its name fits in the crc_electrodes as an 'EEG' 
    % channel (when crc_types = -1).
    if findstr(upper(tempchan.label),'EOG')
        tempchan.type = 'EOG';
    elseif findstr(upper(tempchan.label),'EMG')
        tempchan.type = 'EMG';
    elseif findstr(upper(tempchan.label),'ECG')
        tempchan.type = 'ECG';
    elseif findstr(upper(tempchan.label),'EKG')
        tempchan.type = 'ECG';
    elseif ~isempty(el_ind)
            if CRC_el.crc_types(el_ind)==-1
                tempchan.type = 'EEG';
            end
    end

    % Pick 2D position of channel (for display purpose) in crc_electrodes,
    % if it's available there.
    if ~isempty(el_ind)
        tempchan.X_plot2D = CRC_el.Cpos(el_ind,1);
        tempchan.Y_plot2D = CRC_el.Cpos(el_ind,2);
    else
        tempchan.X_plot2D = NaN;
        tempchan.Y_plot2D = NaN;
    end

    % Get units from header
    tempchan.units = deblank(header.physdime(ii,:));
    % save in channel array
    channels(ii) = tempchan; %#ok<AGROW>
end
% D.channels = channels';
D.channels = channels;

%% Handle the field trial
% Treat marker as seperate channels in edf format.
% Not sure this is always the case, but it's like that with the data at
% hand...

blankevent = struct('type','Unknown','value',[],'time',0, ...
                'duration',0,'offset',0);

% create 1st event = trial/session
events = blankevent;
events.type = 'Trial';
events.value = 'start';
events.duration = D.Nsamples;
events.time = 1/D.Fsample;

% Other real events
types = [1 2];
ch_mkr = strmatch('mkr',lower(strvcat(channels.label)));
N_mkr = length(ch_mkr);
if N_mkr
    d_mkr = data.y(ch_mkr,:); d_mkr = d_mkr/max(abs(d_mkr(:)));
    dd = diff(d_mkr); %figure, plot(dd,'.');
    lpeak_pos = find(dd> .8); % find 'large' positive jumps
    lpeak_neg = find(dd<-.8); % find 'large' negative jumps

    % Sometimes, there are 2 consecutive points with large jumps
    % => Keep only the first of the 2 of points
    lkeep_p = [1 find(diff(lpeak_pos)>1)+1];
    lpeak_pos = lpeak_pos(lkeep_p);
    lkeep_n = [1 find(diff(lpeak_neg)>1)+1];
    lpeak_neg = lpeak_neg(lkeep_n);

    ev_timing = [lpeak_pos lpeak_neg];
    [ev_timing,perm] = sort(ev_timing);
    ev_code = [ones(size(lpeak_pos))*types(1) ones(size(lpeak_neg))*types(2)];
    ev_code = ev_code(perm);
    for ii=1:length(ev_code)
        events(ii+1) = blankevent;
        events(ii+1).type = ev_code(ii);
        events(ii+1).value = ev_code(ii);
        events(ii+1).time = ev_timing(ii)/D.Fsample;
    end
end

trials = struct('label','Undefined','onset',1/D.Fsample, ...
                       'repl',1,'bad',0,'events',events);
D.trials = trials;

%% Deal now with filter, reref, transform and descript field, info
D.fiducials = struct([]);
D.sensors   = struct([]);
D.artifacts = struct([]);
% These fields can be updated later within SPM8, when doing source
% reconstruction or DCM-EEG
D.history   = struct('fun','crc_eeg_rdata_edf', ...
                     'args',struct('D',struct('fname',Fdata)));
D.transform.ID = 'time';    % Time series not frequency spectrogram

% turn into meeg object and add info
D = meeg(D);
D.info = struct('EDFversion',header.version,'PID',header.PID, ...
         'RID',header.RID,'date',header.T0(1:3),'hour',header.T0(4:6));
         % just in case, save the ASCII header, including date and hour of
         % recording => Very useful for sleep studies !
[v,r] = crc_fasst_utils('Ver',mfilename);
D.info.ver = struct('v_nr',v,'rel',r);        
[v,r] = crc_fasst_utils('Ver','crc_eeg_readHeaderEDF');
D.info.ver_readH = struct('v_nr',v,'rel',r);

%% Save results
save(D)

return

























