function D = crc_eeg_rdata_egi(Fdata)
% converts EEG data from EGI-raw to SPM8-format
%
% FORMAT D = crc_eeg_rdata_egi(Fdata)
% 
% Fdata       - filename of raw-file(s)
% D           - SPM8 compatible matlab structure for an M/EEG object.
%
%_______________________________________________________________________
% 
% crc_eeg_rdata_egi reads a continuous *.raw file (exported from EGI), 
% stores everything in struct D and saves struct D to mat-file. The data is
% stored separately in a 'dat-file', created from the raw file.
% There are calls to crc_readHeaderEGI.m, which is adapted from readedf.
%
% If multiple files are passed/selected, then it is assumed that they form
% a series of successive continuous files that were chopped at the
% exportation stage. For example, multi-hours sleep recording data ending 
% up in 5 separate files. Those files will be automatically concatenated
% when reading in the data.
%
% Note:
% - data are assumed continuous.
% - no check for file order, i.e. trust the selection order!
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% But largely inspired from spm_rdata_bdf (from an old version of SPM) 
% which had been written by Stefan Kiebel, thanks mate !
% $Id$

%% Read in file information
if nargin<1
    Fdata = spm_select([1 Inf],'^.*\.[rR][aA][wW]$','select EGI-raw file');
end
Nfd = size(Fdata,1);

%% read data/header, load channel setup file
[fdata,header] = crc_eeg_readHeaderEGI(Fdata);

%Introduce a template with recognized electrodes.
Fchannels = 'CRC_electrodes.mat';
if ~exist(Fchannels,'file'),
    % Create electrode template if not there
    crc_electrodes;
end
CRC_el = load('CRC_electrodes.mat');

%% Create blank spm8 structure
D = struct('type',[],'Nsamples',[],'Fsample',[],'timeOnset',[], ...
    'trials',[], 'channels', struct, 'data', [], 'fname', [], 'path', [], ...
    'sensors', [], 'fiducials', [], 'artifacts', [], 'transform', [], ...
    'other', [], 'history', []);

%% Set key information first
% type, fname, path, Fsample, timeOnset

% Assume data are continuous...
D.type = 'continuous' ;

[pth,fname,ext] = fileparts(deblank(Fdata(1,:)));
if Nfd>1
    D.fname = [crc_get_defaults('cat.prefix2'),fname,'.mat'];
else
    D.fname = [fname,'.mat'];
end
D.path = pth;
D.timeOnset = 0;
D.Fsample = header.Srate;
fnamedat = spm_str_manip(fdata,'t');
Nchannels = header.Nchan;

%% Deal with the data
datatype = header.fodat;
if findstr(datatype,'float');
    datatype_le = [datatype,'-le'];
end
D.Nsamples = header.Nsamples;

%Scaling issues: not sure how it works... to check.
if ~header.Cbits && ~header.Scrange
    scl_slope = []; % No need to provide scaling as floats
                    % used to be: ones(Nchannels,1); % data expressed in µV
else % data in A/D format, assume we use range and bits to convert in µV
    scl_slope = ones(Nchannels,1)*header.Scrange/2^header.Cbits;
end

dat = file_array( ...
    fullfile(D.path, fnamedat), ...     % fname     - filename
    [Nchannels D.Nsamples],...          % dim       - dimensions
    spm_type(datatype), ...             % dtype     - datatype
    0, ...                              % offset    - offset into file (default = 0)
    scl_slope);                         % scl_slope - scalefactor
data = struct(...
    'fnamedat',fnamedat, ...
    'datatype',datatype_le, ... % Not sure if should indicate '-le' at the end, use datatype if not.
    'y',dat);
%     'scale',ones(Nchannels,1), ... % Initialised at 1 as scaling is already stored in file_array object
D.data = data;

%% Deal with channels
blankchan = struct('label',[],'bad', [], 'type', 'unknown', ...
    'X_plot2D', [],'Y_plot2D',[], 'units', 'unknown');

for ii=1:header.Nchan
    tempchan = blankchan;
    % Get channel label from header, i.e. EEG montage
    tempchan.label = ['e',fliplr(deblank(fliplr(header.channames(ii,:))))];
    % add an 'e' before the electrode number to conform to FieldTrip format
    el_ind = strmatch(upper(tempchan.label),upper(CRC_el.names),'exact');
        % Index in crc_electrodes, or empty
    tempchan.bad = 0; % Not a bad channel by default

    % By default type is 'unknown', unless it is clearly some EOG/EMG/ECG
    % channel, or if its name fits in the crc_electrodes as an 'EEG' 
    % channel (when crc_types = -1).
    if findstr(upper(tempchan.label),'EOG')
        tempchan.type = 'EOG';
    elseif findstr(upper(tempchan.label),'EOGV')
        tempchan.type = 'EOG';
    elseif findstr(upper(tempchan.label),'EOGH')
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
    tempchan.units = 'µV';
    % save in channel array
    channels(ii) = tempchan; %#ok<AGROW>
end
D.channels = channels;

%% Handle the field trial

blankevent = struct('type','Unknown','value',[],'time',0, ...
                'duration',0,'offset',0);

% create 1st event = trial/session
events = blankevent;
events.type = 'Trial';
events.value = 'start';
events.duration = D.Nsamples;
events.time = 1/D.Fsample;

% Other real events
Nevent = length(header.event_tb);
if Nevent
    for ii = 1:Nevent
        events(ii+1) = blankevent;
%         % Either give type the event code and value, the event code index
%         events(ii+1).type = header.ecode(header.event_co(ii),:);
%         events(ii+1).value = header.event_co(ii);
        % or define type 'Trigger' and give event code to 'value'
        events(ii+1).type = 'Trigger';
        events(ii+1).value = header.ecode(header.event_co(ii),:);
        if header.event_du(ii)>1
            events(ii+1).duration = header.event_du(ii)/D.Fsample ;
        end
        events(ii+1).time = header.event_tb(ii)/D.Fsample;
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
D.history   = struct('fun','crc_eeg_rdata_egi', ...
                     'args',struct('D',struct('fname',Fdata)));
D.transform.ID = 'time';    % Time series not frequency spectrogram

%% turn into meeg object and add info
D = meeg(D);
header.time(3) = header.time(3) + header.msec*1e-3; % add up ms to time
D.info = struct('EGIversion',header.version, ...
         'date',header.date','hour',header.time');
         % date and hour of recording => Very useful for sleep studies !
[v,r] = crc_fasst_utils('Ver',mfilename);
D.info.ver = struct('v_nr',v,'rel',r);
[v,r] = crc_fasst_utils('Ver','crc_eeg_readHeaderEGI');
D.info.ver_readH = struct('v_nr',v,'rel',r);

%% Save results
save(D)

return