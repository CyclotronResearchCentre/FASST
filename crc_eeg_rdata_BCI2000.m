function D = crc_eeg_rdata_BCI2000(Fdata)
% converts EEG data from BCI2000 to SPM8-format
%
% FORMAT D = crc_eeg_rdata_BCI2000(Fdata)
% 
% Fdata       - filename of data-file
% D           - SPM8 compatible matlab structure for an M/EEG object.
%_______________________________________________________________________
% 
% crc_eeg_rdata_BCI2000 reads a continuous *.dat file in BCI2000 format,
% stores everything in struct D and saves struct D to mat-file. The data is
% stored separately in a 'dat-file', created from the raw file.
% There are calls to crc_readHeaderBCI2000.m, which is reads in all the
% info from the header and generates the seperate binary data file.
%
% NOTE:
% There is room for improvements! All cases of the BCI2000 are certainly
% not covered. Give me some feedback if you encounter a problem.
%__________________________________________________________________
% Copyright (C) 2012 Cyclotron Research Centre

% Written by C. Phillips, 2012.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

%% Read in file information
if nargin<1
    Fdata = spm_select(1,'^.*\.dat$','Select BCI2000 data file');
end

%% read data/header, load channel setup file
[fdat,ev,param] = crc_eeg_readHeaderBCI2000(Fdata);

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
D.fname = ['spm8_',fname,'.mat'];
D.path = pth;
D.timeOnset = 0;
D.Fsample = param.SamplingRate.NumericValue;
fnamedat = spm_str_manip(fdat,'t');
Nchannels = param.NSourceCh;

%% Deal with the data
datatype = 'float32';
datatype_le = [datatype,'-le'];
D.Nsamples = param.Nsamples;
scl_slope = []; % No need to provide scaling as floats
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
D.data = data;

%% Deal with channels
% but there is little information provided in BCI2000 format. :-(
blankchan = struct('label',[],'bad', 0, 'type', 'unknown', ...
    'X_plot2D', NaN,'Y_plot2D',NaN, 'units', 'µV');
channels(Nchannels) = blankchan;
for ii=1:Nchannels
    channels(ii) = blankchan;
    % Get channel label from header, i.e. EEG montage
    if isfield(param,'ChannelNames')
        channels(ii).label = param.ChannelNames.Value{ii};
    else
        % if not available, simply add an 'e' before the electrode number
        channels(ii).label = ['e',num2str(ii)];
    end
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
for ii=1:numel(ev)
    events(ii+1) = ev(ii);
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
D.history   = struct('fun','crc_eeg_rdata_BCI2000', ...
                     'args',struct('D',struct('fname',Fdata)));
D.transform.ID = 'time';    % Time series not frequency spectrogram

%% turn into meeg object and add info
D = meeg(D);

[v,r] = crc_fasst_utils('Ver',mfilename);
D.info.ver = struct('v_nr',v,'rel',r);
[v,r] = crc_fasst_utils('Ver','crc_eeg_readHeaderBCI2000');
D.info.ver_readH = struct('v_nr',v,'rel',r);
D.info = struct('parameters',param);
         % keep the whole set of BCI parameters at hand, just in case.

%% Save results
save(D);

return