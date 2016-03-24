function [D,Fdata] = crc_eeg_rdata_brpr(Fdata,opt)
% converts EEG data from BrainProduct- to SPM-format
% FORMAT D = crc_eeg_rdata_brpr(Fdata)
%
% Fdata       - filename of header-file
% opt         - option flag (can be omited to ensure back compatibility)
% D           - SPM8/12 compatible matlab structure for an M/EEG object.
%_______________________________________________________________________
%
% spm_eeg_rdata_brpr reads a continuous BrainProducts recording with
% *.vhdr/*.vmrk/*.eeg (or .dat) files, stores everything in struct D and
% saves a structure D to a mat-file. This structure is used as a @meeg
% object in SPM8
% The data are stored separately in a data-file. This file is NOT copied
% when importing the data into SPM8/12 !
% Please keep a safe copy of your data...
%
% NOTE: the datafile.eeg file can be renamed into datafile.dat to fit SPM's
% usual format (if opt = true).
%_______________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% But largely inspired from spm_rdata_bdf (from an old version of SPM)
% which had been written by Stefan Kiebel, thanks mate !

%% Read in file information
if nargin<2, opt = false; end
if nargin<1 || isempty(Fdata)
    Fdata = spm_select(1,'^.*\.vhdr$','select EEG data header file');
end

%% read header & marker, load channel setup file
[header,marker] = crc_eeg_readHeaderBrainP(Fdata);
Fchannels = 'CRC_electrodes.mat';
if ~exist(Fchannels,'file'),
    % Create electrode template if not there
    crc_electrodes;
end
CRC_el = load('CRC_electrodes.mat');

%% Create blank spm8/12 structure
D = meeg;
D = struct(D);

%% Get key information first
% type, fname, path, Fsample, timeOnset
if strcmp(header.SegmentationType,'NotSegmented')
    % Data have not been segmented -> continuous
    D.type = 'continuous' ;
else
    % This toolbox has NOT been optimised for segmented data...
    error('This toolbox has NOT been optimised for segmented data, use SPM8 machinery instead for this data set.')
    %     D.type = ''; % Not sure what to put in !
end

[pth,fname,ext] = fileparts(Fdata); %#ok<*NASGU>
D.fname = [fname,'.mat'];
D.path = pth;
D.timeOnset = 0;
D.Fsample = 1e6/header.SamplingInterval;

% Other stuff
Nchannels = header.NumberOfChannels;

%% Deal with the data
if opt
    % ensure that the datafile has .dat extension -> need to rename it!
    [kk,fname,ext] = fileparts(header.DataFile); %#ok<*ASGLU>
    fnamedatO = [fname,ext];
    fnamedat = [fname,'.dat'];
    movefile(fullfile(pth,fnamedatO),fullfile(pth,fnamedat));
else
    % or leave unchanged
    [kk,fname,ext] = fileparts(header.DataFile);
    fnamedat = [fname,ext];
end

if ~strcmp(header.DataFormat,'BINARY')
    error('Data are NOT in binary format, sorry can''t deal with that');
end

% Get from orig data file: binary format, #samples
file_info = dir(fullfile(pth,fnamedat));
switch lower(header.binary.BinaryFormat)
    case 'ieee_float_32'
        datatype = 'float32';
        if isfield(header,'DataPoints') && header.DataPoints>0
            Nsamples = header.DataPoints;
        else
            bin_type = 4;
            Nsamples = floor(file_info.bytes/bin_type/Nchannels);
        end
        scl_slope = ones(Nchannels,1); % No scaling for float data
    case 'int_16'
        datatype = 'int16';
        if isfield(header,'DataPoints') && header.DataPoints>0
            Nsamples = header.DataPoints;
        else
            bin_type = 2;
            Nsamples = floor(abs(file_info.bytes/bin_type/Nchannels));
        end
        scl_slope = ones(Nchannels,1); % Initialize at 1
        l_sc = find(header.channels.Resolution)~=0;
        scl_slope(l_sc) = header.channels.Resolution(l_sc) ;
    otherwise
        error('Don''t know this data format...');
end

l_NewSegment = find(strcmp('New Segment',marker.Type));
D.Nsamples = Nsamples;

% NOTE:
% here above, we're thus assuming that the whole binary file contains data
% so header size, and offset(s) are NOT accounted for !
% If bits had to be left out, it would be necessary to remove them from
% file_info.bytes, with something like
%   Nbytes = file_info.bytes - header.binary.SegmentHeaderSize - header.binary.TrailerSize;
%   D.Nsamples = Nbytes/bin_type/Nchannels;

dat = file_array( ...
    fullfile(D.path, fnamedat), ...     % fname     - filename
    [Nchannels D.Nsamples],...          % dim       - dimensions
    spm_type(datatype), ...             % dtype     - datatype
    0, ...                              % offset    - offset into file (default = 0)
    scl_slope);                         % scl_slope - scalefactor (default = 1)
% NOTE:
% Similar note as above. NO offset in the binary file is considered here.
% If that was the case, then that information should be put in here above.
% Beware also of scaling issues !

data = struct(...
    'fnamedat',fnamedat, ...
    'datatype',datatype, ... % Not sure if should indicate '-le' at the end
    'y',dat);

D.data = data;

%% Deal with channels
blankchan = struct('label',[],'bad', [], 'type', 'unknown', ...
    'X_plot2D', [],'Y_plot2D',[], 'units', 'unknown');
channels(header.NumberOfChannels) = blankchan;
for ii=1:header.NumberOfChannels
    tempchan = blankchan;
    % Get chnnel label from header, i.e. EEG montage
    tempchan.label = header.channels.Cnames{ii};
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
    tempchan.units = header.channels.CUnit{ii};
    
    channels(ii) = tempchan;
end
D.channels = channels';

%% Handle the field trial
% Assume we are dealing with continuous recording !
% Then we need to import the markers from the .vmrk file into SPM8.

% Contrary to SPM8/12 and for ease of coding, we want to restrict the events'
% value to numbers only.
% We therefore attribute a 'code value' to each marker description.
Code = zeros(1,marker.nMrk);
[Des,u_ind] = unique(marker.Descript,'first');
NDes = length(Des);
uCode = zeros(1,NDes);
base = 1000;
fprintf('WARNING: changing event type descriptors into numbers!\n')
for kk=1:NDes
    try
        uCode(kk) = str2double(Des{kk}(2:end));
    catch
        uCode(kk) = NaN;
    end
    if ~isnan(uCode(kk))
        if upper(Des{kk}(1))=='R'
            % To differentiate stuff starting with 'R' from those with 'S'
            % 'S' being the usual 1st character for a stimulus.
            uCode(kk) = -uCode(kk);
        end
        Code(strcmpi(marker.Descript,Des{kk}))=uCode(kk);
    else
        uCode(kk) = base;
        Code(strcmpi(marker.Descript,Des{kk})) = uCode(kk);
        base = base+1;
    end
    fprintf('\tTurning marker "%s" : "%s" into code : "%s" \n', ...
        marker.Type{u_ind(kk)}, Des{kk}, num2str(uCode(kk)))
end

blankevent = struct('type','Unknown','value',[],'time',0, ...
    'duration',0,'offset',0);

% create 1st event = trial/session
events = blankevent;
events.type = 'Trial';
events.value = 0;
events.duration = D.Nsamples;
events.time = 1/D.Fsample;

if length(l_NewSegment)==1
    % only one segment, proceed as originally intended, i.e. no trouble
    % Other real events
    for ii=1:marker.nMrk
        events(ii+1) = blankevent;
        events(ii+1).type = marker.Type{ii};
        events(ii+1).value = Code(ii);
        events(ii+1).time = marker.Position(ii)/D.Fsample;
    end
else
    % EITHER
    %-------
    % Deal with unusual case of multiple data segments because of some
    % pause or disconnections (cf. C. Schmidt data)
    % Need to rebuild the whole data set with zeros filling the interval.
    % OR
    %---
    % Consider only the last bit
    for ii= (l_NewSegment(end):marker.nMrk)-l_NewSegment(end)+1
        jj = ii+l_NewSegment(end)-1;
        events(jj+1) = blankevent;
        events(jj+1).type = marker.Type{jj};
        events(jj+1).value = Code(jj);
        events(jj+1).time = marker.Position(jj)/D.Fsample;
    end
end

trials = struct('label','Undefined','onset',1/D.Fsample, ...
    'repl',1,'bad',0,'events',events);

% A priori users of this toolbox use the native continuous recording...
% Looks like that in SPM8/12, onset of trials is set at the time of the 1st
% event...

D.trials = trials;

%% Deal now with filter, reref, transform and descript field, info
D.fiducials = struct([]);
D.sensors   = struct([]);
D.artifacts = struct([]);
% These fields can be updated later within SPM8/12, when doing source
% reconstruction or DCM-EEG
D.history   = struct('fun','crc_eeg_rdata_brpr', ...
    'args',struct('D',struct('fname',Fdata)));
D.transform.ID = 'time';    % Time series not frequency spectrogram

%% Handling scale

% D.data.scale = header.channels.Resolution;
if all(header.channels.Resolution~=0)
    D.data.y.scl_slope = header.channels.Resolution;
end
% D.data.y.fname = header.DataFile;

%% treat date and hour of recording beginning of recording
% Very useful for sleep studies !
% turn into meeg object and add info
D = meeg(D);
D.info = header.info; % just in case, save the ASCII header
try
    text_dt = marker.Date{1};
    year = str2double(text_dt(1:4));
    month = str2double(text_dt(5:6));
    day = str2double(text_dt(7:8));
    hour = str2double(text_dt(9:10));
    minutes = str2double(text_dt(11:12));
    secondes = str2double(text_dt(13:end))/1e6;
    
    D.info.date = [year month day];
    D.info.hour = [hour minutes secondes];
catch  %#ok<*CTCH>
    warning('FASST:readBPdata', ...
        'Could not extract the date/time of recording.')
end
[v,r] = crc_fasst_utils('Ver',mfilename);
D.info.ver = struct('v_nr',v,'rel',r);
[v,r] = crc_fasst_utils('Ver','crc_eeg_readHeaderBrainP');
D.info.ver_readH = struct('v_nr',v,'rel',r);

%% Save results
save(D);

return
