function newD = crc_spm5tospm8(P)

% Routine to quickly convert SPM5 data into SPM8 format.
%
% WARNING: 
% - not working with ERP data but only continuous data !!!
% - may not be entirely up to date...
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$


if ~exist('CRC_electrodes.mat','file')
    crc_electrodes;
end
load('CRC_electrodes.mat')
    
pos=pos';

%% check filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>0 && isstruct(P)
    D=P;
else
    
    try
        P = deblank(P);
    catch
        P = spm_select(1, '\.mat$', 'Select EEG mat file');
    end

    
    %Find path and filename
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ppath = spm_str_manip(P, 'H');
    if strcmp('.', Ppath) | strcmp('..', Ppath)
        Ppath = pwd;
    end
    matname = spm_str_manip(P, 't')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % load file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        load(P);
    catch
        error(sprintf('Trouble reading file %s', P));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% Create blank spm8 structure

newD = struct('type',[],'Nsamples',[],'Fsample',[],'timeOnset',[], ...
    'trials',[], 'channels', struct, 'data', [], 'fname', [], 'path', [], ...
    'sensors', [], 'fiducials', [], 'artifacts', [], 'transform', [], ...
    'other', [], 'history', []);

%% Handle the field trial
Newtrial = struct('label','Unknown','onset',0,'repl',0,'bad',0,'events',[]);

for ii = 1:D.Nevents
    temptrial = Newtrial;
    % Maybe some code should be added here for ERP data
    newD.trials = temptrial;
end

if D.Nevents == 1

    newD.type = 'continuous';
    newD.timeOnset = 0;

    newD.trials.events = struct('type', 'Unknown', ...
        'value', 0, 'time', 0,...
        'duration', eps);

    for ii=1:length(D.events.code)
        newD.trials.events(ii)=struct('type', 'Unknown', ...
            'value', D.events.code(ii), 'time', D.events.time(ii)/D.Radc,...
            'duration', eps); %Duration fixed arbitrarily at a very small value
    end

else

    newD.type = 'evoked';
    newD.timeOnset = -0.1; % Arbitrarily chosen as standard

    %Need to implement how to handle events with ERP

end


%% Deal with data

newD.data.y  = D.data;
newD.data.fnamedat = D.fnamedat;
newD.data.datatype = D.datatype;
if strcmp(newD.data.datatype,'float')
    newD.data.datatype='float32';
end
% newD.data.scale = D.scale; %%% Not sure about this one, scale seems different in spm8

%% Number of samples and sampling frequency

newD.Nsamples = D.Nsamples;
newD.Fsample = D.Radc;

%% Path & filename

try
    newD.path = Ppath;
    newD.fname = matname;
catch
    newD.path=D.path;
    newD.fname=D.fname;
end

%% Deal with channels

blankchan = struct('label',[],'bad', [], 'type', 'unknown', ...
    'X_plot2D', [],'Y_plot2D',[], 'units', 'unknown');

newD.channels = blankchan;

for ii=1:length(D.channels.name)
    tempchan = blankchan;
    tempchan.label = D.channels.name{ii};
    if isfield(D.channels, 'bad')
        tempchan.bad = ismember(ii, D.channels.bad);
    else
        tempchan.bad = 0; % Not a bad channel by default
    end

    if isfield(D.channels, 'eeg')

        if ismember(ii, D.channels.eeg)
            tempchan.type = 'EEG';
        end
    else
        if findstr(D.channels.name{ii},'EOG')
            tempchan.type = 'EOG';
        elseif findstr(D.channels.name{ii},'EMG')
            tempchan.type = 'EMG';
        elseif findstr(D.channels.name{ii},'ECG')
            tempchan.type = 'ECG';
        else
            tempchan.type = 'EEG';
        end

    end

    if D.channels.order(ii)~=0
        tempchan.X_plot2D = pos(2,D.channels.order(ii));
        tempchan.Y_plot2D = pos(1,D.channels.order(ii));
    else
        tempchan.X_plot2D = NaN;
        tempchan.Y_plot2D = NaN;
    end

    newD.channels(ii) = tempchan;

end

newD.channels = newD.channels';

%% Deal now with filter, reref, transform and descript field, info and CRC if existing

newD.fiducials=struct([]);
newD.sensors=struct([]);
newD.artifacts=struct([]);
newD.history=struct([]);


newD.transform.ID = 'time';

try
    if ~isempty(D.filter)
        newD.other.filter = D.filter;
    end
end

try
    if ~isempty(D.reref)
        newD.other.reref = D.reref;
    end
end

try
    if ~isempty(D.info)
        newD.other.info = D.info;
    end
end

try
    if ~isempty(D.CRC)
        newD.other.CRC = D.CRC;
    end
end

%% If exists CRC.art field, copy result in artifact struct array

if isfield(D,'CRC')
    if isfield(D.CRC,'art')
        newD.artifacts = struct('start', 0,'stop',0);
        if ~(size(D.CRC.art,1) == 0)
            for ii = size(D.CRC.art,1)
                newD.artifacts(ii) = struct('start', D.CRC.art(ii,1),...
                    'stop', D.CRC.art(ii,2));
            end
        end
    end
end
