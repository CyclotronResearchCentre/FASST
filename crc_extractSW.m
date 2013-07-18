function [SW, cleanSW]= crc_extractSW(Ds,an,stages,scorer,flag)

% Function to extract SWS stages of whole nights data sets.
%
% INPUT
%   Ds     : scored whole night recording and scores to extract (vector or 
%            double).
%   an     : analysis type, by default 3
%   stages : stages to be used, by default [3 4]
%   scorer : index of scorer to use, by default 1
%   flag   : saving (1, default) or not (0) the SW episode on disk in a
%            seperate data file.
%
% OUTPUT 
%  - SW: new file with prefix SWS_ containing only SWS stages, without
%        artefacts or arousals
%  - cleanSW: vector containing the time point position in the original 
%        file of each selected data point
%  
%__________________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre
%
% Written by J. Schrouff & C. Phillips, 2010.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

%load data and get defaults
if nargin<1   
    file = spm_select(1, 'mat', 'Select SPM M/EEG file');
    Ds=crc_eeg_load(file);
end
if nargin<2
    an = 3;
end
if nargin<3
    stages = [3 4];
end
if nargin<4
    scorer = 1;
end
if nargin<5
    flag = 1;
end
% if file not scored, return
if ~isfield(Ds,'CRC') && an==3
    warning('File not scored: please score sleep using FASST!')
    return
elseif (isfield(Ds,'CRC') && ~isfield(Ds.CRC, 'score') && an<3) || ...
        (~isfield(Ds,'CRC') && an<3)
    SW = Ds;
    cleanSW = 1:Ds.nsamples;
    return
end

%% Get data index
% get length of time window used to score
sc   = Ds.CRC.score;
lw   = sc{3,scorer};
rate = Ds.fsample;

sww=[];
% get windows scored as wished by the user
for i=1:numel(stages)
    sw  = find(sc{1,scorer}==stages(i));
    sww = union(sww,sw);
end

% build vector of corresponding time bins
swb = zeros(lw*rate, numel(sww));
tmp = (1:lw*rate)';
for ii=1:numel(sww)
    swb(:,ii) = tmp +(sww(ii)-1)*lw*rate;
end
swb = swb(:);
stoppt = min(max(swb),nsamples(Ds));
swb = swb(swb<=stoppt);

%% get artefacts and arousals in considered windows
% artefacts
artb=[];
if ~isempty(sc{5,scorer})
    for i=1:size(sc{5,scorer},1)
        art=floor(sc{5,scorer}(i,1)*rate):ceil(sc{5,scorer}(i,2)*rate);
        if ~isempty(intersect(art,swb))
            artb=[artb; art'];
        end
    end
end

% arousals
arou=[];
if ~isempty(sc{6,scorer})
    for i=1:size(sc{6,scorer},1)
        aro=floor(sc{6,scorer}(i,1)*rate):ceil(sc{6,scorer}(i,2)*rate);
        if ~isempty(intersect(aro,swb))
            arou=[arou; aro'];
        end
    end
end

% Combining artefacts & arousal
if ~isempty(artb) && ~isempty (arou)
    disc = union(artb, arou);
elseif ~isempty(artb) && isempty (arou)
    disc = unique(artb);
elseif isempty(artb) && ~isempty (arou)
    disc = unique(arou);
else
    disc = [];
end
    
% discard them from selected windows
cleanSW   = setdiff(swb,disc);
szcleansw = length(cleanSW);
if szcleansw >= Ds.nsamples
    flag=0;
end
clear swb artb arou disc sww

if flag
    %save only clean SWS stages, avoiding 'out of memory' errors
    szd = szcleansw*Ds.nchannels*8;
    memsz  = 1024^2*200; % memory use limit, 200Mb here
    numblocks = ceil(szd/memsz);
    sizblocks = ceil(szcleansw/numblocks);

    % get first block to init new structure
    stopdat = min(Ds.nsamples,sizblocks);
    swdat = Ds(:,cleanSW(1:stopdat));
    prefix = 'SW_';

    SW = crc_save_spm(prefix,Ds,swdat);
    clear swdat

    % append other blocks
    idat = sizblocks+1; % index of beginning of next block to write 
    for i=2:numblocks
        stopdat = min(szcleansw,idat+sizblocks-1);
        swdat = Ds(:,cleanSW(idat:stopdat));
        if ~isempty(swdat)
            SW = crc_append_spm(SW,swdat);
        end
        idat=idat+sizblocks;
        clear swdat 
    end
    % Removing redundant stuff
    SW.CRC.score = [];
    SW.info.date = [];
    SW.info.hour = [];
    
    save(SW);
else
    SW=Ds;
end







