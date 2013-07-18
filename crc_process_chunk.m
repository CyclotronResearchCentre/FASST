function fn_chunk = crc_process_chunk(Dmeg,Begpts,Endpts, flags)
% fn_chunk = crc_process_chunk(Dmeg,Begpts,Endpts, flags)
%
% Chunks MEEG object Dmeg with start at Begpts and stop at Endpts, both of 
% these expressed in "time points".
% 'flags' is a structure with several optional fields:
% flags
%   .overwr   : previous chunks are overwritten [1] or not [0, def.]
%   .numchunk : fix the chunk index, useful when overwriting. Ortherwise 
%               the index is incremented from already existing files on disk
%               The default value is 1
%   .prefix   : filename prefic, 'chk' by default
%   .clockt   : use the clocktime [1] or not [0, def.]
% All default values are defined in crc_defaults.m
%
% The chunked data file is written on disk, and its filename returned.
%__________________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

flags_def = crc_get_defaults('chk');
flags = crc_check_flag(flags_def,flags);

% Identifying file chunk number, overwriting or not
numchk = flags.numchunk;
if ~flags.overwr
    while exist(fullfile( ...
                path(Dmeg),[flags.prefix,num2str(numchk),'_' fname(Dmeg)]),'file')
        numchk=numchk+1;
    end
end
prefix = [flags.prefix,num2str(numchk),'_'];
fn_chunk = fullfile(Dmeg.path,[prefix Dmeg.fname]);

try
    winsize = median([Dmeg.CRC.score{3,:}]);

    firstwin = floor(Begpts/(fsample(Dmeg)*winsize))+1;
    Begpts = (firstwin-1)*winsize*fsample(Dmeg)+1;

    lastwin = ceil(Endpts/(fsample(Dmeg)*winsize));
    Endpts = lastwin*winsize*fsample(Dmeg);
end

Dcopy = clone(Dmeg,fn_chunk,[Dmeg.nchannels Endpts-Begpts+1 1]);

%==========================================================================
% DIRTY CODE :
% getting rid of offset in orginal data file.
% This should actually be part of the clone function or at least handled
% somewhere in the object definition and methods...
% => have to use the structure form
% DO NOT DO THIS AT HOME. :-)
Dcs = struct(Dcopy);
Dcs.data.y.offset = 0;
Dcopy = meeg(Dcs);
%==========================================================================

try

    for nsc=1:size(Dmeg.CRC.score,2)
        Dcopy.CRC.score{1,nsc} = Dmeg.CRC.score{1,nsc}(firstwin:lastwin);
        Dcopy.CRC.score{4,nsc} = [1/fsample(Dmeg) (Endpts-Begpts-1)/fsample(Dmeg)];

        tmp = Dmeg.CRC.score{5,nsc}-Begpts/fsample(Dmeg);
        idxkept = find(~all((tmp<0) + (tmp>(Endpts-Begpts)/fsample(Dmeg)),2));
        tmp = tmp(idxkept,:);
        tmp(tmp<0) = 1/fsample(Dmeg);
        tmp(tmp>(Endpts-Begpts)/fsample(Dmeg)) = (Endpts-Begpts-1)/fsample(Dmeg);
        Dcopy.CRC.score{5,nsc} = tmp;
        
        if size(Dmeg.CRC.score,1)==8 ...
                && ~isempty(Dmeg.CRC.score{8,nsc})
            if ~isempty(idxkept)
                Dcopy.CRC.score{8,nsc} = {Dmeg.CRC.score{8,nsc}{idxkept}};
            else
                Dcopy.CRC.score{8,nsc} = {};
            end
        end

        tmp = Dmeg.CRC.score{6,nsc}-Begpts/fsample(Dmeg);
        idxkept = find(~all((tmp<0) + (tmp>(Endpts-Begpts)/fsample(Dmeg)),2));
        tmp = tmp(idxkept,:);
        tmp(tmp<0) = 1/fsample(Dmeg);
        tmp(tmp>(Endpts-Begpts)/fsample(Dmeg)) = (Endpts-Begpts-1)/fsample(Dmeg);
        Dcopy.CRC.score{6,nsc}=tmp;

        tmp = Dmeg.CRC.score{7,nsc}-Begpts/fsample(Dmeg);
        idxkept = find(~all((tmp<0) + (tmp>(Endpts-Begpts)/fsample(Dmeg)),2));
        tmp = tmp(idxkept,:);
        tmp(tmp<0) = 1/fsample(Dmeg);
        tmp(tmp>(Endpts-Begpts)/fsample(Dmeg)) = (Endpts-Begpts-1)/fsample(Dmeg);
        Dcopy.CRC.score{7,nsc}=tmp;
    end
end

% Limit size of data chunk loaded in memory in one step.
Maxmem = crc_get_defaults('mem.maxmemload'); % Maxmem in byte

% Assuming data are in 32 bits/ 4 bytes
Maxelts = floor(Maxmem/(nchannels(Dmeg)*4));
Nmbchk  = ceil((Endpts-Begpts)/Maxelts);

% Update structure and save file
% handling the starting time information if available
if flags.clockt
    hour = Dcopy.info.hour;
    date = Dcopy.info.date;

    timeshift = ...
        datenum(double([0 0 0 crc_time_converts(Begpts/fsample(Dmeg))]));
    newbeg = datevec(datenum(double([date hour])) + timeshift);

    Dcopy.info.hour = newbeg(4:6);
    Dcopy.info.date = newbeg(1:3);
end

% Not too sure about his about handling the events/trials !!!!
% Should check with a file with loads of events...
Ev = Dcopy.events;
l_keep = []; % Look for events within the new file
for mm=1:numel(Ev)
    Ev(mm).time = Ev(mm).time - Begpts/fsample(Dmeg);
    if Ev(mm).time>=0 && Ev(mm).time<(Endpts-Begpts)/fsample(Dmeg)
        l_keep = [l_keep mm];
    end
end
Dcopy = events(Dcopy,1,Ev(l_keep));

h = waitbar(0,'Please wait...');

for ii=1:Nmbchk
    string = ['Please wait... ' num2str(round(100*((ii-1)/Nmbchk))) ' %'];
    waitbar((ii-1)/Nmbchk,h,string)
    
    ch_data = Dmeg(:,(Begpts+Maxelts*(ii-1)): ...
                        min([Begpts+Maxelts*ii-1 nsamples(Dmeg) Endpts]));
    Dcopy(:,Maxelts*(ii-1)+(1:size(ch_data,2))) = ch_data;
end

string = ['Please wait... ' num2str(100*1) ' %'];
waitbar(1,h,string)
close(h)

save(Dcopy)