function Do = crc_concatenate(prefile, prefix)
%
% FORMAT fname_out = crc_concatenate(prefile, prefix)
%
% Concatenate 2 M/EEG data files in @meeg format, taking into account
% their recording time, if available.
% The period without data between 2 successive files is left blank, if real
% world recording time is available. Otherwise files are simply appended.
% Very convenient when recording was interupted during the night...
%
% INPUT
% - prefile     file name of 2 data files to concatenate. They could be
%               "raw" (not imported yet) or @meeg data files.
% - prefix      prefix prepended to the 1st file name for the output
%               concatenated data ('cc_' by default).
%
% OUTPUT
% - Do          output data file
%
% NOTE:
% - so far the routine is limited to concatenating 2 files but could
%   be extended to any number of files.
% - files are assumed to be continuous data with a single trial
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

if nargin<2
    prefix = crc_get_defaults('cat.prefix1');
end

%% Loading data
if nargin<1
    prefile = spm_select(2, 'any', 'Select imported EEG file','' ,pwd, ...
        '\.[mMvVeE][dDhHaA][fFDdTt]');
end
D1 = crc_eeg_load(deblank(prefile(1,:)));
D2 = crc_eeg_load(deblank(prefile(2,:)));

%% Checking stuff
% And figuring out which file is the first and compute the time between the
% end of the first and the beginning of the second, if possible

if D1.fsample~=D2.fsample || D1.nchannels~=D2.nchannels
    error('Incompatible files!')
end

bconc=0;
try
    diffdate = datenum([D1.info.date D1.info.hour]) -  ...
                datenum([D2.info.date D2.info.hour]);
catch
    diffdate=NaN;
    txt = sprintf(['Warning: Unable to find the time information of the files provided.\n',...
        '\t => Concatenating files back to back']);
    disp(txt)
    bconc = 1;
end

% switch order if needed
if diffdate>0
    Dtmp = D1;
    D1 = D2;
    D2 = Dtmp;
elseif diffdate<0
    diffdate = -diffdate;
end

% check there isn't some overlap, in case of absolute times
% and get gap in time bins
if bconc
    tb_gap = 0;
else
    t_dur1 = D1.nsamples/D1.fsample;
    if (t_dur1-diffdate*24*3600)>1e4
        error('Duration of file1 is bigger than interval between beginning of the 2 files.')
    end
    t_gap = diffdate*24*3600-t_dur1;
    tb_gap = round(t_gap*D1.fsample);
end

Ntb_o = D1.nsamples + D2.nsamples + tb_gap;

fnamedat_o = fullfile(D1.path,[prefix,spm_str_manip(D1.fnamedat,'t')]);
Do = clone(D1, fnamedat_o, [D1.nchannels Ntb_o 1]);

% Copy D1 data, then D2
MaxMem = crc_get_defaults('mem.maxmemload');
Mtb_chan = MaxMem/8/D1.nchannels;
Nchunks = ceil(D1.nsamples/Mtb_chan);
for ii=1:Nchunks
    if ii==Nchunks
        % last chunk
        ind_b = (ii-1)*Mtb_chan+1;
        Do(:,ind_b:D1.nsamples,1) = D1(:,ind_b:D1.nsamples,1);
    else
        % full chunks
        ind_b = (ii-1)*Mtb_chan+1;
        ind_e = ii*Mtb_chan;
        Do(:,ind_b:ind_e,1) = D1(:,ind_b:ind_e,1);
    end
end
Nchunks = ceil(D2.nsamples/Mtb_chan);
tb_shift = D1.nsamples+tb_gap;
for ii=1:Nchunks
    if ii==Nchunks
        % last chunk
        ind_b = (ii-1)*Mtb_chan+1;
        Do(:,tb_shift+(ind_b:D2.nsamples),1) = D2(:,ind_b:D2.nsamples,1);
    else
        % full chunks
        ind_b = (ii-1)*Mtb_chan+1;
        ind_e = ii*Mtb_chan;
        Do(:,tb_shift+(ind_b:ind_e),1) = D2(:,ind_b:ind_e,1);
    end
end

%% handles events time,
ev1 = D1.events; Nv1 = numel(ev1);
ev2 = D2.events; Nv2 = numel(ev2);
t_shift = tb_shift/D1.fsample;
evo = [ev1 ev2];
if Nv2
    for ii=1:Nv2
        evo(Nv1+ii).time = ev2(ii).time + t_shift;
    end
end
Do = events(Do,1,evo);
save(Do)

% fname_out = fullfile(Do.path,Do.fname);

disp('Concatenation done')