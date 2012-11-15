function [fdata,header] = crc_eeg_readHeaderEGI(filename)

% Reads in EEG data from EGI-raw to SPM8-format filearray format, and
% creates a header structure (not SPM8/@meeg object!) with channel and data
% informations.
%
% FORMAT [fdata,header] = crc_readegi(filename)
%
% filename    - filename(s) of EGI-raw file
% fdata       - filename (including path) to binary EEG data file
% header      - header information
%       header.version  version number, even/odd for continuous/segmented
%       header.date     date [year month day] \
%       header.time 	time [hour min sec]   | of recording
%       header.msec     time [ms]             /
%       header.Srate 	sampling rate (Hz)
%       header.Nchan    # channels
%       header.gain 	board gain (1, 2, 4, 8)     (???)
%       header.Cbits    # conversion bits           (???)
%       header.Scrange 	Full-scale range of amp in µV
%       header.Nsamples # samples (time-bins)
%       header.Nuevents # unique event codes (0 is no events)
%       header.ecode 	list of event codes (Nuevent x 4 char array)
%       header.fodat 	data format ('int16','float32','float64')
%       header.event_tb list of event time bin of occurence
%       header.event_co list of corresponding event code index
%       header.event_du list of event duration in time bins
%
% If multiple files are passed/selected, then it is assumed that they form
% a series of successive continuous files that were chopped at the
% exportation stage. For example, multi-hours sleep recording data ending 
% up in 5 separate files. Those files will be automatically concatenated
% when reading in the data.
%
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2009.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

if nargin < 1
    filename = spm_select([1 Inf],'^.*\.[rR][aA][wW]$', 'Select EGI-raw file');
end;

[fp,hdr] = readHdr(filename);
fn1 = deblank(filename(1,:));
Nfn = length(fp);
header = hdr(1);

%% Data file to build
[pth,fname,ext] = fileparts(fn1); %#ok<NASGU>
if Nfn>1
    fdata = fullfile(pth,[crc_get_defaults('cat.prefix2'),fname,'.dat']);
else
    fdata = fullfile(pth,[fname,'.dat']);
end
fo = fopen(fdata,'w','ieee-le');
event_tb = [];
event_co = [];
for jj = 1:Nfn
    for ii = 1:hdr(jj).Nsamples
        tmp_d = fread(fp(jj),hdr(jj).Nchan,hdr(jj).fodat);
        fwrite(fo,tmp_d,hdr(jj).fodat);
        tmp_e = fread(fp(jj),hdr(jj).Nuevents,hdr(jj).fodat);
        if any(tmp_e)
            % account for possibility of multiple stim at same time bin
            for kk = find(tmp_e')
                event_tb = [event_tb ; ii]; %#ok<AGROW>
                event_co = [event_co ; kk]; %#ok<AGROW>
            end
        end
    end
    if jj>1
        % accumulate #samples
        header.Nsamples = header.Nsamples + hdr(jj).Nsamples;
    end
end
%% Deal with continuous stimulations, 
% i.e. event of one type lasting multiple consecutive time-bins.
header.event_tb = [];
header.event_co = [];
header.event_du = [];
if header.Nuevents ~=0
    for ii=1:header.Nuevents
        le_ii = find(event_co==ii);
        tb_e_ii = event_tb(le_ii); %#ok<FNDSB>
        l_tb_diff = diff(tb_e_ii);
        le_ii_rem = [] ; % consecutive evens to remove
        dur = [] ; % duration of events
        dur_jj = 1;
        for jj=1:length(l_tb_diff)
            if l_tb_diff(jj)==1
                % next event same as previous one
                % - next one to be removed
                % - duration of event block increaed by 1
                le_ii_rem = [le_ii_rem ; jj+1]; %#ok<AGROW>
                dur_jj = dur_jj+1;
            else
                % next event different
                % - save current duration index
                % - and reset it to 1
                dur = [dur ; dur_jj]; %#ok<AGROW>
                dur_jj = 1;
            end
        end
        % Keep track of last event duration
        dur = [dur ; dur_jj]; %#ok<AGROW>
        % remove individual 'events' when they form a continuous block
        event_tb(le_ii(le_ii_rem)) = []; %#ok<AGROW>
        event_co(le_ii(le_ii_rem)) = []; %#ok<AGROW>
        % keep track of stim duration
        le_ii_keep = le_ii; le_ii_keep((le_ii_rem)) = [];
        event_du(le_ii_keep) = dur'; %#ok<AGROW>
    end
    
    header.event_tb = event_tb;
    header.event_co = event_co;
    header.event_du = event_du';
end

fclose(fo);
for ii=1:Nfn
    fclose(fp(ii));
end

% EGI doesn't provide channel names but just their number,
% so use this as channel names.
header.channames = num2str((1:header.Nchan)');

return

%==========================================================================
%% SUBFUNCTION
% Reading the 1st part of the raw file, i.e. header.

function [fp,header] = readHdr(filename)

Nf = size(filename,1);

fp = zeros(Nf,1);
hdr = struct( ...
    'version',[],'date',[],'time',[],'msec',[],'Srate',[],'Nchan',[], ...
    'gain',[],'Cbits',[],'Scrange',[],'Nsamples',[],'Nuevents',[], ...
    'ecode','','fodat','');
header(Nf) = hdr;

for ii=1:Nf
    fp(ii) = fopen(deblank(filename(ii,:)),'r','ieee-be');
    % First assume mac big-endian storage format
    if fp(ii) == -1,
        error('File not found ...!');
    end
    
    version = fread(fp(ii),1,'int32');
    
    if version>7 || version<2
        % Problem with the version number, it should be between 2 and 7, in
        % this case, take the bytes in "little-endian" order.
        fclose(fp(ii));
        fp(ii) = fopen(deblank(filename(ii,:)),'r','ieee-le');
        version = fread(fp(ii),1,'int32');
    end
    
    if rem(version,2)
        error('Can''t deal with segmented data');
    end
    
    header(ii).version  = version;
    header(ii).date     = fread(fp(ii),3,'int16');
    header(ii).time     = fread(fp(ii),3,'int16');
    header(ii).msec     = fread(fp(ii),1,'int32');
    header(ii).Srate    = fread(fp(ii),1,'int16');
    header(ii).Nchan    = fread(fp(ii),1,'int16');
    header(ii).gain     = fread(fp(ii),1,'int16');
    header(ii).Cbits    = fread(fp(ii),1,'int16');
    header(ii).Scrange  = fread(fp(ii),1,'int16');
    header(ii).Nsamples = fread(fp(ii),1,'int32');
    header(ii).Nuevents = fread(fp(ii),1,'int16');
    if header(ii).Nuevents ~=0
        for jj=1:header(ii).Nuevents
            ecode(jj,:) = char(fread(fp(ii),4,'uchar')); %#ok<FDEPR,AGROW>
        end
        header(ii).ecode = ecode;
    end
    
    switch version
        case 2, fodat = 'int16';
        case 4, fodat = 'float32';
        case 6, fodat = 'float64';
        otherwise
            error('Don''t know this data format, sorry.');
    end
    header(ii).fodat = fodat;
end

if Nf>1 % Checking if all files are "similar"
    Nch = [header(:).Nchan];
    if any(diff(Nch)),
        error('All selected files should have same number of channels!')
    end
    vers = [header(:).version];
    if any(diff(vers)),
        error('All selected files should have same version flag!')
    end
    Nuev = [header(:).Nuevents];
    if any(diff(Nuev)),
        error('All selected files should have same number of ''unique events''!')
    end
    Sr = [header(:).Srate];
    if any(diff(Sr)),
        error('All selected files should have same sampling rate!')
    end
end

return
