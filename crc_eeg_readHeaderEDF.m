function [fdata,header] = crc_eeg_readHeaderEDF(fname)
% FORMAT [fdata,header] = crc_eeg_read_HeaderEDF(fname)
%
% 'crc_eeg_readHeaderEDF' reads the header and converts EDF data.
% It produces 2 outputs:
% - fdata  : filename of binary multiplexed file
% - header : structure with all the dat information
% (see below)
%
% This routine is only a wrap up function around the readedf.m function
% from Jeng-Ren Duann. Still readedf.m was improved and cleaned up a bit.
% See below for the details.
%_______________________________________________________________________
% Copyright (C) 2012 Cyclotron Research Centre

% Written by C. Phillips, 2012.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

if nargin < 1
    fname = spm_select(1, 'any', 'Select EDF file');
end;

[fdata,header] = readedf(fname);

return

%%
% SUBFUNCTION doing the job
%==========================================================================

% readedf() - read eeg data in EDF format.
%
% Usage:
%    >> [fdata,header] = readedf(filename);
%
% Input:
%    filename - file name of the eeg data
%
% Output:
%    fdata   - filename of binary multiplexed data file (see note)
%    header - structured information about the read eeg data
%      header.length - length of header to jump to the first entry of eeg data
%      header.records - how many frames in the eeg data file
%      header.duration - duration (measured in second) of one frame
%      header.channels - channel number in eeg data file
%      header.channelname - channel name
%      header.transducer - type of eeg electrods used to acquire
%      header.physdime - details
%      header.physmin - details
%      header.physmax - details
%      header.digimin - details
%      header.digimax - details
%      header.prefilt - pre-filterization spec
%      header.samplerate - sampling rate
%
% Author: Jeng-Ren Duann, CNL/Salk Inst., 2001-12-21

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Jeng-Ren Duann, CNL/Salk Inst., 2001-12-21
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: readedf.m,v $
% Revision 1.2  2002/08/12 19:00:57  arno
% errordlg->error
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 03-21-02 editing header, add help -ad

% 06-12-20 editing header, read the patient info.
% Christophe Phillips

% 09-09-22
% - correcting header structure, re. nr samples/sampling rate.
% - reading data and saving them into multiplexed data file directly
% Christophe Phillips

% 11-05-23
% - allowing data with different sampling rate per channel
% Christophe Phillips

% NOTE:
% - Data are assumed continuous
% - data are read in then written on disk in a more straightforward format
%   (cfr. improvements on 09-09-22 & 11-05-23)

function [fdata,header] = readedf(filename)

if nargin < 1
    filename = spm_select(1, 'any', 'Select EDF file');
end;

fp = fopen(filename,'r','ieee-le');
if fp == -1,
    error('File not found ...!');
    return;
end

hdr = setstr(fread(fp,256,'uchar')');

% Reading in patient and recording information !
header.version=hdr(1:8);           % 8 Byte  Versionsnummer
header.PID = deblank(hdr(9:88));   % 80 Byte local patient identification
header.RID = deblank(hdr(89:168)); % 80 Byte local recording identification
header.T0=[str2num(hdr(168+[7 8])) str2num(hdr(168+[4 5])) ...
    str2num(hdr(168+[1 2])) str2num(hdr(168+[9 10])) ...
    str2num(hdr(168+[12 13])) str2num(hdr(168+[15 16])) ];

% Y2K compatibility until year 2090
if header.version(1)=='0'
    if header.T0(1) < 91
        header.T0(1)=2000+header.T0(1);
    else
        header.T0(1)=1900+header.T0(1);
    end;
else
    % in a future version, this is hopefully not needed
end;

header.length = str2num(hdr(185:192));   %offset
header.records = str2num(hdr(237:244));  %number of data records
header.duration = str2num(hdr(245:252)); %duration in seconds
header.channels = str2num(hdr(253:256)); %number of channels
header.channelname = setstr(fread(fp,[16,header.channels],'char')');
header.transducer = setstr(fread(fp,[80,header.channels],'char')');
header.physdime = setstr(fread(fp,[8,header.channels],'char')');
header.physmin = str2num(setstr(fread(fp,[8,header.channels],'char')'));
header.physmax = str2num(setstr(fread(fp,[8,header.channels],'char')'));
header.digimin = str2num(setstr(fread(fp,[8,header.channels],'char')'));
header.digimax = str2num(setstr(fread(fp,[8,header.channels],'char')'));
header.prefilt = setstr(fread(fp,[80,header.channels],'char')');
header.nrsamples = str2num(setstr(fread(fp,[8,header.channels],'char')'));
if any(diff(header.nrsamples)),
    var_SR = 1;
    %     error('Not all channels with same nr of samples, can''t deal with this !!!')
    % EDF format allows different sampling rates for different channels,
    % this is NOT supported by SPM8 format
    % So in that case, it is necessary to upsample the channels with a
    % lower sampling rate, wrt. the highest sampling rate in the data.
else
    var_SR = 0;
end

if var_SR
    Mnrsamples = max(header.nrsamples);
    lMnrs = (0:Mnrsamples-1)/Mnrsamples;
    header.samplerate = Mnrsamples/header.duration;
    % Save data in multiplexed format, putting records together after
    % upsampling with nearest-neighbour
    [pth,fname,ext] = fileparts(filename);
    fdata = fullfile(pth,[fname,'.dat']);
    fo = fopen(fdata,'w');

    fseek(fp,header.length,-1);
    for ii=1:header.records
        data_r = zeros(header.channels,Mnrsamples);
        for jj=1:header.channels
            tmp = fread(fp,header.nrsamples(jj),'int16');
            if header.nrsamples(jj)==Mnrsamples
                data_r(jj,:) = tmp;
            else
                lresamp = floor(lMnrs*header.nrsamples(jj))+1;
                data_r(jj,:) = tmp(lresamp);
            end
        end
        fwrite(fo,data_r(:),'int16');
    end
else
    header.samplerate = header.nrsamples(1)/header.duration;
    % Save data in multiplexed format, putting records together
    [pth,fname,ext] = fileparts(filename);
    fdata = fullfile(pth,[fname,'.dat']);
    fo = fopen(fdata,'w');
    
    fseek(fp,header.length,-1);
    for ii=1:header.records
        data_r = fread(fp,header.nrsamples(1)*header.channels,'int16');
        data_r = reshape(data_r,header.nrsamples(1),header.channels)';
        fwrite(fo,data_r(:),'int16');
    end
end

% close all files
fclose(fp);
fclose(fo);

return