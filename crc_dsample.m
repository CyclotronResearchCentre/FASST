function D = crc_dsample(prefile,fs)

% FORMAT D = crc_dsample(prefile,fs)
% Routine to subsample a data file
% Input:
%   prefile - EEG data file to subsample, it can either be an imported SPM
%             format data set, or the raw BrainProducts data file.
%   fs      - new sampling rate.
%
% Output:
%   D       - downsampled data file
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liège, Belgium
% $Id$

crcdef = crc_get_defaults('ds');

if nargin<2
    fs = crcdef.fs;
end

% Definition of the file to treat
if nargin<1
    D = crc_eeg_load;
else
    D = crc_eeg_load(prefile);
end
file = fullfile(D.path,D.fname);

if ~exist(file,'file'),
    warning('No valid file selected')
    return
end

% Use SPM function:
S.D = file;
S.fsample_new = fs;
S.prefix = crcdef.prefix;
D = spm_eeg_downsample(S);

