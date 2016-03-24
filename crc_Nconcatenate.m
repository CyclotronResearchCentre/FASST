function Do = crc_Nconcatenate(Pfile, prefix)
%
% FORMAT fname_out = crc_Nconcatenate(Pfile, prefix)
%
% Concatenate 'N' M/EEG data files in @meeg format, taking into account
% their recording time, if available.
% When the starting time of each file is availabel, the interval without 
% data between 2 successive files is filed with 0's. Otherwise files are 
% simply appended.
% Very convenient when recording was interupted during the night...
%
% INPUT
% - Pfile   file name of the N data files to concatenate. They could be
%           "raw" (not imported yet) or @meeg data files.
% - prefix  prefix prepended to the 1st file name for the output
%           concatenated data ('cc_' by default).
%
% OUTPUT
% - Do          output data file
%
% NOTE:
% - files are assumed to be continuous data with a single trial
%__________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by C. Phillips, 2016.
% Cyclotron Research Centre, University of Liege, Belgium
% but derived from the crc_concatenate.m function.

%% Ensuring this works without explicitly launching FASST
persistent fasst_defs
if isempty(fasst_defs) || ~fasst_defs
    fasst_defs = crc_main('SetDefs');
end

%% Deal with options
if nargin<2
    try
        prefix = crc_get_defaults('cat.prefix1');
    catch %#ok<*CTCH>
        prefix = 'cc_';
    end
end

% In the case of BrainProducts data, force the renaming of the data file
% from datafile.eeg into datafile.dat (true) or not (false).
opt = true;

%% Loading data
if nargin<1
    Pfile = spm_select([2 Inf], 'any', 'Select imported EEG file','' ,pwd, ...
        '\.[mMvVeE][dDhHaA][fFDdTt]');
end
Nfile = size(Pfile,1);
Di = cell(Nfile,1);
for ii=1:Nfile
    Di{ii} = crc_eeg_load(deblank(Pfile(ii,:)),opt);
end

%% Checking and sorting stuff
% Basic check in term of samplign rate and nuber of channels available, 
% compared to the 1st file selected.
% And figuring out which file is the first and compute the time between the
% end of the first and the beginning of the second, if possible
%
% TO DO: A more fancied check could be useful...

for ii=2:Nfile
    if Di{1}.fsample~=Di{ii}.fsample || Di{1}.nchannels~=Di{ii}.nchannels
        error('FASST:Nconcat','File %s not compatible with %s.', ...
            Di{1}.fname, Di{ii}.fname)
    end
end

bconc = false; % Assume we have gaps between files

% Get the dates out, if available in all files
Fdate = zeros(Nfile,1);
for ii=1:Nfile
    if isfield(Di{ii},'info') && isfield(Di{ii}.info,'date') && isfield(Di{ii}.info,'hour')
        Fdate(ii) = datenum([Di{ii}.info.date Di{ii}.info.hour]);
    else
        Fdate(ii) = NaN;
    end
end

FdurTb = zeros(Nfile,1);
if any(isnan(Fdate))
    % Some date-time missing, just concatenate with files as they're
    bconc = true;
    Fdate(:) = NaN; %#ok<*NASGU>
    FdurTb(:) = NaN;
    Dis = Di;
    fprintf(['Warning: Unable to find the time information of the files provided.\n',...
        '\t => Concatenating files back to back']);
else
    % Proceed with the date-time provided
    % and reorder
    [Fdate,pFile] = sort(Fdate);
    Dis = cell(Nfile,1);
    for ii=1:Nfile;
        Dis{ii} = Di{pFile(ii)};
        FdurTb(ii) = Dis{ii}.nsamples; % expressed in time bins
    end
    DFdate = diff(Fdate)*24*3600; % Time difference, in sec, between 2 recordings.
end
FdurTs = FdurTb/Dis{1}.fsample; % expressed in sec

% Check there isn't some overlap, in case of date-time
% and get gap in time bins
GapTs = zeros(Nfile-1,1); % Gap is 0 by default
GapTb = zeros(Nfile-1,1); % Gap is 0 by default
if ~bconc % Find the gap between files
    for ii=1:Nfile-1
        GapTs(ii) = DFdate(ii)-FdurTs(ii); % gap expressed in sec
        GapTb(ii) = round(GapTs(ii)*Dis{1}.fsample); % gap expressed in time bins
        if GapTb(ii)<0
            error('FASST:Nconcat','Negative gap between %s and %s.', ...
                Di{ii}.fname, Di{ii+1}.fname)
        end
    end
end

% Prepare whole data file:
NtotTb = sum(FdurTb)+sum(GapTb);
fnamedat_o = fullfile(Dis{1}.path,[prefix,spm_str_manip(Dis{1}.fnamedat,'t')]);
Do = clone(Dis{1}, fnamedat_o, [Dis{1}.nchannels NtotTb 1]);

% Copy files with gap, work in chunks
MaxMem = crc_get_defaults('mem.maxmemload');
Mtb_chan = MaxMem/8/Dis{1}.nchannels;

tb_shift = 0; % Shifting of time bins
for ii=1:Nfile % concatenate all files!
    Nchunks = ceil(Dis{ii}.nsamples/Mtb_chan);
    for jj=1:Nchunks
        if jj==Nchunks
            % last chunk
            ind_b = (jj-1)*Mtb_chan+1;
            Do(:,tb_shift+(ind_b:Dis{ii}.nsamples),1) = Dis{ii}(:,ind_b:Dis{ii}.nsamples,1);
        else
            % full chunks
            ind_b = (jj-1)*Mtb_chan+1;
            ind_e = jj*Mtb_chan;
            Do(:,tb_shift+(ind_b:ind_e),1) = Dis{ii}(:,ind_b:ind_e,1);
        end
    end
    tb_shift = tb_shift + Dis{ii}.nsamples; % account for the i_th data set
    if ii<Nfile % Put 0's in gap but not after the last file
        Nchunks = ceil(GapTb(ii)/Mtb_chan);
        for jj=1:Nchunks
            if jj==Nchunks
                % last chunk
                ind_b = (jj-1)*Mtb_chan+1;
                Do(:,tb_shift+(ind_b:GapTb(ii)),1) = 0;
            else
                % full chunks
                ind_b = (jj-1)*Mtb_chan+1;
                ind_e = jj*Mtb_chan;
                Do(:,tb_shift+(ind_b:ind_e),1) = 0;
            end
        end       
        tb_shift = tb_shift + GapTb(ii); % include gap after the i_th data set
    end
end

%% handles events time,
evo = Dis{1}.events;
ts_shift = 0;
for ii=2:Nfile % add events & adjust times
    ts_shift = ts_shift + Dis{ii-1}.nsamples/Dis{1}.fsample + GapTs(ii-1); 
    ev_ii = Dis{ii}.events;
    for jj=1:numel(ev_ii)
        ev_ii(jj).time = ev_ii(jj).time + ts_shift;
    end
    evo = [evo ev_ii]; %#ok<*AGROW>
end

Do = events(Do,1,evo);
save(Do)
% fname_out = fullfile(Do.path,Do.fname);

disp('Concatenation done')
end


