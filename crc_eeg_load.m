function D = crc_eeg_load(P,opt)
% Load an M/EEG file in SPM format
% FORMAT D = crc_eeg_load(P,opt)
%
% P         - filename of M/EEG file, or any format supported by FASST
% opt       - option flag (can be omited to ensure back compatibility)
% D         - MEEG object
%__________________________________________________________________________
%
% crc_eeg_load loads an M/EEG file using the SPM8 format. You can either
% select a .mat file (the meeg header file) or any other data format
% suported by FASST (.vhdr/.edf/.raw so far), the conversion will happen on
% the fly.
% Importantly, the data array is memory-mapped and the struct is converted
% to MEEG object.
%
% ONLY for BrainProducts data, the datafile.eeg file can be renamed into
% datafile.dat to fit SPM's usual format (if opt = true).
%
% After using spm_eeg_load, some housekeeping is done to handle FASST
% specific add-ons to the object/structure.
% - if power spectrum data are stored in the CRC subfield of the MEEG
%   object then its data array is checked and
%   memory-mapped. This is useful if data have been moved around!
% - checking if the ROI data extracted during SW detection is there, and
%   removing it.
%
% NOTE: if pwrspect data are updated, the updated MEEG object is saved.
%__________________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by C. Phillips, 2011.
% Cyclotron Research Centre, University of Liege, Belgium

%% Check initilization of FASST/SPM
%--------------------------------------------------------------------------
persistent fasst_defs
if isempty(fasst_defs) || ~fasst_defs
    fasst_defs = crc_main('SetDefs');
end

%% Bypass if the input is already an MEEG object
%--------------------------------------------------------------------------
if nargin && isa(P, 'meeg')
    D = P;
    return;
end
%% Select the data and convert if necessary
%--------------------------------------------------------------------------
% Filter for mat, vhdr, raw and edf files
if nargin<1 || isempty(P)
    P = spm_select(1, 'any', 'Select an EEG file','' ,pwd, ...
        '\.[mMvVeErR][dDhHaA][fFDdTtwW]');
end
if isempty(P), return; end

[pth,name,ext] = fileparts(P);
ext = deblank(ext);

% conversion
if strcmpi(ext,'.vhdr') % BrainProducts
    fn_mat = fullfile(pth,[name,'.mat']);
    if ~exist(fn_mat,'file')
        % Conversion from vhdr to mat file
        % and possibly turning the .eeg file into .dat
        if nargin<2, opt = false; end
        crc_eeg_rdata_brpr(P,opt);
        disp(' ')
        disp('.vhdr file converted to .mat file (spm8/12 compatible)')
    else
        disp(' ')
        disp('.mat file already existing, the conversion of the vhdr file is not needed')
    end

elseif strcmpi(ext,'.edf') % EDF
    fn_mat = fullfile(pth,[name,'.mat']);
    if ~exist(fn_mat,'file')
        % Conversion from edf to mat file
        crc_eeg_rdata_edf(P);
        disp(' ')
        disp('.edf file converted to .mat file (spm8 compatible)')
    else
        disp(' ')
        disp('.mat file already existing, the conversion of the edf file is not needed')
    end

elseif strcmpi(ext,'.raw') % raw EGI
    fn_mat = fullfile(pth,[name,'.mat']);
    if ~exist(fn_mat,'file')
        % Conversion from raw to mat file
        disp(' ')
        disp('Converting .raw EGI file to .mat file (spm8 compatible)')
        crc_eeg_rdata_egi(P);
        disp(' ')
        disp('.raw EGI file converted to .mat file (spm8 compatible)')
    else
        disp(' ')
        disp('.mat file already existing, the conversion of the raw-EGI file is not needed')
    end

elseif strcmpi(ext,'.dat') % BCI2000
    fn_mat = fullfile(pth,[name,'.mat']);
    if ~exist(fn_mat,'file')
        % Conversion from dat to mat file
        disp(' ')
        disp('Converting .dat BCI2000 file to .mat file (spm8 compatible)')
        crc_eeg_rdata_BCI2000(P);
        disp(' ')
        disp('.dat BCI2000 file converted to .mat file (spm8 compatible)')
    else
        disp(' ')
        disp('.mat file already existing, the conversion of the dat-BCI2000 file is not needed')
    end

elseif strcmpi(ext,'.mat')
    % mat file, do nothing
    fn_mat = P;
else
    % wrong file selection
    disp(' ')
    disp('Unknown file extension, can''t load anything.')
    D = [];
    return
end

% then load with spm function
D = spm_eeg_load(fn_mat);

if isempty(D), 
    disp('MEEG object problem!');
    return
end

%% Now dealing with the CRC specific informations
% if it's there...
if isfield(D,'CRC')

    %% Check if power spectrum data exist, update if necessary.
    % First move things in the right subdir (for older data)
    frq_fn = fieldnames(D.CRC);
    l_frq_fn = find(strncmp('frq',frq_fn,3))';
    if ~isfield(D.CRC,'pwrspect') && ~isempty(l_frq_fn)
        D.CRC.pwrspect = [];
    end
    if ~isempty(l_frq_fn)
        for ii=l_frq_fn
            tmp = getfield(D.CRC,frq_fn{ii});
            D.CRC.pwrspect.(frq_fn{ii}) = tmp;
            D.CRC = rmfield(D.CRC,frq_fn{ii});
        end
    end
    if isfield(D.CRC,'step')
        tmp = getfield(D.CRC,'step');
        D.CRC.pwrspect.('step') = tmp;
        D.CRC = rmfield(D.CRC,'step');
    end
    if isfield(D.CRC,'duration')
        tmp = getfield(D.CRC,'duration');
        D.CRC.pwrspect.('duration') = tmp;
        D.CRC = rmfield(D.CRC,'duration');
    end

    % Then update the .frq file array
    if isfield(D.CRC,'pwrspect') && ...
            isfield(D.CRC.pwrspect,'frqdata') && ...
            isa(D.CRC.pwrspect.frqdata,'file_array')
        % assume pwrspect data are next to the .mat file!
        fn_pwrs = D.CRC.pwrspect.frqdata.fname;
        if ~ispc, fn_pwrs = strrep(fn_pwrs,'\',filesep); end
        [~, pfname, pext] = fileparts(fn_pwrs);
        if isempty(pfname)
            if isempty(D.CRC.pwrspect.frqname)
                [~, pfname, ~] = fileparts(D.fname);
                pext = '.frq';
            else
                [~, pfname, pext] = fileparts(D.CRC.pwrspect.frqname);
            end
        end
        nfn_pwrs = fullfile(D.path,[pfname,pext]);
        if ~strcmpi(fn_pwrs,nfn_pwrs) % New location -> update fname
            if exist(nfn_pwrs,'file') && ~exist(nfn_pwrs,'dir')
                D.CRC.pwrspect.frqdata.fname = nfn_pwrs;
            else
                warning('Can''t find the powerspectrum file!');
                D.CRC.pwrspect.frqdata.fname = '';
            end
        end
        try
            % can access the pwrspect data ?
            D.CRC.pwrspect.frqdata(1,1);
        catch
            warning('problem accessing the powerspectrum file!');
            D.CRC.pwrspect.frqdata.fname = '';
        end
    end

    %% Deal with ROI data that used to be saved for SW detection
    % -> to be removed from the object
    if isfield(D.CRC,'SW') && ...
            isfield(D.CRC.SW,'DATA4ROI') && ...
            isfield(D.CRC.SW.DATA4ROI,'data')
        D4R = D.CRC.SW.DATA4ROI;
        D4R = rmfield(D4R,'data');
        D.CRC.SW.DATA4ROI = D4R;
    end

    %% Saving updated object!
    save(D);
end
