function [Do] = crc_par(files,flags)

%%
% FORMAT [Do] = crc_par(files,flags)
%
% Function to correct the pulse artefact, also named balisto-cardiogram 
% artefact or cardio-balistic artefact.
%
% INPUT
% - files : (list of) file name(s) of imported data to correct. Keep in
%           mind that if multiple files are selected the same options wil
%           be applied to all of them!
% - flags :  options for the correction
%   .qrsmethod : detection of QRS, (1, default) fmrib solution (2) nothing
%                here yet (could be another method, or simply selecting a
%                'peak' file)
%   .bcgmethod : bcg artifact correction,
%                (1) fmrib PCA
%                (2) fmrib (gaussian mean) AAS
%                (3) multiple application of ICA approach (very slow)
%                       and not offered through the GUI anymore
%                (4) automatic constrained ICA correction method
%                (5) compute cICA correction matrix only
%                (6) combination of fmrib (gaussian mean) AAS & PCA
%   .peaksave  : after detecting the QRS peak save it into the data
%                structure (1, default) or not (0)
%   .ecgchan   : index of ECG channel to use. Set it to 0 (default) to let
%                the routine automatically select the *1st* ECG available.
%   .badchan   : manually pre-specified list of bad channels to skip for 
%                correction (none by default). Channels marked as bad will
%                also be skipped!
%   [.nit]     : For the mICA * only*. Number of repetition of ICA on each 
%                data stretch, by default set at 50
%
% OUTPUT
% - Do : Data object in cell array with all the corrected data
% - written on disk, each data file corrected. Names are prepended with
%   "cpca_", "caas_" or "cICA_" depending on the correction method used.
%_______________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% - Selection menu does not appear in an open figure anymore
% - cICA implemented for BCG artefact rejection/can be used for long files
%__________________________________________________________________

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

crcdef = crc_get_defaults('par');

if ~nargin
    files = spm_select(Inf,'mat', 'Select imported EEG file','' ,pwd);
end

if isobject(files)
    D = files;
    files = fullfile(D.path, D.fname);
end

badchan= cell(0);

flags_def = struct( ...
    'qrsmethod',crcdef.qrsmethod, ...
    'bcgmethod',crcdef.bcgmethod, ...
    'peaksave',crcdef.peaksave, ...
    'nit',crcdef.nit, ...
    'ecgchan',crcdef.ecgchan, ...
    'badchan',crcdef.badchan, ...
    'fqrsdet',crcdef.fqrsdet);
if nargin<2, flags = []; end
flags = crc_check_flag(flags_def,flags);
peaksave = flags.peaksave;
force_qrs_detect = flags.fqrsdet;

if nargin<2 % used through the "old" GUI

    for ifile = 1:size(files,1)
        figz = crc_par_selchan({deblank(files(ifile,:))});
        waitfor(figz,'UserData',1);
        gdat = guidata(figz);
        badchan{ifile}=gdat.output;
        close(figz);
    end

    % Open a new window in order to select the method for the BCG artefact
    % removal and fill in the various processing options
    idxfig=figure;
    set(idxfig,'Position',[10 130 320 150]);
    set(idxfig,'NumberTitle','off');
    set(idxfig,'Name','CRC Pulse Artefact Remover');
    %     qrsmethod =  spm_input('QRS detection',1,'m','fmrib|own');
    qrsmethod =  1; % so far no alternative!
    str_method = ['fmrib (pca)|fmrib (gaussian mean)|Constrained ICA',...
        '(automatic)|Constrained ICA (manual)|"AAS & PCA" combined ',...
        '(experimental method)'];
    str_sel = strvcat('pca','aas','acica','mcica','aaspca');
    %     str_method = ['fmrib (pca)|fmrib (gaussian mean)|iterative ICA|',...
    %             'constrained ICA|PCA & AAS combined (experimental method)|',...
    %             'Computation of cICA correction matrix'];
    %     str_sel = strvcat('pca','aas','iica','acica','mcica','aaspca');
    bcgmethod =  spm_input('BCG rejection','+1','m',str_method,str_sel);
    % accounts for the removed iterative/multiple ICA option!
    if strcmp(bcgmethod,'iica')
        nit   = spm_input('# of repetitive ICA on data stretch?', ...
            '+1','i',50);
    end
    ecgchan = spm_input('Specific EKG channel, or leave 0', ...
        '+1','i',0);
    close(idxfig);
else % used in a script or through batch system, only 1 file at a time!
    qrsmethod = flags.qrsmethod;
    bcgmethod = flags.bcgmethod;
    if strcmp(bcgmethod,'iica')
        nit = flags.nit;
    end
    ecgchan = flags.ecgchan;
    % find badchannels
    badchan{1} = flags.badchan;
end

for ifile = 1:size(files,1)
    Di = crc_eeg_load(deblank(files(ifile,:)));

    disp( '.......................')
    disp(['.... Processing file : ',Di.fname])
    disp( '.......................')

    % make sure that Peaks are saved in the data structure
    % and test if it already exists...
    if ~isfield(Di,'CRC') | ...
            (isfield(Di,'CRC') && ~isfield(Di.CRC,'EKGPeaks')) | ...
            force_qrs_detect
        %%  % -------------- QRS DETECTION -----------------------
        disp('.... Detecting QRS ....')
        disp('.......................')
        switch qrsmethod
            case 1 %'fmrib'
                if ecgchan
                    if length(ecgchan)>1
                        Peaks = fmrib_qrsdetect(Di,ecgchan(1),0);
                        disp(['.... WARNING: Using user specified channel # ',...
                            num2str(ecgchan(1)),...
                            ' for QRS detection'])

                    else
                        Peaks = fmrib_qrsdetect(Di,ecgchan,0);
                        disp(['.... WARNING: Using user specified channel # ',...
                            num2str(ecgchan),...
                            ' for QRS detection'])
                    end

                else
                    ecgchan = ecgchannels(Di);
                    if isempty(ecgchan)
                        error('No ECG channel available for QRS detection!');
                    end
                    Peaks = fmrib_qrsdetect(Di,ecgchan(1),0);
                end
                if length(ecgchan)>1
                    disp(['.... WARNING: Only used channel # ',...
                        num2str(ecgchan(1)),...
                        ' for QRS detection'])
                end
            case 2 % Peaks detected by other means...
                error('Not implemented yet...')
                % Should allow the selection of a file mat/txt with peaks
                % information !!!
        end
        Di.CRC.EKGPeaks = Peaks;
        if peaksave
            save(Di);
        end
    else
        disp('.... Using estimated QRS ....')
        disp('.............................')
        
        Peaks = Di.CRC.EKGPeaks;
    end

    %%    % ---------------- BCG ARTEFACT REMOVAL --------------------------
    disp('................................')
    disp('.... BCG artefact rejection ....')
    disp('................................')

    % Combine manually and file specified bad channel list
    badch = union(badchan{ifile},badchannels(Di));
    % Keep all EEG channels that have not been specifically marked as bad
    l2cor = intersect(setdiff(1:Di.nchannels,badch),meegchannels(Di,'EEG'));
    badch = setdiff(1:Di.nchannels,l2cor);
    Di.cache.l2cor = l2cor; % storing list of channels to correct.
    Di.cache.badch = badch; % storing bad channel info.
    switch bcgmethod

            %% fmrib PCA
        case 'pca' % fmrib pca
            Dclean = crc_par_pca(Di);
            % save version number of routine
            [v,r] = crc_fasst_utils('Ver','crc_par_pca');
            Dclean.info.ver_parpca = struct('v_nr',v,'rel',r);       
            Dclean.cache = [];
            save(Dclean);

            %% fmrib gaussian mean AAS
        case 'aas' % fmrib (gaussian mean AAS)
            Dclean = crc_par_aas(Di);
            % save version number of routine
            [v,r] = crc_fasst_utils('Ver','crc_par_aas');
            Dclean.info.ver_paraas = struct('v_nr',v,'rel',r);        
            Dclean.cache = [];
            save(Dclean);

            %% iICA
        case 'iica' % iterative ICA
            Di.cache.nit = nit; % storing number of iterations.
            Dclean = crc_par_iica(Di);
            % save version number of routine
            [v,r] = crc_fasst_utils('Ver','crc_par_iica');
            Dclean.info.ver_pariica = struct('v_nr',v,'rel',r);        
            Dclean.cache = [];
            save(Dclean);

            %% cICA with automatic selection of the correction matrix
        case 'acica' % constrained ICA based correction method, the best
            % correction matrix is picked using an automatic MI
            % criteria. It relies on crc_par_cICAmx to find these
            % correction matrices and crc_par_cICAqa for the MI
            % selection
            Dclean = crc_par_cICA(Di);
            % save version number of routine
            [v,r] = crc_fasst_utils('Ver','crc_par_cICA');
            Dclean.info.ver_paracica = struct('v_nr',v,'rel',r);        
            Dclean.cache = [];
            save(Dclean);

            %% cICA with manual selection of the correction matrix
        case 'mcica' % constrained ICA based correction method, 
            % but data are NOT corrected! Only the correction matrices are,
            % then up to the user to manually pick the "best" one via the 
            % display GUI.
            Dclean = crc_par_cICAmx(Di);
            % save version number of routine
            [v,r] = crc_fasst_utils('Ver','crc_par_cICAmx');
            Dclean.info.ver_parmcica = struct('v_nr',v,'rel',r);        
            Dclean.cache = [];
            save(Dclean);

            %% fmrib PCA & AAS combined
        case 'aaspca' % combination of PCA & (Gaussian mean) AAS
            Dclean = crc_par_aaspca(Di);
            % save version number of routine
            [v,r] = crc_fasst_utils('Ver','crc_par_aaspca');
            Dclean.info.ver_paraaspca = struct('v_nr',v,'rel',r);        
            Dclean.cache = [];
            save(Dclean);

        otherwise
            disp('..................................................')
            disp('.... Unknown correction method, doing nothing ....')
            disp('..................................................')
            Dclean = Di;

    end % end switch
    Do{ifile} = Dclean; %#ok<AGROW> % Keep all the results

end % ifile

% If I wanted to output object, instead of file array
if ifile==1;
    clear Do
    Do = Dclean;
end

disp('.... Done ....')

return