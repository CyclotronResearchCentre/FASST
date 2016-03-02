function crc_defaults
% Sets the defaults which are used by FASST, fMRI Artefact removal and 
% Sleep Scoring Toolbox.
%
% FORMAT crc_defaults
%_______________________________________________________________________
%
% This file can be customised to any the site/person own setup.
% Individual users can make copies which can be stored on their own
% matlab path. Make sure your 'crc_defaults' is the first one found in the
% path. See matlab documentation for details on setting path.
%
% Care must be taken when modifying this file!
%
% The structure and content of this file are largely inspired by SPM.
%_______________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$


%% NEED to include the scaling of ECG/EMG/EOG

%%
global crc_def

%% %%%%%%%%%%%%%%%%% Artefact Rejection Parameters %%%%%%%%%%%%%%%%%%%%%%%

% Parameters for the gradient artefact rejection, crc_gar.m
%-----------------------------------------------
crc_def.gar.prefix         = 'cga_';
crc_def.gar.output_fs      = 500 ;     % Frequency sampling of the output file
crc_def.gar.Nsc_aver       = 30 ;      % Average computed on ... scans
crc_def.gar.Nsc_skipped_be = 3 ;       % Number of scans skipped before correction
crc_def.gar.UseScanMrk     = 1 ;       % Using volume marker from the scanner (1), or not (0)

% If UseScanMrk = 1;
crc_def.gar.ScanMrk1       = 128 ;     % Marker number coming from the scanner
crc_def.gar.ScanMrk2       = 1002 ;    % Secondary scanner marker.
crc_def.gar.MrkStep        = 1 ;      % Number of marker between 2 successive volumes
% If you have a scanner marker every slice,
% use the number of slices !

% If UseScanMrk = 0;
crc_def.gar.TR             = 2.46; % in sec
crc_def.gar.AutoChk        = 1;    % Autodetection of beginning & end
% If crc_def.gar.AutoChk = 1
crc_def.gar.Threshold      = 350;
crc_def.gar.DetStep        = 1 ;   % Step use for detection (in sec)
crc_def.gar.DetChan        = 1 ;   % Channel use for detection
% If crc_def.gar.AutoChk = 0
crc_def.gar.beg            = 1 ;       % in sec
crc_def.gar.nd             = 2500;     % in sec

% Parameters for cICA pulse artefact calculation,
%-----------------------------------------------
% used in crc_bcgremoval_rt.m

% Filename convention.
crc_def.par.cicapref  = 'cICA_'; %filename prefix added after using cICA
crc_def.par.pcapref   = 'cpca_'; %filename prefix added after using PCA
crc_def.par.aaspref   = 'caas_'; %filename prefix added after using AAS
crc_def.par.iicapref  = 'ciica_'; %filename prefix added after using iICA
crc_def.par.pcaaspref = 'cpcaas_'; %filename prefix added after using PCA/AAS

% Parameters for pulse artefact rejection
%----------------------------------------
crc_def.par.bcgrem.size        = 90 ; % in sec (60 is the minimum)
crc_def.par.bcgrem.step        = 60 ; % about 2/3 size
crc_def.par.bcgrem.scSNR       = 3 ; % nr of signal std to consider signal as noise: abs(signal)>scSNR*max(std)
% Channels used to build the constrain vector.
% They're taken from around the head
crc_def.par.Refnames = ...
    {'AF3' 'AF4' 'FPZ' 'FP1' 'FP2' 'O1' 'O2' 'F1' 'F2' 'C1' 'C2' 'PZ' 'OZ' 'F5' ...
    'F6' 'C5' 'C6' 'P5' 'P6'};
crc_def.par.NitKmeans = 150;

% Parameters for pulse artefact matrix quality assessment
%--------------------------------------------------------
crc_def.par.useinitseg = 1; % Use initial segment to assess the quality of the correction matrix.
crc_def.par.additioseg = 10; % Number of additional random segments.
% Note that if both of these value equals 0, the software use the
% 1st segment anyway.
crc_def.par.length = 60; % Length of each segment in seconds.

% Parameters for PAR Batch (giopia)
%-------------------------
crc_def.par.qrsmethod  = 1;  % QRS detection method, only one so far...
crc_def.par.bcgmethod  = 'acica';  % pulse artefact rejection method
crc_def.par.peaksave   = 1;  % save detected QRS peaks (1) or (0)
crc_def.par.fqrsdet    = 0;  % force QRS detection (1) or not (0), even if peaks already available
crc_def.par.nit        = 50; % #iterations for the "iterative ICA" method
crc_def.par.ecgchan    = 0;  % default ECG channel = first one with name ECG/EKG
crc_def.par.badchan    = []; % e.g. 32 for BrainProducts caps with 64ch setup as at the CRC

%% %%%%%%%%%%%%%%%%% Display parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for display one file,
%-----------------------------------------------
% used in crc_dis_main.m
crc_def.one.winsize    = 20; % in sec
crc_def.one.scale      = 150; % in µV
crc_def.one.filtEEG    = [.1 20]; % in Hz
crc_def.one.filtEMG    = [10 125]; % max is adapted to true sampling rate
crc_def.one.filtEOG    = [.1 4]; % in Hz

% Parameters for scale handling,
%-----------------------------------------------
% used in crc_dis_main.m, crc_SP_detect.m, crc_SW_detect
crc_def.scales.eeg   = 10^-6; % in V
crc_def.scales.eog   = 10^-6; % in V
crc_def.scales.ecg   = 3*10^-5; % in V
crc_def.scales.emg   = 10^-6; % in V
crc_def.scales.other   = 1; % in µV
crc_def.scales.lfp   = 10^-5; % in V
crc_def.scales.megmag   = 5*10^-14; % in µV
crc_def.scales.megplan   = 5*10^-13; % in µV

% Parameters for display power spectrum
%-----------------------------------------------
crc_def.disfrq.scale       = [200 0.3 5]; % Absolute Value / Relative Value / Mongrain Value
crc_def.disfrq.subsmpl     = 0; %subsampling factor computed automatically. The user can set a value instead (enter value)
crc_def.disfrq.maxpix      = 1600;
crc_def.disfrq.window      = 0;

% Parameters for scoring file
%-----------------------------------------------
crc_def.score.scale      = 150; % in µV
crc_def.score.filtEEG    = [.1 125]; % in Hz
crc_def.score.filtEMG    = [10 125]; % max is adapted to true sampling rate
crc_def.score.filtEOG    = [.1 5]; % in Hz
crc_def.score.winsize    = 30; % in sec
crc_def.score.lab_art    = {'blink','body movement','eye movement'};
crc_def.score.nrStage    = 8; % from 0 to 7, no more than 8 though !
crc_def.score.stnames_L          = {...
    'Subject/Patient is awake', ...
    'Sleep Stage 1', ...
    'Sleep Stage 2', ...
    'Sleep Stage 3', ...
    'Sleep Stage 4', ...
    'REM Sleep',...
    'Mouvement time', ...
    'Unscorable'};
crc_def.score.stnames_S          = {...
    'Awake', ...
    'Stage 1', ...
    'Stage 2', ...
    'Stage 3', ...
    'Stage 4', ...
    'REM',...
    'MT', ...
    'Unscor'};
crc_def.score.stnames_sS          = {...
    'Awake', ...
    'St1', ...
    'St2', ...
    'St3', ...
    'St4', ...
    'REM',...
    'MT', ...
    'Unsc'};
crc_def.score.color = [ ...
    [0.2 0.75 0.6] ; ...
    [0 0.8 1] ; ...
    [0.1 0.5 0.9] ; ...
    [0.1 0.2 0.8] ; ...
    [0.1 0.15 0.5] ; ...
    [0.5 0.5 0.9] ; ...
    [0.9 0.4 0.4] ; ...
    [0.9 0.6 0.3] ];
crc_def.score.height = [7 6 5 4 3 2 1 .3];
% NOTE: the hypnogram plotting is hardcoded with 8 levels!!!
crc_def.score.check_disnum = 10; % #empty windows printed out at prompt
crc_def.score.check_jump   = true; % jumping (or not) to next empty window

%% %%%%%%%%%%%%%%%%% Processing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for computing spectrogram or power spectrum, crc_spectcompute.m
%-----------------------------------------------
crc_def.cps.fmax      = 25 ;  % upper frequency \_ bound (Hz)
crc_def.cps.fmin      = .5 ;  % lower frequency /
crc_def.cps.dur       = 4 ;   % Window size in sec
crc_def.cps.step      = 2 ;   % sliding step
crc_def.cps.scorer    = 1 ;   % index of scorer
crc_def.cps.ref       = 1 ;   % index of reference

% Parameters for downsampling
%-----------------------------------------------
crc_def.ds.stepse  = 60; % sec
crc_def.ds.fs      = 250; % Hz
crc_def.ds.prefix  = 'ds_';

% Parameters for chunking
%-----------------------------------------------
crc_def.chk.overwr   = 0; % No overwriting by def
crc_def.chk.numchunk = 1; % Start index at 1
crc_def.chk.prefix   = 'chk';
crc_def.chk.clockt   = 0; % Don't worry about clocktime

% Parameters for concatenating
%-----------------------------------------------
crc_def.cat.prefix1   = 'cc_'; % standard concatenation
crc_def.cat.prefix2   = 'cr_'; % concatenation at reading (raw EGI)


% Parameters for Slow Wave detection
%-----------------------------------------------
crc_def.swsd.dispscale     = 200; % microV
crc_def.swsd.bpfc          = 200; % lowpass cutoff of bandpass filter
crc_def.swsd.ROI_centers   = ... % position of ROI centers on the 2D unit disk
    [0.5    0.6341;  %Fz
    0.3536 0.5013;  %C3
    0.6464 0.5013;  %C4
    0.5    0.3575]; %Pz
crc_def.swsd.SWlength      = [250 1250 1500 200];
% minimum and maximum time durations between down zero crossing and up zero crossing
crc_def.swsd.SWmAmpl       = [-40 -80 75 140];
% Massimini criteria of magnitude for SWS (-80 140) and delta waves (-40 75)
crc_def.swsd.highfc        = 0.2;  % highpass frequency cutoff
crc_def.swsd.lowfc         = 4;    % lowpass frequency cutoff
crc_def.swsd.stagesw       = [3 4];% stages to extract from original scored file
crc_def.swsd.butterorder   = 4;    % order of butterworth filter
crc_def.swsd.prcentile     = 90;    % order of butterworth filter
crc_def.swsd.type          = 'EEG';    % type of channels to detect SWs
crc_def.swsd.usetheor        = 0;% do not use theoretical positions if localisation file available (1 for always use theoretical positions).

% Parameters for memory handling
%-----------------------------------------------
crc_def.mem.cps_maxmemload = 50  * 1024^2; % i.e. 50Mb
crc_def.mem.maxmemload     = 100 * 1024^2; % i.e. 100Mb
crc_def.mem.sp_maxmemload  = 200 * 1024^2; % i.e. 200Mb

% Parameters for Spindle detection
%--------------------------------------------
crc_def.sp.highfc          = 8;
crc_def.sp.lowfc           = 20;
crc_def.sp.elecoi          = {'Fz' 'Cz' 'Pz'};
crc_def.sp.stagesp         = [2 3 4];% stages to extract from original scored file
crc_def.sp.prcthresh       = 95; %percentile used to perform detection of spindles
crc_def.sp.stagethresh     = 2;% stage to compute threshold of detection
crc_def.sp.threshold       = 0;% user chosing the threshold of detection
crc_def.sp.lengthsp        = 0.4;% spindles of 400 ms duration minimum (in s)
crc_def.sp.succsp          = 1;% 1000ms between 2 successive sp (in s)
crc_def.sp.type            = 'EEG';% type of channels to detect spindles
crc_def.sp.usetheor        = 0;% do not use theoretical positions if localisation file available (1 for always use theoretical positions).

return

