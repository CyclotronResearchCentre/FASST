% Automatic Artifact and Arousal Detection for whole sleep recordings.
% The data for which the proposed was developed are assumed to come from a standard sleep study: 
%       * multichannel EEG recording
%       * the electrodes should be placed according to the 10-20 system 
%       * the detection will be finer it it is included:
%             - least one electromyogram, EMG, derivation; 
%             - least one active (i.e. not used as a reference) mastoid, MAS, channel; 
%
% The program proposes to clone your data to avoid losing initial data.
%
% INPUT
%       .file   - data file (.mat files)
%__________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014
% and adapted by F. Rudzik, 2017
% Cyclotron Research Centre, University of Liege, Belgium

% Set pathes
% Pathes to CODE
% addpath 'C:\Users\Franziska\Documents\FRANZISKA\SiRENE\FHR\Cooperation\Liege_Toolbox\ToolboxES\AASM_toolbox\ASDM_Doro\'

% INPUT
pathin  = 'D:\DATA_COF\DATA\';
% OUTPUT: Arousal
outpath  = 'D:\DATA_COF\DATA\';
% OUTPUT: Spectral Analysis
out_spect = 'D:\DATA_COF\DATA\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL for output naming:
nom_start = 1;  % which is the first character for your output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set relevant PARAMETER
cfg.winsize =   30;   % time window for detection; default: 30

%% Actual analysis
dirio = dir(fullfile(pathin, '*.mat'));

for f = 1 %: size(dirio,1)
%     EEGFile = fullfile(pathin, dirio(f).name);
%     D = spm_eeg_load(deblank(EEGFile));
    D = spm_eeg_load;
    % check channels
    chantype(D);
    channels	=   meegchannels(D,'EEG');
    if isempty(channels)
        error(1,'Error ----- No M/EEG channels found ------ \n')
    end
    % % check mastoids channels
    if any(or(or(strncmp(chanlabels(D),'M',1), strncmp(chanlabels(D),'A1',2)), strncmp(chanlabels(D),'A2',2)))
        fprintf(1,'Mastoid channels are not considered in the bad channel detection process \n')
        channels(or(strncmp(chanlabels(D),'M',1), strncmp(chanlabels(D),'A',1))) = [];
    end
    cfg.channels	= channels;

    fprintf(1,'===========================================\n');
    fprintf(1,'TREATING SUBJECT %s \n',fname(D));
    fprintf(1,'===========================================\n');
    
    D       =   DC_preprocessing_aasm(D);     % preprocessing: filtering, cloning  
    Dbadc   =   DC_badchannels(D, cfg);       % detection of bad channels
    Dfeat   =   DC_pwr(Dbadc,cfg);            % power spectrum computation
    Dpop    =   DC_popping(Dfeat);            % Identification of popping artefacts 
    DREM    =   FR_REMclass(Dpop,out_spect);  % REM classifier
    Dbade   =   DC_bademg(D);                 % EMG/EEG shift evaluation
    Darou   =   DC_arousal(Dbade,outpath);    % Main Arousal Output as .csv in indicated output path
    Dprep   =   DC_extra(Darou);              % UNION of all the artifacts: D.CRC.DC.shortartf.artefact
end
