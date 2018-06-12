function DC_artifact_COF
%%% It takes about 13 min/subject
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
% Pathes to CODEcodedir = 'D:\DATA_COF\SCRIPTS\Arousal';
pathdir = 'D:\DATA_COF\SCRIPTS\Arousal';
addpath(pathdir);
origdir = pwd;

% For Cofitage data
% pathtodata= 'D:\GWA';
pathtodata= 'D:\DATA_COF\DATA\BL_Analyse_Berthomier';
data = COF_data(pathtodata);
EEGdir = 'BL';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL for output naming:
nom_start = 1;  % which is the first character for your output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Actual analysis
% spm eeg
% fasst

for isub = 1;%36:47%size(data,2)% Scoring on recording of subject isub=6 is not avalaible
    num = str2num(data(isub).id(5:6));
    cd(strcat(data(isub).dir,'\',EEGdir))
 
    files = dir('*.mat');
   for xx = 1:size(files,1)
                if strfind(files(xx).name,'A4R')
                    SCOFILE = files(xx).name;
                end
   end
   EEGFILES = spm_select('FPList',pwd,'MCOF.*.mat');
        
    D = spm_eeg_load(EEGFILES(1,:));
    AS = load(SCOFILE);
    try D.CRC
    catch D.CRC = [];
        D.CRC.goodevents = [];
        D.CRC.score = [];
    end
%     A = events(D);
%     for zz = 1:size(A,2)
%        if strcmp(A(zz).type,'LIGHTS-OFF') == 1;
%            debut = A(zz).time;
%            break
%        end
%     end
    D.CRC.score{1} = AS.aseega_analysis.scoring.hypno5;
    %%%%% REMETTRE BONS CODES POUR STADES
    D.CRC.score{1}(D.CRC.score{1}==10) = 0;
    D.CRC.score{1}(D.CRC.score{1}==1) = 3;
    D.CRC.score{1}(D.CRC.score{1}==5) = 2;
    D.CRC.score{1}(D.CRC.score{1}==7) = 5;
    D.CRC.score{1}(D.CRC.score{1}==8) = 1;
    D.CRC.score{1}(D.CRC.score{1}==11) = 7;

    %%%%%%%%%%%%%%
    D.CRC.score{2} = 'Aseega';
    D.CRC.score{3} = 30;
    D.CRC.score{4} = [];
    D.CRC.score{5} = [];
    D.CRC.score{6} = [];
    D.CRC.score{7} = [];
    D.CRC.score{8} = [];
    D.CRC.score = D.CRC.score';
    
%     D = chantype(D,14,'EEG');
%     D = chantype(D,15,'EMG');
%     D = chantype(D,16:17,'EOG');
    
    save(D)
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
    cfg.winsize =   30;   % time window for detection; default: 30

    fprintf(1,'===========================================\n');
    fprintf(1,'TREATING SUBJECT %s \n',fname(D));
    fprintf(1,'===========================================\n');
    
    D       =   DC_preprocessing_aasm(D);     % preprocessing: filtering, cloning  %% +- 5min
    Dbadc   =   DC_badchannels(D, cfg);       % detection of bad channels  %% +- 1min
    Dfeat   =   DC_pwr(Dbadc,cfg);            % power spectrum computation   %% +- 5min
    Dpop    =   DC_popping(Dfeat);            % Identification of popping artefacts   %% 10sec
%     DREM    =   FR_REMclass(Dpop,pwd);  % REM classifier
    indREM = find(AS.aseega_analysis.scoring.hypno3 == 7);
    Dpop.CRC.DC.REM = indREM;
    Dbade   =   DC_bademg(Dpop);                 % EMG/EEG shift evaluation   %% +- 1min
    % indSL = find(AS.aseega_analysis.scoring.hypno2 == 1);
    % Only computing arousal from N2,N3 & REM
    indSL = find(AS.aseega_analysis.scoring.hypno5 < 8);
    indSLincl1 = find(AS.aseega_analysis.scoring.hypno5 < 9);
    Dbade.CRC.DC.SL = indSL;
    Dbade.CRC.DC.SL_incl_1 = indSLincl1;
    Darou   =   DC_arousal(Dbade);    % Main Arousal Output as .csv in indicated output path   %% few secs
    Dprep   =   DC_extra(Darou);              % UNION of all the artifacts: D.CRC.DC.shortartf.artefact
%     clear cfg
end
rmpath(pathdir);