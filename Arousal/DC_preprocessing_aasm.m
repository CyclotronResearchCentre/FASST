function Dnew = DC_preprocessing_aasm(D)

% FORMAT DC_preprocessing(args)
% Preprocessing for whole sleep recordings: 
%       * Filtering :
%                  - EEG : 0.5 - 30Hz
%                  - EOG : 0.1 - 5Hz
%                  - EMG : 10 - 125Hz
%         
% The program proposes to clone your data to avoid losing initial data.
% the new file preprocessed starts by a 'P' letter to mention it has been preprocessed.
%
% INPUT
%       .file   - data file (.mat files)
%__________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% load parameters values
fs      =   fsample(D);
if fs/2>2000
    forder = 1;
else
    forder = 3;
end
% clone data 
name_data   =   fnamedat(D);
[pth dataname ext] = fileparts(name_data);
newname    =   ['P' dataname];  
newfilename  =   fullfile(path(D),[newname '.mat']);
Dnew = clone(D,newfilename);  

% dataFT = spm2fieldtrip(data);

% filter EEG channels (0.5-30Hz) relatively to the AASM manual (Iber et
% al., 2007)
% cfg_FT.channel = meegchannels(data);
% cfg_FT.lpfreq = 30;
% cfg_FT.hpfreq = 0.5;
% cfg_FT.lpfilter = 'yes';
% cfg_FT.hpfilter = 'yes';
% cfg.lpfiltord   = forder;    % lowpass  filter order 
% cfg.hpfiltord   = forder;    % highpass filter order 
% 
% prep_data = FT_preprocessing(cfg_FT,dataFT);
% Dnew(cfg_FT.channel,:,1) = prep_data.trial{1};
% fprintf(1,'EEG data filtered between 0.5-30Hz \n');

% filter EMG channels (10-125Hz) relatively to the AASM manual (Iber et
% al., 2007)
% cfg_FT.channel = emgchannels(data);
% cfg_FT.lpfreq = 125;
% cfg_FT.hpfreq = 10;
% 
% prep_data = FT_preprocessing(cfg_FT,dataFT);
% Dnew(cfg_FT.channel,:,1) = prep_data.trial{1};
% fprintf(1,'EMG data filtered between 10-125Hz \n');
% 
%  ************************************************************************
%                       Filtering
%  ************************************************************************
% 
[B_eeg_low,A_eeg_low] = butter(3,30/(fs/2),'low');
[B_eog_low,A_eog_low] = butter(3,5/(fs/2),'low');
[B_emg_low,A_emg_low] = butter(3,min(125,fs/2-1)/(fs/2),'low');

[B_eeg_high,A_eeg_high] = butter(3,0.5/(fs/2),'high');
[B_eog_high,A_eog_high] = butter(3,0.1/(fs/2),'high');
[B_emg_high,A_emg_high] = butter(3,10/(fs/2),'high');

meeg    =   meegchannels(D);
eog     =   eogchannels(D);
emg     =   emgchannels(D);
% channels = 1 : size(D,1);
fprintf(['MEG and EEG channels are filtered ... \n'])
% if isempty(meeg)
%     meeg = 1:9;
% end
for ieeg = 1 : numel(meeg)
    ieeg
    Dnew(meeg(ieeg),:) = filtfilt(B_eeg_low,A_eeg_low,D(meeg(ieeg),:));
    Dnew(meeg(ieeg),:) = filtfilt(B_eeg_high,A_eeg_high,Dnew(meeg(ieeg),:));
end
fprintf(['EOG channels are filtered ... \n'])
if ~isempty(eog)
    Dnew(eog,:) = filtfilt(B_eog_low,A_eog_low,D(eog,:)')';
    Dnew(eog,:) = filtfilt(B_eog_high,A_eog_high,Dnew(eog,:)')';
end
fprintf(['EMG channels are filtered ... \n'])
if ~isempty(emg)
    Dnew(emg,:) = filtfilt(B_emg_low,A_emg_low,D(emg,:)')';
    Dnew(emg,:) = filtfilt(B_emg_high,A_emg_high,Dnew(emg,:)')';
end
% for e=1:length(channels)
%      
%      if any(meeg == e)
%         C(e,:) = filtfilt(B_eeg_low,A_eeg_low,D(e,:)')';
%         C(e,:) = filtfilt(B_eeg_high,A_eeg_high,C(e,:)')';
%      elseif any(eog == e)
%          C(e,:) = filtfilt(B_eog_low,A_eog_low,D(e,:)')';
%          C(e,:) = filtfilt(B_eog_high,A_eog_high,C(e,:)')';
%      elseif any(emg == e)
%          C(e,:) = filtfilt(B_emg_low,A_emg_low,D(e,:)')';
%          C(e,:) = filtfilt(B_emg_high,A_emg_high,C(e,:)')';
%      else 
%      	C(e,:) = D(e,:);    
%      end
%     fprintf([' ' num2str(e) ' channels copied and filtered on ' num2str(length(channels)) '\n'])
% end

% ************************************************************************
%                       Mean correction
% ************************************************************************

% h = waitbar(0,'Progression for mean correction...');
% for w = 1 : NofW
%     win    =   C(meeg,(w-1)*nbr_sc_epoch + 1 : min(w*nbr_sc_epoch,nspl)); 
%     C(meeg,(w-1)*nbr_sc_epoch + 1 : min(w*nbr_sc_epoch,nspl))   =   win - mean(win,2)*ones(1,size(win,2));
%     String  =  ['Progress for mean correction... : ' num2str(w/NofW*100) ' %'];
%     waitbar((w/NofW),h,String);
% end
% close(h);
% sprintf('* Mean correction per scoring window of %dsec duration \n', [winsize])
% C.CRC.DC.correction = '[0.5Hz - 30Hz]';
save(Dnew);

