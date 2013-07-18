function [D, S] = crc_cont2epoch(S)
%
% Converts continous M/EEG-files files scored in FASST into epoched files
% for use with SPM -- Continous files are first filtered then downsampled
% before subsequent epoching
% FORMAT [D, S] = crc_cont2epoch(S)
%
% S                   - input structure (optional)
% (optional) fields of S:
%   S.D               - MEEG object or filename of M/EEG mat-file with
%                       epoched data
%   S.options with entries (all optional):
%     winSize         - Size of desired epochs (in sec), default is 4 sec
%     noverlap        - Overlap between epochs a number between 1 and 99 
%                       (in %) [50, default]
%     savedir         - path of directory where result files will be saved
%                       (default is current directory)
%     lower           - Frequency for highpass filtering [.5, default]
%     higher          - Frequency for lowpass filtering [25, default]
%     choice          - Choice of scorer (default is the last one)
%     preprocess      - Apply preprocessing steps (filtering, downsample,
%                       rereference) [1] or not [0, default]
%
% output:
% S         - can be used to construct script (as in the history-function)
%__________________________________________________________________________
%
% crc_cont2epoch converts continous M/EEG-files files scored in FAST
% into epoched files for use with SPM -- Continous files are first filtered
% then downsampled before subsequent epoching
%
% Output: The converted data are written to files. The header structs, but
% not the data, are returned in D as a cell vector of structs, and the
% struct S is returned to allow for saving the history of function calls.
%__________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by R. Lehembre & C. Phillips, 2011.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

%-Startup
%--------------------------------------------------------------------------

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = crc_eeg_load(D);

% Check options and set defaults

if ~isfield(S, 'winSize')
    S.winSize = 4;
end

if ~isfield(S, 'noverlap')
    S.noverlap = 50;
end
%convert to number between 0 and 1
S.noverlap = 1 - S.noverlap/100;

if ~isfield(S, 'savedir')
    S.savedir = '';
end

if ~isfield(S, 'lower')
    S.lower = 0.5;
end

if ~isfield(S, 'higher')
    S.higher = 25;
end

if ~isfield(S, 'preprocess')
    S.preprocess = false;
end


if ~isfield(S, 'choice')
    S.choice = size(D.CRC.score,2); %Take last scorer (merged?)
    %1 check if there is a score if yes multiple scorers? Choose one ..
    if ~isfield(D.CRC,'score')
        display('Data has not been scored yet !')
        return
    else
        Nscorers = size(D.CRC.score,2);
        if Nscorers>1
            display(['There are ' num2str(Nscorers) ' scorers'])
            display('choose scorer')
            scorerNames = {};
            for i=1:Nscorers
                scorerNames{i} = D.CRC.score{2,i};
            end
            S.choice = menu('choose scorer',scorerNames);
        else
            S.choice = 1;
        end
    end
end

S.choice

%Preprocessing: Filter, Downsample and Rereference
if S.preprocess == 1
    %High pass filter
    S.filter.band = 'high';
    S.filter.PHz = S.lower;
    S.D = spm_eeg_filter(S);
    
    %Low pass filter under 40 Hz
    S.filter.band = 'low';
    S.filter.PHz = S.higher;
    S.D = spm_eeg_filter(S);
    
    %Also downsample data to 200Hz
    S.fsample_new=200;
    S.D = spm_eeg_downsample(S);
    
    %And rereference to average reference for SPM source reconstructions
    S.refchan = 'average';
    S.D = spm_eeg_reref_eeg(S);
    
    D = S.D;
end

% Iterate over each scoring window and get non artefacted overlapping windows
nofSW = length(D.CRC.score{1,S.choice}); % Number of Scored Windows (nofSW)
sw = S.winSize; % size of windows (in sec)
lw = sw*fsample(D)-1; % size of windows in sample
so = S.winSize*S.noverlap; % size of overlap (in sec)
nofT = 1;

for i=1:nofSW   
    w2 = i*D.CRC.score{3,S.choice}*fsample(D); 
        % end (in samples) of i_th scoring window
    w1 = w2-D.CRC.score{3,S.choice}*fsample(D)+1; 
        % start (in samples) of i_th scoring window
    pos = w1;
    nx = w2;
    
    while (pos+lw <= nx)
        % check if window is inside an artefact
        A1 = find(pos>D.CRC.score{5,S.choice}(:,1)*fsample(D) & ...
                D.CRC.score{5,S.choice}(:,2)*fsample(D)>pos+lw);
        % check if any artefact starts inside window
        A2 = find(pos<D.CRC.score{5,S.choice}(:,1)*fsample(D) & ...
                D.CRC.score{5,S.choice}(:,1)*fsample(D)<pos+lw);
        % check if any artefact ends inside window
        A3 = find(pos<D.CRC.score{5,S.choice}(:,2)*fsample(D) & ...
                D.CRC.score{5,S.choice}(:,2)*fsample(D)<pos+lw);
        if isempty(A1) && isempty(A2) && isempty(A3)
            S.epochinfo.trl(nofT,:) = [pos pos+lw];
            S.epochinfo.conditionlabels{nofT} = D.CRC.score{1,S.choice}(i);
            nofT = nofT+1;
            pos = pos+so*fsample(D);
        else
            pos = pos+so*fsample(D);
        end
    end
    
end
D = spm_eeg_epochs(S);

return

