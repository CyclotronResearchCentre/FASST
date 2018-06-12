function D = FR_REMclass(D,out_spect)
% FORMAT FR_REMclass(args)
% It classifies sleep according to NREM/REM parts.
% Based on Doro's 'crc_NREMclassifier' she wrote for the spindle derection
% based on based on a Mixture of Gaussian models.
% REM: All-NREM and aggregated by a duration criterium.
%
% INPUT
%       .file   - data file (.mat files)
%__________________________________________________________________

% Written by F. Rudzik, 2017 based on D. Coppieters 't Wallant, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% from D
fs = fsample(D);
nspl = nsamples(D);
winsize = 30;
NofW = floor(nspl/(fs*winsize));

% Choose the right channels: here: Fz, Cz, and Pz 
chan_x = chanlabels(D);
chan1 = find(ismember(chan_x,'Fz'));
chan2 = find(ismember(chan_x,'Cz'));
chan3 = find(ismember(chan_x,'Pz'));
channels = [chan1 chan2 chan3];
 
ne = numel(channels);
uni = units(D,channels);
scal = zeros(1,size(D,1));
for iu = 1 : length(channels)
    if strcmp(uni(iu),'T/m')
        scal(channels(iu)) = 5e-13;
    elseif strcmp(uni(iu),'T')
        scal(channels(iu)) = 5e-14;
    elseif strcmpi(uni(iu),'V')
        scal(channels(iu)) = 1e-6;
    elseif strcmpi(uni(iu),'UV')
        scal(channels(iu)) = 1;
    else 
        scal(channels(iu)) = 1;
    end
end 

% NREM PS computation for each channel selected at a time 
D = crc_NREMclassifier(D,channels,out_spect); 

% Aggregate all the REM parts to form the REM bouts: all within a REM bot
% (be it 0 or NREM1, NREM2 is then also considered REM)
% Use Electrode Cz for classification, if not available other z-electrodes
if ~isempty(D.CRC.DC.classREM.Cz)
    NREM_est = D.CRC.DC.classREM.Cz;
else
    if ~isempty(D.CRC.DC.classREM.Fz)
        NREM_est = D.CRC.DC.classREM.Fz;
    else
        if ~isempty(D.CRC.DC.classREM.Pz)
            NREM_est = D.CRC.DC.classREM.Pz;
        else
        fprintf('Classification difficult as no z-electrodes in the data set!')
        end
    end
end

% REM epochs
idx = find(diff(NREM_est)>15);
REM = [];
for i = 1:length(idx)
    REM = [REM NREM_est(idx(i))+1:NREM_est(idx(i)+1)-1];
end

% Single epochs before and after are also included
REMAdd = [];
for i = 1:length(idx)
    % BEFORE
    iy = NREM_est(idx(i))+1-20; % start of last REM within the next twenty epochs
    iyy = intersect(iy:iy+20,NREM_est);
    iyyy = find(diff(iyy)>6);  % equals 5 min of REM
    tA = length(iyyy);
    if tA > 0
        REMAdd = [REMAdd iyy(iyyy)+1:iyy(iyyy+1)-1];
    end
    % AFTER
    ix = NREM_est(idx(i)+1)-1+20; % end of last REM within the next twenty epochs
    ixx = intersect(ix-20-1:ix,NREM_est);
    ixxx = find(diff(ixx)>6); % equals 5 min of REM
    tB = length(ixxx);
    if tB > 0
        REMAdd = [REMAdd ixx(ixxx)+1:ixx(ixxx+1)-1];
    end
end
REMAdd = unique(REMAdd);
REM = unique([REM REMAdd]);
D.CRC.DC.REM = REM;
save(D);
end