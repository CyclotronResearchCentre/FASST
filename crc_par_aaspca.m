function Do = crc_par_aaspca(Di)

% Pulse artefact rejection with combined PCA & (Gaussian mean) AAS, based
% on the fmrib toolbox.
%_______________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by C. Phillips, 2011.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Initializaing the cleaned EEG datafile
crcdef = crc_get_defaults('par');
prefPCAAS = crcdef.pcaaspref;
% fn_dat = [ prefPCAAS, fnamedat(Di)];
fn_dat = fullfile(spm_str_manip(fnamedat(Di),'H'), ...
    [ prefAAS, spm_str_manip(fnamedat(Di),'t')]);

% Recover Peaks
Peaks = Di.CRC.EKGPeaks;

if exist(fullfile(Di.path, fn_dat),'file')
    % This loop is because I wanna append data in the new file
    % I have to make sure that it does not exist yet.
    eval(['delete ' fullfile(Di.path, fn_dat)  ])
    disp(['!!! WARNING : EXISTING ' fn_dat ' FILE OVERWRITTEN !!!' ])
end

% create empty object
Do = clone(Di,fn_dat);

Mem_limit = crc_get_defaults('mem.maxmemload'); % or about for 5e7~50Mb for a small system
tmpL_limit = floor(Mem_limit/8/Di.nchannels);% Maximum length of tw to use

Nchunks = ceil(Di.nsamples/tmpL_limit); % Number of chunks to use
t_chunk = zeros(2,Nchunks);
t_chunk(1,:) = (0:Nchunks-1)*tmpL_limit+1;
t_chunk(2,:) = t_chunk(1,:)+tmpL_limit-1;
t_chunk(2,end) = nsamples(Di);
% checking that last chunk is long enough!
if Nchunks>1
    lChunk = diff(t_chunk(:,end));
    if lChunk<tmpL_limit*.2 % less than 20% of chunk length
        t_chunk(2,end-1) = t_chunk(2,end);
        t_chunk(:,end) = [];
        Nchunks = Nchunks-1;
    end
end

fs = Di.fsample;
swfrq = 4; % Cut off frequency, between AAS and PCA correction
Lfrq = .05; % lowest frequency to use
Hfrq = fs/2*.8; % limitting to .8 of Nyquist freq

% Loop over chunks
l_eeg = Di.cache.l2cor;
l_other = setdiff(1:Di.nchannels,l_eeg);
for ii=1:Nchunks
    chunk_ii = t_chunk(:,ii);
    Di.cache.chunk_ii = chunk_ii;
    tmpPeaks = intersect(Peaks,t_chunk(1,ii):t_chunk(2,ii));
    tmpPeaks = tmpPeaks-t_chunk(1,ii)+1;
    tmpPeaks = tmpPeaks(1:end-1);
    tmpPeaks = tmpPeaks-min(tmpPeaks)+1;
    
    [fitted_aas] = fmrib_pas(Di,tmpPeaks,'gmean');
    [fitted_pca] = fmrib_pas(Di,tmpPeaks,'obs',4);

    aas = filterlowhigh( ...
                (Di(l_eeg,t_chunk(1,ii):t_chunk(2,ii))-fitted_aas)', ...
                 fs,[Lfrq swfrq]);
    pca = filterlowhigh(...
                (Di(l_eeg,t_chunk(1,ii):t_chunk(2,ii))-fitted_pca)', ...
                 fs,[swfrq  Hfrq]);

    Do(l_eeg,t_chunk(1,ii):t_chunk(2,ii)) = ...
                            filterlowhigh((aas + pca),fs, [Lfrq Hfrq])';
    Do(l_other,t_chunk(1,ii):t_chunk(2,ii)) = ...
                            Do(l_other,t_chunk(1,ii):t_chunk(2,ii));
    
end

%% Other functions
function Y = filterlowhigh(X,fs,frqcut)

[B,A] = butter(3, [frqcut(1)/(fs/2), frqcut(2)/(fs/2)], 'pass');

% Apply Butterworth filter
Y = filtfilt(B,A,X);

return

%% OLD filter version
% function Y = filterlowhigh(X,fs,frqcut)
% 
% % Plt are the data
% 
% flc = frqcut(1)/fs;
% fhc = frqcut(2)/fs;
% 
% k = .7; % cut-off value
% 
% alphal = (1-k*cos(2*pi*flc)-sqrt(2*k*(1-cos(2*pi*flc))-k^2*sin(2*pi*flc)^2))/(1-k);
% alphah = (1-k*cos(2*pi*fhc)-sqrt(2*k*(1-cos(2*pi*fhc))-k^2*sin(2*pi*fhc)^2))/(1-k);
% 
% 
% % Apply low pass filter
% Y = filtfilt(1-alphal,[1 -alphal],X);
% 
% % Apply high pass filter
% Y = filtfilt(1-alphah,[1 -alphah],X-Y);
% 
% return