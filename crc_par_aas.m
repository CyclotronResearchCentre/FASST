function Do = crc_par_aas(Di)

% Pulse artefact rejection with a AAS method, based on the fmrib toolbox.
%_______________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Written by C. Phillips, 2011.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

crcdef = crc_get_defaults('par');
prefAAS = crcdef.aaspref;
fn_dat = fullfile(spm_str_manip(fnamedat(Di),'H'), ...
    [ prefAAS, spm_str_manip(fnamedat(Di),'t')]);

% Recover Peaks
Peaks = Di.CRC.EKGPeaks;

if exist(fullfile(path(Di), fn_dat),'file')
    % This check is because I want append data in the new file
    % I have to make sure that it does not exist yet.
    delete(fullfile(path(Di), fn_dat))
    disp(['!!! WARNING : EXISTING ',fn_dat,' FILE OVERWRITTEN !!!' ])
end

% create empty object
Do = clone(Di,fn_dat);

Mem_limit = crc_get_defaults('mem.maxmemload'); % or about for 5e7~50Mb for a small system
tmpL_limit = floor(Mem_limit/8/Di.nchannels);% Maximum length of tw to use

Nchunks = ceil(nsamples(Di)/tmpL_limit);% Number of chunks to use
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

% Loop over chunks
l_eeg = Di.cache.l2cor;
l_other = setdiff(1:Di.nchannels,l_eeg);
for ii=1:Nchunks
    chunk_ii = t_chunk(:,ii);
    Di.cache.chunk_ii = chunk_ii;
    tmpPeaks = intersect(Peaks,t_chunk(1,ii):t_chunk(2,ii));
    tmpPeaks = tmpPeaks-t_chunk(1,ii)+1;

    [fitted_art] = fmrib_pas(Di,tmpPeaks,'gmean',4);
    Do(l_eeg,t_chunk(1,ii):t_chunk(2,ii)) = ...
        Di(l_eeg,t_chunk(1,ii):t_chunk(2,ii)) - fitted_art;
    Do(l_other,t_chunk(1,ii):t_chunk(2,ii)) = ...
        Di(l_other,t_chunk(1,ii):t_chunk(2,ii));
end

return
