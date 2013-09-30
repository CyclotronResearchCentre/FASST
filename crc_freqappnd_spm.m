function D = crc_freqappnd_spm(D,data)

% crc_freqappnd_spm(D,data)
%
% Append frequency files.
%
% Additional fields in the object
% D.CRC.pwrspect.frqname     frq filename
% D.CRC.pwrspect.frqNsamples Nsamples of the frqfiles
% D.CRC.pwrspect.frqNbins    Number of frequency bins
% D.CRC.pwrspect.frqdata     File array containing the frequency data
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Check sizes
[Nchan frqNbins frqNsamples] = size(data);
if D.CRC.pwrspect.frqNbins ~= frqNbins || size(D.CRC.pwrspect.frqdata,1) ~= Nchan
    error('Data passed to append have wrong size');
end

if isfield(D.CRC.pwrspect,'frqNsamples')
    D.CRC.pwrspect.frqNsamples = D.CRC.pwrspect.frqNsamples + frqNsamples ;
else
    D.CRC.pwrspect.frqNsamples = frqNsamples;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fpd_clean = fopen(fullfile(D.path, D.CRC.pwrspect.frqname), 'a'); % 'a' append

% append the data in file .frq
fseek(fpd_clean,0,'eof');
fwrite(fpd_clean, data, 'float32');
fclose(fpd_clean);
D.CRC.pwrspect.frqdata     = ...
    file_array(fullfile(D.path, D.CRC.pwrspect.frqname), ... % fname     - filename
               [D.nchannels D.CRC.pwrspect.frqNbins D.CRC.pwrspect.frqNsamples],...  % dim       - dimensions (default = [0 0] )
               spm_type('float32'), ...             % dtype     - datatype   (default = 'uint8-le')
               0, ...                               % offset    - offset into file (default = 0)
               D.CRC.pwrspect.frqscale);                     % scl_slope - scalefactor (default = 1)

if nargout<1
    save(D);
end

return