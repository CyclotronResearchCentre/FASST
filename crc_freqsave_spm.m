function D = crc_freqsave_spm(files,D,data,reference, nscor)

% FORMAT crc_freqsave_spm(files,D,data)
% Save frequency files.
%
% Variable "files" need to be the .frq file.
%
% Additional array in the D object as compared to the "classic" spm object
% D.CRC.pwrspect.frqname       frq filename
% D.CRC.pwrspect.frqNsamples   Nsamples of the frqfiles
% D.CRC.pwrspect.frqNbins      Number of frequency bins
% D.CRC.pwrspect.frqdata       File array containing the frequency data
% D.CRC.pwrspect.frq_reference Re-reference used, if any.
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$
[dumb1 frqNbins frqNsamples] = size(data);
[pth name]=fileparts(files(1,:));
namesscorer=D.CRC.score{2,nscor};
name = strcat(name,'_',namesscorer);
% Adding additional "CRC" array to spm

pathname = D.path;
filename = D.fname;
D.CRC.pwrspect.frqname = strcat(name,'.frq');%[char(n),'.frq'])
D.CRC.pwrspect.frqNsamples = frqNsamples;
D.CRC.pwrspect.frqNbins = frqNbins;
if reference >1
    D.CRC.pwrspect.frq_reference = chanlabels(D,reference-1); % added by PM for AF and LM
elseif reference==1
    D.CRC.pwrspect.frq_reference = 'no re-reference';
else
    D.CRC.pwrspect.frq_reference = reference;
end

%%%%%%%%%%%%%%%%%%% Check file .frq already exists %%%%%%%%%%%%%%%%%%%%
%if exist(fullfile(pathname,filename),'file') 
%  delete(fullfile(D.path,D.CRC.pwrspect.frqname))
% disp(['!!! WARNING : EXISTING ' D.CRC.pwrspect.frqname ' FILE OVERWRITTEN !!!' ])
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fpd_clean = fopen(fullfile(pathname, filename), 'w'); % 'w' write
%write the data in file .frq
fwrite(fpd_clean, data, 'float32');
fclose(fpd_clean);
D.CRC.pwrspect.frqscale = ones(size(data, 1), 1);
D.CRC.pwrspect.frqdata = ...
    file_array(fullfile(pathname, filename), ...    % fname     - filename
               [D.nchannels D.CRC.pwrspect.frqNbins D.CRC.pwrspect.frqNsamples],...  % dim       - dimensions (default = [0 0] )
                spm_type('float32'), ...               % dtype     - datatype   (default = 'uint8-le')
                0, ...                                  % offset    - offset into file (default = 0)
                D.CRC.pwrspect.frqscale);                               % scl_slope - scalefactor (default = 1)
if nargout<1
    save(D);
end

return