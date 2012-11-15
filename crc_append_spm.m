function new_D = crc_append_spm(D,data)
% FORMAT new_D = crc_append_spm(D,data)
%
% Append the data in matrix "data" into the *dat file of the D object
% and updates the object accordingly.
%
% INPUT
%   D      : meeg object to modify
%   data   : data to append
%
% OUTPUT
%   new_D  : meeg object created from D with appended data written on disk.
%
% We are assuming that the data are in float32 format !
% Note also that this routine is using the structure form of the meeg
% object to modify the data part (and some other fields).
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

new_Ds = struct(D);

fpd_clean = fopen(fullfile(new_Ds.path, new_Ds.data.fnamedat), 'a'); % 'a' append

% Write the data in file .dat
fseek(fpd_clean,0,'eof');
fwrite(fpd_clean, data, 'float32');   
fclose(fpd_clean);

new_Ds.Nsamples = new_Ds.Nsamples+size(data,2);

new_Ds.data.y = file_array( ...
    fullfile(new_Ds.path, new_Ds.data.fnamedat), ... % fname - filename
    [length(new_Ds.channels) new_Ds.Nsamples size(new_Ds.trials,1)],...  % dim - dimensions (default = [0 0] )
    spm_type('float32'), ... % dtype - datatype   (default = 'float32')
    0, ... % offset - offset into file (default = 0)
    ones(size(data, 1), 1)); % scl_slope - scalefactor (default = 1)
new_Ds.data.datatype = 'float32';

% save(fullfile(pth,D.fname),'D')
new_D = meeg(new_Ds);
save(new_D)
