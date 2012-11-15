function new_D = crc_save_spm(prefix,D,data)

% FORMAT crc_save_spm(prefix,D,data)
% Save EEG data in SPM8 meeg object format
%
% INPUT:
%   prefix : prefix of filename to be created
%   D      : object from which data come from
%   data   : data to be saved
%
% OUTPUT:
%   new_D  : new meeg object created from D and data.
%            + file written on disk with specified prefix
%
% The routine use the structure format of the meeg object to access
% directly the data part and adjust various fields.
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

new_Ds = struct(D);
[kk1 name kk2] = fileparts(D.fname);
new_Ds.data.fnamedat = [prefix name '.dat'];

% remove scale field
if isfield(new_Ds.data,'scale')
    new_Ds.data = rmfield(new_Ds.data,'scale');
end

if exist(fullfile(new_Ds.path,new_Ds.data.fnamedat),'file')
    spm_unlink(fullfile(new_Ds.path,new_Ds.data.fnamedat))
    disp(['!!! WARNING : EXISTING ' new_Ds.data.fnamedat ' FILE OVERWRITTEN !!!' ])
end

% Create new data file
fpd_new = fopen(fullfile(new_Ds.path,new_Ds.data.fnamedat), 'w'); % 'w' write

%write the data in the new .dat file
fwrite(fpd_new, data, 'float32');
fclose(fpd_new);

new_Ds.fname = [prefix name '.mat'];
new_Ds.Nsamples = size(data,2);
new_Ds.data.y = file_array(...
    fullfile(new_Ds.path, new_Ds.data.fnamedat), ...       % fname     - filename
    [length(new_Ds.channels) new_Ds.Nsamples size(new_Ds.trials,1)],...  % dim       - dimensions (default = [0 0] )
    spm_type('float32'), ...        % dtype     - datatype   (default = 'float')
    0, ...                          % offset    - offset into file (default = 0)
    ones(size(data, 1), 1));        % scl_slope - scalefactor (default = 1)
new_Ds.data.datatype = 'float32';
new_D = meeg(new_Ds);
save(new_D)

