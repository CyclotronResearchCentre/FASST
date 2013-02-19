function flags = crc_check_flag(flags_o,flags)

% FORMAT flags = crc_check_flag(flags_o,flags)
%
% Function to automatically check the content of a "flag" structure, using
% a "default flag structure", adding the missing fields and putting in the 
% default value if none was provided.
%
% INPUT:
% flags_o   default or reference structure
% flags     input flag/option structure that need to be filled for missing
%           fields with default values
%
% OUPUT:
% flags     filled flag/option structure
%__________________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

f_names = fieldnames(flags_o);
% list fields in default structure

Nfields = length(f_names);
for ii=1:Nfields
    if ~isfield(flags,f_names{ii}) || isempty(getfield(flags,f_names{ii}))
        flags = setfield(flags,f_names{ii},getfield(flags_o,f_names{ii}));
    end
end


