function varargout = crc_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defaults = crc_get_defaults
% Return the global "defaults" variable defined in crc_defaults.m.
%
% FORMAT defval = crc_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "crc_def" variable defined in crc_defaults.m.
%
% FORMAT crc_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of FASST. To make
% persistent changes, edit crc_defaults.m.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Copyright (C) 2010 Cyclotron Research Centre

% Volkmar Glauche
% Then modified for use with the FASST toolbox by Christophe Phillips
% $Id$

global crc_def;
if isempty(crc_def)
    crc_defaults;
end

if nargin == 0
    varargout{1} = crc_def;
    return
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(crc_def, subs);
else
    crc_def = subsasgn(crc_def, subs, varargin{1});
end
