function concatdata = crc_cfg_concatenate
% Concatenate a pair of data files
% This builds the whole tree for the crc_concatenate tool.
%_______________________________________________________________________
% Copyright (C) 2012 Cyclotron Research Centre

% Christophe Phillips
% $Id$

% ---------------------------------------------------------------------
% data Data file
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'M/EEG data file';
data.help    = {'Select the 2 data files to concatenate.'};
data.filter  = 'mat';
data.ufilter = '.mat';
data.num     = [2 2];
% ---------------------------------------------------------------------
% readegi Function to correct the pulse artefact
% ---------------------------------------------------------------------
concatdata         = cfg_exbranch;
concatdata.tag     = 'concatdata';
concatdata.name    = 'Concatenate continuous data files';
concatdata.val     = {data};
concatdata.help    = {'Concatenate a pair of data files, accounting '...
                        'for their acquisition time a duration, '...
                        ' if that information is available.'};
concatdata.prog    = @crc_run_concat;
concatdata.vout    = @vout_concat;

%% SUBFUNCTIONS

%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
function out = crc_run_concat(job)

[out.D] = crc_concatenate(job.data{1});
out.Dfname = {fullfile(out.D.path, out.D.fname)};
%------------------------------------------------------------------------
% dependency
%------------------------------------------------------------------------
function dep = vout_concat(job)
% Output is always in field "D", no matter how job is structured

% reference field "D" from output
% this can be entered into any evaluated input
dep(1)            = cfg_dep;
dep(1).sname      = 'Data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

% reference field "Dfname" from output
% this can be entered into any file selector
dep(2)            = cfg_dep;
dep(2).sname      = 'Data filename';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
