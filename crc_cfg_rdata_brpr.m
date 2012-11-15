function readbrpr = crc_cfg_rdata_brpr
% Read Brain Product data configuration file
% This builds the whole tree for the crc_eeg_rdata_brpr tool.
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Giovanni Piantoni & Christophe Phillips
% $Id$

% ---------------------------------------------------------------------
% data Data file
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'M/EEG data file';
data.help    = {'Brain Products dataset, .vhdr header file, to read in.'};
data.filter  = 'vhdr';
data.ufilter = '.vhdr';
data.num     = [1 1];
% ---------------------------------------------------------------------
% readbrpr Function to correct the pulse artefact
% ---------------------------------------------------------------------
readbrpr         = cfg_exbranch;
readbrpr.tag     = 'readbrpr';
readbrpr.name    = 'Read Brain Products data';
readbrpr.val     = {data};
readbrpr.help    = {'Converts EEG data from BrainProduct to SPM8 format.'};
readbrpr.prog    = @crc_run_brpr;
readbrpr.vout    = @vout_brpr;

%% SUBFUNCTIONS

%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
function out = crc_run_brpr(job)

[out.D] = crc_eeg_rdata_brpr(job.data{1});
out.Dfname = {fullfile(out.D.path, out.D.fname)};
%------------------------------------------------------------------------
% dependency
%------------------------------------------------------------------------
function dep = vout_brpr(job)
% Output is always in field "D", no matter how job is structured

% reference field "D" from output
% this can be entered into any evaluated input
dep(1)            = cfg_dep;
dep(1).sname      = 'BrPr Data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

% reference field "Dfname" from output
% this can be entered into any file selector
dep(2)            = cfg_dep;
dep(2).sname      = 'BrPr Datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
