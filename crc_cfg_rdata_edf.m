function readedf = crc_cfg_rdata_edf
% Read EDF data configuration file
% This builds the whole tree for the crc_eeg_rdata_edf tool.
%_______________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Giovanni Piantoni & Christophe Phillips
% $Id$

% ---------------------------------------------------------------------
% data Data file
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'M/EEG data file';
data.help    = {'EDF dataset, .edf file, to read in.'};
data.filter  = 'edf';
data.ufilter = '.edf';
data.num     = [1 1];
% ---------------------------------------------------------------------
% readedf Function to correct the pulse artefact
% ---------------------------------------------------------------------
readedf         = cfg_exbranch;
readedf.tag     = 'readedf';
readedf.name    = 'Read EDF data';
readedf.val     = {data};
readedf.help    = {'Converts EEG data from EDF format to SPM8 format.'};
readedf.prog    = @crc_run_edf;
readedf.vout    = @vout_edf;

%% SUBFUNCTIONS

%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
function out = crc_run_edf(job)

[out.D] = crc_eeg_rdata_edf(job.data{1});
out.Dfname = {fullfile(out.D.path, out.D.fname)};
%------------------------------------------------------------------------
% dependency
%------------------------------------------------------------------------
function dep = vout_edf(job)
% Output is always in field "D", no matter how job is structured

% reference field "D" from output
% this can be entered into any evaluated input
dep(1)            = cfg_dep;
dep(1).sname      = 'EDF Data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

% reference field "Dfname" from output
% this can be entered into any file selector
dep(2)            = cfg_dep;
dep(2).sname      = 'EDF Datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
