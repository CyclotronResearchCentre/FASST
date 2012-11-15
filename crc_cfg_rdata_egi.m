function readegi = crc_cfg_rdata_egi
% Read raw EGI data configuration file
% This builds the whole tree for the crc_eeg_rdata_egi tool.
%_______________________________________________________________________
% Copyright (C) 2011 Cyclotron Research Centre

% Christophe Phillips
% $Id$

% ---------------------------------------------------------------------
% data Data file
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'M/EEG data file';
data.help    = {'raw EGI dataset, .raw file, to read in.' ...,
                'If multiple files are selected then they will be' ...
                '*CONCATENATED* when reading in the data!'};
data.filter  = 'raw';
data.ufilter = '.raw';
data.num     = [1 Inf];
% ---------------------------------------------------------------------
% readegi Function to correct the pulse artefact
% ---------------------------------------------------------------------
readegi         = cfg_exbranch;
readegi.tag     = 'readegi';
readegi.name    = 'Read raw EGI data';
readegi.val     = {data};
readegi.help    = {'Converts EEG data from raw EGI format to SPM8 format.'};
readegi.prog    = @crc_run_egi;
readegi.vout    = @vout_egi;

%% SUBFUNCTIONS

%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
function out = crc_run_egi(job)

[out.D] = crc_eeg_rdata_egi(job.data{1});
out.Dfname = {fullfile(out.D.path, out.D.fname)};
%------------------------------------------------------------------------
% dependency
%------------------------------------------------------------------------
function dep = vout_egi(job)
% Output is always in field "D", no matter how job is structured

% reference field "D" from output
% this can be entered into any evaluated input
dep(1)            = cfg_dep;
dep(1).sname      = 'EGI Data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

% reference field "Dfname" from output
% this can be entered into any file selector
dep(2)            = cfg_dep;
dep(2).sname      = 'EGI Datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
