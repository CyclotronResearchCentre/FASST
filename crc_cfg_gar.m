function gar = crc_cfg_gar
% Removes the gradient artefact induced by fMRI sacnning.
% This builds the whole tree for the crc_gar tool.
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Giovanni Piantoni & Christophe Phillips
% $Id$

% ---------------------------------------------------------------------
% data Data file
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'EEG data file';
data.help    = {'Select the EEG dataset to correct'};
data.filter  = 'mat';
data.ufilter = '.mat';
data.num     = [1 1];
% ---------------------------------------------------------------------
% gar Function to correct the pulse artefact
% ---------------------------------------------------------------------
gar         = cfg_exbranch;
gar.tag     = 'gar';
gar.name    = 'Gradient artifact rejection';
gar.val     = {data};
gar.help    = {'Function to correct the gradient artifact from EEG data acquired during fMRI scanning.' ...
               'Select only one file at a time. To correct multiple files in a row, use the classic GUI or your own script.'};
gar.prog    = @crc_run_gar;
gar.vout    = @vout_gar;

%% SUBFUNCTIONS
%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
function out = crc_run_gar(job)

[out.Dfname] = {crc_gar(job.data{1})};
out.D = crc_eeg_load(out.Dfname{1});
%------------------------------------------------------------------------
% dependency
%------------------------------------------------------------------------
function dep = vout_gar(job)
% Output is always in field "D", no matter how job is structured

% reference field "D" from output
% this can be entered into any evaluated input
dep(1)            = cfg_dep;
dep(1).sname      = 'GA corrected data';
dep(1).src_output = substruct('.','D');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

% reference field "Dfname" from output
% this can be entered into any file selector
dep(2) = cfg_dep;
dep(2).sname      = 'GA corrected datafile';
dep(2).src_output = substruct('.','Dfname');
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
