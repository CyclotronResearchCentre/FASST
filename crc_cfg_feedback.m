function fbck = crc_cfg_feedback
% Pulse artifact (BCG) configuration file
% This builds the whole tree for the feedback tool.
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Giovanni Piantoni
% $Id$

% ---------------------------------------------------------------------
% data Data file
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'M/EEG data file';
data.help    = {'This is the data file with either spindles or slow wave'};
data.filter  = 'mat';
data.ufilter = '.*';
data.num     = [1 1];

% ---------------------------------------------------------------------
% wavtyp Which type of wave to plot
% ---------------------------------------------------------------------
wavtyp         = cfg_menu;
wavtyp.tag     = 'wavtyp';
wavtyp.name    = 'Wave type';
wavtyp.help    = {'Which wave type to use as events'};
wavtyp.labels  = {'Slow Waves'
                  'Spindles'}';
wavtyp.values  = {1 2};

% ---------------------------------------------------------------------
% int Interactive plotting (will have more options in the future)
% ---------------------------------------------------------------------
int           = cfg_const;
int.tag       = 'int';
int.name      = 'Interactive plotting';
int.val       = {true}; 

% ---------------------------------------------------------------------
% comm Comment
% ---------------------------------------------------------------------
comm         = cfg_entry;
comm.tag     = 'comm';
comm.name    = 'comment';
comm.help    = {'Append this comment at the end of the file'
                'Can be empty'
                'Default name of the file is ''fbck_'' filename wavetype'};
comm.strtype = 's';
comm.num     = [0 Inf];
comm.val     = {''};

% ---------------------------------------------------------------------
% pri List of option for pulseartfct
% ---------------------------------------------------------------------
pri         = cfg_branch;
pri.tag     = 'pri';
pri.name    = 'Print to file';
pri.val     = {comm};

% ---------------------------------------------------------------------
% met Methods for plotting
% ---------------------------------------------------------------------
met         = cfg_choice;
met.tag     = 'met';
met.name    = 'Plot method';
met.help    = {'Two plotting methods. ''interactive'' plots on the fly,'
  'and you need to click to scroll through the detection'
  '''print'' prints a png, with fbck.col number of columns (you can change the value in the defaults file)'};
met.values  = {int pri};
met.val     = {int}; %@(val)crc_get_defaults('fbck.met', val{:});

% ---------------------------------------------------------------------
% pad Padding before and after event, in s
% ---------------------------------------------------------------------
pad         = cfg_entry;
pad.tag     = 'pad';
pad.name    = 'Padding before and after event';
pad.help    = {'Padding before and after event (in s)'};
pad.strtype = 'e';
pad.num     = [1 1];
pad.def     = @(val)crc_get_defaults('fbck.pad', val{:});

% ---------------------------------------------------------------------
% scl Scaling for y-axis
% ---------------------------------------------------------------------
scl         = cfg_entry;
scl.tag     = 'scl';
scl.name    = 'Scaling for y-axis';
scl.help    = {'Scaling for y-axis (in uV)'};
scl.strtype = 'r';
scl.num     = [1 1];
scl.def     = @(val)crc_get_defaults('fbck.scl', val{:});

% ---------------------------------------------------------------------
% options List of option for pulseartfct
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Feedback options';
options.val     = {wavtyp met pad scl};
options.help    = {'Options for feedback'};

% ---------------------------------------------------------------------
% fbck Give fbck on events detected 
% ---------------------------------------------------------------------
fbck         = cfg_exbranch;
fbck.tag     = 'fbck';
fbck.name    = 'Detection Feedback';
fbck.val     = {data options};
fbck.help    = {'It plots the events (slow waves or spindles) on screen or disk.'};
fbck.prog    = @crc_run_fbck;

function out = crc_run_fbck(job)
% convert batch to opt for crc_run_feedback

opt.pad = job.options.pad;
opt.scl = job.options.scl;

if job.options.wavtyp == 1
  opt.wavtyp = 'slowwaves';
else
  opt.wavtyp = 'spindles';
end

if isfield(job.options.met, 'int')
  opt.met = 'interactive';
elseif isfield(job.options.met, 'pri')
  opt.met = 'print';
  opt.comm = job.options.met.pri.comm;
end


crc_run_feedback(job.data{1}, opt)

