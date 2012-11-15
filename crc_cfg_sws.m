function sws = crc_cfg_sws
% SWS detection configuration file
% This builds the whole tree for the SWS detection tool.
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
data.help    = {'The data file with slow wave activity'};
data.filter  = 'mat';
data.ufilter = '.*';
data.num     = [1 1];

% ---------------------------------------------------------------------
% allf Use the whole file
% ---------------------------------------------------------------------
allf         = cfg_const;
allf.tag     = 'allf';
allf.name    = 'The whole file';
allf.val     = {true}; 

% ---------------------------------------------------------------------
% score Use score parts
% ---------------------------------------------------------------------
score        = cfg_entry;
score.tag    = 'score';
score.name   = 'Based on sleep stages';
score.help   = {'Select one or more sleep stage:'
                '1 -> NREM 1'
                '2 -> NREM 2'
                '3 -> NREM 3'
                '4 -> NREM 4'
                '5 -> REM'};
score.strtype= 'w';
score.num    = [Inf 1];
score.def    = @(val) crc_get_defaults('swsd.stagesw', val{:});
 
% ---------------------------------------------------------------------
% time Define time window
% ---------------------------------------------------------------------
time         = cfg_entry;
time.tag     = 'time';
time.name    = 'A specified time window';
time.help    = {'Specify time window in s'};
time.strtype = 'r';
time.num     = [2 1];

% ---------------------------------------------------------------------
% sel Seletected data to analyze
% ---------------------------------------------------------------------
sel          = cfg_choice;
sel.tag      = 'sel';
sel.name     = 'Select time window to analyze';
sel.values   = {allf time score};
sel.val      = {allf};

% ---------------------------------------------------------------------
% reref Rereferenced
% ---------------------------------------------------------------------
reref         = cfg_menu;
reref.tag     = 'reref';
reref.name    = 'Re-referenced?';
reref.help    = {'Is the dataset already re-referenced?'};
reref.labels  = {'Yes', 'No'};
reref.values  = {true, false};
reref.val     = {false};

% ---------------------------------------------------------------------
% hpfilt Low-pass filter settings
% ---------------------------------------------------------------------
hpfilt        = cfg_entry;
hpfilt.tag    = 'hpfilt';
hpfilt.name   = 'High-pass cutoff (Hz)';
hpfilt.help   = {'What are the high-pass filter settings?'};
hpfilt.num    = [1 1];
hpfilt.strtype= 'r';
hpfilt.def    = @(val) crc_get_defaults('swsd.highfc', val{:});

% ---------------------------------------------------------------------
% lpfilt Low-pass filter settings
% ---------------------------------------------------------------------
lpfilt        = cfg_entry;
lpfilt.tag    = 'lpfilt';
lpfilt.name   = 'Low-pass cutoff (Hz)';
lpfilt.help   = {'What are the low-pass filter settings?'};
lpfilt.num    = [1 1];
lpfilt.strtype= 'r';
lpfilt.def    = @(val) crc_get_defaults('swsd.lowfc', val{:});

% ---------------------------------------------------------------------
% filt Filter settings
% ---------------------------------------------------------------------
filt          = cfg_branch;
filt.tag      = 'filt';
filt.name     = 'Filter settings';
filt.help     = {'What are the high-pass and low-pass filter settings?'};
filt.val      = {hpfilt lpfilt};

% ---------------------------------------------------------------------
% auto Automatic selection of ROIs
% ---------------------------------------------------------------------
auto         = cfg_const;
auto.tag     = 'auto';
auto.name    = 'Automatic selection of ROIs';
auto.val     = {true};

% ---------------------------------------------------------------------
% roi Seletected ROI
% ---------------------------------------------------------------------
roi          = cfg_choice;
roi.tag      = 'roi';
roi.name     = 'Select electrode ROIs';
roi.help     = {'Automatic or manual ROI definition'
                'Manual ROI definition will be added later'};
roi.values   = {auto};
roi.val      = {auto};

% ---------------------------------------------------------------------
% sensauto Automatic sensor selectioon
% ---------------------------------------------------------------------
sensauto      = cfg_const;
sensauto.tag  = 'sensauto';
sensauto.name = 'Automatic sensors';
sensauto.val  = {true};

% ---------------------------------------------------------------------
% sensmanu Manual sensor selection
% ---------------------------------------------------------------------
sensmanu         = cfg_files; % maybe cfg_entry with strtype = 'e' is better
sensmanu.tag     = 'sensmanu';
sensmanu.name    = 'Sensor file';
sensmanu.help    = {'File with sensor definitions'};
sensmanu.filter  = 'mat';
sensmanu.ufilter = '.*';
sensmanu.num     = [1 1];

% ---------------------------------------------------------------------
% sens Auto or man sensor definition
% ---------------------------------------------------------------------
sens        = cfg_choice;
sens.tag    = 'sens';
sens.name   = 'Define sensors';
sens.values = {sensauto sensmanu};
sens.val    = {sensauto};

% ---------------------------------------------------------------------
% sens Auto or man sensor definition
% ---------------------------------------------------------------------
maps        = cfg_menu;
maps.tag    = 'maps';
maps.name   = 'Maps of...';
maps.help   = {'Maps of potentials have not been implemented'};
maps.labels = {'delays', 'potentials'};
maps.values = {1 2};
maps.val    = {1};

% ---------------------------------------------------------------------
% yesrev Review Slow Waves
% ---------------------------------------------------------------------
yesrev        = cfg_branch;
yesrev.tag    = 'yesrev';
yesrev.name   = 'Yes: review slow waves';
yesrev.val    = {sens maps};

% ---------------------------------------------------------------------
% norev Review Slow Waves
% ---------------------------------------------------------------------
norev         = cfg_const;
norev.tag     = 'norev';
norev.name    = 'Do not review slow waves';
norev.val     = {true};

% ---------------------------------------------------------------------
% review Review Slow Waves
% ---------------------------------------------------------------------
review        = cfg_choice;
review.tag    = 'review';
review.name   = 'Review slow waves';
review.values = {yesrev norev};
review.val    = {norev};

% ---------------------------------------------------------------------
% tr TR of fMRI
% ---------------------------------------------------------------------
tr            = cfg_entry;
tr.tag        = 'tr';
tr.name       = 'fMRI TR (in s)';
tr.num        = [1 1];
tr.strtype    = 'r';

% ---------------------------------------------------------------------
% mrk Marker of fMRI
% ---------------------------------------------------------------------
mrk           = cfg_entry;
mrk.tag       = 'mrk';
mrk.name      = 'Marker';
mrk.help      = {'Marker sent from the fMRI to the EEG to signal the start of a TR'};
mrk.num       = [1 1];
mrk.strtype   = 'r';

% ---------------------------------------------------------------------
% yfmri Yes Joint EEG-fMRI
% ---------------------------------------------------------------------
yfmri         = cfg_branch;
yfmri.tag     = 'yfmri';
yfmri.name    = 'Yes';
yfmri.val     = {tr mrk};

% ---------------------------------------------------------------------
% nfmri No Joint EEG-fMRI
% ---------------------------------------------------------------------
nfmri         = cfg_const;
nfmri.tag     = 'nfmri';
nfmri.name    = 'No';
nfmri.val     = {true};

% ---------------------------------------------------------------------
% fmri Joint EEG-fMRI
% ---------------------------------------------------------------------
fmri          = cfg_choice;
fmri.tag      = 'fmri';
fmri.name     = 'Joint EEG-fMRI';
fmri.values   = {yfmri nfmri};
fmri.val      = {nfmri};

% ---------------------------------------------------------------------
% sws Function to detect slow waves
% ---------------------------------------------------------------------
sws         = cfg_exbranch;
sws.tag     = 'sws';
sws.name    = 'Detect slow waves';
sws.val     = {data sel reref filt roi review fmri};
sws.help    = {'Function to detect slow waves.'};
sws.prog    = @crc_run_sws;
sws.vout    = @vout_sws;

%------------------------------------------------------------------------
% dependency
%------------------------------------------------------------------------
function dep = vout_sws(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'SWS-detected Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'SWS-detected Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
function out = crc_run_sws(job)

opt.fname = job.data{1};

if isfield(job.sel, 'time')
  opt.analyse = 1;
  opt.Begpts  = job.sel.time(1);
  opt.Endpts  = job.sel.time(2);
  
elseif isfield(job.sel, 'allf')
  opt.analyse = 2;
  
elseif isfield(job.sel, 'score')
  opt.analyse = 3;
  opt.stagesw = job.sel.score;
  
end

opt.reref  = job.reref;
opt.highfc = job.filt.hpfilt;
opt.lowfc  = job.filt.lpfilt;

opt.roisel = job.roi.auto;

if isfield(job.review, 'yesrev')
  opt.review = 1;
  
  if isfield(job.review.yesrev.sens, 'sensauto')
    opt.sensauto = 1;
  
  elseif isfield(job.review.yesrev.sens, 'sensmanu')
    opt.sensauto = 0;
    opt.sensfname = job.review.yesrev.sens.sensmanu{1};
  
  end
  
elseif isfield(job.review, 'norev')
  opt.review = 0;
  
end

if isfield(job.fmri, 'yfmri')
  opt.fmri = 1;
  opt.TR   = job.fmri.yfmri.tr;
  opt.marker = job.fmri.yfmri.mrk;
  
elseif isfield(job.fmri, 'nfmri')
  opt.fmri = 0;
  
end

D = crc_SWS_detect(opt);
out.D = meeg(D);

out.Dfname = {fullfile(out.D.path, out.D.fname)};














