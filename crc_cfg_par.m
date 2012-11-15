function pulseartfct = crc_cfg_par
% Pulse artifact (BCG) configuration file
% This builds the whole tree for the pulse artifact tool.
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Written by Giovanni Piantoni & Christophe Phillips
% $Id$

% ---------------------------------------------------------------------
% data Data file
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'M/EEG data file';
data.help    = {'This is the data file to correct BCG on'};
data.filter  = 'mat';
data.ufilter = '.*';
data.num     = [1 1];

% ---------------------------------------------------------------------
% qrsmethod detection of QRS
% ---------------------------------------------------------------------
qrsmethod         = cfg_menu;
qrsmethod.tag     = 'qrsmethod';
qrsmethod.name    = 'Detection of QRS';
qrsmethod.help    = {'which method to detect QRS'};
qrsmethod.labels  = {'fmrib'
                  'other'}'; % nothing here yet (could be another method, or simply selecting a 'peak' file)
qrsmethod.values  = {1 2};
qrsmethod.def     = @(val)crc_get_defaults('par.qrsmethod', val{:});

% ---------------------------------------------------------------------
% bcgmethod bcg artifact correction
% ---------------------------------------------------------------------
bcgmethod         = cfg_menu;
bcgmethod.tag     = 'bcgmethod';
bcgmethod.name    = 'Correction of BCG';
bcgmethod.help    = {'which method to correct BCG'};
bcgmethod.labels  = {'fmrib (PCA)'
                     'fmrib (gaussian mean)' 
                     'constrained ICA (automatic)'
                     'constrained ICA (manual)'
                     'AAS & PCA combined method (experimental)'};
% bcgmethod.values  = {1 2 4 5 6}; % 3 is a deprecated option
bcgmethod.values  = {'pca' 'aas' 'acica' 'mcica' 'aaspca'};
                        % 'iica' is a deprecated option
bcgmethod.def     = @(val)crc_get_defaults('par.bcgmethod', val{:});

% ---------------------------------------------------------------------
% peaksave Save the peaks
% ---------------------------------------------------------------------
peaksave         = cfg_menu;
peaksave.tag     = 'peaksave';
peaksave.name    = 'Save QRS peaks';
peaksave.help    = {'After detecting the QRS peak save it into the data structure, or not.'};
peaksave.labels  = {'yes' 'no'}';
peaksave.values  = {1 2};
peaksave.def     = @(val)crc_get_defaults('par.peaksave', val{:});

% % ---------------------------------------------------------------------
% % nit Number of repetitions of ICA on each data stretch
% % ---------------------------------------------------------------------
% nit         = cfg_entry;
% nit.tag     = 'nit';
% nit.name    = 'number of ICA repetions';
% nit.help    = {'Number of repetitions of ICA on each data stretch'};
% nit.strtype = 'e';
% nit.num     = [1 1];
% nit.def     = @(val)crc_get_defaults('par.nit', val{:});

% ---------------------------------------------------------------------
% ecgchan Index of ECG channel to use
% ---------------------------------------------------------------------
ecgchan         = cfg_entry;
ecgchan.tag     = 'ecgchan';
ecgchan.name    = 'ECG channel';
ecgchan.help    = {'Which channel is the ECG (0 for default)'};
ecgchan.strtype = 'e';
ecgchan.num     = [1 1];
ecgchan.def     = @(val)crc_get_defaults('par.ecgchan', val{:});

% ---------------------------------------------------------------------
% fqrsdet Force detection of QRS, or not
% ---------------------------------------------------------------------
fqrsdet         = cfg_menu;
fqrsdet.tag     = 'fqrsdet';
fqrsdet.name    = 'Force QRS detection';
fqrsdet.help    = {'Force the detection of QRS events, or not. '...
                    'If forcing, then previsouly detected QRS are '...
                    'simply skipped and "new" peaks are detected.'};
fqrsdet.labels  = {'no' 'yes'}';
fqrsdet.values  = {0 1};
fqrsdet.def     = @(val)crc_get_defaults('par.fqrsdet', val{:});

% ---------------------------------------------------------------------
% badchan Index of bad channels to skip
% ---------------------------------------------------------------------
badchan         = cfg_entry;
badchan.tag     = 'badchan';
badchan.name    = 'Bad channels';
badchan.help    = {'Which channels are ''bad'' (EKG for sure)'};
badchan.strtype = 'e';
badchan.num     = [0 Inf];
badchan.def     = @(val)crc_get_defaults('par.badchan', val{:});

% ---------------------------------------------------------------------
% options List of option for pulseartfct
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'BCG option';
% options.val     = {qrsmethod bcgmethod peaksave nit ecgchan fqrsdet badchan};
options.val     = {qrsmethod bcgmethod peaksave ecgchan fqrsdet badchan};
options.help    = {'Options for pulse correction'};

% ---------------------------------------------------------------------
% pulseartfct Function to correct the pulse artefact
% ---------------------------------------------------------------------
pulseartfct         = cfg_exbranch;
pulseartfct.tag     = 'pulseartfct';
pulseartfct.name    = 'BCG artefact correction';
pulseartfct.val     = {data options};
pulseartfct.help    = {'Function to correct the pulse artefact, also named balisto-cardiogram artefact or cardio-balistic artefact, from EEG data acquired inside an MR scanner.'};
pulseartfct.prog    = @crc_run_par;
pulseartfct.vout    = @vout_par;

%% SUBFUNCTIONS
%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
function out = crc_run_par(job)
[out.D] = crc_par(job.data{1}, job.options);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

%------------------------------------------------------------------------
% dependency
%------------------------------------------------------------------------
function dep = vout_par(job)
% Output is always in field "D", no matter how job is structured
dep(1)            = cfg_dep;
dep(1).sname      = 'BCG-corrected Data';
% reference field "D" from output
dep(1).src_output = substruct('.','D');
% this can be entered into any evaluated input
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2)            = cfg_dep;
dep(2).sname      = 'BCG-corrected Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});