function crc_batch
% FASST function to prepare and launch the batch system.
% This builds the whole tree for the various tools and their GUI at the
% first call to this script.
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Christophe Phillips
% $Id$

persistent batch_initialize

if (isempty(batch_initialize) || ~batch_initialize) 
    % Whole initialisation of SPM batch system
    spm('defaults','eeg')
    spm_jobman('initcfg')
    batch_initialize = 1;
end
if batch_initialize~=2 || isempty(cfg_util('tag2cfg_id', 'fasst'))
    % FASST config tree setup
    fasst = crc_cfg_fasst;
    % Adding FASST config tree to the SPM tools
    cfg_util('addapp', fasst)
    % No need to do it again for this session
    batch_initialize = 2;
end

% Launching the batch system
cfg_ui

return

