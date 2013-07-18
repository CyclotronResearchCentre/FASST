function fasst = crc_cfg_fasst
% FASST Configuration file
% This builds the whole tree for the various tools and their GUI.
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Christophe Phillips
% $Id$

% if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','FASST')); end
% Only needed if
% - FASST is installed in SPM8 toolboxes directory, for example SPM8\toolbox\FASST
% - and this directory isn't saved on Matlab path

% ---------------------------------------------------------------------
% wavedetect FASST wave detection tools
% ---------------------------------------------------------------------
wavedetect        = cfg_choice;
wavedetect.tag    = 'wavedetect';
wavedetect.name   = 'Wave detection tools';
wavedetect.help   = {'Use FASST toolbox to automatically detect "waves".'}';
wavedetect.values = {crc_cfg_sws };

% ---------------------------------------------------------------------
% readdata FASST data importation tools
% ---------------------------------------------------------------------
readdata        = cfg_choice;
readdata.tag    = 'readdata';
readdata.name   = 'Read data into SPM8';
readdata.help   = {'Use FASST toolbox to import data'
                   'It reads the recording date and time and it maps better into memory'
                   }';
readdata.values = {crc_cfg_rdata_brpr crc_cfg_rdata_edf crc_cfg_rdata_egi};

% ---------------------------------------------------------------------
% fasst FASST series of tools
% ---------------------------------------------------------------------
fasst         = cfg_choice;
fasst.tag     = 'fasst';
fasst.name    = 'FASST Tools';
fasst.help    = {
                  'This is the batch interface for FASST, i.e. fMRI Artefact rejection and Sleep Scoring Toolbox,'
                  'providing a GUI for the data manipulation tools.'
                  }';
fasst.values  = {readdata;
                crc_cfg_chunk;
                crc_cfg_concatenate;
                crc_cfg_gar;
                crc_cfg_par;
                wavedetect};
%                 crc_cfg_sws;
%                 crc_cfg_feedback};
%------------------------------------------------------------------------


