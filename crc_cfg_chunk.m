function chunking = crc_cfg_chunk
% Chunking configuration file
% This builds the whole tree for the chunking too.
%_______________________________________________________________________
% Copyright (C) 2010 Cyclotron Research Centre

% Christophe Phillips
% $Id$

% ---------------------------------------------------------------------
% data Data file
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'M/EEG data file';
data.help    = {'This is the data file to chunk'};
data.filter  = 'mat';
data.ufilter = '.*';
data.num     = [1 1];

% ---------------------------------------------------------------------
% t_abs Time, beginning or end, specified as absolute time
% ---------------------------------------------------------------------
t_abs         = cfg_entry;
t_abs.tag     = 't_abs';
t_abs.name    = 'Absolute time';
t_abs.val     = {[]};
t_abs.help    = {'Time, beginning or end, specified as absolute time, i.e. real clock time.', ...
                 'Specify the time with 3 numbers: hour minute second'};
t_abs.strtype = 'e';
t_abs.num     = [1 3];

% ---------------------------------------------------------------------
% t_rel Time, beginning or end, specified as relative time
% ---------------------------------------------------------------------
t_rel         = cfg_entry;
t_rel.tag     = 't_rel';
t_rel.name    = 'Relative time';
t_rel.val     = {[]};
t_rel.help    = {'Time, beginning or end, specified as relative time, i.e. time elapsed since beginning of dataset.', ...
                 'Specify the time with 3 numbers: hour minute second'};
t_rel.strtype = 'e';
t_rel.num     = [1 3];

% ---------------------------------------------------------------------
% m_type Marker type to specify beginning or end of chunk
% ---------------------------------------------------------------------
m_type         = cfg_entry;
m_type.tag     = 'm_type';
m_type.name    = 'Marker type';
m_type.val     = {[]};
m_type.help    = {'Marker type to specify beginning or end of chunk.'};
m_type.strtype = 'e';
m_type.num     = [1 1];

% ---------------------------------------------------------------------
% m_ind Marker index to specify beginning or end of chunk
% ---------------------------------------------------------------------
m_ind         = cfg_entry;
m_ind.tag     = 'm_ind';
m_ind.name    = 'Marker index';
m_ind.val     = {[]};
m_ind.help    = {'Marker index to specify beginning or end of chunk.'};
m_ind.strtype = 'e';
m_ind.num     = [1 1];

% ---------------------------------------------------------------------
% chunk Chunk
% ---------------------------------------------------------------------
t_mark         = cfg_branch;
t_mark.tag     = 't_mark';
t_mark.name    = 'Marker time';
t_mark.val     = {m_type m_ind };
t_mark.help    = {'Specify beginning and end of chunk by a marker.'};


% ---------------------------------------------------------------------
% chunk_beg Beginning of chunk
% ---------------------------------------------------------------------
chunk_beg         = cfg_choice;
chunk_beg.tag     = 'chunk_beg';
chunk_beg.name    = 'Beginning of chunk';
chunk_beg.val     = {t_abs};
chunk_beg.help    = {'Specify how the beginning of the chunk is defined: specific marker, relative time, or absolute time'};
chunk_beg.values  = {t_mark t_rel t_abs};

% ---------------------------------------------------------------------
% chunk_end End of chunk
% ---------------------------------------------------------------------
chunk_end         = cfg_choice;
chunk_end.tag     = 'chunk_end';
chunk_end.name    = 'End of chunk';
chunk_end.val     = {t_abs};
chunk_end.help    = {'Specify how the end of the chunk is defined: specific event, relative time, or absolute time.'};
chunk_end.values  = {t_mark t_rel t_abs};

% ---------------------------------------------------------------------
% chunk Chunk
% ---------------------------------------------------------------------
chunk         = cfg_branch;
chunk.tag     = 'chunk';
chunk.name    = 'Chunk';
chunk.val     = {chunk_beg chunk_end };
chunk.help    = {'Specify beginning and end of chunk.'};

% ---------------------------------------------------------------------
% generic Chunks
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Chunks';
generic.help    = {'Add chunks to be extracted.'};
generic.values  = {chunk };
generic.num     = [1 Inf];

% ---------------------------------------------------------------------
% fn_prefix Prefix used for chunked file name
% ---------------------------------------------------------------------
fn_prefix         = cfg_entry;
fn_prefix.tag     = 'fn_prefix';
fn_prefix.name    = 'Filename prefix';
fn_prefix.help    = {'Prefix used for the chunked file, and index is further added to differentiate different chunks.'};
fn_prefix.strtype = 's';
fn_prefix.num     = [1 Inf];
fn_prefix.def     = @(val)crc_get_defaults('chk.prefix', val{:});

% ---------------------------------------------------------------------
% overwr Overwrite or not the chunked file
% ---------------------------------------------------------------------
overwr         = cfg_menu;
overwr.tag     = 'overwr';
overwr.name    = 'Overwrite';
overwr.help    = {'Overwriting or not previously extracted chunks'};
overwr.labels  = {'Yes'
                  'No'}';
overwr.values  = {1 0};
overwr.def     = @(val)crc_get_defaults('chk.overwr', val{:});

% ---------------------------------------------------------------------
% numchunk Marker type to specify beginning or end of chunk
% ---------------------------------------------------------------------
numchunk         = cfg_entry;
numchunk.tag     = 'numchunk';
numchunk.name    = 'First chunk index';
numchunk.def     = @(val)crc_get_defaults('chk.numchunk', val{:});
numchunk.help    = {'First chunk index.'};
numchunk.strtype = 'e';
numchunk.num     = [1 1];

% ---------------------------------------------------------------------
% options List of option for chunking
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Chunk writing option';
options.val     = {overwr fn_prefix numchunk};
options.help    = {'Options for chunk writing.'};

% ---------------------------------------------------------------------
% chunking Chunking one file in one or more chunks
% ---------------------------------------------------------------------
chunking         = cfg_exbranch;
chunking.tag     = 'chunking';
chunking.name    = 'Chunking one file';
chunking.val     = {data generic options};
chunking.help    = {'Specify a M/EEG data file, then the chunks to extract.'};
chunking.prog = @crc_run_chunking;
chunking.vout = @vout_chunking;

%------------------------------------------------------------------------
function dep = vout_chunking(job)
for kk=1:numel(job.chunk)
    cdep            = cfg_dep;
    cdep.sname      = sprintf('Chunked file, #%d', kk);
    cdep.src_output = substruct('.','chunk', '()',{kk}, '.','cfiles');
    cdep.tgt_spec   = cfg_findspec({{ ...
                    'filter',['^',job.options.fn_prefix,'.*\.mat$'], ...
                    'strtype','e'}});
    if kk == 1
        dep = cdep;
    else
        dep = [dep cdep];
    end
end
%------------------------------------------------------------------------
