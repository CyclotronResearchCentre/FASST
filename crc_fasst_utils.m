function varargout = crc_fasst_utils(varargin)
% Function used for a series of 'utilities'.
%
% FORMAT mdir = crc_fasst_utils('Dir',Mfile)
% Returns the current directory of FASST, if Mfile is not defined, or the
% path to Mfile.
%
% FORMAT [v,r] = crc_fasst_utils('Ver',Mfile,ReDo)
% Returns the current version (v) and release (r) of file Mfile. This
% corresponds to the last changed Revision number extracted from the
% Subversion Id tag.
% If Mfile is absent or empty then it returns the current FASST version (v)
% and release (r), extracted from the file Contents.m in the FASST
% directory (these information are cached in a persistent variable to
% enable repeat use without recomputation).
% If Redo [default false] is true, then the cached current FASST
% information are not used but recomputed (and recached).
%_______________________________________________________________________
% Copyright (C) 2012 Cyclotron Research Centre

% Written by C. Phillips, 2012.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% NOTE:
% This was largely inspired by the routine spm.m from the SPM software.
% SPM is developed by members and collaborators of the
% Wellcome Trust Centre for Neuroimaging
%-----------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging


%-Format arguments
%-----------------------------------------------------------------------
if nargin == 0,
    Action = 'XXXX';
else
    Action = varargin{1};
end

%=======================================================================
switch lower(Action)

    case 'dir'                        %-Identify specific (FASST) directory
        %==================================================================
        % crc_fasst_utils('Dir',Mfile)
        %------------------------------------------------------------------
        if nargin<2,
            Mfile = 'crc_main';
        else
            Mfile = varargin{2};
        end
        FASSTdir = which(Mfile);
        if isempty(FASSTdir)             %-Not found or full pathname given
            if exist(Mfile,'file')==2    %-Full pathname
                FASSTdir = Mfile;
            else
                error(['Can''t find ',Mfile,' on MATLABPATH']);
            end
        end
        FASSTdir  = fileparts(FASSTdir);
        varargout = {FASSTdir};
        
    case 'ver'                                              %-FASST version
        %==================================================================
        % [ver, rel] = crc_fasst_utils('Ver',Mfile,ReDo)
        %------------------------------------------------------------------
        if nargin ~= 3,
            ReDo = false;
        else
            ReDo = logical(varargin{3});
        end
        if nargin == 1 || (nargin > 1 && isempty(varargin{2}))
            Mfile = '';
        else
            Mfile = which(varargin{2});
            if isempty(Mfile)
                error('Can''t find %s on MATLABPATH.',varargin{2});
            end
        end
        
        v = crc_version(ReDo);
        
        if isempty(Mfile)
            varargout = {v.Release v.Version};
        else
            unknown = struct('file',Mfile,'id','???','date','','author','');
            fp  = fopen(Mfile,'rt');
            if fp == -1, error('Can''t read %s.',Mfile); end
            str = fread(fp,Inf,'*uchar');
            fclose(fp);
            str = char(str(:)');
            r = regexp(str,['\$Id: (?<file>\S+) (?<id>[0-9]+) (?<date>\S+) ' ...
                '(\S+Z) (?<author>\S+) \$'],'names','once');
            if isempty(r), r = unknown; end
            varargout = {r(1).id v.Release};
        end
        
        %==================================================================
    otherwise % unknown action
        %==================================================================
        return
end

return

%=======================================================================
%=======================================================================
% SUBFUNCTIONS
%=======================================================================
%=======================================================================

%=======================================================================
function v = crc_version(ReDo)                  %-Retrieve FASST version
%=======================================================================
persistent FASST_VER;
v = FASST_VER;
if isempty(FASST_VER) || (nargin > 0 && ReDo)
    v = struct('Name','','Version','','Release','','Date','');
    try
        vfile = fullfile(crc_fasst_utils('Dir'),'Contents.m');
        fid = fopen(vfile,'rt');
        if fid == -1, error(str); end
        l1 = fgetl(fid); l2 = fgetl(fid);
        fclose(fid);
        l1 = strtrim(l1(2:end)); l2 = strtrim(l2(2:end));
        t  = textscan(l2,'%s','delimiter',' '); t = t{1};
        v.Name = l1; v.Date = t{4};
        v.Version = t{2}; v.Release = t{3}(2:end-1);
    catch %#ok<CTCH>
        error('FASST:getversion', ...
            'Can''t obtain FASST Revision information.');
    end
    FASST_VER = v;
end
