function Do = crc_convert20to30s(D, flags)
% Do = crc_convert20to30s(D)
%
% Looks for sleep scores in windows of 20 seconds and converts them into
% windows of 30 seconds. Quite simply it takes the windows [1 3 4 6 7 9...]
% from the list of scores.
% A "new" scorer is created for each 20s-window score, with '_c30s'
% appended at the end of the scorer's name.
% The 
%
% INPUT
%   D     : (array of) filename(s) or (cell array of) SPM/meeg object(s) 
%   flags : structure with options
%       .save    : save or not the output updatred structure [1, def]
%
% OUTPUT
%   Do       : updated meeg object(s) (also saved on disk)
%__________________________________________________________________________
% Copyright (C) 2014 Cyclotron Research Centre

% Written by C. Phillips, 2014.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

if ~nargin
    % no argument
    D = spm_select(Inf,'mat','Select MEEG file(s)');
end
Do = D;

def_flags = struct('save', true);
if nargin<2
    flags = def_flags;
else
    flags = crc_check_flag(def_flags,flags);
end

if ischar(D)
    P = D;
    Ndat = size(P,1);
    for ii=1:Ndat
        dat{ii} = crc_eeg_load(deblank(P(ii,:))); %#ok<*AGROW>
    end
elseif isa(D,'meeg')
    Ndat = 1;
    dat{1} = D;
elseif iscell(D)
    Ndat = numel(D);
    dat = D;
else
    error('fasst:conv20to30sc','Unknown data format');
end

clear D;
for ii=1:Ndat
    % Loop over all data sets
    D = dat{ii};
    if isfield(D,'CRC') && isfield(D.CRC,'score')
        sc = D.CRC.score;
        Nsc = size(sc,2);
        Nconv = 0; 
        for jj=1:Nsc
            % Loop around scorer
            if sc{3,jj}==20
                % check windows size & convert if 20s
                Nconv = Nconv+1; % conversion counter
                sc(:,Nsc+Nconv) = sc(:,jj); % copy at the end
                lsc = 1:numel(sc{1,Nsc+Nconv});
                lsc(2:3:end) = [];
                sc{1,Nsc+Nconv} = sc{1,Nsc+Nconv}(lsc);
                sc{2,Nsc+Nconv} = [sc{2,Nsc+Nconv},'_c30s'];
                sc{3,Nsc+Nconv} = 30;
            end
        end
        fprintf('\n%d scores converted in file %s\n',Nconv,D.fname);
        D.CRC.score = sc;
        if flags.save
            save(D)
        end
    else
        fprintf('\nNo score in file %s\n',D.fname);
    end
end
Do = D;
